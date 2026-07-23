#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/shapegroup.h>
#include <mitsuba/render/kdtree.h>
#include <drjit/tensor.h>

#if defined(MI_ENABLE_EMBREE)
#  include <embree3/rtcore.h>
#endif

NAMESPACE_BEGIN(mitsuba)

/* InstanceList below instantiates TShapeKDTree with its own private
   InstanceKDTree as the CRTP Derived type. TShapeKDTree::BuildTask's
   thread_local build context is only ever defined (in kdtree.cpp) for the
   ShapeKDTree specialization, so a distinct definition is needed here for
   our own specialization to link. */
template <typename B, typename I, typename C, typename D>
thread_local typename TShapeKDTree<B, I, C, D>::LocalBuildContext
    TShapeKDTree<B, I, C, D>::BuildTask::m_local = {};

/**!

.. _shape-instancelist:

InstanceList (:monosp:`instancelist`)
-------------------------------------------------

.. pluginparameters::

 * - (Nested plugin)
   - :paramtype:`shapegroup`
   - A reference to a shape group that should be instantiated.

 * - transforms
   - |tensor|
   - Specifies a buffer of object-to-world transformations of shape [N,4,4] with
     N being the number of instances.

This is a utility plugin that allows the creation of multiple instances of a
shape group from a buffer of object-to-world transforms. Unlike the
:ref:`shape-instance` plugin used N times, this plugin stores a single copy
of the transform buffer and implements ray tracing directly against all N
instances, which avoids the memory overhead of creating N separate instance
objects. This makes it suitable for scenes with millions of instances (e.g.
forests). For details on how to create instances, refer to the
``shapegroup`` and ``instance`` plugins from Mitsuba's core plugin library.

.. note::

    On the GPU (CUDA) variants, this plugin currently falls back to
    expanding into N :ref:`shape-instance` objects (as OptiX support for the
    memory-efficient path has not been implemented yet).
 */

template <typename Float, typename Spectrum>
class InstanceList final : public Shape<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Shape, id)
    MI_IMPORT_TYPES()

    using ShapeGroup_ = ShapeGroup<Float, Spectrum>;
    using typename Base::ScalarSize;
    using typename Base::ScalarIndex;

    using JitIndex = dr::uint32_array_t<Float>;

    /* The kd-tree/embree leaf test can only report two indices per hit
       (`shape_index`, `prim_index`), but we need three: which instance,
       which child shape within the shapegroup, and which primitive of that
       child shape. We pack (instance_index, child_shape_index) into a
       single `shape_index` field and leave `prim_index` untouched. */
    static constexpr uint32_t ShapeIndexBits  = 8;
    static constexpr uint32_t ShapeIndexMask  = (1u << ShapeIndexBits) - 1;
    static constexpr uint32_t MaxInstanceIndex =
        (1u << (32 - ShapeIndexBits)) - 1;

    /// Result of testing a ray against a single instance
    struct Hit {
        ScalarFloat t              = dr::Infinity<ScalarFloat>;
        ScalarPoint2f uv           = ScalarPoint2f(0.f);
        ScalarUInt32 instance_index = (ScalarUInt32) -1;
        ScalarUInt32 shape_index    = (ScalarUInt32) -1;
        ScalarUInt32 prim_index     = (ScalarUInt32) -1;

        bool is_valid() const { return t != dr::Infinity<ScalarFloat>; }
    };

    /**
     * \brief Private acceleration structure used to cull the N instances of
     * this shape. Reuses the generic SAH kd-tree builder from \ref
     * TShapeKDTree, but with a custom leaf test that delegates to the
     * shared shapegroup instead of expecting a list of \ref Shape pointers.
     */
    class InstanceKDTree
        : public TShapeKDTree<ScalarBoundingBox3f, uint32_t,
                              SurfaceAreaHeuristic3<ScalarFloat>, InstanceKDTree> {
    public:
        using Base = TShapeKDTree<ScalarBoundingBox3f, uint32_t,
                                  SurfaceAreaHeuristic3<ScalarFloat>, InstanceKDTree>;
        using typename Base::Size;
        using typename Base::Index;
        using typename Base::KDNode;
        using Base::m_bbox;
        using Base::m_nodes;
        using Base::m_indices;
        using ScalarRay3f = Ray<ScalarPoint3f, Spectrum>;

        InstanceKDTree(const InstanceList *owner)
            : Base(SurfaceAreaHeuristic3<ScalarFloat>(20.f, 15.f, 0.9f)),
              m_owner(owner) {}

        void build_tree() {
            for (ScalarIndex i = 0; i < m_owner->m_instance_count; ++i)
                m_bbox.expand(bbox(i));
            Base::build();
        }

        Size primitive_count() const { return m_owner->m_instance_count; }

        ScalarBoundingBox3f bbox(Index i) const {
            return m_owner->instance_bbox(i);
        }

        ScalarBoundingBox3f bbox(Index i, const ScalarBoundingBox3f &clip) const {
            ScalarBoundingBox3f result = bbox(i);
            result.clip(clip);
            return result;
        }

        template <bool ShadowRay>
        Hit ray_intersect_scalar(ScalarRay3f ray) const {
            struct KDStackEntry {
                ScalarFloat mint, maxt;
                const KDNode *node;
            };
            KDStackEntry stack[MI_KD_MAXDEPTH];
            int32_t stack_index = 0;

            Hit result;

            auto bbox_result = m_bbox.ray_intersect(ray);
            ScalarFloat mint = std::max(ScalarFloat(0), std::get<1>(bbox_result)),
                        maxt = std::min(ray.maxt, std::get<2>(bbox_result));

            ScalarVector3f d_rcp = dr::rcp(ray.d);

            const KDNode *node = m_nodes.get();
            while (mint <= maxt) {
                if (likely(!node->leaf())) { // Inner node
                    const ScalarFloat split = node->split();
                    const uint32_t axis     = node->axis();

                    ScalarFloat t_plane = (split - ray.o[axis]) * d_rcp[axis];

                    bool left_first  = (ray.o[axis] < split) ||
                                       (ray.o[axis] == split && ray.d[axis] >= 0.f),
                         start_after = t_plane<mint, end_before = t_plane> maxt ||
                                       t_plane < 0.f || !dr::isfinite(t_plane),
                         single_node = start_after || end_before;

                    if (likely(single_node)) {
                        bool visit_left = end_before == left_first;
                        node = node->left() + (visit_left ? 0 : 1);
                        continue;
                    }

                    Index node_offset = left_first ? 0 : 1;
                    const KDNode *left   = node->left(),
                                 *n_cur  = left + node_offset,
                                 *n_next = left + (1 - node_offset);

                    KDStackEntry& entry = stack[stack_index++];
                    entry.mint = t_plane;
                    entry.maxt = maxt;
                    entry.node = n_next;

                    node = n_cur;
                    maxt = t_plane;
                    continue;
                } else if (node->primitive_count() > 0) { // Leaf node
                    Index prim_start = node->primitive_offset();
                    Index prim_end = prim_start + node->primitive_count();
                    for (Index i = prim_start; i < prim_end; i++) {
                        Index instance_index = m_indices[i];

                        Hit hit = m_owner->template test_instance<ShadowRay>(
                            instance_index, ray);

                        if (unlikely(hit.is_valid())) {
                            if constexpr (ShadowRay)
                                return hit;
                            result = hit;
                            ray.maxt = result.t;
                        }
                    }
                }

                if (likely(stack_index > 0)) {
                    --stack_index;
                    KDStackEntry& entry = stack[stack_index];
                    mint = entry.mint;
                    maxt = std::min(entry.maxt, ray.maxt);
                    node = entry.node;
                } else {
                    break;
                }
            }

            return result;
        }

        MI_DECLARE_CLASS(InstanceKDTree)
    private:
        const InstanceList *m_owner;
    };

    InstanceList(const Properties &props) : Base(props) {
        for (auto &prop : props.objects()) {
            if (Base *shape = prop.try_get<Base>()) {
                if (shape->is_shape_group()) {
                    if (m_shapegroup)
                        Throw("Only a single shapegroup can be specified per instance list.");
                    m_shapegroup = (ShapeGroup_ *) shape;
                }
            }
        }

        if (!m_shapegroup)
            Throw("A reference to a 'shapegroup' must be specified!");

        m_tensor = props.get_any<TensorXf>("transforms");
        if (m_tensor.ndim() != 3)
            Throw("transforms tensor has %zu dimensions, expected 3", m_tensor.ndim());
        if (m_tensor.shape(1) != 4 || m_tensor.shape(2) != 4)
            Throw("transforms tensor must have shape [N, 4, 4], got [%zu, %zu, %zu]",
                  m_tensor.shape(0), m_tensor.shape(1), m_tensor.shape(2));

        m_instance_count = (uint32_t) m_tensor.shape(0);
        if (m_instance_count > MaxInstanceIndex)
            Throw("Too many instances (%u); instancelist supports at most %u.",
                  m_instance_count, MaxInstanceIndex);

        for (ScalarIndex i = 0; i < m_instance_count; ++i)
            m_bbox.expand(instance_bbox(i));

        if constexpr (!dr::is_cuda_v<Float>) {
            m_tree = std::make_unique<InstanceKDTree>(this);
            m_tree->build_tree();
        }
    }

    std::vector<ref<Object>> expand() const override {
        std::vector<ref<Object>> shapes;

        // The memory-efficient path below only supports CPU (scalar/LLVM)
        // variants. On CUDA, fall back to the original N-object expansion
        // so existing OptiX support keeps working unmodified.
        if constexpr (!dr::is_cuda_v<Float>)
            return shapes;

        auto pmgr = PluginManager::instance();

        std::string id_prefix = std::string(id());
        ScalarFloat ilog10    = dr::rcp(dr::log(10.f));
        int32_t max_zeros =
            (int32_t) (dr::log(ScalarFloat(m_instance_count)) * ilog10) + 1;

        for (uint32_t i = 0; i < m_instance_count; ++i) {
            ScalarAffineTransform4f to_world(
                dr::gather<ScalarMatrix4f>(m_tensor.array(), i));
            Properties props_instance("instance");
            props_instance.set("shapegroup", (Object *) m_shapegroup.get());
            props_instance.set("to_world", to_world);

            ref<Object> shape =
                (Object *) pmgr->create_object<Base>(props_instance);

            if (id_prefix != "") {
                std::string number = std::to_string(i);
                number.insert(number.begin(), max_zeros - number.length(), '0');
                std::string id = id_prefix + "_" + number;
                shape->set_id(id);
            }

            shapes.push_back(shape);
        }
        return shapes;
    }

    ScalarBoundingBox3f bbox() const override { return m_bbox; }

    ScalarSize primitive_count() const override { return 1; }

    ScalarSize effective_primitive_count() const override {
        return m_instance_count * m_shapegroup->primitive_count();
    }

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Ray tracing routines
    // =============================================================

    template <typename FloatP, typename Ray3fP>
    std::tuple<FloatP, Point<FloatP, 2>, dr::uint32_array_t<FloatP>,
               dr::uint32_array_t<FloatP>>
    ray_intersect_preliminary_impl(const Ray3fP &ray, ScalarIndex /*prim_index*/,
                                   dr::mask_t<FloatP> active) const {
        MI_MASK_ARGUMENT(active);

        if constexpr (dr::is_cuda_v<Float>) {
            DRJIT_MARK_USED(ray);
            DRJIT_MARK_USED(active);
            Throw("ray_intersect_preliminary_impl() should not be reached on CUDA "
                  "variants (instancelist falls back to expand() there).");
        } else if constexpr (!dr::is_array_v<FloatP>) {
            Hit hit = m_tree->template ray_intersect_scalar<false>(ray);
            uint32_t packed = hit.is_valid()
                ? (uint32_t) ((hit.instance_index << ShapeIndexBits) | pack_shape_index(hit.shape_index))
                : (uint32_t) -1;
            return { hit.t, hit.uv, packed, hit.prim_index };
        } else if constexpr (dr::is_jit_v<FloatP>) {
            // Only reached via the generic (non-scalar, non-packet) entry
            // point, which the Scene never actually calls for this shape
            // (native backend drives per-ray via ray_intersect_preliminary_scalar,
            // embree drives fixed-width packets via ray_intersect_preliminary_packet).
            DRJIT_MARK_USED(ray);
            DRJIT_MARK_USED(active);
            Throw("ray_intersect_preliminary_impl(): dynamic-width queries are "
                  "not supported by instancelist; use scalar or packet queries.");
        } else {
            constexpr size_t N = FloatP::Size;
            using UInt32P = dr::uint32_array_t<FloatP>;

            FloatP t           = dr::Infinity<Float>;
            Point<FloatP, 2> uv = dr::zeros<Point<FloatP, 2>>();
            UInt32P shape_idx   = dr::full<UInt32P>((uint32_t) -1);
            UInt32P prim_idx    = dr::full<UInt32P>((uint32_t) -1);

            for (size_t i = 0; i < N; ++i) {
                if (!active.entry(i))
                    continue;

                ScalarRay3f ray_i(ScalarPoint3f(ray.o.x().entry(i), ray.o.y().entry(i), ray.o.z().entry(i)),
                                  ScalarVector3f(ray.d.x().entry(i), ray.d.y().entry(i), ray.d.z().entry(i)),
                                  (ScalarFloat) ray.maxt.entry(i), (ScalarFloat) ray.time.entry(i),
                                  dr::zeros<Wavelength>());

                Hit hit = m_tree->template ray_intersect_scalar<false>(ray_i);
                if (hit.is_valid()) {
                    t.entry(i)         = hit.t;
                    uv.x().entry(i)    = hit.uv.x();
                    uv.y().entry(i)    = hit.uv.y();
                    shape_idx.entry(i) = (uint32_t) ((hit.instance_index << ShapeIndexBits) |
                                                     pack_shape_index(hit.shape_index));
                    prim_idx.entry(i)  = hit.prim_index;
                }
            }

            return { t, uv, shape_idx, prim_idx };
        }
    }

    template <typename FloatP, typename Ray3fP>
    dr::mask_t<FloatP> ray_test_impl(const Ray3fP &ray, ScalarIndex /*prim_index*/,
                                     dr::mask_t<FloatP> active) const {
        MI_MASK_ARGUMENT(active);

        if constexpr (dr::is_cuda_v<Float>) {
            DRJIT_MARK_USED(ray);
            DRJIT_MARK_USED(active);
            Throw("ray_test_impl() should not be reached on CUDA variants "
                  "(instancelist falls back to expand() there).");
        } else if constexpr (!dr::is_array_v<FloatP>) {
            return m_tree->template ray_intersect_scalar<true>(ray).is_valid();
        } else if constexpr (dr::is_jit_v<FloatP>) {
            DRJIT_MARK_USED(ray);
            DRJIT_MARK_USED(active);
            Throw("ray_test_impl(): dynamic-width queries are not supported by "
                  "instancelist; use scalar or packet queries.");
        } else {
            constexpr size_t N = FloatP::Size;
            dr::mask_t<FloatP> result = false;

            for (size_t i = 0; i < N; ++i) {
                if (!active.entry(i))
                    continue;

                ScalarRay3f ray_i(ScalarPoint3f(ray.o.x().entry(i), ray.o.y().entry(i), ray.o.z().entry(i)),
                                  ScalarVector3f(ray.d.x().entry(i), ray.d.y().entry(i), ray.d.z().entry(i)),
                                  (ScalarFloat) ray.maxt.entry(i), (ScalarFloat) ray.time.entry(i),
                                  dr::zeros<Wavelength>());

                result.entry(i) = m_tree->template ray_intersect_scalar<true>(ray_i).is_valid();
            }

            return result;
        }
    }

    MI_SHAPE_DEFINE_RAY_INTERSECT_METHODS()

    SurfaceInteraction3f compute_surface_interaction(const Ray3f &ray,
                                                     const PreliminaryIntersection3f &pi,
                                                     uint32_t ray_flags,
                                                     uint32_t recursion_depth,
                                                     Mask active) const override {
        MI_MASK_ARGUMENT(active);

        // Nested instancing is not supported
        if (recursion_depth > 0)
            return dr::zeros<SurfaceInteraction3f>();

        JitIndex instance_index = pi.shape_index >> ShapeIndexBits;
        UInt32 child_shape_index = pi.shape_index & ShapeIndexMask;

        Matrix4f coefficients = dr::gather<Matrix4f>(m_tensor.array(), instance_index, active);
        AffineTransform4f to_world(coefficients);
        AffineTransform4f to_object = to_world.inverse();

        PreliminaryIntersection3f child_pi = pi;
        child_pi.shape_index = child_shape_index;

        SurfaceInteraction3f si = m_shapegroup->compute_surface_interaction(
            to_object * ray, child_pi, ray_flags, recursion_depth, active);

        si.p = to_world * si.p;
        si.n = dr::normalize(to_world * si.n);
        if (likely(has_flag(ray_flags, RayFlags::ShadingFrame)))
            si.sh_frame.n = dr::normalize(to_world * si.sh_frame.n);

        if (likely(has_flag(ray_flags, RayFlags::ShadingFrame)))
            si.initialize_sh_frame();

        if (likely(has_flag(ray_flags, RayFlags::dPdUV))) {
            si.dp_du = to_world * si.dp_du;
            si.dp_dv = to_world * si.dp_dv;
        }

        if (has_flag(ray_flags, RayFlags::dNGdUV) || has_flag(ray_flags, RayFlags::dNSdUV)) {
            Normal3f n = has_flag(ray_flags, RayFlags::dNGdUV) ? si.n : si.sh_frame.n;

            Normal3f tn = to_world * dr::normalize(to_object * n);
            Float inv_len = dr::rcp(dr::norm(tn));
            tn *= inv_len;

            si.dn_du = to_world * Normal3f(si.dn_du) * inv_len;
            si.dn_dv = to_world * Normal3f(si.dn_dv) * inv_len;

            si.dn_du -= tn * dr::dot(tn, si.dn_du);
            si.dn_dv -= tn * dr::dot(tn, si.dn_dv);
        }

        si.prim_index = pi.prim_index;
        si.instance = this;

        return si;
    }

    //! @}
    // =============================================================

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "InstanceList[" << std::endl
            << "  shapegroup = " << string::indent(m_shapegroup) << std::endl
            << "  instance_count = " << m_instance_count << "," << std::endl
            << "]";
        return oss.str();
    }

#if defined(MI_ENABLE_EMBREE)
    RTCGeometry embree_geometry(RTCDevice device) override {
        if constexpr (!dr::is_cuda_v<Float>)
            m_embree_child_scene = m_shapegroup->embree_scene(device);
        return Base::embree_geometry(device);
    }
#endif

    MI_DECLARE_CLASS(InstanceList)
private:
    /// Bounding box of a single instance, in world space
    ScalarBoundingBox3f instance_bbox(ScalarIndex i) const {
        const ScalarBoundingBox3f &group_bbox = m_shapegroup->bbox();
        if (!group_bbox.valid())
            return group_bbox;

        ScalarAffineTransform4f to_world(dr::gather<ScalarMatrix4f>(m_tensor.array(), i));
        ScalarBoundingBox3f result;
        for (int c = 0; c < 8; ++c)
            result.expand(to_world * group_bbox.corner(c));
        return result;
    }

    /// Local shapegroup shape_index values must fit in ShapeIndexBits
    static uint32_t pack_shape_index(uint32_t shape_index) {
        if (unlikely(shape_index > ShapeIndexMask))
            Throw("instancelist: shapegroup has too many child shapes (%u) to "
                  "pack into the instance hit encoding (max %u).",
                  shape_index, ShapeIndexMask);
        return shape_index;
    }

    /// Test a single instance for intersection/occlusion with a scalar ray
    template <bool ShadowRay>
    Hit test_instance(ScalarIndex instance_index, const ScalarRay3f &ray) const {
        ScalarAffineTransform4f to_world(
            dr::gather<ScalarMatrix4f>(m_tensor.array(), instance_index));
        ScalarRay3f local_ray = to_world.inverse() * ray;

        Hit hit;

#if defined(MI_ENABLE_EMBREE)
        if constexpr (!dr::is_cuda_v<Float>) {
            RTCIntersectContext context;
            rtcInitIntersectContext(&context);

            if constexpr (ShadowRay) {
                RTCRay rtc_ray;
                rtc_ray.org_x = (float) local_ray.o.x();
                rtc_ray.org_y = (float) local_ray.o.y();
                rtc_ray.org_z = (float) local_ray.o.z();
                rtc_ray.dir_x = (float) local_ray.d.x();
                rtc_ray.dir_y = (float) local_ray.d.y();
                rtc_ray.dir_z = (float) local_ray.d.z();
                rtc_ray.tnear = 0.f;
                rtc_ray.tfar  = (float) local_ray.maxt;
                rtc_ray.time  = (float) local_ray.time;
                rtc_ray.mask  = (unsigned int) -1;
                rtc_ray.id    = 0;
                rtc_ray.flags = 0;

                rtcOccluded1(m_embree_child_scene, &context, &rtc_ray);
                if (rtc_ray.tfar != (float) local_ray.maxt)
                    hit.t = 0.f;
            } else {
                RTCRayHit rh;
                rh.ray.org_x = (float) local_ray.o.x();
                rh.ray.org_y = (float) local_ray.o.y();
                rh.ray.org_z = (float) local_ray.o.z();
                rh.ray.dir_x = (float) local_ray.d.x();
                rh.ray.dir_y = (float) local_ray.d.y();
                rh.ray.dir_z = (float) local_ray.d.z();
                rh.ray.tnear = 0.f;
                rh.ray.tfar  = (float) local_ray.maxt;
                rh.ray.time  = (float) local_ray.time;
                rh.ray.mask  = (unsigned int) -1;
                rh.ray.id    = 0;
                rh.ray.flags = 0;
                rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;

                rtcIntersect1(m_embree_child_scene, &context, &rh);

                if (rh.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                    hit.t              = rh.ray.tfar;
                    hit.uv             = ScalarPoint2f(rh.hit.u, rh.hit.v);
                    hit.instance_index = instance_index;
                    hit.shape_index    = rh.hit.geomID;
                    hit.prim_index     = rh.hit.primID;
                }
            }
        }
#else
        if constexpr (ShadowRay) {
            if (m_shapegroup->ray_test_scalar(local_ray))
                hit.t = 0.f;
        } else {
            auto [t, uv, shape_index, prim_index] =
                m_shapegroup->ray_intersect_preliminary_scalar(local_ray);
            if (t != dr::Infinity<ScalarFloat>) {
                hit.t              = t;
                hit.uv              = uv;
                hit.instance_index  = instance_index;
                hit.shape_index     = shape_index;
                hit.prim_index      = prim_index;
            }
        }
#endif
        return hit;
    }

    TensorXf m_tensor;
    uint32_t m_instance_count = 0;
    ref<ShapeGroup_> m_shapegroup;
    ScalarBoundingBox3f m_bbox;
    std::unique_ptr<InstanceKDTree> m_tree;

#if defined(MI_ENABLE_EMBREE)
    RTCScene m_embree_child_scene = nullptr;
#endif
};

MI_EXPORT_PLUGIN(InstanceList)
NAMESPACE_END(mitsuba)
