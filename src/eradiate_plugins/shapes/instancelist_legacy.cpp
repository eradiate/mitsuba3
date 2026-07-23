#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/shapegroup.h>
#include <mitsuba/core/transform.h>
#include <drjit/tensor.h>

NAMESPACE_BEGIN(mitsuba)
/**!

.. _shape-instancelist_legacy:

InstanceList legacy baseline (:monosp:`instancelist_legacy`)
--------------------------------------------------------------

.. pluginparameters::

 * - (Nested plugin)
   - :paramtype:`shapegroup`
   - A reference to a shape group that should be instantiated.

 * - transforms
   - |tensor|
   - Specifies a buffer of object-to-world transformations of shape [N,4,4] with
     N being the number of instances.
   - |exposed|

This is the original :monosp:`expand()`-based implementation of
:ref:`shape-instancelist`, kept around under a separate plugin name purely as
a baseline for memory/performance A/B testing against the rewritten,
non-expanding :monosp:`instancelist` plugin. It expands the transform tensor
into N separate :ref:`shape-instance` objects at load time, which is exactly
the behavior the rewrite was meant to avoid. Prefer :monosp:`instancelist`
for actual scenes; this plugin has no other reason to exist.
 */

template <typename Float, typename Spectrum>
class InstanceListLegacy final : public Shape<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Shape, id)
    MI_IMPORT_TYPES()

    using ShapeGroup_ = ShapeGroup<Float, Spectrum>;

    InstanceListLegacy(const Properties &props) : Base(props) {
        m_shapegroup = nullptr;
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

        m_tensor =
            const_cast<TensorXf *>(&props.get_any<TensorXf>("transforms"));
        if (m_tensor->ndim() != 3)
            Throw("Tensor->has %ul dimensions. Expected 3", m_tensor->ndim());

        m_instance_count = (uint32_t) m_tensor->shape(0);
    }

    std::vector<ref<Object>> expand() const override {
        std::vector<ref<Object>> shapes;
        auto pmgr = PluginManager::instance();

        std::string id_prefix = std::string(id());
        ScalarFloat ilog10    = dr::rcp(dr::log(10.f));
        int32_t max_zeros =
            (int32_t) (dr::log(ScalarFloat(m_instance_count)) * ilog10) + 1;

        for (uint32_t i = 0; i < m_instance_count; ++i) {
            ScalarAffineTransform4f to_world(
                dr::gather<ScalarMatrix4f>(m_tensor->array(), i));
            Properties props_instance("instance");
            props_instance.set("shapegroup", (Object *) m_shapegroup);
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

    ScalarBoundingBox3f bbox() const override {
        return dr::zeros<ScalarBoundingBox3f>();
    }

    MI_DECLARE_CLASS(InstanceListLegacy)
private:
    uint32_t m_instance_count;
    TensorXf *m_tensor;
    ShapeGroup_ *m_shapegroup;
};

MI_EXPORT_PLUGIN(InstanceListLegacy)
NAMESPACE_END(mitsuba)
