#include <mitsuba/core/properties.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/eradiate/extremum.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/volumegrid.h>
#include <nanothread/nanothread.h>

NAMESPACE_BEGIN(mitsuba)

enum class SphericalTraversalType { RadialOnly, Full3D };

template <typename Float, typename Spectrum, SphericalTraversalType TraversalType>
class ExtremumSphericalImpl;

/**!
.. _extremum-extremum_spherical:

Extremum spherical structure (:monosp:`extremum_spherical`)
-----------------------------------------------------------

.. pluginparameters::
 * - resolution
   - |vector|
   - Grid resolution as :math:`(r, \theta, \phi)`. Grids with variations on only
     the radial resolution have optimized traversal. Default: [1,1,1]

This plugin creates a spherical extremum structure storing local extremum values
for efficient delta tracking in spherical media. The grid is constructed by
querying the underlying volume's extrema over each spherical cell.

At runtime, concentric shell traversal provides tight-fitting local extremum for
radially-varying media such as planetary atmospheres.

.. warning:: Azimuth is ill-defined for a ray whose line meets the polar
   axis at an angle (``cross(o, d).z == 0`` without the ray being on the
   axis).
*/

template <typename Float, typename Spectrum>
class ExtremumSpherical final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox, m_scale)
    MI_IMPORT_TYPES(Volume)

    ExtremumSpherical(const Properties &props) : Base(props), m_props(props) {
        ScalarVector3i resolution = props.get<ScalarVector3i>("resolution", ScalarVector3i(1, 1, 1));

        if (resolution.x() < 1 || resolution.y() < 1 || resolution.z() < 1)
            Throw("All resolution components must be >= 1!");

        // Determine traversal type from resolution
        if (resolution.y() == 1 && resolution.z() == 1) {
            m_traversal_type = SphericalTraversalType::RadialOnly;
        } else {
            m_traversal_type = SphericalTraversalType::Full3D;
        }

        // Mark all properties as queried so they don't warn in expand()
        props.mark_queried("resolution");
    }

    template <SphericalTraversalType TT>
    using Impl = ExtremumSphericalImpl<Float, Spectrum, TT>;

    std::vector<ref<Object>> expand() const override {
        ref<Object> result;
        switch (m_traversal_type) {
            case SphericalTraversalType::RadialOnly:
                result = (Object *) new Impl<SphericalTraversalType::RadialOnly>(m_props);
                break;
            case SphericalTraversalType::Full3D:
                result = (Object *) new Impl<SphericalTraversalType::Full3D>(m_props);
                break;
            default:
                Throw("Unsupported spherical traversal type!");
        }
        return { result };
    }

    void build(const Volume *) override {
        NotImplementedError("build");
    }

    // Stub overrides — never called, expand() replaces this object
    std::tuple<Float, Float> eval_1(const Interaction3f &,
                                    Mask) const override {
        NotImplementedError("eval_1");
    }

    DDAState dda_init(const Ray3f &, Float, Float, Mask) const override {
        NotImplementedError("dda_init");
    }

    std::pair<ExtremumSegment, DDAState> dda_next(const DDAState &,
                                                  Mask) const override {
        NotImplementedError("dda_next");
    }

    MI_DECLARE_CLASS(ExtremumSpherical)

protected:
    Properties m_props;
    SphericalTraversalType m_traversal_type;
};


// ---------------------------------------------------------------------------
// Implementation class
// ---------------------------------------------------------------------------

template <typename Float, typename Spectrum, SphericalTraversalType TraversalType>
class ExtremumSphericalImpl final : public ExtremumStructure<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ExtremumStructure, m_bbox, m_scale)
    MI_IMPORT_TYPES(Volume)

    using FloatStorage         = DynamicBuffer<Float>;

    static constexpr size_t Dim =
        TraversalType == SphericalTraversalType::RadialOnly ? 1 : 3;
    using CoordF = dr::Array<Float, Dim>;
    using CoordI = dr::Array<Int32, Dim>;


    ExtremumSphericalImpl(const Properties &props) : Base(props) {
        m_resolution =
            props.get<ScalarVector3i>("resolution", ScalarVector3i(1, 1, 1));
    }

    void build(const Volume *volume) override {

        VolumeParametrization<ScalarFloat> volume_param = volume->parametrization();

        if (volume_param.flag != VolumeCoordFlag::Spherical)
            Throw("ExtremumSpherical only compatible with volumes that have a spherical parametrization.");

        m_center   = volume_param.to_world.translation();
        m_to_local = volume_param.to_world.inverse();
        m_rmin     = volume_param.uv_range.min[0];
        m_rmax     = volume_param.uv_range.max[0];

        if (m_rmin >= m_rmax)
            Throw("rmin must be less than rmax!");

        m_dr = (m_rmax - m_rmin) / m_resolution.x();
        m_idr = dr::rcp(m_dr);
        m_dtheta = dr::Pi<ScalarFloat> / m_resolution.y();
        m_dphi   = dr::TwoPi<ScalarFloat> / m_resolution.z();
        m_eps =
            ScalarVector3f(m_resolution) * math::RayEpsilon<ScalarFloat>;

        build_grid(volume);
        build_angle_tables();
    }

    DDAState dda_init(const Ray3f &ray, Float mint, Float maxt,
                      Mask active) const override {
        auto [hit, d0, d1] = m_bbox.ray_intersect(ray);

        DDAState state;
        state.o = m_to_local * ray.o;
        state.d = m_to_local * ray.d;

        state.mint = dr::select(hit, dr::maximum(mint, d0), mint);
        state.maxt =
            dr::select(hit && active, dr::minimum(maxt, d1), state.mint);

        const Point3f  &o = state.o;
        const Vector3f &d = state.d;
        RayCoeffs rc = ray_coeffs(o, d);

        Point3f p = dr::fmadd(d, state.mint, o);
        Float   r = dr::norm(p);

        state.pi = dr::zeros<Vector3i>();
        state.pi.x() = radial_idx(
            r, idx_bias(m_eps.x(), dr::fmadd(state.mint, rc.a, rc.b)));

        if constexpr (TraversalType == SphericalTraversalType::Full3D) {
            Mask on_axis = on_axis_ray(o, d);
            state.pi.y() = theta_idx(
                p, r, on_axis,
                idx_bias(m_eps.y(), -dr::fmadd(state.mint, rc.n1, rc.n0)));
            state.pi.z() = phi_idx(p, idx_bias(m_eps.z(), rc.phi_dot));
        }

        return state;
    }

    std::pair<ExtremumSegment, DDAState>
    dda_next(const DDAState &state, Mask active) const override {
        const Point3f  &o = state.o;
        const Vector3f &d = state.d;
        RayCoeffs rc = ray_coeffs(o, d);

        Vector3f t_turn     = -dr::Infinity<Float>;
        Vector3i step_after = dr::zeros<Vector3i>();

        t_turn.x()     = -rc.b * rc.inv_a;
        step_after.x() = 1;

        Mask on_axis = false;

        if constexpr (TraversalType == SphericalTraversalType::Full3D) {
            on_axis = on_axis_ray(o, d);

            Mask no_turn = rc.n1 == 0.f;
            t_turn.y() =
                dr::select(no_turn, -dr::Infinity<Float>, -rc.n0 / rc.n1);
            dr::masked(t_turn.y(), on_axis) = -o.z() * dr::rcp(d.z());

            step_after.y() =
                dr::select(dr::select(no_turn, rc.n0, rc.n1) < 0.f, 1, -1);
            step_after.z() = dr::select(rc.phi_dot > 0.f, 1, -1);
        }

        Vector3f dt = dr::Infinity<Float>;
        const Float threshold = state.mint + dr::Epsilon<Float> * 2.f;

        Vector3i step =
            dr::select(state.mint < t_turn, -step_after, step_after);
        Vector3i test_idx =
            state.pi + dr::select(step > 0, Vector3i(1), Vector3i(0));

        dt.x() = sphere_crossing(
            rc, shell_radius(dr::clip(test_idx.x(), 0, m_resolution.x())),
            threshold);

        if constexpr (TraversalType == SphericalTraversalType::Full3D) {
            dr::masked(dt.y(), !on_axis) =
                cone_crossing(rc, o, d, cos_theta(test_idx.y()), threshold);
            dt.z() = plane_crossing(o, d, phi_normal(test_idx.z()), threshold,
                                    !on_axis);
        }

        Vector3f t_turn_fwd =
            dr::select(t_turn > threshold, t_turn, dr::Infinity<Float>);

        Float t_next =
            dr::minimum(dr::min(dr::minimum(dt, t_turn_fwd)), state.maxt);
        dr::masked(t_next, !active) = state.maxt;

        UInt32 idx = UInt32(dr::clip(state.pi.x(), 0, m_resolution.x() - 1));
        if constexpr (TraversalType == SphericalTraversalType::Full3D)
            idx += UInt32(state.pi.y() * m_resolution.x() +
                          state.pi.z() * (m_resolution.x() * m_resolution.y()));

        Vector2f extremum = dr::gather<Vector2f>(m_extremum_grid, idx, active);
        dr::masked(extremum, state.pi.x() < 0) = Vector2f(m_fillmin);
        dr::masked(extremum, state.pi.x() >= m_resolution.x()) =
            Vector2f(m_fillmax);

        DDAState next = state;
        next.mint = t_next;

        const Float next_threshold = t_next + dr::Epsilon<Float> * 2.f;

        auto crossed =
            (dt <= next_threshold) && !(t_turn_fwd <= next_threshold);

        if constexpr (TraversalType == SphericalTraversalType::RadialOnly) {
            dr::masked(next.pi.x(), crossed.x() && active) += step.x();
        } else {
            dr::masked(next.pi, crossed && active) += step;

            next.pi.y() = dr::clip(next.pi.y(), 0, m_resolution.y() - 1);

            dr::masked(next.pi.y(), on_axis) = dr::select(
                dr::fmadd(t_next, d.z(), o.z()) > 0.f, 0, m_resolution.y() - 1);

            dr::masked(next.pi.z(), next.pi.z() < 0) = m_resolution.z() - 1;
            dr::masked(next.pi.z(), next.pi.z() >= m_resolution.z()) = 0;
        }

        return { ExtremumSegment(state.mint, t_next, m_scale * extremum),
                 next };
    }

    std::tuple<Float, Float> eval_1(
        const Interaction3f &it,
        Mask active
    ) const override {
        Point3f po = m_to_local * it.p;
        Float r = dr::norm(po);

        Float fillval = -1.f;
        dr::masked(fillval, r < m_rmin) = m_fillmin;
        dr::masked(fillval, r > m_rmax) = m_fillmax;
        Mask fill = fillval >= 0.f;
        Vector2f extremum = dr::zeros<Vector2f>();

        UInt32 ir = UInt32(dr::clip(radial_idx(r), 0, m_resolution.x() - 1));

        if constexpr (TraversalType == SphericalTraversalType::RadialOnly) {
            extremum = dr::gather<Vector2f>(m_extremum_grid, ir, active && !fill);
        } else if constexpr (TraversalType == SphericalTraversalType::Full3D) {
            UInt32 itheta = UInt32(theta_idx(po, r, Mask(false)));
            UInt32 iphi   = UInt32(phi_idx(po));
            UInt32 idx = ir + itheta * UInt32(m_resolution.x()) +
                        iphi * UInt32(m_resolution.x() * m_resolution.y());
            extremum = dr::gather<Vector2f>(m_extremum_grid, idx, active && !fill);
        }

        extremum = dr::select(
                fill, Vector2f(fillval, fillval), extremum
        );

        return { m_scale * extremum.x(), m_scale * extremum.y() };
    }

    void traverse(TraversalCallback *cb) override {
        cb->put("data", m_extremum_grid, ParamFlags::NonDifferentiable);
        cb->put("resolution", m_resolution, ParamFlags::NonDifferentiable);
        cb->put("scale", m_scale, ParamFlags::NonDifferentiable);
        Base::traverse(cb);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ExtremumSpherical[" << std::endl
            << "  traversal = "
            << (TraversalType == SphericalTraversalType::RadialOnly
                ? "RadialOnly" : "Full3D")
            << "," << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  to_world = "   << m_to_local.inverse() << "," << std::endl
            << "  rmin = "       << m_rmin << "," << std::endl
            << "  rmax = "       << m_rmax << std::endl
            << "  fillmin = "    << m_fillmin << "," << std::endl
            << "  fillmax = "    << m_fillmax << std::endl
            << "  scale = " << m_scale << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS(ExtremumSphericalImpl)

private:

    // ------------------------------------------------------------------
    // Grid construction
    // ------------------------------------------------------------------

    void build_grid(const Volume *volume) {
        // Cell size in normalized [0,1]^3 space
        const ScalarVector3f cell_size = dr::rcp(ScalarVector3f(m_resolution));
        size_t n = dr::prod(m_resolution);

        ScalarVector2f safety_factor(
            1.f - dr::Epsilon<Float>,
            1.f + dr::Epsilon<Float>
        );

        size_t n_threads = pool_size() + 1;
        size_t grain_size = std::max(n / (4 * n_threads), (size_t) 1);

        m_extremum_grid = dr::empty<FloatStorage>(n * 2);

        if constexpr (!dr::is_jit_v<Float>) {
            auto guard = volume->pin();

            dr::parallel_for(
                dr::blocked_range<size_t>(0, n, grain_size),
                [&](const dr::blocked_range<size_t> &range) {
                    for (auto idx = range.begin(); idx != range.end(); ++idx) {
                        // r-fastest indexing: idx = ir + itheta * res_r + iphi * res_r * res_theta
                        int32_t ir     = idx % m_resolution.x();
                        int32_t itheta = (idx / m_resolution.x()) % m_resolution.y();
                        int32_t iphi   = idx / (m_resolution.x() * m_resolution.y());

                        ScalarPoint3f cell_min =
                            ScalarVector3f(ir, itheta, iphi) * cell_size;
                        ScalarPoint3f cell_max = cell_min + cell_size;
                        ScalarBoundingBox3f cell_bounds(
                            cell_min + math::RayEpsilon<Float>,
                            cell_max - math::RayEpsilon<Float>);

                        auto [min, maj] = volume->extremum(cell_bounds);

                        dr::scatter(m_extremum_grid,
                                    Vector2f(min, maj) * safety_factor,
                                    UInt32(idx));
                    }
                }
            );
        } else {
            UInt32 idx = dr::arange<UInt32>((uint32_t) n);

            UInt32 ir     = idx % m_resolution.x();
            UInt32 itheta = (idx / m_resolution.x()) % m_resolution.y();
            UInt32 iphi   = idx / (m_resolution.x() * m_resolution.y());

            Point3f cell_min = Vector3f(ir, itheta, iphi) * cell_size;
            Point3f cell_max = cell_min + cell_size;
            BoundingBox3f cell_bounds(
                cell_min + math::RayEpsilon<Float>,
                cell_max - math::RayEpsilon<Float>
            );

            auto [min, maj] = volume->extremum(cell_bounds);

            dr::scatter(m_extremum_grid, min * safety_factor.x(), idx * 2);
            dr::scatter(m_extremum_grid, maj * safety_factor.y(), idx * 2 + 1);
            dr::sync_thread();
        }

        // Retrieve fillmin and fillmax
        Interaction3f it = dr::zeros<Interaction3f>();

        it.p          = m_center;
        Float fillmin = volume->eval_1(it, true);

        it.p          = m_to_local.inverse() * ScalarPoint3f(0.f, 0.f, m_rmax + 1.f);
        Float fillmax = volume->eval_1(it, true);

        if constexpr (dr::is_jit_v<Float>) {
            m_fillmin = fillmin[0];
            m_fillmax = fillmax[0];
        } else {
            m_fillmin = fillmin;
            m_fillmax = fillmax;
        }

        Log(Info, "Extremum spherical grid constructed successfully");
    }

    /// Precomputes every cos(theta) for every theta cone boudary and 2D
    /// normal (sin(phi), -cos(phi)) for every azimuth boundary.
    void build_angle_tables() {
        std::vector<ScalarFloat> cosines(m_resolution.y() + 1);
        for (int32_t i = 0; i <= m_resolution.y(); ++i) {
            cosines[i] = dr::cos(ScalarFloat(i) * m_dtheta);
            // For even resolution, the middle boundary degenartes into a plane
            // Store an exact zero to ensure numerical stability.
            dr::masked(cosines[i], 2 * i == m_resolution.y()) = ScalarFloat(0);
        }

        std::vector<ScalarFloat> normals(2 * (m_resolution.z() + 1));
        for (int32_t i = 0; i <= m_resolution.z(); ++i) {
            // Azimuth boundary at phi = -pi + i * dphi
            auto [s, c] = dr::sincos(dr::fmadd(ScalarFloat(i), m_dphi, -dr::Pi<ScalarFloat>));
            normals[2 * i]     = s;
            normals[2 * i + 1] = -c;
        }

        m_cos_theta  = dr::load<FloatStorage>(cosines.data(), cosines.size());
        m_phi_normal = dr::load<FloatStorage>(normals.data(), normals.size());
    }

    // ------------------------------------------------------------------
    // Boundary-crossing helpers
    // ------------------------------------------------------------------

    /// Ray-only terms shared by the crossing tests, independent of boundaries.
    struct RayCoeffs {
        Float a, b, inv_a, o_sqr;
        /// dot(o,d)^2 - |d|^2 |o|^2, i.e. -|o x d|^2
        Float disc_base;
        /// Cone terms: d_z^2, o_z^2, d_z o_z, and k = |d_z o - o_z d|^2
        Float dz2, oz2, dzoz, k;
        /// d/dt cos(theta) has the sign of n0 + n1 t; d/dt phi that of phi_dot
        Float n0, n1, phi_dot;
    };

    static RayCoeffs ray_coeffs(const Point3f &o, const Vector3f &d) {
        RayCoeffs rc {};
        rc.a         = dr::squared_norm(d);
        rc.b         = dr::dot(o, d);
        rc.o_sqr     = dr::squared_norm(o);
        rc.inv_a     = dr::rcp(rc.a);

        Vector3f c   = dr::cross(o, d);
        rc.disc_base = -dr::squared_norm(c);

        if constexpr (TraversalType == SphericalTraversalType::Full3D) {
            rc.dz2       = dr::square(d.z());
            rc.oz2       = dr::square(o.z());
            rc.dzoz      = d.z() * o.z();
            rc.k         = dr::squared_norm(Vector2f(c.x(), c.y()));
            rc.n0        = dr::fmsub(d.z(), rc.o_sqr, o.z() * rc.b);
            rc.n1        = dr::fmsub(rc.b, d.z(), rc.a * o.z());
            rc.phi_dot   = c.z();
        }
        return rc;
    }

    /// Signed tolerance sending a query that sits on a cell boundary into the
    /// cell it is heading into. `rate` is the derivative of the cell coordinate.
    static Float idx_bias(ScalarFloat tol, Float rate) {
        return dr::select(rate < 0.f, Float(-tol), Float(tol));
    }

    static Float forward(Float t, Mask valid, Float threshold) {
        return dr::select(valid && (t > threshold), t, dr::Infinity<Float>);
    }

    /// Nearest crossing, strictly ahead of `threshold`, of ray (o, d) with
    /// the origin-centered sphere of radius `r_test`; +inf if there is none.
    static Float sphere_crossing(const RayCoeffs &rc, Float r_test,
                                 Float threshold, Mask valid = true) {
        Float disc = dr::fmadd(rc.a, dr::square(r_test), rc.disc_base);
        Float sq   = dr::sqrt(dr::maximum(disc, 0.f));
        Float near = (-rc.b - sq) * rc.inv_a,
              far  = (-rc.b + sq) * rc.inv_a;
        // TODO: now disc should never be invalid because we check if we are in the midpoint ahead of time.
        return forward(dr::select(near > threshold, near, far),
                       valid && (disc >= 0.f), threshold);
    }

    /// Ray (o,d)/Cone intersection. Cone defined by cos(theta) (c) w.r.t the
    /// vertical axis in local frame. +inf for no intersection.
    static Float cone_crossing(const RayCoeffs &rc, const Point3f &o,
                               const Vector3f &d, Float c, Float threshold,
                               Mask valid = true) {
        Float c2 = dr::square(c);

        // Factored discriminant form so that disc = 0 when c = 0.
        Float disc = c2 * dr::fmadd(c2, rc.disc_base, rc.k);
        Float sq   = dr::sqrt(dr::maximum(disc, 0.f));

        Float qa = dr::fnmadd(c2, rc.a,     rc.dz2),
              qb = dr::fnmadd(c2, rc.b,     rc.dzoz),
              qc = dr::fnmadd(c2, rc.o_sqr, rc.oz2);

        // Citardauq form of the quadratic formulation, generic expression
        // suffers from catastrophic cancellation at discriminant = 0.
        // Stable root form (t1 = c/a*t0), avoid double counting equatorial plane.
        Float tmp = -(qb + dr::copysign(sq, qb));
        Float t0  = tmp / qa;
        Float t1  = dr::select(disc > 0.f, qc / tmp, t0);

        valid &= disc >= 0.f;

        // discard intersection with opposite cones.
        auto valid_side = [&](Float t) DRJIT_INLINE_LAMBDA {
            return valid && (dr::fmadd(t, d.z(), o.z()) * c >= 0.f);
        };
        return dr::minimum(forward(t0, valid_side(t0), threshold),
                           forward(t1, valid_side(t1), threshold));
    }

    /// Ray(o,d)/Half-plane intersection test.
    /// The plane goes through the vertical axis and has normal n.
    static Float plane_crossing(const Point3f &o, const Vector3f &d,
                                const Vector2f &n, Float threshold,
                                Mask valid = true) {
        return forward(-(o.x() * n.x() + o.y() * n.y()) /
                        (d.x() * n.x() + d.y() * n.y()), valid, threshold);
    }

    /// Radius of shell boundary `idx`.
    Float shell_radius(Int32 idx) const {
        return dr::fmadd(Float(idx), m_dr, m_rmin);
    }

    /// Tabulated cos(theta) of colatitude boundary `idx`.
    Float cos_theta(Int32 idx) const {
        return dr::gather<Float>(m_cos_theta, UInt32(idx));
    }

    /// Tabulated 2D normal of the plane at azimuth boundary `idx`.
    Vector2f phi_normal(Int32 idx) const {
        return dr::gather<Vector2f>(m_phi_normal, UInt32(idx));
    }

    /// Whether the local ray runs along the vertical axis
    Mask on_axis_ray(const Point3f &o, const Vector3f &d) const {
        const ScalarFloat tol = math::RayEpsilon<ScalarFloat>;
        return dr::squared_norm(Vector2f(o.x(), o.y())) <= tol * dr::squared_norm(o)
            && dr::squared_norm(Vector2f(d.x(), d.y())) <= tol * dr::squared_norm(d);
    }

    /// Shell index for radius `r`. Using -1 and resolution.x() as sentinels
    /// for the fillmin and fillmax exterior regions. `bias` shifts the cell
    /// coordinate, see \ref m_eps.
    Int32 radial_idx(Float r, Float bias = 0.f) const {
        return dr::clip(
            dr::floor2int<Int32>(dr::fmadd(r - m_rmin, m_idr, bias)),
            -1, m_resolution.x());
    }

    /// Zenithal cell index w.r.t vertical axis. for position `p` at radius
    /// `r`. `on_axis` re-derives the pole cell from the sign of z instead.
    Int32 theta_idx(const Point3f &p, Float r, Mask on_axis, Float bias = 0.f) const {
        Int32 result( 0 );

        Float theta =
            dr::acos(p.z() * dr::rcp(dr::maximum(r, dr::Smallest<Float>) ) )
            * dr::InvPi<Float>;
        dr::masked(result, !on_axis) =
            dr::clip(
                dr::floor2int<Int32>( dr::fmadd(theta, Float(m_resolution.y()), bias) ),
                0, m_resolution.y() - 1);

        if (unlikely(dr::any_or<true>(on_axis)))
            dr::masked(result, on_axis)  =
                dr::select(p.z() > 0.f, 0, m_resolution.y() - 1);

        return result;
    }

    /// Azimuth cell index for position `p`, wrapped into [0, resolution.z()).
    Int32 phi_idx(const Point3f &p, Float bias = 0.f) const {
        Int32 i = dr::floor2int<Int32>(
            dr::fmadd(dr::atan2(p.y(), p.x()) * dr::InvTwoPi<Float> + 0.5f,
                      Float(m_resolution.z()), bias));
        dr::masked(i, i < 0) += m_resolution.z();
        dr::masked(i, i >= m_resolution.z()) -= m_resolution.z();
        return i;
    }

private:
    FloatStorage m_extremum_grid;
    FloatStorage m_cos_theta, m_phi_normal;
    ScalarVector3i m_resolution;
    ScalarFloat m_rmin, m_rmax;
    ScalarFloat m_fillmin, m_fillmax;
    ScalarPoint3f m_center;
    ScalarFloat m_dr, m_idr;
    ScalarFloat m_dtheta, m_dphi;

    /// Per-dimension tolerance on a cell dimension
    ScalarVector3f m_eps;

    ScalarAffineTransform4f m_to_local;
};

// ---------------------------------------------------------------------------
// Class name helpers (for expand pattern)
// ---------------------------------------------------------------------------

MI_EXPORT_PLUGIN(ExtremumSpherical)

NAMESPACE_BEGIN(detail)
template <SphericalTraversalType TT>
constexpr const char *extremum_spherical_class_name() {
    if constexpr (TT == SphericalTraversalType::RadialOnly) {
        return "ExtremumSpherical_RadialOnly";
    } else if constexpr (TT == SphericalTraversalType::Full3D) {
        return "ExtremumSpherical_Full3D";
    }
}
NAMESPACE_END(detail)

NAMESPACE_END(mitsuba)
