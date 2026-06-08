#pragma once

#include <mitsuba/core/object.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/sampler.h>
#include <drjit/call.h>

NAMESPACE_BEGIN(mitsuba)

/**
 * \brief This enumeration is used to classify phase functions into different types,
 * i.e. into isotropic, anisotropic and microflake phase functions.
 *
 * This can be used to optimize implementations to for example have less overhead
 * if the phase function is not a microflake phase function.
 */
enum class PhaseFunctionFlags : uint32_t {
    Empty       = 0x00,
    Isotropic   = 0x01,
    Anisotropic = 0x02,
    Microflake  = 0x04
};

MI_DECLARE_ENUM_OPERATORS(PhaseFunctionFlags)

/**
 * \brief Context data structure for phase function evaluation and sampling
 *
 * Phase function models in Mitsuba can be queried and sampled using a variety of
 * different modes. Using this data structure, a rendering algorithm can indicate whether
 * radiance or importance is being transported.
 *
 * The context further holds a pointer to a sampler object, in case
 * the evaluation or sampling functions need additional random numbers.
 *
 */
MI_VARIANT
struct MI_EXPORT_LIB PhaseFunctionContext {
    MI_IMPORT_TYPES(Sampler);

    // =============================================================
    //! @{ \name Fields
    // =============================================================

    /// Transported mode (radiance or importance)
    TransportMode mode = TransportMode::Radiance;

    /// Sampler object
    Sampler *sampler = nullptr;

    /*
     * Bit mask for requested phase function component types to be
     * sampled/evaluated.
     * The default value (equal to \ref PhaseFunctionFlags::All) enables all
     * components.
     */
    uint32_t type_mask = (uint32_t) 0x7u;

    /*
     * Integer value of requested phase function component index to be
     * sampled/evaluated.
     */
    uint32_t component = (uint32_t) -1;

    //! @}
    // =============================================================

    PhaseFunctionContext() = default;

    PhaseFunctionContext(Sampler *sampler,
                         TransportMode mode = TransportMode::Radiance)
        : mode(mode), sampler(sampler) { }

    PhaseFunctionContext(Sampler *sampler, TransportMode mode,
                         uint32_t type_mask, uint32_t component)
        : mode(mode), sampler(sampler), type_mask(type_mask),
          component(component) { }

    /**
     * \brief Reverse the direction of light transport in the record
     *
     * This updates the transport mode (radiance to importance and vice versa).
     */
    void reverse() { mode = (TransportMode)(1 - (int) mode); }

    /**
     * Checks whether a given phase function component type and index are
     * enabled in this context.
     */
    bool is_enabled(PhaseFunctionFlags type_, uint32_t component_ = 0) const {
        uint32_t type = (uint32_t) type_;
        return (type_mask == (uint32_t) -1 || (type_mask & type) == type) &&
               (component == (uint32_t) -1 || component == component_);
    }
};

/**
 * \brief Abstract phase function base-class.
 *
 * This class provides an abstract interface to all Phase function plugins in
 * Mitsuba. It exposes functions for evaluating and sampling the model.
 */

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB PhaseFunction : public JitObject<PhaseFunction<Float, Spectrum>> {
public:
    MI_IMPORT_TYPES(PhaseFunctionContext);

// #ERADIATE_CHANGE_BEGIN: DDIS
    using FloatStorage = mitsuba::DynamicBuffer<Float>;
// #ERADIATE_CHANGE_END
    /**
     * \brief Importance sample the phase function model
     *
     * The function returns a sampled direction.
     *
     * \param ctx
     *     A phase function sampling context, contains information
     *     about the transport mode
     *
     * \param mi
     *     A medium interaction data structure describing the underlying
     *     medium position. The incident direction is obtained from
     *     the field <tt>mi.wi</tt>.
     *
     * \param sample1
     *     A uniformly distributed sample on \f$[0,1]\f$. It is used
     *     to select the phase function component in multi-component models.
     *
     * \param sample2
     *     A uniformly distributed sample on \f$[0,1]^2\f$. It is
     *     used to generate the sampled direction.
     *
     * \return A sampled direction wo and its corresponding weight and PDF
     */
    virtual std::tuple<Vector3f, Spectrum, Float> sample(const PhaseFunctionContext &ctx,
                                                         const MediumInteraction3f &mi,
                                                         Float sample1, const Point2f &sample2,
                                                         Mask active = true) const = 0;
    /**
     * \brief Evaluates the phase function model value and PDF
     *
     * The function returns the value (which often equals the PDF) of the phase
     * function in the query direction.
     *
     * \param ctx
     *     A phase function sampling context, contains information
     *     about the transport mode
     *
     * \param mi
     *     A medium interaction data structure describing the underlying
     *     medium position. The incident direction is obtained from
     *     the field <tt>mi.wi</tt>.
     *
     * \param wo
     *     An outgoing direction to evaluate.
     *
     * \return The value and the sampling PDF of the phase function in direction wo
     */
    virtual std::pair<Spectrum, Float> eval_pdf(const PhaseFunctionContext &ctx,
                                                const MediumInteraction3f &mi,
                                                const Vector3f &wo,
                                                Mask active = true) const = 0;

    /**
     * \brief Returns the microflake projected area
     *
     * The function returns the projected area of the microflake distribution defining the phase
     * function. For non-microflake phase functions, e.g. isotropic or Henyey-Greenstein, this
     * should return a value of 1.
     *
     * \param mi
     *     A medium interaction data structure describing the underlying
     *     medium position. The incident direction is obtained from
     *     the field <tt>mi.wi</tt>.
     *
     * \return The projected area in direction <tt>mi.wi</tt> at position <tt>mi.p</tt>
     */
    virtual Float projected_area(const MediumInteraction3f & /* mi */, Mask /* active */ = true) const {
        return 1.f;
    }

    /// Return the maximum projected area of the microflake distribution
    virtual Float max_projected_area() const { return 1.f; }

    /// Flags for this phase function.
    uint32_t flags(Mask /*active*/ = true) const { return m_flags; }

    /// Flags for a specific component of this phase function.
    uint32_t flags(size_t i, Mask /*active*/ = true) const {
        Assert(i < m_components.size());
        return m_components[i];
    }

    /// Number of components this phase function is comprised of.
    size_t component_count(Mask /*active*/ = true) const {
        return m_components.size();
    }

    /// Return a human-readable representation of the phase function
    std::string to_string() const override = 0;

    /// Return type of phase function
    uint32_t get_flags() const { return m_flags; }

    /// Set type of phase function
    void set_flags(uint32_t flags) { m_flags = flags; }
// #ERADIATE_CHANGE_BEGIN: DDIS
    /// Return whether this phase function's parameters have changed since the
    /// last call to set_dirty(false). Set by parameters_changed(), cleared by
    /// the scene after all dependent media have been updated.
    bool dirty() const { return m_dirty; }

    /// Modify the phase function's dirty flag
    void set_dirty(bool dirty) { m_dirty = dirty; }

    void parameters_changed(const std::vector<std::string> &keys = {}) override;

    /**
     * \brief Populate a set of cos_theta nodes suitable for representing this
     * phase function.
     *
     * The nodes are used to build the piecewise-linear envelope required by
     * the DDIS importance sampling scheme. Subclasses may override this method
     * to supply irregularly spaced nodes that better resolve sharp features
     * (e.g. a strong forward-scattering peak). The default implementation
     * places \c m_node_count nodes uniformly in [-1, 1].
     *
     * \note cos_theta follows the physics convention: a value of +1 corresponds
     *     to aligned (forward-scattering) incoming and outgoing directions, and
     *     -1 corresponds to exact backscatter.
     *
     * \return A sorted buffer of cos_theta values at which the phase function
     *     should be evaluated.
     */
    virtual FloatStorage get_nodes() const;


     /**
     * \brief Evaluate the phase function at the given cos_theta nodes and
     * accumulate the result into \c values by taking the elementwise maximum.
     *
     * For each node \f$\mu_i\f$ in \c nodes this method evaluates the phase
     * function value \f$p(\mu_i)\f$ and updates
     * \f$\texttt{values}[i] \leftarrow \max(\texttt{values}[i],\, p(\mu_i))\f$.
     *
     * Delegating the comparison to the callee rather than the caller enables
     * natural recursion through composite phase functions (e.g.
     * \ref BlendPhaseFunction, \ref MultiPhaseFunction): a composite
     * implementation simply calls \c eval_max on each child with the same
     * buffer, and each child accumulates its contribution independently. The
     * resulting buffer holds the pointwise supremum over the entire
     * phase-function tree without the caller needing to know its structure.
     *
     * \note cos_theta follows the physics convention (see \ref get_nodes).
     *
     * \param nodes
     *     cos_theta values at which to evaluate the phase function, as
     *     returned by \ref get_nodes.
     * \param values
     *     In/out buffer. Must have the same length as \c nodes and be
     *     zero-initialised before the first comparison. On return, each entry
     *     holds the maximum of its previous value and the phase function
     *     evaluated at the corresponding node.
     */
    virtual void eval_max(const FloatStorage &nodes, FloatStorage &values) const;
// #ERADIATE_CHANGE_END
    //! @}
    // -----------------------------------------------------------------------

    MI_DECLARE_PLUGIN_BASE_CLASS(PhaseFunction)

protected:
    PhaseFunction(const Properties &props);

protected:
    /// Type of phase function (e.g. anisotropic)
    uint32_t m_flags;

    /// Flags for each component of this phase function.
    std::vector<uint32_t> m_components;
    // #ERADIATE_CHANGE_BEGIN: DDIS
    /// Number of nodes to regularly discretize the cos theta dimension in \ref eval_max
    size_t m_node_count;
    /// True if the phase function's parameters have changed since the last scene update
    bool m_dirty = false;
    // #ERADIATE_CHANGE_END
};

MI_VARIANT
std::ostream &operator<<(std::ostream &os, const PhaseFunctionContext<Float, Spectrum>& ctx) {
    os << "PhaseFunctionContext[" << std::endl
       << "  mode = " << ctx.mode << "," << std::endl
       << "  sampler = " << ctx.sampler << "," << std::endl
       << "  component = ";
    if (ctx.component == (uint32_t) -1)
        os << "all";
    else
        os << ctx.component;
    os << std::endl << "]";
    return os;
}

//! @}
// -----------------------------------------------------------------------

MI_EXTERN_CLASS(PhaseFunction)
NAMESPACE_END(mitsuba)

// -----------------------------------------------------------------------
//! @{ \name Enables vectorized calls on Dr.Jit arrays of phase functions
// -----------------------------------------------------------------------

DRJIT_CALL_TEMPLATE_BEGIN(mitsuba::PhaseFunction)
    DRJIT_CALL_METHOD(sample)
    DRJIT_CALL_METHOD(eval_pdf)
    DRJIT_CALL_METHOD(projected_area)
    DRJIT_CALL_METHOD(max_projected_area)
    DRJIT_CALL_GETTER(flags)
    DRJIT_CALL_GETTER(component_count)
// #ERADIATE_CHANGE_BEGIN: DDIS
    DRJIT_CALL_METHOD(get_nodes)
    DRJIT_CALL_METHOD(eval_max)
// #ERADIATE_CHANGE_END
DRJIT_CALL_END()

//! @}
// -----------------------------------------------------------------------
