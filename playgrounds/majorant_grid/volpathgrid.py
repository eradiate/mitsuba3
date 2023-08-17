import sys
import typing as t

import drjit as dr
import mitsuba as mi

MI_PATH = "/home/leroyv/Documents/src/thirdparty/mitsuba3/build/python"

if MI_PATH not in sys.path:
    sys.path.insert(0, MI_PATH)


class VolPathGridIntegratorPy(mi.SamplingIntegrator):
    """
    Known limitations:

    * Refractive index set to 1
    * Single medium with null BSDF
    """

    def __init__(self, props):
        super().__init__(props)
        self.max_depth = props.get("max_depth", 999)
        self.rr_depth = props.get("rr_depth", 5)

    def sample(
        self,
        scene: mi.Scene,
        sampler: mi.Sampler,
        ray_: mi.Ray3f,
        initial_medium: mi.Medium,
        active: mi.Bool = None,
    ) -> t.Tuple[mi.Spectrum, mi.Bool, mi.Spectrum]:
        """
        This is basically Merlin's VolPathSimple.sample() implementation
        with the AD parts removed (at least this is the intention).

        TODO: Add surface interaction handling
        TODO: Fix basic issues
        """

        ray = mi.Ray3f(ray_)
        wavefront_size = dr.width(ray.d)
        result = mi.Spectrum(0.0)
        throughput = mi.Spectrum(1.0)
        medium = initial_medium
        mei = dr.zeros(mi.MediumInteraction3f)
        depth = dr.zeros(mi.Int32, wavefront_size)
        active = mi.Mask(True)

        si, escaped = self.reach_medium(
            scene, ray, active
        )  # TODO: this assumes ray starts outside of medium
        needs_intersection = mi.Mask(True)
        last_scatter_it = dr.zeros(mi.Interaction3f)
        last_scatter_direction_pdf = mi.Float(1.0)
        has_scattered = mi.Mask(False)

        loop = mi.Loop(
            "VolPathGridIntegratorPy::Sample loop",
            lambda: (
                active,
                depth,
                ray,
                throughput,
                si,
                result,
                sampler,
                last_scatter_it,
                last_scatter_direction_pdf,
            ),
        )

        while loop(active):
            # Russian Roulette (taking eta = 1.0)
            # RR threshold 0.95 in volpath, 0.99 in volpathsimple
            q = dr.minimum(dr.max(throughput), 0.95)
            perform_rr = depth > self.rr_depth
            active &= dr.any(dr.neq(throughput, 0.0)) & (
                (~perform_rr | (sampler.next_1d(active) < q))
            )
            throughput[perform_rr] = throughput * dr.rcp(dr.detach(q))

            # Handle medium sampling and potential medium escape
            mei, mei_weight = self.sample_real_interaction(medium, ray, sampler, active)
            throughput[active] = throughput * mei_weight

            did_escape = active & (~mei.is_valid())
            still_in_medium = active & mei.is_valid()

            # Handle null and real scatter events
            did_scatter = mi.Mask(still_in_medium)

            # Rays that have still not escaped but reached max depth
            # are killed inside of the medium (zero contribution)
            depth[did_scatter] = depth + 1
            active &= still_in_medium & (depth < self.max_depth)

            # --- Emitter sampling
            phase_ctx = mi.PhaseFunctionContext(sampler)
            phase = mei.medium.phase_function()
            phase[~did_scatter] = dr.zeros(mi.PhaseFunctionPtr, 1)
            active_e = did_scatter & active
            nee_contrib = self.sample_emitter_for_nee(
                mei, scene, sampler, medium, phase_ctx, phase, throughput, active_e
            )
            result[active_e] += nee_contrib
            del nee_contrib
            # ----------

            # --- Phase function sampling
            # Note: assuming phase_pdf = 1 (perfect importance sampling)
            wo, phase_pdf = phase.sample(
                phase_ctx,
                mei,
                sampler.next_1d(did_scatter),
                sampler.next_2d(did_scatter),
                did_scatter,
            )
            new_ray = mei.spawn_ray(wo)
            ray[did_scatter] = new_ray
            # ----------

            # --- Handle escaped rays: cross the null boundary
            ray[did_escape] = si.spawn_ray(
                ray.d
            )  # Continue on the other side of the boundary
            escaped |= did_escape
            # ----------

        # --- Envmap contribution
        si_update_needed = escaped & si.is_valid()
        si[si_update_needed] = scene.ray_intersect(ray, si_update_needed)
        # All escaped rays can now query the envmap
        emitter = si.emitter(scene)
        active_e = (
            escaped & dr.neq(emitter, None) & ~((depth <= 0) & self.hide_emitters)
        )

        if self.use_nee:
            assert last_scatter_it is not None
            ds = mi.DirectionSample3f(scene, si, last_scatter_it)
            emitter_pdf = emitter.pdf_direction(last_scatter_it, ds, active_e)
            # MIS should be disabled (i.e. MIS weight = 1) if there wasn't even
            # a valid interaction from which the emitter could have been sampled,
            # e.g. in the case a ray escaped directly.
            emitter_pdf = dr.select(has_scattered, emitter_pdf, 0.0)
            hit_mis_weight = mi.ad.common.mis_weight(
                last_scatter_direction_pdf, emitter_pdf
            )
        else:
            emitter_pdf = None
            hit_mis_weight = 1.0

        # TODO: envmap gradients
        contrib = emitter.eval(si, active_e)
        result[active_e] += throughput * hit_mis_weight * contrib

        del si_update_needed, emitter, active_e, contrib

        return result, active, result

    def reach_medium(self, scene, ray, active):
        """
        In this simplified setting, rays either hit the medium's bbox and
        go in or escape directly to infinity.
        Warning: this function mutates its inputs.
        """
        si = scene.ray_intersect(ray, active)
        escaped = active & (~si.is_valid())
        active &= si.is_valid()
        # By convention, crossing the medium's bbox does *not*
        # count as an interaction (depth++) when it's a null BSDF.

        # Continue on the other side of the boundary and find
        # the opposite side (exit point from the medium).
        ray[active] = si.spawn_ray(ray.d)
        si_new = scene.ray_intersect(ray, active)
        # We might have hit a corner case and escaped despite
        # originally hitting the medium bbox.
        active &= si_new.is_valid()

        # If we lifted the restriction on the number of media,
        # we would need to get the correct pointers now.
        # medium[active] = si.target_medium(ray.d)

        ray.maxt[active] = dr.select(
            dr.isfinite(si_new.t), si_new.t, dr.largest(mi.Float)
        )
        si[active] = si_new

        return si, escaped

    def sample_real_interaction(self, medium, ray, sampler, _active):
        """
        `Medium::sample_interaction` returns an interaction that could be a null interaction.
        Here, we loop until a real interaction is sampled.

        The given ray's `maxt` value must correspond to the closest surface
        interaction (e.g. medium bounding box) in the direction of the ray.
        """
        # TODO: could make this faster for the homogeneous special case
        # TODO: could make this faster when there's a majorant supergrid
        #       by performing both DDA and "real interaction" sampling in
        #       the same loop.
        # We will keep updating the origin of the ray during traversal.
        running_ray = dr.detach(type(ray)(ray))
        # So we also keep track of the offset w.r.t. the original ray's origin.
        running_t = mi.Float(0.0)

        active = mi.Mask(_active)
        weight = mi.Spectrum(1.0)
        mei = dr.zeros(mi.MediumInteraction3f, dr.width(ray))
        mei.t = dr.select(active, dr.nan, dr.inf)

        loop = mi.Loop(
            "medium_sample_interaction_real",
            lambda: (active, weight, mei, running_ray, running_t, sampler),
        )
        while loop(active):
            mei_next = medium.sample_interaction(
                running_ray, sampler.next_1d(active), channel, active
            )
            mei[active] = mei_next
            mei.t[active] = mei.t + running_t

            majorant = mei_next.combined_extinction[channel]
            r = dr.select(dr.neq(majorant, 0), mei_next.sigma_t[channel] / majorant, 0)

            # Some lanes escaped the medium. Others will continue sampling
            # until they find a real interaction.
            active &= mei_next.is_valid()
            did_null_scatter = active & (sampler.next_1d(active) >= r)

            active &= did_null_scatter
            # Update ray to only sample points further than the
            # current null interaction.
            next_t = dr.detach(mei_next.t)
            running_ray.o[active] = running_ray.o + next_t * running_ray.d
            running_ray.maxt[active] = running_ray.maxt - next_t
            running_t[active] = running_t + next_t

        did_sample = _active & mei.is_valid()
        mei.p = dr.select(did_sample, ray(mei.t), dr.nan)
        mei.mint = mi.Float(
            dr.nan
        )  # Value was probably wrong, so we make sure it's unused
        with dr.resume_grad(when=not is_primal):
            mei.sigma_s, mei.sigma_n, mei.sigma_t = medium.get_scattering_coefficients(
                mei, did_sample
            )

        return mei, weight

    def sample_emitter(self, ref_interaction, scene, sampler, medium, channel, active):
        """
        Starting from the given `ref_interaction` inside of a medium, samples a direction
        toward an emitter and estimates transmittance with ratio tracking.

        This simplified implementation does not support:
        - presence of surfaces within the medium
        - propagating adjoint radiance (adjoint pass)
        """
        active = mi.Mask(active)

        dir_sample = sampler.next_2d(active)
        ds, emitter_val = scene.sample_emitter_direction(
            ref_interaction, dir_sample, False, active
        )
        sampling_worked = dr.neq(ds.pdf, 0.0)
        emitter_val &= sampling_worked
        active &= sampling_worked

        # Trace a ray toward the emitter and find the medium's bbox
        # boundary in that direction.
        ray = ref_interaction.spawn_ray(ds.d)
        si = scene.ray_intersect(ray, active)
        ray.maxt = si.t
        transmittance = self.estimate_transmittance(
            ray,
            0,
            si.t,
            medium,
            sampler,
            channel,
            active & si.is_valid(),
        )

        return emitter_val * transmittance, ds

    def estimate_transmittance(
        self, ray_full, tmin, tmax, medium, sampler, channel, active
    ):
        """Estimate the transmittance between two points along a ray.

        This simplified implementation does not support:
        - presence of surfaces within the medium
        - propagating adjoint radiance (adjoint pass)
        """

        # Support tmax < tmin, but not negative tmin or tmax
        needs_swap = tmax < tmin
        tmp = tmin
        tmin = dr.select(needs_swap, tmax, tmin)
        tmax = dr.select(needs_swap, tmp, tmax)
        del needs_swap, tmp

        active = mi.Mask(active)
        ray = type(ray_full)(ray_full)
        ray.o = ray_full(tmin)
        tmax = tmax - tmin
        ray.maxt = tmax
        del ray_full, tmin

        transmittance = mi.Spectrum(dr.select(active, 1.0, 0.0))

        # --- Estimate transmittance with Ratio Tracking
        # Simplified assuming that we start from within a medium, there's a single
        # medium in the scene and no surfaces.
        loop = mi.Loop(
            "VolpathSimpleNEELoop", lambda: (active, ray, tmax, transmittance, sampler)
        )
        while loop(active):
            # TODO: support majorant supergrid in-line to avoid restarting DDA traversal each time
            # Handle medium interactions / transmittance
            mei = medium.sample_interaction(
                ray, sampler.next_1d(active), channel, active
            )
            # Ratio tracking for transmittance estimation:
            # update throughput estimate with probability of sampling a null-scattering event.
            tr_contribution = dr.select(
                dr.neq(mei.combined_extinction, 0),
                mei.sigma_n / mei.combined_extinction,
                mei.sigma_n,
            )

            # If interaction falls out of bounds, we don't have anything
            # valid to accumulate.
            mei.t[active & (mei.t > tmax)] = dr.inf
            active &= mei.is_valid()

            # Apply effect of this interaction
            transmittance[active] *= tr_contribution
            # Adopt newly sampled position in the medium
            ray.o[active] = mei.p
            tmax[active] = tmax - mei.t
            ray.maxt[active] = tmax

            # Continue walking through medium
            active &= dr.any(dr.neq(transmittance, 0.0))

        return transmittance


mi.register_integrator("volpathgridpy", lambda props: VolPathGridIntegratorPy(props))
