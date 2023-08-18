import sys
import typing as t

import drjit as dr
import mitsuba as mi


MI_PATH = "/home/leroyv/Documents/src/thirdparty/mitsuba3/build/python"

if MI_PATH not in sys.path:
    sys.path.insert(0, MI_PATH)


def mis_weight(pdf_a, pdf_b):
    """
    Compute the Multiple Importance Sampling (MIS) weight given the densities
    of two sampling strategies according to the power heuristic.
    """
    a2 = dr.sqr(pdf_a)
    b2 = dr.sqr(pdf_b)
    w = a2 / (a2 + b2)
    return dr.detach(dr.select(dr.isfinite(w), w, 0))


def index_spectrum(spec, idx):
    m = spec[0]
    if mi.is_rgb:
        m[dr.eq(idx, 1)] = spec[1]
        m[dr.eq(idx, 2)] = spec[2]
    return m


class VolPathGridIntegratorPy(mi.SamplingIntegrator):
    """
    Known limitations:

    * Refractive index set to 1
    * Single medium with null BSDF
    """

    def __init__(self, props):
        super().__init__(props)
        self.max_depth = props.get("max_depth", -1)
        self.rr_depth = props.get("rr_depth", 5)
        self.hide_emitters = props.get("hide_emitters", False)

        self.use_nee = False
        self.nee_handle_homogeneous = False
        self.handle_null_scattering = False
        self.is_prepared = False

    def prepare_scene(self, scene):
        if self.is_prepared:
            return

        for shape in scene.shapes():
            for medium in [shape.interior_medium(), shape.exterior_medium()]:
                if medium:
                    # Enable NEE if a medium specifically asks for it
                    self.use_nee = self.use_nee or medium.use_emitter_sampling()
                    self.nee_handle_homogeneous = (
                        self.nee_handle_homogeneous or medium.is_homogeneous()
                    )
                    self.handle_null_scattering = self.handle_null_scattering or (
                        not medium.is_homogeneous()
                    )
        self.is_prepared = True
        # By default enable always NEE in case there are surfaces
        self.use_nee = True

    def sample(
        self,
        scene: mi.Scene,
        sampler: mi.Sampler,
        ray_: mi.Ray3f,
        initial_medium: mi.Medium,
        active: mi.Bool = None,
    ) -> t.Tuple[mi.Spectrum, mi.Bool, mi.Spectrum]:
        """
        TODO: Add surface interaction handling
        TODO: Fix basic issues
        """

        ray = mi.Ray3f(ray_)
        depth = mi.UInt32(0)
        result = mi.Spectrum(0.0)  # Radiance accumulator
        throughput = mi.Spectrum(1.0)  # Path throughput weight
        eta = mi.Float(1.0)  # TODO: add refractive index tracking
        active = mi.Mask(True)

        si = dr.zeros(mi.SurfaceInteraction3f)
        needs_intersection = mi.Bool(True)
        last_scatter_it = dr.zeros(mi.Interaction3f)
        last_scatter_direction_pdf = mi.Float(1.0)

        # # TODO: Support sensors inside media
        # medium = mi.MediumPtr(medium)
        medium = dr.zeros(mi.MediumPtr)

        channel = 0
        valid_ray = mi.Bool(False)

        if mi.is_rgb:  # Sample a color channel to sample free-flight distances
            n_channels = dr.size_v(mi.Spectrum)
            channel = dr.minimum(n_channels * sampler.next_1d(active), n_channels - 1)

        loop = mi.Loop(
            name="Volumetric path tracing with local majorant grid",
            state=lambda: (
                sampler,
                active,
                depth,
                ray,
                medium,
                si,
                throughput,
                result,
                needs_intersection,
                last_scatter_it,
                last_scatter_direction_pdf,
                valid_ray,
            ),
        )

        while loop(active):
            active &= dr.any(dr.neq(throughput, 0.0))
            q = dr.minimum(dr.max(throughput) * dr.sqr(eta), 0.99)
            perform_rr = depth > self.rr_depth
            active &= (sampler.next_1d(active) < q) | ~perform_rr
            throughput[perform_rr] = throughput * dr.rcp(q)

            active_medium = active & dr.neq(medium, None)  # TODO this is not necessary
            active_surface = active & ~active_medium

            # Handle medium sampling and potential medium escape
            u = sampler.next_1d(active_medium)
            mei = medium.sample_interaction(ray, u, channel, active_medium)
            mei.t = dr.detach(mei.t)

            ray.maxt[active_medium & medium.is_homogeneous() & mei.is_valid()] = mei.t
            intersect = needs_intersection & active_medium
            si_new = scene.ray_intersect(ray, intersect)
            si[intersect] = si_new

            needs_intersection &= ~active_medium
            mei.t[active_medium & (si.t < mei.t)] = dr.inf

            # Evaluate ratio of transmittance and free-flight PDF
            tr, free_flight_pdf = medium.eval_tr_and_pdf(mei, si, active_medium)
            tr_pdf = index_spectrum(free_flight_pdf, channel)
            weight = mi.Spectrum(1.0)
            weight[active_medium] *= dr.select(
                tr_pdf > 0.0, tr / dr.detach(tr_pdf), 0.0
            )

            escaped_medium = active_medium & ~mei.is_valid()
            active_medium &= mei.is_valid()

            # Handle null and real scatter events
            if self.handle_null_scattering:
                scatter_prob = index_spectrum(mei.sigma_t, channel) / index_spectrum(
                    mei.combined_extinction, channel
                )
                act_null_scatter = (
                    sampler.next_1d(active_medium) >= scatter_prob
                ) & active_medium
                act_medium_scatter = ~act_null_scatter & active_medium
                weight[act_null_scatter] *= mei.sigma_n / dr.detach(1 - scatter_prob)
            else:
                scatter_prob = mi.Float(1.0)
                act_medium_scatter = active_medium

            depth[act_medium_scatter] += 1
            last_scatter_it[act_medium_scatter] = dr.detach(mei)

            # Don't estimate lighting if we exceeded number of bounces
            active &= depth < self.max_depth
            act_medium_scatter &= active
            if self.handle_null_scattering:
                ray.o[act_null_scatter] = dr.detach(mei.p)
                si.t[act_null_scatter] = si.t - dr.detach(mei.t)

            weight[act_medium_scatter] *= mei.sigma_s / dr.detach(scatter_prob)
            throughput[active_medium] *= dr.detach(weight)

            mei = dr.detach(mei)
            phase_ctx = mi.PhaseFunctionContext(sampler)
            phase = mei.medium.phase_function()
            phase[~act_medium_scatter] = dr.zeros(mi.PhaseFunctionPtr)

            valid_ray |= act_medium_scatter
            with dr.suspend_grad():
                wo, phase_pdf = phase.sample(
                    phase_ctx,
                    mei,
                    sampler.next_1d(act_medium_scatter),
                    sampler.next_2d(act_medium_scatter),
                    act_medium_scatter,
                )
            act_medium_scatter &= phase_pdf > 0.0
            new_ray = mei.spawn_ray(wo)
            ray[act_medium_scatter] = new_ray
            needs_intersection |= act_medium_scatter
            last_scatter_direction_pdf[act_medium_scatter] = phase_pdf

            # --------------------- Surface Interactions ---------------------
            active_surface |= escaped_medium
            intersect = active_surface & needs_intersection
            si[intersect] = scene.ray_intersect(ray, intersect)

            # ---------------- Intersection with emitters ----------------
            ray_from_camera = active_surface & dr.eq(depth, 0)
            count_direct = ray_from_camera  # | specular_chain
            emitter = si.emitter(scene)
            active_e = (
                active_surface
                & dr.neq(emitter, None)
                & ~(dr.eq(depth, 0) & self.hide_emitters)
            )

            # Get the PDF of sampling this emitter using next event estimation
            ds = mi.DirectionSample3f(scene, si, last_scatter_it)
            if self.use_nee:
                emitter_pdf = scene.pdf_emitter_direction(last_scatter_it, ds, active_e)
            else:
                emitter_pdf = 0.0
            emitted = emitter.eval(si, active_e)
            contrib = dr.select(
                count_direct,
                throughput * emitted,
                throughput
                * mis_weight(last_scatter_direction_pdf, emitter_pdf)
                * emitted,
            )
            result[active_e] += dr.detach(contrib)

            active_surface &= si.is_valid()
            ctx = mi.BSDFContext()
            bsdf = si.bsdf(ray)

            # --------------------- Emitter sampling ---------------------
            if self.use_nee:
                active_e_surface = (
                    active_surface
                    & mi.has_flag(bsdf.flags(), mi.BSDFFlags.Smooth)
                    & (depth + 1 < self.max_depth)
                )
                sample_emitters = mei.medium.use_emitter_sampling()
                # specular_chain &= ~act_medium_scatter
                # specular_chain |= act_medium_scatter & ~sample_emitters
                active_e_medium = act_medium_scatter & sample_emitters
                active_e = active_e_surface | active_e_medium
                ref_interaction = dr.zeros(mi.Interaction3f)
                ref_interaction[act_medium_scatter] = mei
                ref_interaction[active_surface] = si
                nee_sampler = sampler  # if is_primal else sampler.clone()
                emitted, ds = self.sample_emitter(
                    ref_interaction,
                    scene,
                    nee_sampler,
                    medium,
                    channel,
                    active_e,
                )
                # Query the BSDF for that emitter-sampled direction
                bsdf_val, bsdf_pdf = bsdf.eval_pdf(
                    ctx, si, si.to_local(ds.d), active_e_surface
                )
                phase_val = phase.eval(phase_ctx, mei, ds.d, active_e_medium)
                nee_weight = dr.select(active_e_surface, bsdf_val, phase_val)
                nee_directional_pdf = dr.select(
                    ds.delta, 0.0, dr.select(active_e_surface, bsdf_pdf, phase_val)
                )

                contrib = (
                    throughput
                    * nee_weight
                    * mis_weight(ds.pdf, nee_directional_pdf)
                    * emitted
                )
                result[active_e] += dr.detach(contrib)

            # ----------------------- BSDF sampling ----------------------
            with dr.suspend_grad():
                bs, bsdf_weight = bsdf.sample(
                    ctx,
                    si,
                    sampler.next_1d(active_surface),
                    sampler.next_2d(active_surface),
                    active_surface,
                )
                active_surface &= bs.pdf > 0

            bsdf_eval = bsdf.eval(ctx, si, bs.wo, active_surface)

            throughput[active_surface] *= bsdf_weight
            eta[active_surface] *= bs.eta
            bsdf_ray = si.spawn_ray(si.to_world(bs.wo))
            ray[active_surface] = bsdf_ray

            needs_intersection |= active_surface
            non_null_bsdf = active_surface & ~mi.has_flag(
                bs.sampled_type, mi.BSDFFlags.Null
            )
            depth[non_null_bsdf] += 1

            # update the last scatter PDF event if we encountered a non-null scatter event
            last_scatter_it[non_null_bsdf] = si
            last_scatter_direction_pdf[non_null_bsdf] = bs.pdf

            valid_ray |= non_null_bsdf
            # specular_chain |= non_null_bsdf & mi.has_flag(
            #     bs.sampled_type, mi.BSDFFlags.Delta
            # )
            # specular_chain &= ~(
            #     active_surface & mi.has_flag(bs.sampled_type, mi.BSDFFlags.Smooth)
            # )
            has_medium_trans = active_surface & si.is_medium_transition()
            medium[has_medium_trans] = si.target_medium(ray.d)
            active &= active_surface | active_medium

        return result, valid_ray, result

    def sample_emitter(
        self,
        ref_interaction,
        scene,
        sampler,
        medium,
        channel,
        active,
    ):
        active = mi.Bool(active)
        medium = dr.select(active, medium, dr.zeros(mi.MediumPtr))

        ds, emitter_val = scene.sample_emitter_direction(
            ref_interaction, sampler.next_2d(active), False, active
        )
        ds = dr.detach(ds)
        invalid = dr.eq(ds.pdf, 0.0)
        emitter_val[invalid] = 0.0
        active &= ~invalid

        ray = ref_interaction.spawn_ray(ds.d)
        total_dist = mi.Float(0.0)
        si = dr.zeros(mi.SurfaceInteraction3f)
        needs_intersection = mi.Bool(True)
        transmittance = mi.Spectrum(1.0)
        loop = mi.Loop(
            name=f"VolPathGrid Next Event Estimation",
            state=lambda: (
                sampler,
                active,
                medium,
                ray,
                total_dist,
                needs_intersection,
                si,
                transmittance,
            ),
        )
        while loop(active):
            remaining_dist = ds.dist * (1.0 - mi.math.ShadowEpsilon) - total_dist
            ray.maxt = dr.detach(remaining_dist)
            active &= remaining_dist > 0.0

            # This ray will not intersect if it reached the end of the segment
            needs_intersection &= active
            si[needs_intersection] = scene.ray_intersect(ray, needs_intersection)
            needs_intersection &= False

            active_medium = active & dr.neq(medium, None)
            active_surface = active & ~active_medium

            # Handle medium interactions / transmittance
            mei = medium.sample_interaction(
                ray, sampler.next_1d(active_medium), channel, active_medium
            )
            mei.t[active_medium & (si.t < mei.t)] = dr.inf
            mei.t = dr.detach(mei.t)

            tr_multiplier = mi.Spectrum(1.0)

            # Special case for homogeneous media: directly advance to the next
            # surface / end of the segment
            if self.nee_handle_homogeneous:
                active_homogeneous = active_medium & medium.is_homogeneous()
                mei.t[active_homogeneous] = dr.minimum(remaining_dist, si.t)
                tr_multiplier[active_homogeneous] = medium.eval_tr_and_pdf(
                    mei, si, active_homogeneous
                )[0]
                mei.t[active_homogeneous] = dr.inf

            escaped_medium = active_medium & ~mei.is_valid()

            # Ratio tracking transmittance computation
            active_medium &= mei.is_valid()
            ray.o[active_medium] = dr.detach(mei.p)
            si.t[active_medium] = dr.detach(si.t - mei.t)
            tr_multiplier[active_medium] *= mei.sigma_n / mei.combined_extinction

            # Handle interactions with surfaces
            active_surface |= escaped_medium
            active_surface &= si.is_valid() & ~active_medium
            bsdf = si.bsdf(ray)
            bsdf_val = bsdf.eval_null_transmission(si, active_surface)
            tr_multiplier[active_surface] = tr_multiplier * bsdf_val

            transmittance *= dr.detach(tr_multiplier)

            # Update the ray with new origin & t parameter
            new_ray = si.spawn_ray(mi.Vector3f(ray.d))
            ray[active_surface] = dr.detach(new_ray)
            ray.maxt = dr.detach(remaining_dist)
            needs_intersection |= active_surface

            # Continue tracing through scene if non-zero weights exist
            active &= (active_medium | active_surface) & dr.any(
                dr.neq(transmittance, 0.0)
            )
            total_dist[active] += dr.select(active_medium, mei.t, si.t)

            # If a medium transition is taking place: Update the medium pointer
            has_medium_trans = active_surface & si.is_medium_transition()
            medium[has_medium_trans] = si.target_medium(ray.d)

        return emitter_val * transmittance, ds

    def to_string(self):
        return "\n".join(
            [
                "VolPathGridIntegratorPy[",
                f"  max_depth = {self.max_depth},",
                f"  rr_depth = {self.rr_depth},",
                f"  hide_emitters = {self.hide_emitters}",
                "]",
            ]
        )


mi.register_integrator("volpathgridpy", lambda props: VolPathGridIntegratorPy(props))
