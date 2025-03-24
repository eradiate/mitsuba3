import pytest
import drjit as dr
import mitsuba as mi


def test_create(variant_scalar_mono_polarized):
    p = mi.load_dict({"type": "rayleigh_polarized", "depolarization": 0.05})
    assert p is not None


def test_eval(variant_scalar_mono_polarized):
    import numpy as np

    def sph_to_dir(theta, phi):
        """Map spherical to Euclidean coordinates."""
        st = np.sin(theta)
        ct = np.cos(theta)
        sp = np.sin(phi)
        cp = np.cos(phi)
        return np.array([cp * st, sp * st, ct])

    # Reference implementation
    def rayleigh_scatter(cos_theta, rho=0.0):
        # Based on Hansen & Travis (1974) eq. (2.15)
        delta = (1.0 - rho) / (1.0 + 0.5 * rho)
        delta_prime = (1.0 - 2.0 * rho) / (1.0 - rho)

        cos_theta2 = cos_theta**2
        a = cos_theta2 + 1.0
        b = cos_theta2 - 1.0
        c = 2.0 * cos_theta

        norm_factor = 3.0 / (16.0 * np.pi)
        return norm_factor * np.array(
            [
                [a * delta + (1.0 - delta), b * delta, 0, 0],
                [b * delta, a * delta, 0, 0],
                [0, 0, c * delta, 0],
                [0, 0, 0, c * delta * delta_prime],
            ]
        )

    rayleigh_polarized = mi.load_dict(
        {
            "type": "rayleigh_polarized",
        }
    )

    # Create a (dummy) surface interaction to use for the evaluation of the BSDF
    mei = dr.zeros(mi.MediumInteraction3f)

    # Specify an incident direction with 45 degrees elevation
    theta_i = 0.0
    phi_i = 0.0
    mei.wi = sph_to_dir(dr.deg2rad(theta_i), dr.deg2rad(phi_i))
    mei.sh_frame = mi.Frame3f(mei.wi)

    # Create grid in spherical coordinates and map it onto the sphere
    test_res = 15
    theta_os = np.linspace(0.00001, dr.pi, test_res)
    phi_os = np.linspace(0.0, 2.0 * dr.pi, 2 * test_res)

    data_reference = np.zeros((test_res, 2 * test_res, 4, 4))
    data_rayleigh = np.zeros((test_res, 2 * test_res, 4, 4))

    for i, theta_o in enumerate(theta_os):
        for j, phi_o in enumerate(phi_os):
            wo = sph_to_dir(theta_o, phi_o)

            data_rayleigh[i, j, :, :] = rayleigh_polarized.eval_pdf(
                mi.PhaseFunctionContext(None), mei, wo
            )[0].numpy().squeeze()

            phase_ref = rayleigh_scatter(np.dot(wo, -mei.wi))

            wo = mi.Vector3f(wo)

            # rotate reference frame to stokes basis
            x_hat = dr.normalize(dr.cross(-wo, mei.wi))
            p_axis_in = dr.normalize(dr.cross(x_hat, -wo))
            p_axis_out = dr.normalize(dr.cross(x_hat, mei.wi))
            phase_ref = mi.mueller.rotate_mueller_basis(
                phase_ref,
                -wo,
                p_axis_in,
                mi.mueller.stokes_basis(-wo),
                mei.wi,
                p_axis_out,
                mi.mueller.stokes_basis(mei.wi),
            )
            data_reference[i, j, :, :] = phase_ref
    assert np.allclose(data_rayleigh, data_reference, rtol=1e-4, atol=1e-4)


def test_chi2(variant_llvm_mono_polarized):
    sample_func, pdf_func = mi.chi2.PhaseFunctionAdapter("rayleigh_polarized", "")

    chi2 = mi.chi2.ChiSquareTest(
        domain=mi.chi2.SphericalDomain(),
        sample_func=sample_func,
        pdf_func=pdf_func,
        sample_dim=3,
    )

    assert chi2.run()

