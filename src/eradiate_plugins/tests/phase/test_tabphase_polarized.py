import pytest
import drjit as dr
import mitsuba as mi


def test_create_m11_only(variant_scalar_mono_polarized):
    p = mi.load_dict(
        {
            "type": "tabphase_polarized",
            "nodes": "-1.0, 0.0, 1.0",
            "m11": "0.5, 1.0, 1.5",
        }
    )
    assert p is not None


def test_create_all_coefs(variant_scalar_mono_polarized):
    p = mi.load_dict(
        {
            "type": "tabphase_polarized",
            "nodes": "-1.0, 0.0, 1.0",
            "m11": "0.5, 1.0, 1.5",
            "m12": "0.5, 1.0, 1.5",
            "m33": "0.5, 1.0, 1.5",
            "m34": "0.5, 1.0, 1.5",
        }
    )
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

    res = 100
    # regular grid
    # cos_thetas = np.linspace(-1., 1., res)
    # irregular grid
    cos_thetas = np.concatenate(
        (
            np.linspace(-1, 0, int(res * 0.3), False),
            np.linspace(0.0, 1.0, int(res * 0.7)),
        )
    )

    reference_val = np.zeros((res, 4, 4))
    rho = 0.0

    for i, cos_theta in enumerate(cos_thetas):
        reference_val[i, :, :] = rayleigh_scatter(cos_theta, rho)

    m_ref = np.zeros((4, res))
    m_ref[0, :] = reference_val[:, 0, 0]
    m_ref[1, :] = reference_val[:, 0, 1]
    m_ref[2, :] = reference_val[:, 2, 2]
    m_ref[3, :] = reference_val[:, 2, 3]

    cos_theta_str = ", ".join([str(x) for x in cos_thetas])

    m_ref_str = []
    for i in range(4):
        m_ref_str.append(", ".join([str(x) for x in m_ref[i, :]]))

    tabphase_polarized = mi.load_dict(
        {
            "type": "tabphase_polarized",
            "nodes": cos_theta_str,
            "m11": m_ref_str[0],
            "m12": m_ref_str[1],
            "m22": m_ref_str[0],
            "m33": m_ref_str[2],
            "m34": m_ref_str[3],
            "m44": m_ref_str[2],
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
    data_tabphase = np.zeros((test_res, 2 * test_res, 4, 4))

    print(f"mei.wi : {mei.wi}")

    for i, theta_o in enumerate(theta_os):
        for j, phi_o in enumerate(phi_os):
            wo = sph_to_dir(theta_o, phi_o)

            data_tabphase[i, j, :, :] = tabphase_polarized.eval_pdf(
                mi.PhaseFunctionContext(None), mei, wo
            )[0]

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

    assert np.allclose(data_tabphase, data_reference, rtol=1e-4, atol=1e-4)


def test_sample(variant_scalar_mono_polarized):
    """
    Check if the sampling routine uses consistent incoming-outgoing orientation
    conventions.
    """

    tab = mi.load_dict(
        {
            "type": "tabphase_polarized",
            "nodes": "-1.0, 0.0, 1.0",
            "m11": "0.0, 0.5, 1.0",
        }
    )
    ctx = mi.PhaseFunctionContext(None)
    mei = mi.MediumInteraction3f()
    mei.t = 0.1
    mei.p = [0, 0, 0]
    mei.sh_frame = mi.Frame3f([0, 0, 1])
    mei.wi = [0, 0, 1]

    # The passed sample corresponds to forward scattering
    wo, w, pdf = tab.sample(ctx, mei, 0, (1, 0))

    # The sampled direction indicates forward scattering in the "graphics"
    # convention
    assert dr.allclose(wo, [0, 0, -1])

    # The expected value was derived manually from the PDF expression.
    # An incorrect convention (i.e. using -cos θ to fetch the PDF value) will
    # yield 0 thanks to the values used to initialize the distribution and will
    # make the test fail.
    assert dr.allclose(pdf, 0.5 / dr.pi)


def test_chi2(variant_llvm_mono_polarized):
    sample_func, pdf_func = mi.chi2.PhaseFunctionAdapter(
        "tabphase_polarized",
        "<string name='nodes' value='-1.0, 0.0, 1.0'/>"
        "<string name='m11' value='0.5, 1.0, 1.5'/>",
    )

    chi2 = mi.chi2.ChiSquareTest(
        domain=mi.chi2.SphericalDomain(),
        sample_func=sample_func,
        pdf_func=pdf_func,
        sample_dim=3,
    )

    result = chi2.run()
    # chi2._dump_tables()
    assert result


def test_traverse(variant_scalar_mono_polarized):
    # Phase function table definition
    import numpy as np

    ref_y = np.array([0.5, 1.0, 1.5])
    ref_x = np.linspace(-1, 1, len(ref_y))
    ref_integral = np.trapz(ref_y, ref_x)

    # Initialise as isotropic and update with parameters
    phase = mi.load_dict(
        {
            "type": "tabphase_polarized",
            "m11": "1, 1, 1",
            "m12": "1, 1, 1",
            "m22": "1, 1, 1",
            "m33": "1, 1, 1",
            "m34": "1, 1, 1",
            "m44": "1, 1, 1",
            "nodes": "-1, 0, 1",
        }
    )
    params = mi.traverse(phase)
    params["m11"] = [0.5, 1.0, 1.5]
    params["m12"] = [0.5, 1.0, 1.5]
    params["m22"] = [0.5, 1.0, 1.5]
    params["m33"] = [0.5, 1.0, 1.5]
    params["m34"] = [0.5, 1.0, 1.5]
    params["m44"] = [0.5, 1.0, 1.5]
    params["nodes"] = [-1.0, 0.5, 1.0]
    params.update()

    # Distribution parameters are updated
    params = mi.traverse(phase)
    assert dr.allclose(params["m11"], [0.5, 1.0, 1.5])
    assert dr.allclose(params["m12"], [0.5, 1.0, 1.5])
    assert dr.allclose(params["m33"], [0.5, 1.0, 1.5])
    assert dr.allclose(params["m34"], [0.5, 1.0, 1.5])
    assert dr.allclose(params["nodes"], [-1.0, 0.5, 1.0])

    # The plugin itself evaluates consistently
    ctx = mi.PhaseFunctionContext(None)
    mei = mi.MediumInteraction3f()
    mei.wi = np.array([0, 0.0000000001, -0.999999999])
    wo = [0, 0, 1]
    print(phase.eval_pdf(ctx, mei, wo)[0])
    a = dr.inv_two_pi * 1.5 / ref_integral
    ref = (
        mi.Spectrum(
            [
                [1.0, -1.0, 0.0, 0.0],
                [-1.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, -1.0],
                [0.0, 0.0, 1.0, 1.0],
            ]
        )
        * a
    )
    print(ref)

    assert dr.allclose(phase.eval_pdf(ctx, mei, wo)[0], ref, rtol=1e-5, atol=1e-6)

