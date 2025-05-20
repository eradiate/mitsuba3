import drjit as dr
import mitsuba as mi
import numpy as np
import pytest

_bsdf_dict = {
    "type": "maignan",
    "rho_0": 0.1,
    "k": 0.1,
    "g": 0.0,
    "ndvi": 0.0,
    "ref_re": 1.5,
    "ref_im": 0.0
}


@pytest.mark.slow
def test_chi2_maignan(variants_vec_backends_once_rgb):
    """
    Test the consistency of the oceanic BSDF using the chi2 test.
    """
    sample_func, pdf_func = mi.chi2.BSDFAdapter("ocean_maignan", _bsdf_dict)

    chi2 = mi.chi2.ChiSquareTest(
        domain=mi.chi2.SphericalDomain(),
        sample_func=sample_func,
        pdf_func=pdf_func,
        sample_dim=3,
        ires=16,
        res=201,
    )

    assert chi2.run()


def test_create_oceanic(variants_vec_backends_once_rgb):
    # Test constructor of oceanic BSDF
    brdf = mi.load_dict(_bsdf_dict)
    gloss = brdf.flags()

    # Obtain binary Mitsuba flags
    gloss_flag = mi.BSDFFlags.GlossyReflection | mi.BSDFFlags.FrontSide

    assert isinstance(brdf, mi.BSDF)
    assert gloss == gloss_flag


def test_traverse_oceanic(variant_scalar_mono_polarized):
    # Mishchenko reference :
    # wavelenth 550, windspeed 2, eta 1.33, vza 15, vaa 0, sza 15, saa 180
    ref_550_2 = [
        [0.125155, -0.0132689, 0.0, 0.0],
        [-0.0132689, 0.125155, 0.0, -0.0],
        [0.0, 0.0, -0.124450, 0.0],
        [0.0, 0.0, -0.0, -0.124450],
    ]

    # wavelenth 900, windspeed 10, eta 1.39, vza 60, vaa 0, sza 40, saa 170
    ref2_900_10 = [
        [0.733924e-01, -0.713385e-01, 0.000000e00, -0.000000e00],
        [-0.713385e-01, 0.733924e-01, 0.000000e00, -0.000000e00],
        [0.000000e00, 0.000000e00, -0.172412e-01, 0.000000e00],
        [0.000000e00, 0.000000e00, -0.000000e00, -0.172412e-01],
    ]

    def sph_to_eucl(theta, phi):
        """angles in radians, return Vector3f"""
        x = dr.sin(theta) * dr.cos(phi)
        y = dr.sin(theta) * dr.sin(phi)
        z = dr.cos(theta)
        return mi.Vector3f(x, y, z)

    vza = 15
    vaa = 0.0
    sza = 15
    saa = 180.0

    wi = sph_to_eucl(dr.deg2rad(vza), dr.deg2rad(vaa))
    wo = sph_to_eucl(dr.deg2rad(sza), dr.deg2rad(saa))
    n = dr.zeros(mi.Vector3f, dr.width(wi))
    n.z = 1.0
    si = dr.zeros(mi.SurfaceInteraction3f, dr.width(wi))
    si.wi = wi
    si.n = n

    ctx = mi.BSDFContext(mi.TransportMode.Radiance)

    brdf = mi.load_dict(_bsdf_dict)

    brdf_dr = brdf.eval(ctx, si, wo)
    brdf_np = brdf_dr.numpy()[0]

    assert dr.allclose(brdf_np, ref_550_2, 0.0001, 0.00001)

    params = mi.traverse(brdf)

    params["ndvi"] = 0.2
    params["ref_re"] = 1.7
    params.update()

    vza = 60
    vaa = 0.0
    sza = 40
    saa = 180.0

    wo = sph_to_eucl(dr.deg2rad(sza), dr.deg2rad(saa))
    si.wi = sph_to_eucl(dr.deg2rad(vza), dr.deg2rad(vaa))

    brdf_dr = brdf.eval(ctx, si, wo)
    brdf_np = brdf_dr.numpy()[0]

    assert dr.allclose(brdf_np, ref2_900_10, 0.001, 0.0001)
