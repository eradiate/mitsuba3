import drjit as dr
import mitsuba as mi
import pytest

_bsdf_dict = {
    "type": "ocean_legacy",
    "component": 0,
    "wavelength": 1500,
    "wind_speed": 1.0,
    "wind_direction": 90.,
    "chlorinity": 19,
    "pigmentation": 0.3,
    "shadowing": False,
}


@pytest.mark.slow
def test_chi2_oceanic(variants_vec_backends_once_rgb):
    """
    Test the consistency of the oceanic BSDF using the chi2 test.
    """
    sample_func, pdf_func = mi.chi2.BSDFAdapter("ocean_legacy", {**_bsdf_dict, "shadowing": True})

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
    diff, gloss, comb = brdf.flags(0), brdf.flags(1), brdf.flags()

    # Obtain binary Mitsuba flags
    diff_flag = mi.BSDFFlags.DiffuseReflection | mi.BSDFFlags.FrontSide
    gloss_flag = mi.BSDFFlags.GlossyReflection | mi.BSDFFlags.FrontSide
    comb_flag = diff_flag | gloss_flag

    assert isinstance(brdf, mi.BSDF)
    assert diff == diff_flag
    assert gloss == gloss_flag
    assert comb == comb_flag


def test_traverse_oceanic(variants_vec_backends_once_rgb):
    # 6SV reference :
    # wavelenth 1500, windspeed 1
    ref_1500_1 = [
        1.91408132e-03,
        -5.44804487e-08,
        -5.45187861e-08,
        -5.45187861e-08,
        -5.45187861e-08,
    ]
    # wavelenth 550, windspeed 30
    ref2_550_30 = [0.11300096, 0.10733355, 0.10722339, 0.1055309, 0.11909937]

    def sph_to_eucl(theta, phi):
        """angles in radians, return Vector3f"""
        x = dr.sin(theta) * dr.cos(phi)
        y = dr.sin(theta) * dr.sin(phi)
        z = dr.cos(theta)
        return mi.Vector3f(x, y, z)

    vza = mi.Float(0.0, 22.475, 44.95, 67.425, 89.9)
    vaa = 0.0
    sza = 22.475
    saa = 0.0

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
    brdf_np = brdf_dr.numpy()[:, 0]
    assert dr.allclose(brdf_np * dr.pi, ref_1500_1, 0.001, 0.0001)

    params = mi.traverse(brdf)

    params["wavelength"] = 550.0
    params["wind_speed"] = 30.0
    params.update()

    brdf_dr = brdf.eval(ctx, si, wo)
    brdf_np = brdf_dr.numpy()[:, 0]
    assert dr.allclose(brdf_np*dr.pi, ref2_550_30, 0.001, 0.0001)
