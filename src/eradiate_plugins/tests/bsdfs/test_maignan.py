import drjit as dr
import mitsuba as mi
import pytest

_bsdf_dict = {
    "type": "maignan",
    "C": 5.0,
    "ndvi": 0.8,
    "refr_re": 1.5,
    "refr_im": 0.0,
    "ext_ior": 1.0,
}


def test_maignan_construct(variant_scalar_rgb):
    # Test constructor of Maignan BSDF
    bsdf = mi.load_dict(_bsdf_dict)
    gloss = bsdf.flags()

    # Obtain binary Mitsuba flags
    gloss_flag = mi.BSDFFlags.GlossyReflection | mi.BSDFFlags.FrontSide

    assert isinstance(bsdf, mi.BSDF)
    assert gloss == gloss_flag


def test_maignan_eval(variant_scalar_mono_polarized):
    # wavelenth 550, C 4.98, ndvi 0.8, vza 40, vaa 0, sza 40, saa 5
    ref_evergreen_needleleaf = [
        [ 1.42013086e-02,  1.48094732e-05, -1.60061711e-06,  0.00000000e+00],
        [ 1.48624285e-05,  1.41894910e-02, -5.79237472e-04,  0.00000000e+00],
        [ 9.95336109e-07, -5.79236308e-04, -1.41894845e-02,  0.00000000e+00],
        [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -1.42013021e-02],
    ]

    def sph_to_dir(theta, phi):
        """angles in radians, return Vector3f"""
        x = dr.sin(theta) * dr.cos(phi)
        y = dr.sin(theta) * dr.sin(phi)
        z = dr.cos(theta)
        return mi.Vector3f(x, y, z)

    vza = 40
    vaa = 0
    sza = 40
    saa = 5.0

    bsdf = mi.load_dict({**_bsdf_dict, "ndvi": 0.8, "C": 4.98})

    wi = sph_to_dir(dr.deg2rad(vza), dr.deg2rad(vaa))
    wo = sph_to_dir(dr.deg2rad(sza), dr.deg2rad(saa))
    n = dr.zeros(mi.Vector3f, dr.width(wi))
    n.z = 1.0
    si = dr.zeros(mi.SurfaceInteraction3f, dr.width(wi))
    si.wi = wi
    si.n = n

    ctx = mi.BSDFContext(mi.TransportMode.Radiance)
    result = bsdf.eval(ctx, si, wo)
    assert dr.allclose(result, ref_evergreen_needleleaf)

    # wavelenth 900, C 6.66, ndvi 0.3, vza 60, vaa 0, sza 40, saa 170
    ref_savanna = [
        [ 1.9543538e-02, -3.1128337e-03, -1.1320484e-03,  0.0000000e+00],
        [ 3.1127632e-03, -1.9510504e-02, -8.9618654e-05,  0.0000000e+00],
        [-1.1322410e-03,  9.2017806e-05,  1.9293837e-02,  0.0000000e+00],
        [ 0.0000000e+00,  0.0000000e+00,  0.0000000e+00, -1.9260805e-02],
    ]

    vza = 0
    vaa = 0
    sza = 40
    saa = 170

    bsdf = mi.load_dict({**_bsdf_dict, "ndvi": 0.3, "C": 6.66})

    wi = sph_to_dir(dr.deg2rad(vza), dr.deg2rad(vaa))
    wo = sph_to_dir(dr.deg2rad(sza), dr.deg2rad(saa))
    n = dr.zeros(mi.Vector3f, dr.width(wi))
    n.z = 1.0
    si = dr.zeros(mi.SurfaceInteraction3f, dr.width(wi))
    si.wi = wi
    si.n = n

    ctx = mi.BSDFContext(mi.TransportMode.Radiance)

    result = bsdf.eval(ctx, si, wo)
    assert dr.allclose(result, ref_savanna)


@pytest.mark.slow
def test_maignan_chi2(variants_vec_backends_once_rgb):
    """
    Test the consistency of the Maignan BSDF using the chi2 test.
    """
    sample_func, pdf_func = mi.chi2.BSDFAdapter("maignan", _bsdf_dict)

    chi2 = mi.chi2.ChiSquareTest(
        domain=mi.chi2.SphericalDomain(),
        sample_func=sample_func,
        pdf_func=pdf_func,
        sample_dim=3,
        ires=16,
        res=201,
    )

    assert chi2.run()
