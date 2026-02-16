import drjit as dr
import mitsuba as mi


def test_01_construct(variant_scalar_spectral):
    assert mi.load_dict({"type": "delta", "wavelength": 550.0})
    assert mi.load_dict({"type": "delta", "wavelength": 550.0, "value": 0.5})


def test_02_eval(variant_scalar_spectral):
    s = mi.load_dict({"type": "delta", "wavelength": 550.0})
    si = dr.zeros(mi.SurfaceInteraction3f)
    si.wavelengths = mi.Spectrum([400, 500, 600, 700])
    assert dr.allclose(s.eval(si), 0.0)
    assert dr.allclose(s.eval_1(si), 0.0)
    assert dr.allclose(s.eval_3(si), 0.0)


def test_03_sample(variant_scalar_spectral):
    s = mi.load_dict({"type": "delta", "wavelength": 550.0, "value": 0.5})
    si = dr.zeros(mi.SurfaceInteraction3f)
    sample = mi.Spectrum([0.2, 0.4, 0.6, 0.8])
    wavelengths, weights = s.sample_spectrum(si, sample)
    assert dr.allclose(wavelengths, 550.0)
    assert dr.allclose(weights, 0.5)
