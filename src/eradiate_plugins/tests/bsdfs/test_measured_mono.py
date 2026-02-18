import os

import drjit as dr
import mitsuba as mi
import pytest
from mitsuba.scalar_rgb.test.util import fresolver_append_path


@fresolver_append_path
def _lookup_sample_data_path():
    if not mi.variant():
        mi.set_variant("scalar_rgb")

    with mi.variant_context("scalar_rgb"):
        fresolver = mi.file_resolver()
        fname = str(
            fresolver.resolve("src/eradiate_plugins/tests/bsdfs/cm_white_spec.bsdf")
        )

    if not os.path.isfile(fname):
        return None
    else:
        return fname


CM_WHITE_PATH = _lookup_sample_data_path()
SKIPIF_REASON = (
    "Sample data file 'src/eradiate_plugins/tests/bsdfs/cm_white_spec.bsdf' "
    "not found; download it from the RGL material database and place it at the "
    "right location."
)
# DB link: https://rgl.epfl.ch/materials
# File link: https://d38rqfq1h7iukm.cloudfront.net/media/materials/cm_white/cm_white_spec.bsdf


@pytest.mark.skipif(CM_WHITE_PATH is None, reason=SKIPIF_REASON)
def test01_construct(variant_scalar_mono):
    # Check if the plugin can be instantiate and the wavelength set to the
    # requested value
    bsdf = mi.load_dict(
        {"type": "measured_mono", "filename": CM_WHITE_PATH, "wavelength": 500.0}
    )
    assert "wavelength = 500" in str(bsdf)


@pytest.mark.skipif(CM_WHITE_PATH is None, reason=SKIPIF_REASON)
def test02_traverse(variant_scalar_mono):
    # Check if wavelength updates are applied as expected
    bsdf = mi.load_dict(
        {"type": "measured_mono", "filename": CM_WHITE_PATH, "wavelength": 550.0}
    )

    for w in [330.0, 440.0, 550.0, 660.0]:
        params = mi.traverse(bsdf)
        params.update({"wavelength": w})
        assert f"wavelength = {w:.0f}" in str(bsdf)


@pytest.mark.skipif(CM_WHITE_PATH is None, reason=SKIPIF_REASON)
def test03_mono_vs_spectral():
    # Check if mono and spectral variants produce consistent eval() and sample()
    # outputs
    wavelengths = [330.0, 440.0, 550.0, 660.0]
    bsdf_dict = {"type": "measured_mono", "filename": CM_WHITE_PATH}

    with mi.variant_context("scalar_mono"):
        bsdf = mi.load_dict(bsdf_dict)
        ctx = mi.BSDFContext()
        si = dr.zeros(mi.SurfaceInteraction3f)
        si.wi = mi.ScalarVector3f(0, 0, 1)
        wo = mi.ScalarVector3f(0, 0, 1)
        sample1 = 0.5
        sample2 = mi.Point2f(0.25, 0.75)

        eval_mono = []
        sample_weight_mono = []
        sample_wo_mono = None
        sample_pdf_mono = None

        for w in wavelengths:
            params = mi.traverse(bsdf)
            params.update({"wavelength": w})
            eval_mono.append(bsdf.eval(ctx, si, wo)[0])
            bs, sample_weight = bsdf.sample(ctx, si, sample1, sample2)

            sample_weight_mono.append(sample_weight[0])

            if sample_wo_mono is None:
                sample_wo_mono = bs.wo
            else:
                assert dr.allclose(bs.wo, sample_wo_mono), f"{w = }"

            if sample_pdf_mono is None:
                sample_pdf_mono = bs.pdf
            else:
                assert dr.allclose(bs.pdf, sample_pdf_mono), f"{w = }"

    with mi.variant_context("scalar_spectral"):
        bsdf = mi.load_dict(bsdf_dict)
        ctx = mi.BSDFContext()
        si = dr.zeros(mi.SurfaceInteraction3f)
        si.wi = mi.ScalarVector3f(0, 0, 1)
        wo = mi.ScalarVector3f(0, 0, 1)
        si.wavelengths = mi.Spectrum(wavelengths)
        sample1 = 0.5
        sample2 = mi.Point2f(0.25, 0.75)

        eval_spectral = bsdf.eval(ctx, si, wo)
        bs, sample_weight_spectral = bsdf.sample(ctx, si, sample1, sample2)
        sample_wo_spectral = bs.wo
        sample_pdf_spectral = bs.pdf

    assert dr.allclose(eval_mono, eval_spectral)
    assert dr.allclose(sample_wo_mono, sample_wo_spectral)
    assert dr.allclose(sample_pdf_mono, sample_pdf_spectral)
    assert dr.allclose(sample_weight_mono, sample_weight_spectral)
