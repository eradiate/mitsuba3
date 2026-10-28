import drjit as dr
import mitsuba as mi
import pytest


def _component(g):
    """A component carrying a phase function that is recognisable by its value."""
    return {
        "type": "eoheterogeneous",
        "sigma_t": {"type": "constvolume", "value": 1.0},
        "albedo": 0.8,
        "phase": {"type": "hg", "g": g},
        "extremum": {"type": "globalextremum"},
    }


def _eval(phase, wo):
    ctx = mi.PhaseFunctionContext(None)
    mei = dr.zeros(mi.MediumInteraction3f)
    mei.wi = mi.Vector3f(0, 0, 1)
    return phase.eval_pdf(ctx, mei, wo)[1]


@pytest.mark.parametrize("g", [(0.4, -0.4), (0.0, 0.7)])
def test_phase_function_per_component(variants_vec_rgb, g):
    """``phase_function(mei, i)`` returns component ``i``'s phase function.

    The medium resolves them from a buffer built at load time rather than
    through a nested vcall on the components: a nested call is where Dr.Jit
    substitutes the enclosing call's symbolic ``self`` for any instance whose
    registry id collides numerically, which used to make this dispatch fail.
    """
    medium = mi.load_dict(
        {"type": "overlap", "a": _component(g[0]), "b": _component(g[1])}
    )

    wo = mi.Vector3f(0.0, 0.6, -0.8)
    phase, _ = mi.MediumPtr(medium).phase_function(
        dr.zeros(mi.MediumInteraction3f, 2), mi.UInt32([0, 1])
    )
    got = _eval(phase, wo)
    expected = mi.Float(
        [_eval(mi.load_dict({"type": "hg", "g": gi}), wo)[0] for gi in g]
    )

    assert dr.allclose(got, expected)


def test_overlap_rejects_nested_overlap(variant_scalar_rgb):
    inner = {"type": "overlap", "a": _component(0.0), "b": _component(0.0)}
    repeat = {
        "type": "repeat",
        "inner": inner,
        "aabb_min": mi.ScalarPoint3f(-20.0),
        "aabb_max": mi.ScalarPoint3f(20.0),
    }
    for nested in (inner, repeat):
        with pytest.raises(RuntimeError, match="nested overlaps"):
            mi.load_dict({"type": "overlap", "a": nested, "b": _component(0.0)})
