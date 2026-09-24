import drjit as dr
import mitsuba as mi
import numpy as np
import pytest


def _component(g):
    """A component carrying a phase function that is recognisable by its value."""
    return {
        "type": "eoheterogeneous",
        "sigma_t": {"type": "constvolume", "value": 1.0},
        "albedo": 0.8,
        "phase": {"type": "hg", "g": g},
        "extremum_structure": {"type": "extremum_global"},
    }


def _eval(phase, wo):
    ctx = mi.PhaseFunctionContext(None)
    mei = dr.zeros(mi.MediumInteraction3f)
    mei.wi = mi.Vector3f(0, 0, 1)
    return phase.eval_pdf(ctx, mei, wo)[1]


@pytest.mark.parametrize("g", [(0.4, -0.4), (0.0, 0.7)])
def test_phase_function_per_component(variants_vec_rgb, g):
    """``phase_function(i)`` returns component ``i``'s phase function.

    The medium resolves them from a buffer built at load time rather than
    through a nested vcall on the components: a nested call is where Dr.Jit
    substitutes the enclosing call's symbolic ``self`` for any instance whose
    registry id collides numerically, which used to make this dispatch fail.
    """
    medium = mi.load_dict(
        {"type": "overlap", "a": _component(g[0]), "b": _component(g[1])}
    )

    wo = mi.Vector3f(0.0, 0.6, -0.8)
    got = _eval(mi.MediumPtr(medium).phase_function(mi.UInt32([0, 1]), True), wo)
    expected = np.array(
        [_eval(mi.load_dict({"type": "hg", "g": gi}), wo)[0] for gi in g]
    )

    assert dr.allclose(got, mi.Float(expected))
