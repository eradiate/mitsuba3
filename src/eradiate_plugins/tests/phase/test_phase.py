import drjit as dr
import mitsuba as mi
import numpy as np
import pytest


@pytest.fixture
def hg():
    return mi.load_dict({"type": "hg", "g": 0.6})


@pytest.fixture
def tabphase_nodes():
    return [-1.0, -0.3, 0.4, 1.0]


@pytest.fixture
def tabphase_irregular(tabphase_nodes):
    values = [0.5, 0.8, 1.2, 1.5]
    return mi.load_dict(
        {
            "type": "tabphase_irregular",
            "values": ",".join(map(str, values)),
            "nodes": ",".join(map(str, tabphase_nodes)),
        }
    )


def reference_eval_max(phase, nodes):
    """Replicate the C++ ``PhaseFunction::eval_max`` evaluation for a single
    (non-composite) phase function.

    For each node ``mu`` (cos_theta in physics convention, +1 = forward
    scattering), the envelope evaluates the phase value at the reconstructed
    outgoing direction and accumulates the elementwise maximum into a
    zero-initialised buffer. We mirror exactly the direction reconstruction
    used in ``src/render/phase.cpp`` so the comparison is meaningful.
    """
    ctx = mi.PhaseFunctionContext(None)
    mei = mi.MediumInteraction3f()
    mei.wi = mi.Vector3f(0, 0, 1)

    nodes_np = np.array(nodes)
    out = np.zeros(nodes_np.shape[0])
    for i, mu in enumerate(nodes_np):
        sin_theta = np.sqrt(max(0.0, 1.0 - mu * mu))
        wo = mi.Vector3f(sin_theta, 0.0, -mu)
        # value == pdf for the isotropic/hg/tabphase plugins exercised here
        out[i] = phase.eval_pdf(ctx, mei, wo)[1]
    return out


def test_get_nodes_default(variant_scalar_rgb):
    # A leaf phase function without irregular nodes falls back to the base
    # implementation: 256 nodes spread uniformly over [-1, 1].
    phase = mi.load_dict({"type": "isotropic"})
    nodes = np.array(phase.get_nodes())

    assert nodes.shape[0] == 256
    # The grid is built in single precision; compare with an absolute tolerance.
    assert np.allclose(nodes, np.linspace(-1, 1, 256), atol=1e-6)
    # Strictly increasing, bounded by [-1, 1].
    assert np.all(np.diff(nodes) > 0)


def test_get_nodes_tabphase_irregular(
    variant_scalar_rgb, tabphase_nodes, tabphase_irregular
):
    # tabphase_irregular overrides get_nodes() to expose its own grid.
    assert np.allclose(np.array(tabphase_irregular.get_nodes()), tabphase_nodes)


def test_get_nodes_blendphase_merged(
    variant_scalar_rgb, tabphase_nodes, tabphase_irregular
):
    # A composite phase function merges (sorts + de-duplicates) the node grids
    # of its children. We mix an irregular grid with the default uniform grid
    # so the union contains nodes contributed by both children.

    blend = mi.load_dict(
        {
            "type": "blendphase",
            "weight": 0.5,
            "phase_0": {"type": "isotropic"},
            "phase_1": tabphase_irregular,
        }
    )

    nodes = np.array(blend.get_nodes())
    expected = np.unique(np.concatenate([np.linspace(-1, 1, 256), tabphase_nodes]))

    # Single-precision grid: compare with an absolute tolerance.
    assert np.allclose(nodes, expected, atol=1e-6)
    # The shared endpoints (-1, 1) collapse: 256 uniform + 2 interior nodes.
    assert nodes.shape[0] == 258
    assert np.all(np.diff(nodes) > 0)


def test_eval_max_single(variant_scalar_rgb, hg):
    # The envelope of a single phase function is just the phase function
    # evaluated at every node.
    nodes = hg.get_nodes()
    values = dr.zeros(type(nodes), dr.width(nodes))

    hg.eval_max(nodes, values)
    assert np.allclose(np.array(values), reference_eval_max(hg, nodes))


def test_eval_max_blendphase(variant_scalar_rgb, hg):
    # A composite envelope is the pointwise maximum over its (unweighted)
    # children, regardless of the blend weight.
    iso = mi.load_dict({"type": "isotropic"})
    blend = mi.load_dict(
        {"type": "blendphase", "weight": 0.3, "phase_0": iso, "phase_1": hg}
    )

    nodes = blend.get_nodes()
    values = dr.zeros(type(nodes), dr.width(nodes))
    blend.eval_max(nodes, values)

    ref = np.maximum(
        reference_eval_max(iso, nodes),
        reference_eval_max(hg, nodes),
    )
    assert np.allclose(np.array(values), ref)


def test_eval_max_multiphase(variant_scalar_rgb, hg, tabphase_irregular):
    # Same envelope semantics for the N-component mixture.
    iso = mi.load_dict({"type": "isotropic"})
    multi = mi.load_dict(
        {
            "type": "multiphase",
            "phase0": iso,
            "weight0": 1.0,
            "phase1": hg,
            "weight1": 1.0,
            "phase2": tabphase_irregular,
            "weight2": 1.0,
        }
    )

    nodes = multi.get_nodes()
    values = dr.zeros(type(nodes), dr.width(nodes))
    multi.eval_max(nodes, values)

    ref = np.maximum.reduce(
        [
            reference_eval_max(iso, nodes),
            reference_eval_max(hg, nodes),
            reference_eval_max(tabphase_irregular, nodes),
        ]
    )
    assert np.allclose(np.array(values), ref)
