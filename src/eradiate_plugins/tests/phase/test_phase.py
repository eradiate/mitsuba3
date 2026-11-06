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


def reference_accumulate_envelope(phase, nodes):
    """Replicate the C++ ``PhaseFunction::accumulate_envelope`` evaluation for a single
    (non-composite) phase function.

    For each node ``mu`` (cos_theta in physics convention, +1 = forward
    scattering), the envelope evaluates the phase value at the reconstructed
    outgoing direction and accumulates the elementwise maximum into a
    zero-initialized buffer. We mirror exactly the direction reconstruction
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
        # value == pdf for the isotropic/hg/tabphase plugins exercized here
        out[i] = phase.eval_pdf(ctx, mei, wo)[1]
    return out


def test_get_envelope_nodes_default(variant_scalar_rgb):
    # A leaf phase function without irregular nodes falls back to the base
    # implementation: 256 nodes spread uniformly over [-1, 1].
    phase = mi.load_dict({"type": "isotropic"})
    nodes = np.array(phase.get_envelope_nodes())

    assert nodes.shape[0] == 256
    # The grid is built in single precision; compare with an absolute tolerance.
    assert np.allclose(nodes, np.linspace(-1, 1, 256), atol=1e-6)
    # Strictly increasing, bounded by [-1, 1].
    assert np.all(np.diff(nodes) > 0)


def test_get_envelope_nodes_tabphase_irregular(
    variant_scalar_rgb, tabphase_nodes, tabphase_irregular
):
    # tabphase_irregular overrides get_envelope_nodes() to expose its own grid.
    assert np.allclose(
        np.array(tabphase_irregular.get_envelope_nodes()), tabphase_nodes
    )


def test_get_envelope_nodes_blendphase_merged(
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

    nodes = np.array(blend.get_envelope_nodes())
    expected = np.unique(np.concatenate([np.linspace(-1, 1, 256), tabphase_nodes]))

    # Single-precision grid: compare with an absolute tolerance.
    assert np.allclose(nodes, expected, atol=1e-6)
    # The shared endpoints (-1, 1) collapse: 256 uniform + 2 interior nodes.
    assert nodes.shape[0] == 258
    assert np.all(np.diff(nodes) > 0)


def test_accumulate_envelope_single(variant_scalar_rgb, hg):
    # The envelope of a single phase function is just the phase function
    # evaluated at every node.
    nodes = hg.get_envelope_nodes()
    values = dr.zeros(type(nodes), dr.width(nodes))

    hg.accumulate_envelope(nodes, values)
    assert np.allclose(np.array(values), reference_accumulate_envelope(hg, nodes))


def test_accumulate_envelope_blendphase(variant_scalar_rgb, hg):
    # A composite envelope is the pointwise maximum over its (unweighted)
    # children, regardless of the blend weight.
    iso = mi.load_dict({"type": "isotropic"})
    blend = mi.load_dict(
        {"type": "blendphase", "weight": 0.3, "phase_0": iso, "phase_1": hg}
    )

    nodes = blend.get_envelope_nodes()
    values = dr.zeros(type(nodes), dr.width(nodes))
    blend.accumulate_envelope(nodes, values)

    ref = np.maximum(
        reference_accumulate_envelope(iso, nodes),
        reference_accumulate_envelope(hg, nodes),
    )
    assert np.allclose(np.array(values), ref)


def test_accumulate_envelope_multiphase(variant_scalar_rgb, hg, tabphase_irregular):
    # Same envelope semantics for the N-component mixture.
    iso = mi.load_dict({"type": "isotropic"})
    multi = mi.load_dict(
        {
            "type": "multiphase",
            "phase_0": iso,
            "weight_0": 1.0,
            "phase_1": hg,
            "weight_1": 1.0,
            "phase_2": tabphase_irregular,
            "weight_2": 1.0,
        }
    )

    nodes = multi.get_envelope_nodes()
    values = dr.zeros(type(nodes), dr.width(nodes))
    multi.accumulate_envelope(nodes, values)

    ref = np.maximum.reduce(
        [
            reference_accumulate_envelope(iso, nodes),
            reference_accumulate_envelope(hg, nodes),
            reference_accumulate_envelope(tabphase_irregular, nodes),
        ]
    )
    assert np.allclose(np.array(values), ref)


def _make_ddis_medium(medium_type, phase, ddis_threshold=0.5):
    return mi.load_dict(
        {
            "type": medium_type,
            "sigma_t": 0.5,
            "albedo": 0.8,
            "phase": phase,
            "ddis_threshold": ddis_threshold,
        }
    )


def _ddis_phase_tree(combo, shared):
    """Embed ``shared`` as the medium's phase function, either directly or as
    one branch of a composite phase function."""
    if combo == "hg":
        return shared
    if combo == "blendphase":
        return {
            "type": "blendphase",
            "weight": 0.5,
            "phase_0": shared,
            "phase_1": {"type": "isotropic"},
        }
    if combo == "multiphase":
        return {
            "type": "multiphase",
            "phase_0": shared,
            "weight_0": 1.0,
            "phase_1": {"type": "isotropic"},
            "weight_1": 1.0,
        }
    raise ValueError(combo)


@pytest.mark.parametrize(
    "medium_type", ["homogeneous", "heterogeneous", "piecewise", "eoheterogeneous"]
)
@pytest.mark.parametrize("phase_combo", ["hg", "blendphase", "multiphase"])
def test_ddis_phase_function_update(
    variant_scalar_mono_double, medium_type, phase_combo
):
    # Two media share a single leaf phase function: one uses it as its direct
    # phase function, the other reaches it only as a child of a composite
    # phase function (except for the "hg" combo, where both use it
    # directly). The layout of the scene relies on the more robust traversal
    # function of `eradiate.traverse`.
    shared = mi.load_dict({"type": "hg", "g": 0.2})

    media = {
        "medium_direct": _make_ddis_medium(medium_type, shared),
        "medium_composite": _make_ddis_medium(
            medium_type, _ddis_phase_tree(phase_combo, shared)
        ),
    }
    scene = mi.load_dict(
        {
            "type": "scene",
            **{
                name: {"type": "cube", "interior": medium}
                for name, medium in media.items()
            },
        }
    )

    params = mi.eradiate.traverse(scene)
    for name, medium in media.items():
        _assert_ddis_table(params, f"{name}.interior_medium", [medium.phase_function()])

    g_key = next(k for k in params.keys() if k.endswith(".g"))  # noqa: SIM118
    params[g_key] = 0.8
    params.update()

    for name, medium in media.items():
        # Recompute the expected DDIS envelope directly from the medium's
        # (now updated) phase function, and compare against the DDIS phase
        # function's tabulated values exposed as scene parameters.
        _assert_ddis_table(params, f"{name}.interior_medium", [medium.phase_function()])


def _scene_with_medium(medium, attach="interior"):
    """Load a scene that reaches ``medium`` through a shape's interior or
    exterior, or through a sensor. Return the scene parameters and the key
    prefix of the medium."""
    if attach == "sensor":
        scene = {
            "type": "scene",
            "sensor": {
                "type": "radiancemeter",
                "origin": [0, 0, 0],
                "direction": [0, 0, 1],
                "film": {"type": "hdrfilm", "width": 1, "height": 1},
                "medium": medium,
            },
        }
        prefix = "sensor.medium"
    else:
        scene = {"type": "scene", "box": {"type": "cube", attach: medium}}
        prefix = f"box.{attach}_medium"
    return mi.eradiate.traverse(mi.load_dict(scene)), prefix


def _ddis_table(params, prefix):
    return (
        params[f"{prefix}.ddis_phase_function.nodes"],
        params[f"{prefix}.ddis_phase_function.values"],
    )


def _assert_ddis_table(params, prefix, phases):
    """The DDIS table at ``prefix`` is the envelope of ``phases``: their merged
    nodes and, at each node, the maximum of their values."""
    storage = type(phases[0].get_envelope_nodes())
    nodes = storage(
        np.unique(np.concatenate([np.array(p.get_envelope_nodes()) for p in phases]))
    )
    values = dr.zeros(storage, dr.width(nodes))
    for phase in phases:
        phase.accumulate_envelope(nodes, values)

    got_nodes, got_values = _ddis_table(params, prefix)
    assert dr.width(got_nodes) == dr.width(nodes)
    assert dr.allclose(got_nodes, nodes)
    assert dr.allclose(got_values, values)


@pytest.mark.parametrize("attach", ["interior", "exterior", "sensor"])
def test_ddis_phase_function_scene_build(variant_scalar_mono_double, attach):
    """The scene builds the DDIS phase function of a medium it reaches through a
    shape's interior, a shape's exterior or a sensor."""
    phase = mi.load_dict({"type": "hg", "g": 0.7})
    medium = {"type": "homogeneous", "phase": phase, "ddis_threshold": 0.5}
    params, prefix = _scene_with_medium(medium, attach)
    _assert_ddis_table(params, prefix, [phase])


def test_ddis_phase_function_disabled(variant_scalar_mono_double):
    """With the default ``ddis_threshold``, no DDIS phase function is built."""
    medium = {"type": "homogeneous"}
    params, prefix = _scene_with_medium(medium)
    assert f"{prefix}.ddis_phase_function.nodes" not in params


@pytest.mark.parametrize(
    "component_thresholds, overlap_threshold, enabled",
    [
        ((-1.0, -1.0), None, False),
        ((-1.0, 0.3), None, True),  # maximum over the components
        ((0.3, 0.3), -1.0, False),  # an explicit value wins
        ((-1.0, -1.0), 0.2, True),
    ],
)
def test_ddis_threshold_overlap(
    variant_scalar_mono_double, component_thresholds, overlap_threshold, enabled
):
    """An ``overlap`` uses its own ``ddis_threshold`` when set, the largest of its
    components' otherwise, and builds a DDIS phase function only if it is
    positive."""
    overlap = {"type": "overlap"}
    for i, threshold in enumerate(component_thresholds):
        overlap[f"comp_{i}"] = {
            "type": "piecewise",
            "sigma_t": 0.5,
            "albedo": 0.8,
            "ddis_threshold": threshold,
        }
    if overlap_threshold is not None:
        overlap["ddis_threshold"] = overlap_threshold

    params, prefix = _scene_with_medium(overlap)
    assert (f"{prefix}.ddis_phase_function.nodes" in params) == enabled


def test_ddis_phase_function_overlap_matches_multiphase(variant_scalar_mono_double):
    # `overlap`'s DDIS envelope is built by merging each component's own
    # phase-function envelope. An overlap medium's aggregated DDIS
    # envelope must match a single ordinary medium whose phase function is a
    # `multiphase` combining the same per-component phase functions..
    phase0 = mi.load_dict({"type": "hg", "g": 0.3})
    phase1 = mi.load_dict({"type": "isotropic"})
    phase2 = mi.load_dict({"type": "hg", "g": -0.6})

    overlap = {
        "type": "overlap",
        "comp_a": {
            "type": "piecewise",
            "sigma_t": 0.5,
            "albedo": 0.8,
            "phase": phase0,
            "ddis_threshold": 0.5,
        },
        "comp_b": {
            "type": "piecewise",
            "sigma_t": 0.3,
            "albedo": 0.6,
            "phase": phase1,
            "ddis_threshold": 0.5,
        },
        "comp_c": {
            "type": "piecewise",
            "sigma_t": 0.7,
            "albedo": 0.9,
            "phase": phase2,
            "ddis_threshold": 0.5,
        },
    }
    multiphase_ref = {
        "type": "piecewise",
        "sigma_t": 0.5,
        "albedo": 0.8,
        "phase": {
            "type": "multiphase",
            "phase0": phase0,
            "weight_0": 1.0,
            "phase1": phase1,
            "weight_1": 2.5,
            "phase2": phase2,
            "weight_2": 0.4,
        },
        "ddis_threshold": 0.5,
    }

    overlap_params, overlap_prefix = _scene_with_medium(overlap)
    ref_params, ref_prefix = _scene_with_medium(multiphase_ref)

    got_nodes, got_values = _ddis_table(overlap_params, overlap_prefix)
    expected_nodes, expected_values = _ddis_table(ref_params, ref_prefix)

    assert dr.width(got_nodes) == dr.width(expected_nodes)
    assert dr.allclose(got_nodes, expected_nodes)
    assert dr.allclose(got_values, expected_values)

    # The scene only reaches the components through the overlap, so they get
    # no DDIS phase function of their own.
    for i in range(3):
        key = f"{overlap_prefix}.component{i}.ddis_phase_function.nodes"
        assert key not in overlap_params


def test_ddis_phase_function_overlap_update(variant_scalar_mono_double):
    """Changing a component's phase function refreshes the ``overlap``'s DDIS
    phase function, and does not build one for the component."""
    phases = [mi.load_dict({"type": "hg", "g": g}) for g in (0.3, -0.6)]
    overlap = {"type": "overlap"}
    for i, phase in enumerate(phases):
        overlap[f"comp_{i}"] = {
            "type": "piecewise",
            "sigma_t": 0.5,
            "albedo": 0.8,
            "phase": phase,
            "ddis_threshold": 0.5,
        }

    params, prefix = _scene_with_medium(overlap)
    params[f"{prefix}.component1.phase_function.g"] = 0.9
    params.update()
    _assert_ddis_table(params, prefix, phases)
    for i in range(2):
        assert f"{prefix}.component{i}.ddis_phase_function.nodes" not in params


def test_ddis_phase_function_repeat(variant_scalar_mono_double):
    """A ``repeat`` of an ``overlap`` builds its DDIS phase function from every
    component of the inner medium and refreshes it when one of them changes; the
    inner ``overlap`` builds none."""
    phases = [mi.load_dict({"type": "hg", "g": g}) for g in (0.3, -0.6)]
    inner = {"type": "overlap"}
    for i, phase in enumerate(phases):
        inner[f"comp_{i}"] = {
            "type": "eoheterogeneous",
            "sigma_t": {"type": "constvolume", "value": 1.0},
            "albedo": 0.8,
            "aabb_min": mi.ScalarPoint3f(0.0),
            "aabb_max": mi.ScalarPoint3f(1.0),
            "phase": phase,
            "ddis_threshold": 0.5,
        }
    repeat = {
        "type": "repeat",
        "inner": inner,
        "aabb_min": mi.ScalarPoint3f(-20.0),
        "aabb_max": mi.ScalarPoint3f(20.0),
    }

    params, prefix = _scene_with_medium(repeat)
    _assert_ddis_table(params, prefix, phases)

    params[f"{prefix}.inner.component1.phase_function.g"] = 0.9
    params.update()
    _assert_ddis_table(params, prefix, phases)
    assert f"{prefix}.inner.ddis_phase_function.nodes" not in params
