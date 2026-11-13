import drjit as dr
import mitsuba as mi
import numpy as np
import pytest

MAPPINGS = ["guided", "search"]
CART_BACKENDS = ["volume", "flat"]


def _build_extremum(volume, resolution, extremum_type="gridextremum"):
    """Construct an extremum structure and wire it to volume's majorant
    grid. Construction only takes resolution; the volume and its validity
    bbox are supplied afterwards via update_extremum()."""
    extremum_struct = mi.load_dict({"type": extremum_type, "resolution": resolution})
    extremum_struct.update_extremum(volume.bbox(), volume, None)
    return extremum_struct


def _volume_backend(values):
    data = np.array(values, dtype=float).reshape(len(values), 1, 1)
    return mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(data),
            "filter_type": "nearest",
            "accel": False,
            "raw": True,
        }
    )


def _make_volume(boundaries_z, values, mapping="guided", backend="volume"):
    d = {
        "type": "mapcoordsvolume",
        "mapping": mapping,
        "backend": backend,
        "boundaries_z": mi.TensorXf(np.array(boundaries_z, dtype=float)),
    }
    if backend == "volume":
        d["volume"] = _volume_backend(values)
    elif backend == "flat":
        d["data"] = mi.TensorXf(np.array(values, dtype=float).reshape(len(values), 1, 1))
    return mi.load_dict(d)


@pytest.mark.parametrize("backend", CART_BACKENDS)
@pytest.mark.parametrize("mapping", MAPPINGS)
def test_construct(variant_scalar_rgb, mapping, backend):
    """Boundary arrays avoid length 3, which the dict parser mistakes for an RGB color."""
    volume = _make_volume([0.0, 0.3, 0.7, 1.0], [1.0, 2.0, 3.0], mapping, backend)
    assert volume is not None

    with pytest.raises(Exception):
        _make_volume([1.0], [1.0], mapping, backend)

    with pytest.raises(Exception):
        _make_volume([0.0, 0.5, 0.5, 1.0], [1.0, 2.0, 3.0], mapping, backend)

    with pytest.raises(Exception):
        _make_volume([0.1, 0.4, 0.7, 1.0], [1.0, 2.0, 3.0], mapping, backend)
    with pytest.raises(Exception):
        _make_volume([0.0, 0.3, 0.6, 0.9], [1.0, 2.0, 3.0], mapping, backend)


def test_construct_unknown_mapping(variant_scalar_rgb):
    with pytest.raises(Exception):
        _make_volume([0.0, 1.0], [1.0], mapping="nonsense", backend="volume")


def test_construct_unknown_backend(variant_scalar_rgb):
    with pytest.raises(Exception):
        _make_volume([0.0, 1.0], [1.0], mapping="guided", backend="nonsense")


@pytest.mark.parametrize("mapping", MAPPINGS)
def test_construct_stale_property(variant_scalar_rgb, mapping):
    """data/volume belong to one backend each; expand() leaves the other one unreferenced."""
    with pytest.raises(Exception):
        mi.load_dict(
            {
                "type": "mapcoordsvolume",
                "mapping": mapping,
                "backend": "volume",
                "boundaries_z": mi.TensorXf([0.0, 1.0]),
                "volume": _volume_backend([1.0]),
                "data": mi.TensorXf(np.array([1.0]).reshape(1, 1, 1)),
            }
        )
    with pytest.raises(Exception):
        mi.load_dict(
            {
                "type": "mapcoordsvolume",
                "mapping": mapping,
                "backend": "flat",
                "boundaries_z": mi.TensorXf([0.0, 1.0]),
                "data": mi.TensorXf(np.array([1.0]).reshape(1, 1, 1)),
                "volume": _volume_backend([1.0]),
            }
        )


@pytest.mark.parametrize("backend", CART_BACKENDS)
@pytest.mark.parametrize("mapping", MAPPINGS)
def test_traverse_update(variant_scalar_rgb, mapping, backend):
    """mi.traverse() must expose the backend's data as a live, updatable
    parameter"""
    boundaries_z = [0.0, 0.3, 0.7, 1.0]
    values = [1.0, 2.0, 3.0]
    volume = _make_volume(boundaries_z, values, mapping, backend)

    key = "data" if backend == "flat" else "volume.data"
    params = mi.traverse(volume)
    assert key in params

    it = mi.Interaction3f()
    it.p = [0.5, 0.5, 0.15]  # cell 0 midpoint
    assert dr.allclose(volume.eval_1(it), 1.0)

    new_values = np.array([10.0, 20.0, 30.0], dtype=float).reshape(3, 1, 1)
    params[key] = mi.TensorXf(new_values)
    params.update()

    assert dr.allclose(volume.eval_1(it), 10.0)


@pytest.mark.parametrize("backend", CART_BACKENDS)
@pytest.mark.parametrize("mapping", MAPPINGS)
def test_single_cell(variant_scalar_rgb, mapping, backend):
    """Degenerate one-cell profile: boundaries_z=[0, 1], a single value.
    Binary search and the bit-pattern guide table must both handle a
    single cell without a dedicated special case."""
    volume = _make_volume([0.0, 1.0], [1.0], mapping, backend)

    for z in [0.0, 0.3, 0.5, 0.7, 1.0]:
        it = mi.Interaction3f()
        it.p = [0.5, 0.5, z]
        assert dr.allclose(volume.eval_1(it), 1.0)

    resolution = mi.ScalarVector3i(1, 1, 2)
    extremum_struct = _build_extremum(volume, resolution)
    params = mi.traverse(extremum_struct)
    actual_resolution = params["resolution"]
    n_cells = int(actual_resolution[0]) * int(actual_resolution[1]) * int(actual_resolution[2])
    extremum_grid = params["data"].numpy().reshape(n_cells, 2)
    assert np.allclose(extremum_grid, [[1.0, 1.0]] * n_cells)


@pytest.mark.parametrize("backend", CART_BACKENDS)
@pytest.mark.parametrize("mapping", MAPPINGS)
def test_eval_matches_cell(variant_scalar_rgb, mapping, backend):
    boundaries_z = [0.0, 0.1, 0.2, 0.3, 0.4, 0.7, 1.0]
    values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    volume = _make_volume(boundaries_z, values, mapping, backend)

    cell_midpoints = [
        0.5 * (boundaries_z[i] + boundaries_z[i + 1]) for i in range(len(values))
    ]
    for z, expected in zip(cell_midpoints, values):
        it = mi.Interaction3f()
        it.p = [0.5, 0.5, z]
        assert dr.allclose(volume.eval_1(it), expected)


def test_mapping_backend_agree(variant_scalar_rgb):
    """All four (mapping, backend) combinations compute the exact same inverse mapping."""
    boundaries_z = np.concatenate([[0.0], np.sort(np.random.default_rng(0).random(9)), [1.0]])
    boundaries_z = np.unique(boundaries_z)
    values = np.random.default_rng(1).uniform(0.0, 5.0, len(boundaries_z) - 1)
    zs = np.random.default_rng(2).uniform(0.0, 1.0, 50)

    reference = None
    for mapping in MAPPINGS:
        for backend in CART_BACKENDS:
            volume = _make_volume(list(boundaries_z), list(values), mapping, backend)
            results = []
            for z in zs:
                it = mi.Interaction3f()
                it.p = [0.5, 0.5, float(z)]
                results.append(float(volume.eval_1(it)))
            results = np.array(results)
            if reference is None:
                reference = results
            else:
                assert np.allclose(results, reference), (mapping, backend)


@pytest.mark.parametrize("backend", CART_BACKENDS)
@pytest.mark.parametrize("mapping", MAPPINGS)
def test_extremum_irregular_profile(variant_scalar_rgb, mapping, backend):
    """Regions must narrow through the mapping before bounding the backend."""
    boundaries_z = [0.0, 0.1, 0.2, 0.3, 0.4, 0.7, 1.0]
    values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    volume = _make_volume(boundaries_z, values, mapping, backend)

    resolution = mi.ScalarVector3i(1, 1, 2)
    extremum_struct = _build_extremum(volume, resolution)
    extremum_grid = mi.traverse(extremum_struct)["data"].numpy().reshape(2, 2)

    expected = np.array([[1.0, 5.0], [5.0, 6.0]])
    assert np.allclose(extremum_grid, expected)


@pytest.mark.parametrize("backend", CART_BACKENDS)
@pytest.mark.parametrize("mapping", MAPPINGS)
@pytest.mark.parametrize("seed", range(3))
def test_extremum_safety_random(variant_scalar_rgb, mapping, backend, seed):
    """Every sampled value lies within its region's reported minorant/majorant."""
    rng = np.random.default_rng(seed)

    n_cells = int(rng.integers(3, 12))
    raw = rng.random(n_cells - 1)
    boundaries_z = np.concatenate([[0.0], np.sort(raw), [1.0]])
    boundaries_z = np.unique(boundaries_z)
    if len(boundaries_z) < 3:
        pytest.skip("degenerate boundary draw")
    n_cells = len(boundaries_z) - 1

    values = rng.uniform(0.0, 5.0, n_cells)
    volume = _make_volume(list(boundaries_z), list(values), mapping, backend)

    n_extremum_cells = int(rng.integers(2, 6))
    resolution = mi.ScalarVector3i(1, 1, n_extremum_cells)
    extremum_struct = _build_extremum(volume, resolution)
    params = mi.traverse(extremum_struct)
    # gridextremum clamps resolution to the underlying volume's own
    # resolution, so the actual cell count may be smaller than requested.
    actual_resolution = params["resolution"]
    n_extremum_cells = int(actual_resolution[2])
    extremum_grid = params["data"].numpy().reshape(n_extremum_cells, 2)

    for i in range(n_extremum_cells):
        lo = i / n_extremum_cells
        hi = (i + 1) / n_extremum_cells
        minorant, majorant = extremum_grid[i]

        samples_z = rng.uniform(lo, hi, 200)
        for z in samples_z:
            it = mi.Interaction3f()
            it.p = [0.5, 0.5, float(z)]
            value = float(volume.eval_1(it))
            assert minorant - 1e-5 <= value <= majorant + 1e-5, (
                f"seed={seed}, mapping={mapping}, backend={backend}, "
                f"cell={i}, z={z}: value {value} outside [{minorant}, {majorant}]"
            )


def _make_wrap_volume(values_x, wrap_mode=None):
    data = np.array(values_x, dtype=float).reshape(1, 1, len(values_x))
    d = {
        "type": "mapcoordsvolume",
        "mapping": "guided",
        "backend": "flat",
        "boundaries_z": mi.TensorXf(np.array([0.0, 1.0])),
        "data": mi.TensorXf(data),
    }
    if wrap_mode is not None:
        d["wrap_mode"] = wrap_mode
    return mi.load_dict(d)


def test_wrap_mode_default_is_clamp(variant_scalar_rgb):
    volume = _make_wrap_volume([10.0, 20.0, 30.0, 40.0])
    it = mi.Interaction3f()
    it.p = [1.4, 0.5, 0.5]
    assert dr.allclose(volume.eval_1(it), 40.0)


@pytest.mark.parametrize(
    "wrap_mode,x,expected",
    [
        ("clamp", 1.4, 40.0),
        ("repeat", 1.4, 20.0),
        ("mirror", 1.4, 30.0),
        ("clamp", -0.4, 10.0),
        ("repeat", -0.4, 30.0),
        ("mirror", -0.4, 20.0),
    ],
)
def test_wrap_mode_out_of_range(variant_scalar_rgb, wrap_mode, x, expected):
    """Out-of-[0,1] X queries resolve per wrap_mode; gridvolume follows the same convention."""
    volume = _make_wrap_volume([10.0, 20.0, 30.0, 40.0], wrap_mode=wrap_mode)
    it = mi.Interaction3f()
    it.p = [x, 0.5, 0.5]
    assert dr.allclose(volume.eval_1(it), expected)


def _make_spherical_volume(
    boundaries_x, values, rmin=0.0, rmax=1.0, fillmin=0.0, fillmax=0.0, to_world=None
):
    d = {
        "type": "mapcoordsvolume",
        "mapping": "spherical",
        "backend": "flat",
        "boundaries_x": mi.TensorXf(np.array(boundaries_x, dtype=float)),
        "data": mi.TensorXf(np.array(values, dtype=float).reshape(1, 1, len(values))),
        "rmin": rmin,
        "rmax": rmax,
        "fillmin": fillmin,
        "fillmax": fillmax,
    }
    if to_world is not None:
        d["to_world"] = to_world
    return mi.load_dict(d)


def test_spherical_matches_sphericalcoordsvolume(variant_scalar_rgb):
    """mapcoordsvolume(mapping="spherical") with a trivial radial boundary
    array must match sphericalcoordsvolume over its full domain: inside
    the shell, below rmin and above rmax."""
    values = [1.0, 5.0, 2.0, 8.0, 3.0, 9.0, 4.0, 6.0]
    data_flat = np.array(values).reshape(1, 1, len(values))

    rmin, rmax = 0.1, 0.4
    fillmin, fillmax = -7.0, -8.0

    ours = _make_spherical_volume([0.0, 1.0], values, rmin, rmax, fillmin, fillmax)

    ref_volume = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(data_flat),
            "filter_type": "nearest",
            "accel": False,
            "raw": True,
        }
    )
    ref = mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": ref_volume,
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": fillmin,
            "fillmax": fillmax,
        }
    )

    points = [
        (0.0, 0.0, 0.0),      # r=0, below rmin
        (0.05, 0.0, 0.0),     # below rmin
        (0.1, 0.0, 0.0),      # exactly rmin
        (0.25, 0.0, 0.0),     # mid-shell
        (0.4, 0.0, 0.0),      # exactly rmax
        (0.6, 0.0, 0.0),      # above rmax
        (0.0, 0.25, 0.0),     # mid-shell, off x-axis
        (0.0, 0.0, 0.25),     # mid-shell, +Z pole (theta=0)
        (0.0, 0.0, -0.25),    # mid-shell, -Z pole (theta=pi)
        (-0.2, -0.2, 0.1),    # mid-shell, negative coordinates
        (0.15, -0.1, 0.2),    # mid-shell, general off-axis point
        (0.5, 0.5, 0.5),      # above rmax
        (-0.5, -0.5, -0.5),   # above rmax, negative coordinates
    ]
    for p in points:
        it = mi.Interaction3f()
        it.p = list(p)
        assert dr.allclose(ours.eval_1(it), ref.eval_1(it)), p


def test_spherical_irregular_radius_shifts_result(variant_scalar_rgb):
    """Irregular boundaries_x with a spherical mapping must remap the radial axis"""
    values = [1.0, 2.0, 3.0]
    trivial = _make_spherical_volume([0.0, 1.0], values)
    irregular = _make_spherical_volume([0.0, 0.1, 0.9, 1.0], values)

    it = mi.Interaction3f()
    it.p = [0.3, 0.0, 0.0]  # r = 0.3, mid-shell either way

    assert not dr.allclose(trivial.eval_1(it), irregular.eval_1(it))


def test_spherical_extremum_matches_sphericalcoordsvolume(variant_scalar_rgb):
    """sphericalcoordsextremum built over mapcoordsvolume(mapping="spherical")
    with a trivial radial boundary array must match the same structure
    built over sphericalcoordsvolume, shell by shell."""
    values = [1.0, 5.0, 2.0, 8.0]
    data_flat = np.array(values).reshape(1, 1, len(values))
    rmin, rmax = 0.0, 1.0
    fillmin, fillmax = -5.0, -6.0

    ours = _make_spherical_volume([0.0, 1.0], values, rmin, rmax, fillmin, fillmax)

    ref_volume = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(data_flat),
            "filter_type": "nearest",
            "accel": False,
            "raw": True,
        }
    )
    ref = mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": ref_volume,
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": fillmin,
            "fillmax": fillmax,
        }
    )

    resolution = mi.ScalarVector3i(len(values), 1, 1)
    ours_extremum = _build_extremum(ours, resolution, extremum_type="sphericalcoordsextremum")
    ref_extremum = _build_extremum(ref, resolution, extremum_type="sphericalcoordsextremum")

    ours_grid = mi.traverse(ours_extremum)["data"].numpy()
    ref_grid = mi.traverse(ref_extremum)["data"].numpy()
    assert np.allclose(ours_grid, ref_grid)


@pytest.mark.slow
@pytest.mark.parametrize("backend", CART_BACKENDS)
@pytest.mark.parametrize("mapping", MAPPINGS)
def test_render_beer_lambert(mapping, backend):
    """Beer-Lambert closed-form transmittance check through the full render path, scalar and LLVM."""
    boundaries_z = np.array([0.0, 0.1, 0.25, 0.4, 0.6, 0.8, 1.0])
    sigma = np.array([0.2, 1.2, 0.3, 2.5, 0.6, 1.8])
    tau = float(np.sum(sigma * np.diff(boundaries_z)))
    expected = np.exp(-tau)

    for variant in ["scalar_rgb", "llvm_ad_rgb"]:
        if variant.startswith("llvm") and not dr.has_backend(dr.JitBackend.LLVM):
            continue
        mi.set_variant(variant)

        sigma_t = _make_volume(boundaries_z.tolist(), sigma.tolist(), mapping, backend)
        scene = mi.load_dict(
            {
                "type": "scene",
                "integrator": {"type": "volpath", "max_depth": 2},
                "sensor": {
                    "type": "orthographic",
                    "to_world": mi.ScalarTransform4f.look_at(
                        origin=[0.5, 0.5, 3.0], target=[0.5, 0.5, 0.5], up=[0, 1, 0]
                    )
                    @ mi.ScalarTransform4f.scale(0.3),
                    "film": {
                        "type": "hdrfilm",
                        "width": 8,
                        "height": 8,
                        "pixel_format": "luminance",
                        "rfilter": {"type": "box"},
                    },
                    "sampler": {"type": "independent", "sample_count": 4096},
                },
                "light": {"type": "constant", "radiance": {"type": "rgb", "value": 1.0}},
                "medium_shape": {
                    "type": "cube",
                    "to_world": mi.ScalarTransform4f.translate([0.5, 0.5, 0.5])
                    @ mi.ScalarTransform4f.scale(0.5),
                    "bsdf": {"type": "null"},
                    "interior": {
                        "type": "heterogeneous",
                        "sigma_t": sigma_t,
                        "albedo": 0.0,
                    },
                },
            }
        )
        img = np.array(mi.render(scene, seed=0))
        assert np.isclose(img.mean(), expected, rtol=0.05), (
            f"variant={variant}, mapping={mapping}, backend={backend}: "
            f"{img.mean()} vs {expected}"
        )
