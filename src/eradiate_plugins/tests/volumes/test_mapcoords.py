import drjit as dr
import mitsuba as mi
import numpy as np
import pytest

BACKENDS = ["volume", "flat"]


def _build_extremum(volume, resolution, extremum_type="gridextremum"):
    """Construct an extremum structure and wire it to volume's majorant
    grid. Construction only takes resolution; the volume and its validity
    bbox are supplied afterwards via update_extremum()."""
    extremum_struct = mi.load_dict({"type": extremum_type, "resolution": resolution})
    extremum_struct.update_extremum(volume.bbox(), volume, None)
    return extremum_struct


def _volume_backend(values):
    return mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(np.reshape(values, (-1, 1, 1))),
            "filter_type": "nearest",
            "accel": False,
            "raw": True,
        }
    )


def _make_volume(boundaries_z, values, backend="volume", mapping="search"):
    d = {
        "type": "mapcoordsvolume",
        "mapping": mapping,
        "backend": backend,
        "boundaries_z": mi.TensorXf(np.asarray(boundaries_z)),
    }
    if backend == "volume":
        d["volume"] = _volume_backend(values)
    elif backend == "flat":
        d["data"] = mi.TensorXf(np.reshape(values, (-1, 1, 1)))
    return mi.load_dict(d)


@pytest.mark.parametrize("backend", BACKENDS)
def test_construct(variant_scalar_mono, backend):
    """Boundary arrays avoid length 3, which the dict parser mistakes for an RGB color."""
    volume =_make_volume([0.0, 0.3, 0.7, 1.0], [1.0, 2.0, 3.0], backend)
    assert volume is not None

    with pytest.raises(RuntimeError, match="at least 2 points"):
        _make_volume([1.0], [1.0], backend)

    with pytest.raises(RuntimeError, match="strictly increasing"):
        _make_volume([0.0, 0.5, 0.5, 1.0], [1.0, 2.0, 3.0], backend)

    with pytest.raises(RuntimeError, match="must start at 0 and end at 1"):
        _make_volume([0.1, 0.4, 0.7, 1.0], [1.0, 2.0, 3.0], backend)
    with pytest.raises(RuntimeError, match="must start at 0 and end at 1"):
        _make_volume([0.0, 0.3, 0.6, 0.9], [1.0, 2.0, 3.0], backend)


def test_construct_unknown_mapping(variant_scalar_mono):
    with pytest.raises(RuntimeError, match='unknown "mapping"'):
        _make_volume([0.0, 1.0], [1.0], mapping="nonsense", backend="volume")


def test_construct_unknown_backend(variant_scalar_mono):
    with pytest.raises(RuntimeError, match='unknown "backend"'):
        _make_volume([0.0, 1.0], [1.0], backend="nonsense")


def test_construct_stale_property(variant_scalar_mono):
    """data/volume belong to one backend each; expand() leaves the other one unreferenced."""
    with pytest.raises(RuntimeError, match="unreferenced property"):
        mi.load_dict(
            {
                "type": "mapcoordsvolume",
                "mapping": "search",
                "backend": "volume",
                "boundaries_z": mi.TensorXf([0.0, 1.0]),
                "volume": _volume_backend([1.0]),
                "data": mi.TensorXf(np.reshape([1.0], (1, 1, 1))),
            }
        )
    with pytest.raises(RuntimeError, match="unreferenced property"):
        mi.load_dict(
            {
                "type": "mapcoordsvolume",
                "mapping": "search",
                "backend": "flat",
                "boundaries_z": mi.TensorXf([0.0, 1.0]),
                "data": mi.TensorXf(np.reshape([1.0], (1, 1, 1))),
                "volume": _volume_backend([1.0]),
            }
        )


def test_construct_flat_multichannel_variant(variant_scalar_rgb):
    with pytest.raises(RuntimeError, match="single spectral channel"):
        _make_volume([0.0, 1.0], [1.0], backend="flat")


@pytest.mark.parametrize("backend", BACKENDS)
def test_traverse_update(variant_scalar_mono, backend):
    """mi.traverse() must expose the backend's data as a live, updatable
    parameter"""
    volume = _make_volume([0.0, 0.3, 0.7, 1.0], [1.0, 2.0, 3.0], backend)

    key = "data" if backend == "flat" else "volume.data"
    params = mi.traverse(volume)
    assert key in params

    it = mi.Interaction3f()
    it.p = [0.5, 0.5, 0.15]  # cell 0 midpoint
    assert dr.allclose(volume.eval_1(it), 1.0)

    params[key] = mi.TensorXf(np.reshape([10.0, 20.0, 30.0], (3, 1, 1)))
    params.update()

    assert dr.allclose(volume.eval_1(it), 10.0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_single_cell(variant_scalar_mono, backend):
    """Degenerate one-cell profile: boundaries_z=[0, 1], a single value."""
    volume = _make_volume([0.0, 1.0], [1.0], backend)

    for z in [0.0, 0.3, 0.5, 0.7, 1.0]:
        it = mi.Interaction3f()
        it.p = [0.5, 0.5, z]
        assert dr.allclose(volume.eval_1(it), 1.0)

    extremum_struct = _build_extremum(volume, mi.ScalarVector3i(1, 1, 2))
    extremum_grid = mi.traverse(extremum_struct)["data"].numpy()
    assert np.allclose(extremum_grid, 1.0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_eval_matches_cell(variant_scalar_mono, backend):
    boundaries_z = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.7, 1.0])
    values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    volume = _make_volume(boundaries_z, values, backend)

    cell_midpoints = 0.5 * (boundaries_z[:-1] + boundaries_z[1:])
    for z, expected in zip(cell_midpoints, values):
        it = mi.Interaction3f()
        it.p = [0.5, 0.5, z]
        assert dr.allclose(volume.eval_1(it), expected)


def test_backends_agree(variant_scalar_mono):
    """Both backends compute the exact same inverse mapping."""
    boundaries_z = np.unique(
        np.concatenate([[0.0], np.random.default_rng(0).random(9), [1.0]])
    )
    values = np.random.default_rng(1).uniform(0.0, 5.0, len(boundaries_z) - 1)
    zs = np.random.default_rng(2).uniform(0.0, 1.0, 50)

    results = {}
    for backend in BACKENDS:
        volume = _make_volume(boundaries_z, values, backend)
        results[backend] = []
        for z in zs:
            it = mi.Interaction3f()
            it.p = [0.5, 0.5, z]
            results[backend].append(volume.eval_1(it))

    assert np.allclose(results["volume"], results["flat"])


@pytest.mark.parametrize("backend", BACKENDS)
def test_extremum_irregular_profile(variant_scalar_mono, backend):
    """Regions must narrow through the mapping before bounding the backend."""
    boundaries_z = [0.0, 0.1, 0.2, 0.3, 0.4, 0.7, 1.0]
    values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    volume = _make_volume(boundaries_z, values, backend)

    extremum_struct = _build_extremum(volume, mi.ScalarVector3i(1, 1, 2))
    extremum_grid = mi.traverse(extremum_struct)["data"].numpy().reshape(2, 2)

    assert np.allclose(extremum_grid, [[1.0, 5.0], [5.0, 6.0]])


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("seed", range(3))
def test_extremum_safety_random(variant_scalar_mono, backend, seed):
    """Every sampled value lies within its region's reported minorant/majorant."""
    rng = np.random.default_rng(seed)

    n_cells = rng.integers(3, 12)
    boundaries_z = np.unique(np.concatenate([[0.0], rng.random(n_cells - 1), [1.0]]))
    if len(boundaries_z) < 3:
        pytest.skip("degenerate boundary draw")

    values = rng.uniform(0.0, 5.0, len(boundaries_z) - 1)
    volume = _make_volume(boundaries_z, values, backend)

    resolution = mi.ScalarVector3i(1, 1, int(rng.integers(2, 6)))
    extremum_struct = _build_extremum(volume, resolution)
    params = mi.traverse(extremum_struct)
    # gridextremum clamps resolution to the underlying volume's
    # resolution, so the actual cell count may be smaller than requested.
    n_extremum_cells = params["resolution"][2]
    extremum_grid = params["data"].numpy().reshape(n_extremum_cells, 2)

    for i in range(n_extremum_cells):
        minorant, majorant = extremum_grid[i]
        samples_z = rng.uniform(i / n_extremum_cells, (i + 1) / n_extremum_cells, 200)
        for z in samples_z:
            it = mi.Interaction3f()
            it.p = [0.5, 0.5, z]
            value = volume.eval_1(it)
            assert minorant - 1e-5 <= value <= majorant + 1e-5, (
                f"seed={seed}, backend={backend}, "
                f"cell={i}, z={z}: value {value} outside [{minorant}, {majorant}]"
            )


def _make_wrap_volume(values_x, wrap_mode=None):
    d = {
        "type": "mapcoordsvolume",
        "mapping": "search",
        "backend": "flat",
        "boundaries_z": mi.TensorXf([0.0, 1.0]),
        "data": mi.TensorXf(np.reshape(values_x, (1, 1, -1))),
    }
    if wrap_mode is not None:
        d["wrap_mode"] = wrap_mode
    return mi.load_dict(d)


def test_wrap_mode_default_is_clamp(variant_scalar_mono):
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
def test_wrap_mode_out_of_range(variant_scalar_mono, wrap_mode, x, expected):
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
        "boundaries_x": mi.TensorXf(np.asarray(boundaries_x)),
        "data": mi.TensorXf(np.reshape(values, (1, 1, -1))),
        "rmin": rmin,
        "rmax": rmax,
        "fillmin": fillmin,
        "fillmax": fillmax,
    }
    if to_world is not None:
        d["to_world"] = to_world
    return mi.load_dict(d)


def test_spherical_irregular_radius_shifts_result(variant_scalar_mono):
    """Irregular boundaries_x with a spherical mapping must remap the radial axis"""
    values = [1.0, 2.0, 3.0]
    trivial = _make_spherical_volume([0.0, 1.0], values)
    irregular = _make_spherical_volume([0.0, 0.1, 0.9, 1.0], values)

    it = mi.Interaction3f()
    it.p = [0.3, 0.0, 0.0]  # r = 0.3, mid-shell either way

    assert not dr.allclose(trivial.eval_1(it), irregular.eval_1(it))


def test_spherical_partial_sphere(variant_scalar_mono):
    """2 x 2 patch covering colatitudes [60, 64] deg and azimuths [0, 40] deg:

                    colatitude 61   colatitude 63
        azimuth 10        1               2
        azimuth 30        3               4
    """
    uv_to_spherical = (
        mi.ScalarTransform4f()
        .translate([0.0, 60.0 / 180.0, 0.5])
        .scale([1.0, 4.0 / 180.0, 40.0 / 360.0])
    )
    volume = mi.load_dict(
        {
            "type": "mapcoordsvolume",
            "mapping": "spherical",
            "backend": "flat",
            "data": mi.TensorXf(np.reshape([[1.0, 2.0], [3.0, 4.0]], (2, 2, 1))),
            "boundaries_x": mi.TensorXf([0.0, 1.0]),
            "uv_to_spherical": uv_to_spherical,
        }
    )

    theta, phi = np.deg2rad(63.0), np.deg2rad(30.0)
    it = mi.Interaction3f()
    it.p = [0.5 * np.sin(theta) * np.cos(phi), 0.5 * np.sin(theta) * np.sin(phi), 0.5 * np.cos(theta)]
    assert dr.allclose(volume.eval_1(it), 4.0)

    # Colatitude cells [0, 90] and [90, 180] deg. The first contains the
    # patch; the second only sees the clamped colatitude-63 row.
    extremum = _build_extremum(
        volume, mi.ScalarVector3i(1, 2, 1), extremum_type="sphericalcoordsextremum"
    )
    min_max = mi.traverse(extremum)["data"].numpy().reshape(-1, 2)
    assert np.allclose(min_max, [[1, 4], [2, 4]], rtol=1e-5)


@pytest.mark.slow
@pytest.mark.parametrize("backend", BACKENDS)
def test_render_beer_lambert(variant_scalar_mono, backend):
    """Beer-Lambert closed-form transmittance check through the full render path."""
    boundaries_z = np.array([0.0, 0.1, 0.25, 0.4, 0.6, 0.8, 1.0])
    sigma = np.array([0.2, 1.2, 0.3, 2.5, 0.6, 1.8])
    expected = np.exp(-np.sum(sigma * np.diff(boundaries_z)))

    sigma_t = _make_volume(boundaries_z, sigma, backend)
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
            "light": {"type": "constant", "radiance": 1.0},
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
        f"backend={backend}: {img.mean()} vs {expected}"
    )
