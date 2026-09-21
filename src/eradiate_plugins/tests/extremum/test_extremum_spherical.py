import drjit as dr
import mitsuba as mi
import numpy as np
import pytest


def generate_extremum_spherical(
    volume_grid,
    extremum_res,
    filter_type,
    rmin=0.5,
    rmax=1.0,
    fillmin=1.0,
    fillmax=0.0,
    transform=None,
):
    if transform is None:
        transform = mi.ScalarAffineTransform4f()

    volume = mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": {
                "type": "gridvolume",
                "grid": volume_grid,
                "filter_type": filter_type,
                "accel": False,
            },
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": fillmin,
            "fillmax": fillmax,
            "to_world": transform,
        }
    )
    extremum_struct = mi.load_dict(
        {"type": "extremum_spherical", "resolution": extremum_res}
    )
    extremum_struct.update_extremum(volume.bbox(), volume)

    extremum_grid = mi.traverse(extremum_struct)["data"].numpy()
    extremum_grid = extremum_grid.reshape(
        extremum_res.z, extremum_res.y, extremum_res.x, 2
    )
    extremum_grid = extremum_grid.transpose(2, 1, 0, 3)
    return extremum_struct, extremum_grid


def test_radial_build_high_res(variant_scalar_rgb):

    n_x = 4
    n_y = 1
    n_z = 1
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(n_x, n_y, n_z)
    _, extremum_grid = generate_extremum_spherical(
        volume_grid, extremum_resolution, "nearest"
    )

    assert np.allclose(data, extremum_grid[:, :, :, 0])
    assert np.allclose(data, extremum_grid[:, :, :, 1])


def test_radial_build_half_res(variant_scalar_rgb):

    n_x = 4
    n_y = 1
    n_z = 1
    n_prod = n_x * n_y * n_z
    data = np.linspace(1, n_prod, n_prod).reshape(n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(2, 1, 1)
    _, extremum_grid = generate_extremum_spherical(
        volume_grid, extremum_resolution, "nearest"
    )

    assert np.allclose(data[::2, ::2, ::2], extremum_grid[:, :, :, 0])
    assert np.allclose(data[1::2, 1::2, 1::2], extremum_grid[:, :, :, 1])


def test_full3d_build_matching_res(variant_scalar_rgb):
    n_r, n_theta, n_phi = 2, 3, 4
    data = np.linspace(1, n_r * n_theta * n_phi, n_r * n_theta * n_phi).reshape(
        n_r, n_theta, n_phi
    )
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(n_r, n_theta, n_phi)
    _, extremum_grid = generate_extremum_spherical(
        volume_grid, extremum_resolution, "nearest"
    )

    assert np.allclose(data, extremum_grid[:, :, :, 0])
    assert np.allclose(data, extremum_grid[:, :, :, 1])


def test_full3d_build_half_res(variant_scalar_rgb):
    data = np.random.default_rng(0).random((4, 2, 4))
    volume_grid = mi.VolumeGrid(data.transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(2, 1, 2)
    _, extremum_grid = generate_extremum_spherical(
        volume_grid, extremum_resolution, "nearest"
    )

    blocks = data.reshape(2, 2, 1, 2, 2, 2)
    assert np.allclose(blocks.min(axis=(1, 3, 5)), extremum_grid[:, :, :, 0])
    assert np.allclose(blocks.max(axis=(1, 3, 5)), extremum_grid[:, :, :, 1])


def _make_spherical_volume(values, rmin, rmax, fillmin=1.0, fillmax=0.0):
    # `values` is indexed [r] or [r, theta, phi]; the grid's axis order is the
    # reverse, radius being the fastest (last) axis.
    data = np.asarray(values, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1, 1)
    data = data.transpose(2, 1, 0)
    grid = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(data),
            "filter_type": "nearest",
            "accel": False,
        }
    )
    return mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": grid,
            "rmin": rmin,
            "rmax": rmax,
            "fillmin": fillmin,
            "fillmax": fillmax,
        }
    )


def _make_extremum(volume, n, n_y=1, n_z=1):
    extremum = mi.load_dict(
        {
            "type": "extremum_spherical",
            "resolution": mi.ScalarVector3i(n, n_y, n_z),
        }
    )
    extremum.update_extremum(volume.bbox(), volume)
    return extremum


def _make_medium(medium_type, volume, extremum):
    return mi.load_dict(
        {
            "type": medium_type,
            "sigma_t": volume,
            "albedo": 0.5,
            "extremum": extremum,
        }
    )


@pytest.mark.parametrize("medium_type", ["heterogeneous", "eoheterogeneous"])
def test_update_on_sigma_t_change(variant_scalar_mono, medium_type):
    n = 4
    rmin, rmax = 0.5, 1.0
    resolution = mi.ScalarVector3i(n, 1, 1)
    before = [1.0, 2.0, 3.0, 4.0]
    after = [5.0, 6.0, 7.0, 8.0]

    volume = _make_spherical_volume(before, rmin, rmax)
    extremum = mi.load_dict({"type": "extremum_spherical", "resolution": resolution})
    medium = _make_medium(medium_type, volume, extremum)

    # Ground truth: an extremum structure built directly from the "after"
    # data, independently of the update mechanism under test.
    ref_volume = _make_spherical_volume(after, rmin, rmax)
    ref_extremum = mi.load_dict(
        {"type": "extremum_spherical", "resolution": resolution}
    )
    ref_extremum.update_extremum(ref_volume.bbox(), ref_volume)
    expected = np.array(mi.traverse(ref_extremum)["data"])

    params = mi.eradiate.traverse(medium)
    params["sigma_t.volume.data"] = mi.TensorXf(np.array(after).reshape(n, 1, 1))
    params.update()

    got = np.array(mi.traverse(medium.extremum_structure())["data"])
    assert np.allclose(got, expected)


def test_sample_multi_layer(variant_scalar_mono):
    # Radial ray entering exactly at rmax, heading to the center. 4 shells of
    # width 0.2, values [1, 2, 3, 4] from rmin=0.2 outward. Sampled inside the
    # 3rd layer (COT per layer: 0.8, 1.4, 1.8) at mid-distance.
    volume = _make_spherical_volume([1.0, 2.0, 3.0, 4.0], rmin=0.2, rmax=1.0)
    extremum = _make_extremum(volume, 4)

    ray = mi.Ray3f(o=[1, 0, 0], d=[-1, 0, 0])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 2.0, target_ot=1.6)

    assert np.allclose(distance, 0.5)
    assert np.allclose(leftover_ot, 0.2)


def test_sample_escapes(variant_scalar_mono):
    # Same setup as `test_sample_mutli_layer`, but the ray escapes the structure.
    volume = _make_spherical_volume([1.0, 2.0, 3.0, 4.0], rmin=0.2, rmax=1.0)
    extremum = _make_extremum(volume, 4)

    ray = mi.Ray3f(o=[1, 0, 0], d=[-1, 0, 0])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 3.0, target_ot=6.0)

    assert np.isinf(distance)
    assert np.allclose(leftover_ot, 1.6)


def test_sample_tangent_shell(variant_scalar_mono):
    # Single shell [0.6, 1.0] (value 2.0). Ray tangent to rmin. Samples past
    # tangent point.
    volume = _make_spherical_volume([2.0], rmin=0.6, rmax=1.0)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[0.6, 0, -3], d=[0, 0, 1])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 10.0, target_ot=2.4)

    assert np.allclose(distance, 3.4, rtol=0.001)
    assert np.allclose(leftover_ot, 0.8)


def test_sample_through_origin(variant_scalar_mono):
    # Ray straight through the origin is tangent to the degenerate r=0
    # boundary. Sampled past the origin.
    volume = _make_spherical_volume([5.0], rmin=0.0, rmax=1.0)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[0, 0, -3], d=[0, 0, 1])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 10.0, target_ot=7.5)

    assert np.allclose(distance, 3.5)
    # `leftover_ot = 2.5` instead of 7.5 means that the origin, in this case the
    # closest point to the ray is treated as a separate segment.
    assert np.allclose(leftover_ot, 2.5)


def test_sample_fillmax_sampled(variant_scalar_mono):
    # fillmax=0.5: the region outside rmax is not vacuum, so an interaction
    # can be sampled there before the ray ever reaches the shell.
    volume = _make_spherical_volume([2.0], rmin=0.5, rmax=1.0, fillmax=0.5)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[3, 0, 0], d=[-1, 0, 0])
    distance, leftover_ot = extremum.sample_test(ray, 0.0, 10.0, target_ot=0.6)

    assert np.allclose(distance, 1.2)
    assert np.allclose(leftover_ot, 0.6)


def test_sample_traverse_outside_rmax(variant_scalar_mono):
    # Test that a ray correctly traverse the extremum outside rmax.
    rmin, rmax, sigma_t, target_ot = 0.8, 1.5, 2.0, 0.6
    volume = _make_spherical_volume(
        [sigma_t], rmin=rmin, rmax=rmax, fillmin=0.0, fillmax=0.0
    )
    extremum = _make_extremum(volume, 1, n_y=2, n_z=4)

    ray = mi.Ray3f(
        o=mi.Vector3f([3.0, 1.0, 0.0]), d=dr.normalize(mi.Vector3f([-1.0, -0.8, 0.0]))
    )

    o_sqr = dr.dot(ray.o, ray.o)
    r_sqr = dr.square(rmax)
    b = dr.dot(ray.o, ray.d)
    t_in = -b - dr.sqrt(b * b - (o_sqr - r_sqr))

    # The ray crosses the azimuth boundary (y=0) before entering the shell.
    assert -ray.o.y / ray.d.y < t_in

    distance, _ = extremum.sample_test(ray, 0.0, 10.0, target_ot=target_ot)

    assert np.allclose(distance, t_in + target_ot / sigma_t)


def test_sample_turning_point_vs_boundary(variant_scalar_mono):
    # The traversal algorithm should favour turning point against boundary
    # events. This test makes them coincide by casting rays from the +X axis.
    # Such rays will have zenith turning points coinciding with the azimuth
    # boundary on the zy plane. Test that the correct event is chosen
    sigma_t, target_ot = 2.0, 0.6

    # 4 azimuth sectors: the ray enters through an empty one and reaches the
    # absorbing one only by crossing the boundary at x = 0.
    grid = mi.load_dict(
        {
            "type": "gridvolume",
            "grid": mi.VolumeGrid(
                np.array([0.0, 0.0, 0.0, sigma_t], dtype=np.float32).reshape(4, 1, 1)
            ),
            "filter_type": "nearest",
            "accel": False,
        }
    )
    volume = mi.load_dict(
        {
            "type": "sphericalcoordsvolume",
            "volume": grid,
            "rmin": 0.0,
            "rmax": 1.0,
            "fillmin": 0.0,
            "fillmax": 0.0,
        }
    )
    extremum = _make_extremum(volume, 1, n_z=4)

    o = mi.Point3f(2.0, 0.0, 0.0)
    for a in np.linspace(0.02, 0.3, 15):
        for b in np.linspace(-0.3, 0.3, 16):  # avoids b = 0, no turning point
            d = dr.normalize(mi.Vector3f(-1.0, a, b))

            t_plane = -o.x / d.x
            expected = t_plane + target_ot / sigma_t
            if dr.norm(o + expected * d) > 0.95:
                continue  # interaction would fall outside the shell

            # The turning point -n0/n1 coincides with the x = 0 crossing
            n0 = d.z * dr.squared_norm(o) - o.z * dr.dot(o, d)
            n1 = dr.dot(o, d) * d.z - dr.squared_norm(d) * o.z
            assert np.allclose(-n0 / n1, t_plane)

            distance, _ = extremum.sample_test(
                mi.Ray3f(o=o, d=d), 0.0, 10.0, target_ot=target_ot
            )

            assert np.allclose(distance, expected, rtol=1e-4), (a, b)


def test_sample_across_zenith_boundary(variant_scalar_mono):
    # Full3D, 2 shells x 2 zenith cells. Only the inner southern cell
    # absorbs, so the interaction is reachable only after the ray has crossed
    # both a shell boundary and the theta = pi/2 cone -- a wrong linear index
    # over (r, theta) picks up a vacuum cell instead.
    sigma_t, target_ot = 2.0, 0.6
    values = np.zeros((2, 2, 1))
    values[0, 1, 0] = sigma_t
    volume = _make_spherical_volume(values, rmin=0.0, rmax=1.0)
    extremum = _make_extremum(volume, 2, n_y=2)

    # Runs parallel to -Z at x = 0.2: crosses r = 0.5 at z = sqrt(0.21), then
    # the z = 0 cone at t = 2, and leaves the absorbing cell at z = -sqrt(0.21).
    ray = mi.Ray3f(o=[0.2, 0, 2], d=[0, 0, -1])
    distance, _ = extremum.sample_test(ray, 0.0, 10.0, target_ot=target_ot)

    assert np.allclose(distance, 2.0 + target_ot / sigma_t)


def test_sample_on_z_axis(variant_scalar_mono):
    # Same split, with the ray on the polar axis where the azimuth is undefined
    # and the zenith index flips discontinuously at the origin rather than
    # by a cone crossing. Only the southern half absorbs.
    sigma_t, target_ot = 2.0, 0.6
    values = np.zeros((1, 2, 4))
    values[0, 1, :] = sigma_t
    volume = _make_spherical_volume(values, rmin=0.0, rmax=1.0)
    extremum = _make_extremum(volume, 1, n_y=2, n_z=4)

    ray = mi.Ray3f(o=[0, 0, 2], d=[0, 0, -1])
    distance, _ = extremum.sample_test(ray, 0.0, 10.0, target_ot=target_ot)

    assert np.allclose(distance, 2.0 + target_ot / sigma_t)


def test_next_segment_miss(variant_scalar_mono):
    # Ray entirely outside the [-1, 1]^3 domain bbox, heading away from it:
    # the whole ray from t is reported as one empty segment.
    volume = _make_spherical_volume([2.0], rmin=0.2, rmax=0.8)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[5, 5, 5], d=[1, 0, 0])
    segment = extremum.next_segment(ray, 0.0)

    assert np.isinf(segment.maxt)
    assert np.allclose(segment.minorant(), 0.0)
    assert np.allclose(segment.majorant(), 0.0)


def test_next_segment_starts_inside_shell(variant_scalar_mono):
    # 4 shells of width 0.2 from rmin=0.2 outward, values [1, 2, 3, 4].
    # Starting at r=0.9 lands in the 4th (outermost) shell, exiting at its
    # inner boundary r=0.8.
    volume = _make_spherical_volume([1.0, 2.0, 3.0, 4.0], rmin=0.2, rmax=1.0)
    extremum = _make_extremum(volume, 4)

    ray = mi.Ray3f(o=[0.9, 0, 0], d=[-1, 0, 0])
    segment = extremum.next_segment(ray, 0.0)

    assert np.allclose(segment.mint, 0.0)
    assert np.allclose(segment.maxt, 0.1)
    assert np.allclose(segment.minorant(), 4.0)
    assert np.allclose(segment.majorant(), 4.0)


def test_next_segment_starts_outside_rmax(variant_scalar_mono):
    # Single shell [0.2, 0.8] (value 2.0), fillmax=0.5. Starting at r=0.95,
    # inside the domain bbox but outside rmax, lands in the fillmax region
    # and exits at the rmax boundary r=0.8.
    volume = _make_spherical_volume([2.0], rmin=0.2, rmax=0.8, fillmax=0.5)
    extremum = _make_extremum(volume, 1)

    ray = mi.Ray3f(o=[0.95, 0, 0], d=[-1, 0, 0])
    segment = extremum.next_segment(ray, 0.0)

    assert np.allclose(segment.mint, 0.0)
    assert np.allclose(segment.maxt, 0.15)
    assert np.allclose(segment.minorant(), 0.5)
    assert np.allclose(segment.majorant(), 0.5)


def test_next_segment_ends_at_zenith_boundary(variant_scalar_mono):
    # Full3D, 2 zenith cells. Starting just above the z = 0 cone and
    # heading down, the nearest boundary of any coordinate is that cone.
    values = np.array([3.0, 7.0]).reshape(1, 2, 1)
    volume = _make_spherical_volume(values, rmin=0.0, rmax=1.0)
    extremum = _make_extremum(volume, 1, n_y=2)

    ray = mi.Ray3f(o=[0.5, 0, 0.1], d=[0, 0, -1])
    segment = extremum.next_segment(ray, 0.0)

    assert np.allclose(segment.mint, 0.0)
    assert np.allclose(segment.maxt, 0.1)
    assert np.allclose(segment.minorant(), 3.0)
    assert np.allclose(segment.majorant(), 3.0)


def test_next_segment_ends_at_azimuth_boundary(variant_scalar_mono):
    # Full3D, 4 azimuth sectors with distinct values. The ray sits in sector 2
    # and reaches the y = 0 half-plane before the outer shell.
    values = np.array([1.0, 2.0, 3.0, 4.0]).reshape(1, 1, 4)
    volume = _make_spherical_volume(values, rmin=0.0, rmax=1.0)
    extremum = _make_extremum(volume, 1, n_z=4)

    ray = mi.Ray3f(o=[0.3, 0.2, 0], d=[0, -1, 0])
    segment = extremum.next_segment(ray, 0.0)

    assert np.allclose(segment.maxt, 0.2)
    assert np.allclose(segment.minorant(), 3.0)
    assert np.allclose(segment.majorant(), 3.0)


def test_next_segment_on_axis_ends_at_pole(variant_scalar_mono):
    # On the polar axis no cone or half-plane is crossed; the segment must
    # still be cut at the origin, where the zenith index jumps.
    values = np.array([3.0, 7.0]).reshape(1, 2, 1)
    volume = _make_spherical_volume(values, rmin=0.0, rmax=1.0)
    extremum = _make_extremum(volume, 1, n_y=2)

    ray = mi.Ray3f(o=[0, 0, 0.5], d=[0, 0, -1])
    segment = extremum.next_segment(ray, 0.0)

    assert np.allclose(segment.maxt, 0.5)
    assert np.allclose(segment.minorant(), 3.0)
    assert np.allclose(segment.majorant(), 3.0)


def test_next_segment_ends_all_dimensions(variant_scalar_mono):
    # 2 shells x 2 zenith cells x 4 azimuth sectors, every cell holding a
    # distinct value. The ray crosses one boundary of each kind in turn, so
    # every segment must be cut by the nearest event across all three
    # dimensions -- a dimension left out of the minimum overshoots into the
    # next cell and reports the wrong extremum.
    values = np.arange(1, 17, dtype=float).reshape(2, 2, 4)
    volume = _make_spherical_volume(values, rmin=0.0, rmax=1.0)
    extremum = _make_extremum(volume, 2, n_y=2, n_z=4)

    # Starts at r = 0.78 with x, y, z > 0, then dives towards -x, -z at
    # constant y: it dips inside the inner shell, crosses the equatorial cone
    # and the x = 0 half-plane, and comes back out. y never changes sign, so
    # the y = 0 half-plane is never reached.
    o = mi.Point3f(0.6, 0.3, 0.4)
    d = dr.normalize(mi.Vector3f(-1.0, 0.0, -1.0))
    ray = mi.Ray3f(o=o, d=d)

    # With rmin = 0 the only interior shell is r = 0.5; the only cone is
    # theta = pi/2, i.e. z = 0; the only half-plane the ray meets is x = 0.
    b = dr.dot(o, d)
    half_chord = dr.sqrt(b * b - (dr.squared_norm(o) - 0.25))
    expected = [
        (-b - half_chord, values[1, 0, 2]),  # into the inner shell
        (-o.z / d.z, values[0, 0, 2]),  # across the equatorial cone
        (-o.x / d.x, values[0, 1, 2]),  # across the x = 0 half-plane
        (-b + half_chord, values[0, 1, 3]),  # back out of the inner shell
    ]

    t = 0.0
    for maxt, value in expected:
        segment = extremum.next_segment(ray, t)
        assert np.allclose(segment.mint, t)
        assert np.allclose(segment.maxt, maxt)
        assert np.allclose(segment.minorant(), value)
        assert np.allclose(segment.majorant(), value)
        t = segment.maxt

    # The last cell runs to rmax, past which fillmax takes over.
    segment = extremum.next_segment(ray, t)
    assert np.allclose(segment.maxt, -b + dr.sqrt(b * b - (dr.squared_norm(o) - 1.0)))
    assert np.allclose(segment.majorant(), values[1, 1, 3])
