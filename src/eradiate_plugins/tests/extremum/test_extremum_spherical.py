import drjit as dr
import mitsuba as mi
import numpy as np


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
        {
            "type": "extremum_spherical",
            "volume": volume,
            "resolution": extremum_res,
            "rmin": rmin,
            "rmax": rmax,
        }
    )

    extremum_grid = mi.traverse(extremum_struct)["data"].numpy()
    extremum_grid = extremum_grid.reshape(
        extremum_res.z, extremum_res.y, extremum_res.x, 2
    )
    extremum_grid = extremum_grid.transpose(2, 1, 0, 3)
    return extremum_struct, extremum_grid


def assert_compare_segment(ref, other):
    assert np.allclose(ref.tmin, other.tmin)
    assert np.allclose(ref.tmax, other.tmax)
    assert np.allclose(ref.minorant, other.minorant)
    assert np.allclose(ref.majorant, other.majorant)
    assert np.allclose(ref.tau_acc, other.tau_acc)


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


def test_radial_sample_vertical_heterogeneous_1(variants_any_scalar):
    """
    Test the sampling routine of the RadialOnly variant.
    Direct a ray in the perfect downward vertical direction and test different
    tau values.
    """
    n_x = 4
    n_y = 1
    n_z = 1
    mult = 0.1

    data = np.linspace(1, n_x, n_x).reshape(-1, 1, 1) * mult
    data = np.ones((n_x, n_y, n_z)) * data
    volume_grid = mi.VolumeGrid(data[::-1, :, :].transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(2, 1, 1)
    extremum_struct, _ = generate_extremum_spherical(
        volume_grid,
        extremum_resolution,
        "nearest",
        rmin=0.5,
        rmax=1.0,
        fillmin=1.0,
        fillmax=0.0,
    )

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.0, 0.0, 1.0),
        d=mi.ScalarVector3f(0.0, 0.0, -1.0),
    )
    mint = 0.0
    maxt = 10.0
    active = True

    # Test sampling segment before rmin.
    desired_tau = 0.11
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=0.25,
        tmax=0.5,
        majorant=0.4,
        minorant=0.3,
        tau_acc=0.05,
    )
    assert_compare_segment(ref_segment, res)

    # Test sampling segment in rmin.
    desired_tau = 0.7
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=0.5,
        tmax=1.5,
        majorant=1.0,
        minorant=1.0,
        tau_acc=0.15,
    )
    assert_compare_segment(ref_segment, res)

    # Test sampling segment passed rmin.
    desired_tau = 1.27
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=1.75,
        tmax=2.0,
        majorant=0.2,
        minorant=0.1,
        tau_acc=1.25,
    )
    assert_compare_segment(ref_segment, res)

    # Test sampling passed volume.
    desired_tau = 1.5
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=dr.inf,
        tmax=-dr.inf,
        majorant=0.0,
        minorant=0.0,
        tau_acc=0.0,
    )
    assert_compare_segment(ref_segment, res)


def test_radial_sample_vertical_heterogeneous_2(variants_any_scalar):
    """
    Test the sampling routine of the RadialOnly variant.
    Direct a ray in the perfect downward vertical direction and test different
    ray origin values.
    """
    n_x = 4
    n_y = 1
    n_z = 1
    mult = 0.1

    data = np.linspace(1, n_x, n_x).reshape(-1, 1, 1) * mult
    data = np.ones((n_x, n_y, n_z)) * data
    volume_grid = mi.VolumeGrid(data[::-1, :, :].transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(2, 1, 1)
    extremum_struct, _ = generate_extremum_spherical(
        volume_grid,
        extremum_resolution,
        "nearest",
        rmin=0.5,
        rmax=1.0,
        fillmin=1.0,
        fillmax=0.0,
    )

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.0, 0.0, 1.0),
        d=mi.ScalarVector3f(0.0, 0.0, -1.0),
    )
    mint = 0.0
    maxt = 10.0
    active = True
    desired_tau = 0.11

    # Test ray starting outside of the volume.
    ray.o = mi.ScalarVector3f(0.0, 0.0, 1.2)
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=0.45,
        tmax=0.7,
        majorant=0.4,
        minorant=0.3,
        tau_acc=0.05,
    )
    assert_compare_segment(ref_segment, res)

    # Test ray starting in a layer.
    ray.o = mi.ScalarVector3f(0.0, 0.0, 0.8)
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=0.05,
        tmax=0.3,
        majorant=0.4,
        minorant=0.3,
        tau_acc=0.01,
    )
    assert_compare_segment(ref_segment, res)

    # Test sampling in rmin, passed the midpoint.
    ray.o = mi.ScalarVector3f(0.0, 0.0, -0.25)
    desired_tau = 0.3
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=0.25,
        tmax=0.5,
        majorant=0.4,
        minorant=0.3,
        tau_acc=0.25,
    )
    assert_compare_segment(ref_segment, res)

    # Test sampling passed the volume.
    ray.o = mi.ScalarVector3f(0.0, 0.0, -1.25)
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=dr.inf,
        tmax=-dr.inf,
        majorant=0.0,
        minorant=0.0,
        tau_acc=0.0,
    )
    assert_compare_segment(ref_segment, res)


def test_radial_sample_tangent_heterogeneous(variants_any_scalar):
    """
    Test the sampling routine of the RadialOnly variant.
    Direct a ray in the perfect downward vertical direction tangent to one of
    the radii.
    """
    n_x = 4
    n_y = 1
    n_z = 1
    mult = 0.1

    data = np.linspace(1, n_x, n_x).reshape(-1, 1, 1) * mult
    data = np.ones((n_x, n_y, n_z)) * data
    volume_grid = mi.VolumeGrid(data[::-1, :, :].transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(2, 1, 1)
    extremum_struct, _ = generate_extremum_spherical(
        volume_grid,
        extremum_resolution,
        "nearest",
        rmin=0.5,
        rmax=1.0,
        fillmin=1.0,
        fillmax=0.0,
    )

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.75, 0.0, 1.0),
        d=mi.ScalarVector3f(0.0, 0.0, -1.0),
    )
    mint = 0.0
    maxt = 10.0
    active = True
    desired_tau = 0.2

    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=0.338562,
        tmax=1.66144,
        majorant=0.2,
        minorant=0.1,
        tau_acc=0.0,
    )
    assert_compare_segment(ref_segment, res)


def test_radial_sample_rmin_0_heterogeneous(variants_any_scalar):
    """
    Test the sampling routine of the RadialOnly variant.
    Direct a ray in the perfect downward vertical direction through rmin=0. 
    This should act like a tangent case!
    """
    n_x = 4
    n_y = 1
    n_z = 1
    mult = 0.1

    data = np.linspace(1, n_x, n_x).reshape(-1, 1, 1) * mult
    data = np.ones((n_x, n_y, n_z)) * data
    volume_grid = mi.VolumeGrid(data[::-1, :, :].transpose(2, 1, 0))

    extremum_resolution = mi.ScalarVector3i(2, 1, 1)
    extremum_struct, _ = generate_extremum_spherical(
        volume_grid,
        extremum_resolution,
        "nearest",
        rmin=0.,
        rmax=1.0,
        fillmin=1.0,
        fillmax=0.0,
    )

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.0, 0.0, 1.0),
        d=mi.ScalarVector3f(0.0, 0.0, -1.0),
    )
    mint = 0.0
    maxt = 10.0
    active = True
    desired_tau = 0.3

    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    # The segment tmin and tmax indicate that rmin is dealt as a tangent.
    ref_segment = mi.ExtremumSegment(
        tmin=0.5,
        tmax=1.5,
        majorant=0.4,
        minorant=0.3,
        tau_acc=0.1,
    )
    assert_compare_segment(ref_segment, res)