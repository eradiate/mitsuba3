import mitsuba as mi
import numpy as np

  
def generate_extremum_grid(
        volume_grid, 
        extremum_res,
        filter_type,
        transform=None,
):
    if transform is None:
        transform = mi.ScalarAffineTransform4f()

    volume = mi.load_dict({
        "type":"gridvolume", 
        "grid":volume_grid,
        "filter_type":filter_type,
        "accel":False,
        "to_world":transform,
    })
    extremum_struct = mi.load_dict({
        "type": "extremum_grid",
        "volume": volume,
        "resolution": extremum_res
    })

    extremum_grid =  mi.traverse(extremum_struct)["extremum_grid"].numpy()
    extremum_grid = extremum_grid.reshape(extremum_res.z, extremum_res.y, extremum_res.x, 2)
    extremum_grid = extremum_grid.transpose(2, 1, 0, 3)
    return extremum_struct, extremum_grid
    

def test_build_high_res(variant_scalar_rgb):
    
    n_x = 4
    n_y = 8
    n_z = 3
    n_prod = n_x*n_y*n_z
    data = np.linspace(1, n_prod, n_prod).reshape( n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(n_x, n_y, n_z)    
    _, extremum_grid = generate_extremum_grid(volume_grid, extremum_resolution, "nearest")
    
    assert np.allclose( data, extremum_grid[ :, :, :, 0] )
    assert np.allclose( data, extremum_grid[ :, :, :, 1] )


def test_build_half_res(variant_scalar_rgb):
    
    n_x = 4
    n_y = 8
    n_z = 4
    n_prod = n_x*n_y*n_z
    data = np.linspace(1, n_prod, n_prod).reshape( n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(2, 4, 2)
    _, extremum_grid = generate_extremum_grid(volume_grid, extremum_resolution, "nearest")
    
    assert np.allclose( data[::2,::2,::2], extremum_grid[ :, :, :, 0] )
    assert np.allclose( data[1::2,1::2,1::2], extremum_grid[ :, :, :, 1] )


def test_build_not_multiple(variant_scalar_rgb):
    n_x = 4
    n_y = 9
    n_z = 1
    n_prod = n_x*n_y*n_z
    data = np.linspace(1, n_prod, n_prod).reshape( n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(3, 4, 1)
    _, extremum_grid = generate_extremum_grid(volume_grid, extremum_resolution, "nearest")

    minorant_reference = np.array([1.,3.,5.,7.,10.,12.,14.,16.,19.,21.,23.,25.]).reshape(3,4,1)
    majorant_reference = np.array([12.0,14.0,16.0,18.0,21.0,23.0,25.0,27.0,30.0,32.0,34.0,36.0]).reshape(3,4,1)

    assert np.allclose(minorant_reference, extremum_grid[:,:,:,0])
    assert np.allclose(majorant_reference, extremum_grid[:,:,:,1])

def test_build_trilinear(variant_scalar_rgb):
    n_x = 3
    n_y = 6
    n_z = 1
    n_prod = n_x*n_y*n_z
    data = np.linspace(1, n_prod, n_prod).reshape( n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(3, 6, 1)
    _, extremum_grid = generate_extremum_grid(volume_grid, extremum_resolution, "trilinear")

    minorant_reference = np.array([1.,1.,2.,3.,4.,5.,1.,1.,2.,3.,4.,5.,7.,7.,8.,9.,10.,11.]).reshape(3,6,1)
    majorant_reference = np.array([8.,9.,10.,11.,12.,12.,14.,15.,16.,17.,18.,18.,14.,15.,16.,17.,18.,18.]).reshape(3,6,1)

    assert np.allclose(minorant_reference, extremum_grid[:,:,:,0])
    assert np.allclose(majorant_reference, extremum_grid[:,:,:,1])

def test_build_rotated(variant_scalar_rgb):
    n_x = 3
    n_y = 4
    n_z = 1
    n_prod = n_x*n_y*n_z
    data = np.linspace(1, n_prod, n_prod).reshape( n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(3, 4, 1)
    _, extremum_grid = generate_extremum_grid(
        volume_grid, 
        extremum_resolution, 
        "nearest",
        transform=mi.ScalarAffineTransform4f.rotate([0,0,1], 180)
    )

    assert np.allclose(data.squeeze(), extremum_grid[:,:,:,1].squeeze())

def test_build_scaled(variant_scalar_rgb):
    pass

def assert_compare_segment(ref, other):
    assert np.allclose(ref.tmin, other.tmin)
    assert np.allclose(ref.tmax, other.tmax)
    assert np.allclose(ref.minorant, other.minorant)
    assert np.allclose(ref.majorant, other.majorant)
    assert np.allclose(ref.tau_acc, other.tau_acc)

def test_sample_horizontal_homogeneous(variants_any_scalar, variants_any_llvm):
    n_x = 4
    n_y = 3
    n_z = 1
    mult = 0.5
    data = np.ones((n_x, n_y, n_z))*mult
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(4, 3, 1)
    extremum_struct, _ = generate_extremum_grid(volume_grid, extremum_resolution, "nearest")

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.,0.5,0.5),
        d=mi.ScalarVector3f(1.,0.,0.),
    )
    mint = 0.
    maxt = 4.
    desired_tau = 0.2
    active = True

    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=0.25, tmax=0.5, majorant=0.5, minorant=0.5, tau_acc=0.125,
    )
    # Test ray starting at segment start
    assert_compare_segment(ref_segment, res)

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(-1.,0.5,0.5),
        d=mi.ScalarVector3f(1.,0.,0.),
    )
    ref_segment = mi.ExtremumSegment(
        tmin=1.25, tmax=1.5, majorant=0.5, minorant=0.5, tau_acc=0.125,
    )
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    # Test ray starting outside of the grid
    assert_compare_segment(ref_segment, res)

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.125,0.5,0.5),
        d=mi.ScalarVector3f(1.,0.,0.),
    )
    ref_segment = mi.ExtremumSegment(
        tmin=0.375, tmax=0.625, majorant=0.5, minorant=0.5, tau_acc=0.1875,
    )
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    # Test ray starting outside of the grid
    assert_compare_segment(ref_segment, res)

    desired_tau = 0.8
    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    # Test ray exiting grid
    assert not res.valid()

def test_sample_horizontal_heterogeneous(variants_any_scalar, variants_any_llvm):
    n_x = 8
    n_y = 6
    n_z = 1
    mult = 0.1

    data = np.linspace(1,n_x,n_x).reshape(-1,1,1)*mult
    data = np.ones((n_x, n_y, n_z)) * data
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(4, 3, 1)
    extremum_struct, _ = generate_extremum_grid(volume_grid, extremum_resolution, "nearest")

    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.,0.5,0.5),
        d=mi.ScalarVector3f(1.,0.,0.),
    )
    mint = 0.
    maxt = 10.
    desired_tau = 0.2
    active = True

    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=0.5, tmax=0.75, majorant=0.6, minorant=0.5, tau_acc=0.15,
    )
    # Test ray starting at segment start
    assert_compare_segment(ref_segment, res)

def test_sample_diagonal(variants_any_scalar, variants_any_llvm):
    data_res = mi.ScalarVector3i(2,2,2)
    extremum_res = mi.ScalarVector3i(2, 2, 2)
    mult = 0.5

    data = np.linspace(1,data_res.x,data_res.x).reshape(-1,1,1)*mult
    data = np.ones((data_res.x, data_res.y, data_res.y)) * data
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_struct, _ = generate_extremum_grid(volume_grid, extremum_res, "nearest")

    # perfect diagonal direction
    ray = mi.Ray3f(
        o=mi.ScalarVector3f(0.,0.,0.),
        d=mi.ScalarVector3f(0.5,0.5,0.5),
    )
    mint = 0.
    maxt = 10.
    desired_tau = 0.6
    active = True

    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=1., tmax=2., majorant=1., minorant=1., tau_acc=0.5,
    )
    # Test ray starting at segment start
    assert_compare_segment(ref_segment, res)



# Start adding tests for extremum and medium integration...
# Case with valid sample, case with invalid segment...