import mitsuba as mi
import numpy as np

  
def generate_extremum_grid(
        volume_grid, 
        extremum_res,
        filter_type,
):
    volume = mi.load_dict({
        "type":"gridvolume", 
        "grid":volume_grid,
        "filter_type":filter_type,
        "accel":False,
    })
    extremum_struct = mi.load_dict({
        "type": "extremum_grid",
        "volume": volume,
        "resolution": extremum_res,
    })

    extremum_grid =  mi.traverse(extremum_struct)["extremum_grid"].numpy()
    extremum_grid = extremum_grid.reshape(extremum_res.x, extremum_res.y, extremum_res.z, 2)
    return extremum_struct, extremum_grid
    

def test_extremum_grid_high_res(variant_scalar_rgb):
    
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


def test_extremum_grid_half_res(variant_scalar_rgb):
    
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


def test_extremum_grid_not_multiple(variant_scalar_rgb):
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

def test_extremum_grid_trilinear(variant_scalar_rgb):
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

def test_extremum_grid_segment_horizontal(variants_any_scalar, variants_any_llvm):
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
    maxt = 2.
    desired_tau = 0.2
    active = True

    res = extremum_struct.sample_segment(ray, mint, maxt, desired_tau, active)
    ref_segment = mi.ExtremumSegment(
        tmin=mi.Float(0.25),
        tmax=mi.Float(0.5),
        sigma_maj=mi.Float(0.5),
        sigma_min=mi.Float(0.5),
        tau_acc=mi.Float(0.125),
    )
    print(res)
    assert np.allclose(ref_segment.tmin, res.tmin)
    assert np.allclose(ref_segment.tmax, res.tmax)
    assert np.allclose(ref_segment.sigma_min, res.sigma_min)
    assert np.allclose(ref_segment.sigma_maj, res.sigma_maj)
    assert np.allclose(ref_segment.tau_acc, res.tau_acc)