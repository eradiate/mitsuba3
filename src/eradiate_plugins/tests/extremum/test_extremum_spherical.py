import mitsuba as mi
import numpy as np

def generate_extremum_grid(
        volume_grid, 
        extremum_res,
        filter_type,
        rmin=0.5,
        rmax=1.,
        fillmin=1.,
        fillmax=0.,
        transform=None,
):
    if transform is None:
        transform = mi.ScalarAffineTransform4f()

    volume = mi.load_dict({
        "type": "sphericalcoordsvolume",
        "volume": {
            "type": "gridvolume",
            "grid":volume_grid,
            "filter_type": filter_type,
            "accel": False,
        },
        "rmin":rmin,
        "rmax":rmax,
        "fillmin":fillmin,
        "fillmax":fillmax,
        "to_world": transform,
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


def test_radial_build_high_res(variant_scalar_rgb):
    
    n_x = 4
    n_y = 1
    n_z = 1
    n_prod = n_x*n_y*n_z
    data = np.linspace(1, n_prod, n_prod).reshape( n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(n_x, n_y, n_z)    
    _, extremum_grid = generate_extremum_grid(volume_grid, extremum_resolution, "nearest")
    
    assert np.allclose( data, extremum_grid[ :, :, :, 0] )
    assert np.allclose( data, extremum_grid[ :, :, :, 1] )


def test_radial_build_half_res(variant_scalar_rgb):
    
    n_x = 4
    n_y = 1
    n_z = 1
    n_prod = n_x*n_y*n_z
    data = np.linspace(1, n_prod, n_prod).reshape( n_x, n_y, n_z)
    volume_grid = mi.VolumeGrid(data.transpose(2,1,0))

    extremum_resolution = mi.ScalarVector3i(2, 1, 1)    
    _, extremum_grid = generate_extremum_grid(volume_grid, extremum_resolution, "nearest")
    
    assert np.allclose( data[::2,::2,::2], extremum_grid[ :, :, :, 0] )
    assert np.allclose( data[1::2,1::2,1::2], extremum_grid[ :, :, :, 1] )

