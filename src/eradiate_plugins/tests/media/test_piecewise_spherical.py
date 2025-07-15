import drjit as dr
import mitsuba as mi
import numpy as np

def generate_spherical_volume_grid(height, num_shells, integral, lower_bound, debug=False):
    scale = integral / lower_bound
    z = np.linspace(0.0, height, num_shells, False)
    ext_coefs = scale * np.exp(-z / lower_bound)
    ext_coefs = np.reshape(ext_coefs, (num_shells, 1, 1))

    return ext_coefs

def create_medium():
    '''
    Create a dictionary for a spherical-shell piecewise medium.
    '''
    # Parameters
    medium_height = 100000
    layers = 2
    angular_samples = 1
    transform = mi.ScalarTransform4f().translate([-medium_height / 2, -medium_height / 2, -medium_height]).scale([medium_height, medium_height, 2 * medium_height])

    grid = generate_spherical_volume_grid(medium_height, layers, 1, 8300)
    exp_volume_grid = mi.VolumeGrid(grid)

    return mi.load_dict(
        {
            "type": "piecewise_spherical",
            "albedo": 0.8,
            "sigma_t": {
                "type": "gridvolume",
                "grid": exp_volume_grid,
                "use_grid_bbox": True,
                "filter_type": "nearest",
                "to_world": transform,
            },
            "angular_samples": angular_samples,
        }
    )

def test01_sample_distances_all_layers(variant_scalar_mono_double):
    medium = create_medium()

    si = dr.zeros(mi.SurfaceInteraction3f)

    origin = mi.Point3f([0.0, 0.0, 100000.0])
    direction = mi.Vector3f([-0.1, 0.0, -1.0])
    ray = mi.Ray3f(origin, direction)
    mei, tr, pdf = medium.sample_interaction_real(ray, si, .01, 0, True)

    '''
    gt_dists = [0.0, 74578.93910366, 82117.51365975, 88515.28423884, 98460.51859237]
    gt_trs = [1.0, 0.75, 0.5, 0.25, 0.01]
    gt_pdfs = [
        7.06034444e-09,
        2.43563215e-05,
        5.41709938e-05,
        2.70854969e-05,
        3.61445783e-06,
    ]

    samples = [0.0, 0.25, 0.5, 0.75, 0.99]
    si = dr.zeros(mi.SurfaceInteraction3f)

    origins = mi.Point3f([0.0, 0.0, 100000.0])
    directions = mi.Vector3f([0.0, 0.0, -1.0])
    ray = mi.Ray3f(origins, directions)

    dists = []
    trs = []
    pdfs = []
    for sample in samples:
        mei, tr, pdf = medium.sample_interaction_real(ray, si, sample, 0, True)
        dists.append(mei.t)
        trs.append(tr)
        pdfs.append(pdf)

    assert dr.allclose(dists, gt_dists)
    assert dr.allclose(trs, gt_trs)
    assert dr.allclose(pdfs, gt_pdfs)
    '''