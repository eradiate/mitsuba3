import numpy as np
import mitsuba as mi
import drjit as dr

mi.set_variant("scalar_rgb")

volume = mi.load_dict(
    {
        "type": "gridvolume",
        "accel": False,
        "grid": mi.VolumeGrid(
            np.exp(-10 * np.linspace(0, 1, 10)).reshape((10, 1, 1)),
        ),
    }
)

# Majorant grid generation tests
print(f"{volume = }")
print(f"{volume.max() = }")
print(f"{volume.max_per_channel() = }")
params = mi.traverse(volume)
volume_array = params["data"].numpy()
print(f"{volume_array.shape = }")
print(f"{volume_array.squeeze() = }")
print()

volume_majorants = volume.get_majorant_grid([1, 1, 5])
print(f"{volume_majorants = }")
print(f"{volume_majorants.max() = }")
print(f"{volume_majorants.max_per_channel() = }")
params = mi.traverse(volume_majorants)
volume_majorants_array = params["data"].numpy()
print(f"{volume_majorants_array.shape = }")
print(f"{volume_majorants_array.squeeze() = }")
print(f"{volume_majorants_array.max() = }")
print()

# Majorant grid traversal preparation tests
ray = mi.Ray3f([0.5, 0.5, 0], [0.0, 0.0, 1.0])
traversal_parameters = volume.prepare_majorant_grid_traversal(ray, 0.0, 1e10)
print(f"{traversal_parameters = }")
print()

# Medium majorant supergrid tests
medium = mi.load_dict(
    {
        "type": "heterogeneous",
        "sigma_t": volume,
        "majorant_resolution_factor": [1, 1, 5],
    }
)
print(medium)
params = mi.traverse(medium)
print(params["majorant_grid.data"].numpy())
