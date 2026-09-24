import numpy as np


def write_binary_grid3d(filename, values):
    with open(filename, "wb") as f:
        f.write(b"VOL")
        f.write(np.uint8(3).tobytes())
        f.write(np.int32(1).tobytes())
        f.write(np.int32(values.shape[2]).tobytes())
        f.write(np.int32(values.shape[1]).tobytes())
        f.write(np.int32(values.shape[0]).tobytes())
        f.write(np.int32(1).tobytes())
        f.write(np.array([0, 0, 0, 1, 1, 1], dtype=np.float32).tobytes())
        f.write(values.ravel().astype(np.float32).tobytes())


# Radial profile on [rmin, rmax]: exponential decay plus a thin dense shell.
# Grids are laid out as zyx = (phi, theta, r).
def radial(r):
    return np.exp(-4 * r) + 3 * np.exp(-(((r - 0.6) / 0.05) ** 2))


r = np.linspace(0, 1, 32)
write_binary_grid3d("textures/sigmat_radial.vol", radial(r).reshape(1, 1, -1))

# Same profile modulated in theta and phi
phi, theta, r = np.meshgrid(
    np.linspace(-np.pi, np.pi, 16),
    np.linspace(0, np.pi, 8),
    np.linspace(0, 1, 32),
    indexing="ij",
)
angular = 1 + 0.8 * np.sin(theta) ** 2 * np.cos(3 * phi)
write_binary_grid3d("textures/sigmat_spherical.vol", radial(r) * angular)
