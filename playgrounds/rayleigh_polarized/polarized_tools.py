import numpy as np
import xarray as xr


def load_reference(path):
    # Load the data
    da = xr.load_dataset(path).stokes

    # Add zenith angle coordinates
    da = da.assign_coords(
        {
            "theta": (
                "mu",
                np.rad2deg(np.arccos(da.mu.values)),
                {"long_name": "Viewing zenith angle", "units": "degree"},
            ),
            "theta0": (
                "mu0",
                np.rad2deg(np.arccos(da.mu0.values)),
                {"long_name": "Viewing zenith angle", "units": "degree"},
            ),
        }
    )

    # Split Stokes parameter dimension into variables
    vars = {
        x: da.sel(parameter=x, drop=True).rename(x)
        for x in da.coords["parameter"].values
    }

    # Swap conventions for Q component
    vars["Q"] *= -1

    # Compute degree of polarization
    vars["DLP"] = 100 * np.sqrt(vars["Q"] ** 2 + vars["U"] ** 2) / vars["I"]

    return xr.Dataset(data_vars=vars)


def sph_to_dir(theta, phi):
    """Map spherical to Euclidean coordinates."""
    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)
    return np.array([cp * st, sp * st, ct])


def directions(mus, phis):
    """
    Generate a sequence of directions from cosine zenith and azimuth angle
    sequences.

    Parameters
    ----------
    mus : array-like
        Cosine of the zenith angle.

    phis : array-like
        Relative azimuth angle in degree.

    Returns
    -------
    directions : ndarray
        An array with shape (3, n) holding the (x, y, z) of directions defined
        by the Cartesian product of (mus, phis).
    """
    # Compute the Cartesian product of the zenith and azimuth sequences
    thetas = np.arccos(mus)
    phis = np.deg2rad(phis)
    angles = np.meshgrid(thetas, phis)

    # Turn it into a sequence of directions
    return sph_to_dir(angles[0].ravel(), angles[1].ravel())
