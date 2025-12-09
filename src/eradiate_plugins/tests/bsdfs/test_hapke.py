from os.path import join

import drjit as dr
import matplotlib.pyplot as plt
import mitsuba as mi
import pytest
import xarray as xr
from mitsuba.test.util import find_resource

PARAMS = {
    "default": {},
    "null": {
        "w": 0.0,
        "b": 0.0,
        "c": 0.0,
        "theta": 0.0,
        "h": 0.0,
        "B_0": 0.0,
    },
    "pommerol": {
        # From Pommerol et al. (2013), identical to what is in the
        # "hapke_*_example.nc" datasets below
        "w": 0.526,
        "theta": 13.3,
        "b": 0.187,
        "c": (1.0 + 0.273) / 2.0,
        "h": 0.083,
        "B_0": 1.0,
    },
}


def angles_to_directions(theta, phi):
    return mi.Vector3f(
        dr.sin(theta) * dr.cos(phi), dr.sin(theta) * dr.sin(phi), dr.cos(theta)
    )


def eval_bsdf(bsdf, wi, wo):
    si = mi.SurfaceInteraction3f()
    si.wi = wi
    ctx = mi.BSDFContext()
    return bsdf.eval(ctx, si, wo, True)[0]


@pytest.fixture
def static_pplane():
    references = find_resource("resources/tests/eradiate_plugins/bsdfs")
    hapke_reference_filename = join(references, "hapke_principal_plane_example.nc")
    return xr.load_dataset(hapke_reference_filename)


@pytest.fixture
def static_hemisphere():
    references = find_resource("resources/tests/eradiate_plugins/bsdfs")
    hapke_reference_filename = join(references, "hapke_hemisphere_example.nc")
    return xr.load_dataset(hapke_reference_filename)


@pytest.fixture
def plot_figures():
    return False


def test_create_hapke(variant_scalar_rgb):
    bsdf = mi.load_dict({"type": "hapke", **PARAMS["null"]})

    assert isinstance(bsdf, mi.BSDF)
    assert bsdf.component_count() == 1
    assert bsdf.flags(0) == mi.BSDFFlags.GlossyReflection | mi.BSDFFlags.FrontSide
    assert bsdf.flags() == bsdf.flags(0)

    params = mi.traverse(bsdf)
    assert "w.value" in params
    assert "b.value" in params
    assert "c.value" in params
    assert "theta.value" in params
    assert "B_0.value" in params
    assert "h.value" in params


def test_eval_hotspot(variant_scalar_rgb):
    hapke = mi.load_dict({"type": "hapke", **PARAMS["pommerol"]})

    theta_i = dr.deg2rad(30.0)
    theta_o = dr.deg2rad(30.0)
    phi_i = dr.deg2rad(0.0)
    phi_o = dr.deg2rad(0.0)

    wi = angles_to_directions(theta_i, phi_i)
    wo = angles_to_directions(theta_o, phi_o)

    values = eval_bsdf(hapke, wi, wo) / dr.abs(dr.cos(theta_o))
    assert dr.allclose(values, 0.24746648)


def test_hapke_grazing_outgoing_direction(variant_scalar_rgb):
    c = 0.273
    c = (1 + c) / 2

    hapke = mi.load_dict({"type": "hapke", **PARAMS["pommerol"]})

    theta_i = dr.deg2rad(30.0)
    theta_o = dr.deg2rad(-89.0)
    phi_i = dr.deg2rad(0.0)
    phi_o = dr.deg2rad(0.0)

    wi = angles_to_directions(theta_i, phi_i)
    wo = angles_to_directions(theta_o, phi_o)

    values = eval_bsdf(hapke, wi, wo) / dr.abs(dr.cos(theta_o))
    assert dr.allclose(values, 0.15426355)


def test_eval_backward(variant_scalar_rgb):
    hapke = mi.load_dict({"type": "hapke", **PARAMS["pommerol"]})

    theta_i = dr.deg2rad(30.0)
    theta_o = dr.deg2rad(80.0)
    phi_i = dr.deg2rad(0.0)
    phi_o = dr.deg2rad(0.0)

    wi = angles_to_directions(theta_i, phi_i)
    wo = angles_to_directions(theta_o, phi_o)

    values = eval_bsdf(hapke, wi, wo) / dr.abs(dr.cos(theta_o))
    assert dr.allclose(values, 0.19555340)


def test_hapke_hemisphere(variant_llvm_ad_rgb, static_hemisphere, plot_figures):
    hapke = mi.load_dict({"type": "hapke", **PARAMS["pommerol"]})

    azimuths = dr.deg2rad(mi.Float(static_hemisphere["phi"].values))
    zeniths = dr.deg2rad(mi.Float(static_hemisphere["vza"].values))

    theta_o, phi_o = dr.meshgrid(zeniths, azimuths)
    theta_i = dr.deg2rad(30.0)
    phi_i = 0.0

    wi = angles_to_directions(theta_i, phi_i)
    wo = angles_to_directions(theta_o, phi_o)
    values = eval_bsdf(hapke, wi, wo) / dr.abs(dr.cos(theta_o))

    npref = static_hemisphere["reflectance"].values
    ref = mi.Float(npref.ravel())

    if False:
        import numpy as np
        from matplotlib import colors

        r, theta = np.meshgrid(dr.sin(zeniths).numpy(), azimuths.numpy())
        npvalues = np.reshape(values, npref.shape)

        nrows = 1
        ncols = 3
        fig, axs = plt.subplots(
            1,
            3,
            subplot_kw=dict(projection="polar"),
            figsize=(4 * ncols, 3 * nrows),
            layout="constrained",
        )

        ax = axs[0]
        contour = ax.contourf(theta, r, npvalues, levels=25, cmap="turbo")
        plt.colorbar(contour)
        ax.set_title("plugin")

        ax = axs[1]
        contour = ax.contourf(theta, r, npref, levels=25, cmap="turbo")
        plt.colorbar(contour)
        ax.set_title("reference")

        ax = axs[2]
        contour = ax.contourf(
            theta,
            r,
            npref - npvalues,
            levels=25,
            cmap="RdBu_r",
            norm=colors.CenteredNorm(),
        )
        plt.colorbar(contour)
        ax.set_title("plugin − reference")

        for ax in axs:
            ax.grid(False)
        fig.suptitle("Hapke BSDF plugin")
        plt.savefig("hapke_hemisphere.png", bbox_inches="tight")
        plt.close()

    assert dr.allclose(ref, values)


def test_hapke_static_principal_plane_reference(
    variant_llvm_ad_rgb, static_pplane, plot_figures
):
    hapke = mi.load_dict({"type": "hapke", **PARAMS["pommerol"]})

    theta_o = mi.Float(dr.deg2rad(static_pplane.svza.values))
    phi_o = dr.zeros(mi.Float, dr.shape(theta_o))
    theta_i = dr.deg2rad(30.0) * dr.ones(mi.Float, dr.shape(theta_o))
    phi_i = dr.zeros(mi.Float, dr.shape(theta_o))

    wi = angles_to_directions(theta_i, phi_i)
    wo = angles_to_directions(theta_o, phi_o)

    values = eval_bsdf(hapke, wi, wo) / dr.abs(dr.cos(theta_o))

    if plot_figures:
        nrows = 2
        ncols = 1
        fig, axs = plt.subplots(nrows, ncols, sharex=True)

        ax = axs[0]
        ax.plot(
            dr.rad2deg(theta_o),
            static_pplane.reflectance,
            color="C0",
            label="Reference (Nguyen, Jacquemoud et al.)",
        )
        ax.plot(
            dr.rad2deg(theta_o),
            values,
            color="C0",
            marker=".",
            ls="",
            label="Eradiate/Mitsuba",
        )
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1))
        ax.set_ylabel("BRDF")

        ax = axs[1]
        ax.plot(
            dr.rad2deg(theta_o),
            (values - static_pplane.reflectance) / static_pplane.reflectance,
            color="C0",
            label="Reference (Nguyen, Jacquemoud et al.)",
        )
        ax.set_xlabel(r"$\omega_\mathrm{o}$ [°]")
        ax.set_xticks(dr.arange(-90, 91, 30))
        ax.set_ylabel("(plugin − ref) / ref")

        for ax in axs:
            ax.grid(True)
        plt.savefig("hapke_pplane.png", bbox_inches="tight")
        plt.close()

    ref = static_pplane.reflectance.values

    assert dr.allclose(ref, values)


@pytest.mark.parametrize("params", ["default", "pommerol"])
def test_reciprocal(variant_llvm_ad_rgb, params):
    """This BSDF is expected to be reciprocal, and this test checks for it."""

    bsdf = mi.load_dict({"type": "hapke", **PARAMS[params]})
    num_samples = 1_000_000
    rng = dr.rng()

    theta_i = rng.uniform(mi.Float, (num_samples,)) * dr.pi / 2.0
    theta_o = rng.uniform(mi.Float, (num_samples,)) * dr.pi / 2.0
    phi_i = rng.uniform(mi.Float, (num_samples,)) * dr.two_pi
    phi_o = rng.uniform(mi.Float, (num_samples,)) * dr.two_pi

    wi = angles_to_directions(theta_i, phi_i)
    wo = angles_to_directions(theta_o, phi_o)

    values_a = eval_bsdf(bsdf, wi, wo) / dr.cos(theta_o)
    values_b = eval_bsdf(bsdf, wo, wi) / dr.cos(theta_i)

    assert dr.allclose(values_a, values_b)


def test_chi2_hapke(variants_vec_backends_once_rgb):
    from mitsuba.chi2 import BSDFAdapter, ChiSquareTest, SphericalDomain

    sample_func, pdf_func = BSDFAdapter(
        "hapke",
        """
        <float name="w" value="0.526"/>
        <float name="theta" value="13.3"/>
        <float name="b" value="0.187"/>
        <float name="c" value="0.273"/>
        <float name="h" value="0.083"/>
        <float name="B_0" value="1"/>
    """,
    )

    chi2 = ChiSquareTest(
        domain=SphericalDomain(),
        sample_func=sample_func,
        pdf_func=pdf_func,
        sample_dim=3,
    )

    assert chi2.run()
