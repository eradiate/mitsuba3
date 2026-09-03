import drjit as dr
import mitsuba as mi
import numpy as np
import pytest

LENS_A, LENS_B = 0.7195, -0.0637
LENS_COEFFS = f"{LENS_A}, {LENS_B}"
# Calibration sheet values, normalised (pixels at 6000 x 4000)
CENTER_X, CENTER_Y, RADIUS = 3017 / 6000, 1989 / 4000, 1643 / 6000
# A plausible 220-degree lens, fitted out to theta = 110 deg where rho reaches
WIDE_COEFFS = "0.56, -0.0206"


def dict_sensor(projection_model="equisolid", width=32, height=32, **kwargs):
    return {
        "type": "fisheye",
        "projection_model": projection_model,
        "film": {
            "type": "hdrfilm",
            "width": width,
            "height": height,
            "rfilter": {"type": "box"},
        },
        **kwargs,
    }


def poly_kwargs(coefficients):
    """Full sensor kwargs selecting the polynomial model."""
    return {"projection_model": "polynomial", "lens_coefficients": coefficients}


def sample(sensor, position_sample):
    ray, weight = sensor.sample_ray(0.0, 0.5, position_sample, [0.5, 0.5])
    return ray, np.array(weight).max()


def theta_of(ray):
    d = np.array(ray.d).ravel()
    return float(np.arctan2(np.hypot(d[0], d[1]), d[2]))


def poly_radius(coefficients, theta):
    """rho(theta) = c0*theta + c1*theta^2 + ..., of any order."""
    c = [float(x) for x in coefficients.replace(",", " ").split()]
    return sum(v * theta ** (i + 1) for i, v in enumerate(c))


def test_construct(variant_scalar_rgb):
    for model in [
        "equidistant",
        "equisolid",
        "stereographic",
        "orthographic",
        "equisolid_full",
        "polynomial",
    ]:
        kwargs = {"lens_coefficients": LENS_COEFFS} if model == "polynomial" else {}
        sensor = mi.load_dict(dict_sensor(model, **kwargs))
        assert sensor is not None
        assert sensor.bbox().volume() == 0.0  # Degenerate bounding box
        assert "FisheyeCamera" in str(sensor)
        assert not sensor.needs_aperture_sample()  # Infinitely small aperture

    # Instantiation is impossible with an ill-formed specification
    for kwargs in [
        {"projection_model": "bogus"},  # Unknown projection model
        {"fov": 0.0},  # fov outside (0, 360)
        {"fov": 360.0},
        {"projection_model": "orthographic", "fov": 200.0},  # Not invertible
        {"radius": 0.0},  # Degenerate image circle
        {"center_x": 1.5},  # Optical centre off the film
        {"near_clip": 10.0, "far_clip": 1.0},  # Inverted clip planes
        poly_kwargs("abc"),  # Unparseable coefficients
        poly_kwargs("0.0, 0.0"),  # Zero linear coefficient
        poly_kwargs("0.7, nan"),  # Non-finite coefficient
        poly_kwargs(f"{LENS_A}, {LENS_B}, -0.071"),  # rho is not monotonic
        poly_kwargs(f"{LENS_A}, {LENS_B}, -0.070"),  # Monotonic but uninvertible
        poly_kwargs(", ".join(["0.7"] + ["1e-6"] * 20)),  # Above the order cap
    ]:
        with pytest.raises(RuntimeError):
            mi.load_dict(dict_sensor(**kwargs))


def test_sample_ray_direction(variant_scalar_rgb):
    # Sampled halfway to the rim, each analytic model has a closed form
    for model, expected in [
        ("equidistant", np.pi / 4),
        ("equisolid", 2.0 * np.arcsin(0.5 * np.sin(np.pi / 4))),
        ("stereographic", 2.0 * np.arctan(0.5)),
        ("orthographic", np.arcsin(0.5)),
    ]:
        sensor = mi.load_dict(dict_sensor(model, fov=180.0))

        # The optical centre images the optical axis
        ray, weight = sample(sensor, [0.5, 0.5])
        assert weight > 0.0
        assert dr.allclose(ray.d, [0, 0, 1], atol=1e-6)

        # Halfway to the rim along +x azimuth, which is the image left
        ray, weight = sample(sensor, [0.25, 0.5])
        assert weight > 0.0
        assert dr.allclose(theta_of(ray), expected, atol=1e-5)
        assert ray.d.x > 0.0
        assert dr.allclose(ray.d.y, 0.0, atol=1e-6)

        # The rim reaches theta = fov / 2; the corners lie outside the circle
        ray, weight = sample(sensor, [0.0, 0.5])
        assert weight > 0.0
        assert dr.allclose(theta_of(ray), np.pi / 2, atol=1e-5)
        assert sample(sensor, [0.02, 0.02])[1] == 0.0

    # Rays start on the near clip sphere and end on the far clip sphere
    sensor = mi.load_dict(dict_sensor(near_clip=1.0, far_clip=100.0))
    ray, _ = sample(sensor, [0.5, 0.5])
    assert dr.allclose(ray.o, [0, 0, 1], atol=1e-6)
    assert dr.allclose(ray.maxt, 99.0)


def test_sample_ray_differential(variant_scalar_rgb):
    # The differentials must be the primary ray one pixel away on each axis. An
    # odd/even film makes a swapped or shared step observable.
    width, height = 40, 20
    sensor = mi.load_dict(
        dict_sensor("equidistant", width=width, height=height, fov=180.0)
    )
    ray, weight = sensor.sample_ray_differential(0.0, 0.5, [0.55, 0.5], [0.5, 0.5])
    assert np.array(weight).max() > 0.0
    assert ray.has_differentials

    for offset, d_n in (
        ([0.55 + 1.0 / width, 0.5], ray.d_x),
        ([0.55, 0.5 + 1.0 / height], ray.d_y),
    ):
        expected, _ = sensor.sample_ray(0.0, 0.5, offset, [0.5, 0.5])
        assert dr.allclose(d_n, expected.d, atol=1e-6)

    # Every ray leaves the same point, so the differentials share the
    # primary ray's origin rather than being offset like a perspective camera's
    assert dr.allclose(ray.o_x, ray.o) and dr.allclose(ray.o_y, ray.o)

    # Near fov = 360 the invalid region closes in on the rim, so a valid rim
    # sample's neighbour crosses it. Those differentials must not be NaN.
    sensor = mi.load_dict(dict_sensor("equisolid", fov=359.0, width=64, height=64))
    for r in np.linspace(0.80, 1.00, 21):
        ray, weight = sensor.sample_ray_differential(
            0.0, 0.5, [0.5 - 0.5 * r, 0.5], [0.5, 0.5]
        )
        if np.array(weight).max() > 0.0:
            assert np.all(np.isfinite(np.array(ray.d_x))), f"d_x NaN at r = {r}"
            assert np.all(np.isfinite(np.array(ray.d_y))), f"d_y NaN at r = {r}"


def test_non_square_film_keeps_circle(variant_scalar_rgb):
    # The pixel scale is isotropic: on a 40x20 film the image circle (default:
    # the inscribed disk, 10 px radius) is a circle, not an ellipse
    sensor = mi.load_dict(dict_sensor("equidistant", width=40, height=20, fov=180.0))

    # 10 px right of centre and 10 px above centre both land on the rim
    for position_sample in [[0.75, 0.5], [0.5, 0.0]]:
        ray, weight = sample(sensor, position_sample)
        assert weight > 0.0
        assert dr.allclose(theta_of(ray), np.pi / 2, atol=1e-5)

    # up points to the image top, and 16 px right of centre is outside the circle
    assert sample(sensor, [0.5, 0.0])[0].d.y > 0.0
    assert sample(sensor, [0.9, 0.5])[1] == 0.0


def test_polynomial_projection(variant_scalar_rgb):
    # theta(rho) must invert rho(theta), whatever the order and the declared domain.
    radius = 0.25
    for coefficients, fov in [
        (LENS_COEFFS, 180.0),  # The real calibration
        ("0.70, -0.05, 0.010, -0.002, 0.0005", 180.0),  # Order 5: order is free
        (WIDE_COEFFS, 220.0),  # fov declares the calibrated domain
    ]:
        sensor = mi.load_dict(
            dict_sensor(
                "polynomial",
                width=512,
                height=512,
                fov=fov,
                lens_coefficients=coefficients,
                radius=radius,
            )
        )
        theta_max = np.radians(0.5 * fov)
        rho_max = poly_radius(coefficients, theta_max)

        for theta in np.linspace(0.0, theta_max, 40):
            rho = poly_radius(coefficients, theta)
            # Exactly at the edge, r == rho_max is a float32 boundary tie and is
            # the limit of the valid domain by definition
            if rho >= min(rho_max, 1.0) - 1e-6:
                continue
            ray, weight = sample(sensor, [0.5 - rho * radius, 0.5])
            assert weight > 0.0, f"{coefficients}: theta={theta} unexpectedly invalid"
            assert dr.allclose(theta_of(ray), theta, atol=1e-5)

    # The calibration is expressed in normalised film coordinates, so an
    # off-centre optical axis and a non-square film must both be honoured
    sensor = mi.load_dict(
        dict_sensor(
            "polynomial",
            width=600,
            height=400,
            lens_coefficients=LENS_COEFFS,
            center_x=CENTER_X,
            center_y=CENTER_Y,
            radius=RADIUS,
        )
    )
    # The exact order-2 inverse at rho = 0.5, pinning the table to a closed form
    theta_half = 1.0 / (LENS_A + np.sqrt(LENS_A**2 + 2.0 * LENS_B))

    ray, weight = sample(sensor, [CENTER_X, CENTER_Y])
    assert weight > 0.0
    assert dr.allclose(ray.d, [0, 0, 1], atol=1e-6)  # Offset centre is the axis

    ray, _ = sample(sensor, [CENTER_X + 0.5 * RADIUS, CENTER_Y])
    assert dr.allclose(theta_of(ray), theta_half, atol=1e-5)

    # The same pixel distance below the centre gives the same angle
    ray, _ = sample(sensor, [CENTER_X, CENTER_Y + 0.5 * RADIUS / (400 / 600)])
    assert dr.allclose(theta_of(ray), theta_half, atol=1e-5)
    assert ray.d.y < 0.0

    # rho(pi/2) = 0.973 falls inside the image circle, so the ring beyond the
    # horizon and out to the rim images nothing
    rho_horizon = poly_radius(LENS_COEFFS, 0.5 * np.pi)
    assert rho_horizon < 1.0
    assert sample(sensor, [CENTER_X + 1.02 * rho_horizon * RADIUS, CENTER_Y])[1] == 0.0
    assert sample(sensor, [CENTER_X + 0.98 * rho_horizon * RADIUS, CENTER_Y])[1] > 0.0


def test_out_of_circle_directions_finite(variant_scalar_rgb):
    # Samples outside the image circle are discarded, but must still yield a
    # finite, unit-length direction: a NaN direction survives the zero sample
    # weight as 0 * NaN and corrupts neighbouring ray differentials.
    for model, fov, radius in [
        ("equidistant", 180.0, None),
        ("equisolid", 359.0, None),
        ("stereographic", 300.0, None),
        ("orthographic", 180.0, None),
        ("polynomial", 180.0, None),
        # An image circle smaller than the inscribed disk, the normal geometry
        # for a real circular fisheye, pushes the corners past the domain even
        # at fov = 180
        ("equisolid", 180.0, 0.25),
    ]:
        kwargs = {"lens_coefficients": LENS_COEFFS} if model == "polynomial" else {}
        if radius is not None:
            kwargs["radius"] = radius
        sensor = mi.load_dict(dict_sensor(model, fov=fov, **kwargs))

        # The film corner is the farthest point from the centre, so it is the
        # worst case for the domain
        ray, weight = sample(sensor, [0.0, 0.0])
        d = np.array(ray.d).ravel()
        assert np.all(np.isfinite(d)), f"{model}: non-finite direction {d}"
        assert dr.allclose(np.linalg.norm(d), 1.0, atol=1e-5)
        assert weight == 0.0  # Still outside the image circle


def test_render(variants_all_rgb):
    # A constant environment gives the same radiance for every in-disk pixel,
    # for every projection model.
    for model in ["equidistant", "equisolid", "stereographic", "polynomial"]:
        kwargs = {"lens_coefficients": LENS_COEFFS} if model == "polynomial" else {}
        scene = mi.load_dict(
            {
                "type": "scene",
                "integrator": {"type": "path"},
                "sensor": {
                    **dict_sensor(model, width=16, height=16, **kwargs),
                    "sampler": {"type": "independent", "sample_count": 32},
                },
                "emitter": {
                    "type": "constant",
                    "radiance": {"type": "uniform", "value": 0.5},
                },
            }
        )
        img = np.array(mi.render(scene))[:, :, 0]
        # Stay well inside the disk, and inside the polynomial horizon
        assert np.allclose(img[5:11, 5:11], 0.5, atol=1e-3), model
        assert img[0, 0] == 0.0 and img[-1, -1] == 0.0, model
