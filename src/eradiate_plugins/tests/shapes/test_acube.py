import drjit as dr
import mitsuba as mi
import numpy as np
import pytest


def make_scene(shape_type, to_world, flip_normals=False):
    return mi.load_dict(
        {
            "type": "scene",
            "box": {
                "type": shape_type,
                "to_world": to_world,
                "flip_normals": flip_normals,
                "bsdf": {"type": "null"},
                "interior": {"type": "homogeneous", "sigma_t": 1.0},
            },
        }
    )


def to_world():
    return (
        mi.ScalarTransform4f()
        .translate([0.3, -1.2, 2.0])
        .rotate(list(np.array([1, 2, 3]) / np.sqrt(14)), 27.0)
        .scale([2.0, 0.5, 1.5])
    )


def scalar_rays(n, seed=0):
    """World-space rays starting inside and outside the transformed box."""
    rng = np.random.default_rng(seed)
    o = rng.uniform(-1.6, 1.6, (n, 3))
    d = rng.normal(size=(n, 3))
    m = np.array(to_world().matrix)
    o = o @ m[:3, :3].T + m[:3, 3]
    d = d @ m[:3, :3].T
    d /= np.linalg.norm(d, axis=1)[:, None]
    return [mi.Ray3f(mi.Point3f(*o[i]), mi.Vector3f(*d[i])) for i in range(n)]


def test01_create(variant_scalar_rgb):
    s = mi.load_dict({"type": "acube"})
    assert s is not None
    assert s.primitive_count() == 1
    assert dr.allclose(s.surface_area(), 24.0)


def test02_bbox(variant_scalar_rgb):
    for r in [1, 2, 4]:
        s = mi.load_dict(
            {"type": "acube", "to_world": mi.Transform4f().scale([r, r, r])}
        )

        b = s.bbox()
        assert b.valid()
        assert dr.allclose(b.center(), [0, 0, 0])
        assert dr.allclose(b.min, [-r, -r, -r])
        assert dr.allclose(b.max, [r, r, r])

    # Under a general transform, bounds and area match the mesh
    acube = mi.load_dict({"type": "acube", "to_world": to_world()})
    cube = mi.load_dict({"type": "cube", "to_world": to_world()})
    assert dr.allclose(acube.bbox().min, cube.bbox().min)
    assert dr.allclose(acube.bbox().max, cube.bbox().max)
    assert dr.allclose(acube.surface_area(), cube.surface_area())


def test03_ray_intersect(variant_scalar_rgb):
    for scale in [
        mi.Vector3f([1.0, 1.0, 1.0]),
        mi.Vector3f([2.0, 1.0, 1.0]),
        mi.Vector3f([1.0, 2.0, 1.0]),
    ]:
        for coordsX in [-1.5, -0.9, -0.5, 0, 0.5, 0.9, 1.5]:
            for coordsY in [-1.5, -0.9, -0.5, 0, 0.5, 0.9, 1.5]:
                s = mi.load_dict(
                    {
                        "type": "scene",
                        "foo": {
                            "type": "acube",
                            "to_world": mi.Transform4f().scale(scale),
                        },
                    }
                )

                ray = mi.Ray3f(
                    o=[coordsX, coordsY, -8],
                    d=[0.0, 0.0, 1.0],
                    time=0.0,
                    wavelengths=[],
                )

                si_found = s.ray_test(ray)
                assert si_found == (
                    (abs(coordsX) <= scale.x) and (abs(coordsY) <= scale.y)
                )

                if si_found:
                    si = s.ray_intersect(ray)
                    ray_u = mi.Ray3f(ray)
                    ray_v = mi.Ray3f(ray)
                    eps = 1e-4
                    ray_u.o += si.dp_du * eps
                    ray_v.o += si.dp_dv * eps
                    si_u = s.ray_intersect(ray_u)
                    si_v = s.ray_intersect(ray_v)
                    if si_u.is_valid():
                        du = (si_u.uv - si.uv) / eps
                        assert dr.allclose(du, [1, 0], atol=2e-2)
                    if si_v.is_valid():
                        dv = (si_v.uv - si.uv) / eps
                        assert dr.allclose(dv, [0, 1], atol=2e-2)

                    # Check normal
                    assert dr.all(si.n == mi.Vector3f([0, 0, -1]))

                    # Check UV: on the -z face, u follows x and v is mirrored y
                    oo = (coordsY - (-scale.y)) / ((scale.y) - (-scale.y))
                    tt = (coordsX - (-scale.x)) / (scale.x - (-scale.x))
                    assert dr.allclose(
                        si.uv, mi.Vector2f([tt, 1.0 - oo]), atol=1e-5, rtol=1e-5
                    )


def test04_ray_intersect_vec(variant_scalar_rgb):
    from mitsuba.scalar_rgb.test.util import check_vectorization

    def kernel(o):
        scene = mi.load_dict(
            {
                "type": "scene",
                "foo": {
                    "type": "acube",
                    "to_world": mi.ScalarTransform4f().scale((0.5, 0.5, 0.5)),
                },
            }
        )

        o = 2.0 * o - 1.0
        o.z = -8.0

        t = scene.ray_intersect(mi.Ray3f(o, [0, 0, 1])).t
        dr.eval(t)
        return t

    # acube does not support CUDA variants
    check_vectorization(kernel, arg_dims=[3], atol=1e-5, modes=["llvm", "llvm_ad"])


def test05_check_normals(variant_scalar_rgb):
    s = mi.load_dict({"type": "scene", "foo": {"type": "acube"}})

    for face in [
        (mi.Vector3f(0, 0, -8), mi.Vector3f(0, 0, 1), mi.Vector3f(0, 0, -1)),
        (mi.Vector3f(0, 0, 8), mi.Vector3f(0, 0, -1), mi.Vector3f(0, 0, 1)),
        (mi.Vector3f(-8, 0, 0), mi.Vector3f(1, 0, 0), mi.Vector3f(-1, 0, 0)),
        (mi.Vector3f(8, 0, 0), mi.Vector3f(-1, 0, 0), mi.Vector3f(1, 0, 0)),
        (mi.Vector3f(0, -8, 0), mi.Vector3f(0, 1, 0), mi.Vector3f(0, -1, 0)),
        (mi.Vector3f(0, 8, 0), mi.Vector3f(0, -1, 0), mi.Vector3f(0, 1, 0)),
    ]:
        ray = mi.Ray3f(o=face[0], d=face[1], time=0.0, wavelengths=[])
        si = s.ray_intersect(ray)
        assert dr.all(si.n == face[2])

    ss = mi.load_dict({"type": "scene", "foo": {"type": "acube", "flip_normals": True}})

    for face in [
        (mi.Vector3f(0, 0, -8), mi.Vector3f(0, 0, 1), mi.Vector3f(0, 0, 1)),
        (mi.Vector3f(0, 0, 8), mi.Vector3f(0, 0, -1), mi.Vector3f(0, 0, -1)),
        (mi.Vector3f(-8, 0, 0), mi.Vector3f(1, 0, 0), mi.Vector3f(1, 0, 0)),
        (mi.Vector3f(8, 0, 0), mi.Vector3f(-1, 0, 0), mi.Vector3f(-1, 0, 0)),
        (mi.Vector3f(0, -8, 0), mi.Vector3f(0, 1, 0), mi.Vector3f(0, 1, 0)),
        (mi.Vector3f(0, 8, 0), mi.Vector3f(0, -1, 0), mi.Vector3f(0, -1, 0)),
    ]:
        ray = mi.Ray3f(o=face[0], d=face[1], time=0.0, wavelengths=[])
        si2 = ss.ray_intersect(ray)
        assert dr.all(si2.n == face[2])


def test06_medium_accesors(variants_any_llvm):
    shape = mi.load_dict(
        {
            "interior": {
                "type": "homogeneous",
                "albedo": {"type": "rgb", "value": [0.99, 0.9, 0.96]},
            },
            "exterior": {
                "type": "homogeneous",
                "albedo": {"type": "rgb", "value": [0.1, 0.1, 0.12]},
            },
            "type": "acube",
        }
    )

    shape_ptr = mi.ShapePtr(shape)

    assert type(shape.interior_medium()) is mi.Medium
    assert type(shape_ptr.interior_medium()) is mi.MediumPtr

    assert type(shape.exterior_medium()) is mi.Medium
    assert type(shape_ptr.exterior_medium()) is mi.MediumPtr


def test07_bsdf_accesors(variants_any_llvm):
    shape = mi.load_dict({"type": "acube", "bsdf": {"type": "diffuse"}})
    shape_ptr = mi.ShapePtr(shape)

    assert type(shape.bsdf()) is mi.BSDF
    assert type(shape_ptr.bsdf()) is mi.BSDFPtr


@pytest.mark.parametrize("flip_normals", [False, True])
def test08_matches_cube_mesh(variant_scalar_rgb, flip_normals):
    s_ref = make_scene("cube", to_world(), flip_normals)
    s_new = make_scene("acube", to_world(), flip_normals)
    assert s_new.shapes()[0].has_flipped_normals() == flip_normals
    n_hit = 0
    for ray in scalar_rays(2000):
        a, b = s_ref.ray_intersect(ray), s_new.ray_intersect(ray)
        assert a.is_valid() == b.is_valid()
        if not a.is_valid():
            continue
        n_hit += 1
        assert dr.allclose(a.t, b.t, rtol=1e-4, atol=1e-5)
        assert dr.allclose(a.p, b.p, rtol=1e-4, atol=1e-4)
        assert dr.allclose(a.n, b.n, atol=1e-4)
        # Interior medium is on the side opposite to the normal
        assert (a.target_medium(ray.d) is None) == (b.target_medium(ray.d) is None)
        # As on meshes, dp_du x dp_dv points outward regardless of flip_normals
        assert dr.allclose(
            dr.dot(dr.cross(b.dp_du, b.dp_dv), b.n) > 0, not flip_normals
        )
        assert dr.allclose(b.sh_frame.n, b.n)
        assert dr.all((b.uv >= -1e-5) & (b.uv <= 1 + 1e-5))
    assert n_hit > 500


def test09_spawned_rays(variant_scalar_rgb):
    """Rays spawned from a hit cross the box exactly once more, or leave it."""
    s_ref, s_new = make_scene("cube", to_world()), make_scene("acube", to_world())
    for ray in scalar_rays(500, seed=1):
        si = s_new.ray_intersect(ray)
        if not si.is_valid():
            continue
        a = s_ref.ray_intersect(si.spawn_ray(ray.d))
        b = s_new.ray_intersect(si.spawn_ray(ray.d))
        assert a.is_valid() == b.is_valid()
        if b.is_valid():
            assert dr.allclose(a.t, b.t, rtol=1e-4, atol=1e-5)


def test10_enclosed_geometry(variant_scalar_rgb):
    """With geometry inside the box, kd-tree nodes inside the cube leave it
    out; rays must still find the cube surface exactly like the mesh."""
    rng = np.random.default_rng(3)
    n = 20000
    c = rng.uniform(-0.9, 0.9, (n, 1, 3))
    v = (c + rng.normal(scale=0.02, size=(n, 3, 3))).reshape(-1, 3)
    m = np.array(to_world().matrix)
    v = v @ m[:3, :3].T + m[:3, 3]

    def scene(shape_type):
        leaves = mi.Mesh("leaves", 3 * n, n)
        params = mi.traverse(leaves)
        params["vertex_positions"] = dr.scalar.ArrayXf(v.ravel().astype(np.float32))
        params["faces"] = dr.scalar.ArrayXu(np.arange(3 * n, dtype=np.uint32))
        params.update()
        return mi.load_dict(
            {
                "type": "scene",
                "leaves": leaves,
                "box": {"type": shape_type, "to_world": to_world()},
            }
        )

    s_ref, s_new = scene("cube"), scene("acube")
    n_box = 0
    for ray in scalar_rays(2000, seed=4):
        a, b = s_ref.ray_intersect(ray), s_new.ray_intersect(ray)
        assert a.is_valid() == b.is_valid()
        if a.is_valid():
            assert dr.allclose(a.t, b.t, rtol=1e-4, atol=1e-5)
            n_box += b.shape.id() == "box"
    assert n_box > 100
