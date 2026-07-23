import drjit as dr
import mitsuba as mi
import numpy as np
import pytest


@pytest.fixture
def shapegroup():
    return mi.load_dict({"type": "shapegroup", "cube": {"type": "cube"}})


def transforms(n_inst):
    patch_scale = 10
    patch_offset = 5
    translate = np.random.rand(n_inst, 3)
    translate *= np.array([[patch_offset, patch_offset, 0]])
    translate -= np.array([patch_scale, patch_scale, 0])
    to_worlds = np.expand_dims(np.identity(4), 0)
    to_worlds = np.tile(to_worlds, (n_inst, 1, 1))
    to_worlds[:, 0:3, 3] = translate
    return to_worlds


def translations_along_x(n_inst, spacing):
    """[N,4,4] identity transforms, translated by i * spacing along x."""
    to_worlds = np.tile(np.expand_dims(np.identity(4), 0), (n_inst, 1, 1))
    to_worlds[:, 0, 3] = np.arange(n_inst) * spacing
    return to_worlds


def test01_create(variant_scalar_rgb, shapegroup):
    n_inst = 10
    to_worlds = transforms(n_inst)
    s = mi.load_dict(
        {
            "type": "instancelist",
            "shapegroup": shapegroup,
            "transforms": mi.TensorXf(to_worlds),
        }
    )

    assert s is not None
    # The shape presents itself as a single top-level primitive, but reports
    # the true number of primitives it stands in for.
    assert s.primitive_count() == 1
    assert s.effective_primitive_count() == n_inst * shapegroup.primitive_count()


def test02_ray_intersect_hits_correct_instance(variants_all_backends_once):
    n_inst = 5
    spacing = 10.0
    to_worlds = translations_along_x(n_inst, spacing)

    scene = mi.load_dict(
        {
            "type": "scene",
            "group": {"type": "shapegroup", "sphere": {"type": "sphere"}},
            "instances": {
                "type": "instancelist",
                "shapegroup": {"type": "ref", "id": "group"},
                "transforms": mi.TensorXf(to_worlds),
            },
        }
    )

    for i in range(n_inst):
        x = i * spacing
        ray = mi.Ray3f(o=[x, 0.0, 5.0], d=[0.0, 0.0, -1.0], time=0.0, wavelengths=[])

        assert dr.all(scene.ray_test(ray))

        si = scene.ray_intersect(ray)
        assert dr.all(si.is_valid())
        assert dr.allclose(si.p, mi.Point3f(x, 0.0, 1.0), atol=1e-4)
        assert dr.allclose(si.t, 4.0, atol=1e-4)
        assert dr.allclose(si.n, mi.Normal3f(0.0, 0.0, 1.0), atol=1e-4)


def test03_ray_miss_between_instances(variants_all_backends_once):
    n_inst = 5
    spacing = 10.0
    to_worlds = translations_along_x(n_inst, spacing)

    scene = mi.load_dict(
        {
            "type": "scene",
            "group": {"type": "shapegroup", "sphere": {"type": "sphere"}},
            "instances": {
                "type": "instancelist",
                "shapegroup": {"type": "ref", "id": "group"},
                "transforms": mi.TensorXf(to_worlds),
            },
        }
    )

    # Halfway between instance 0 and instance 1: well clear of both spheres.
    ray = mi.Ray3f(
        o=[0.5 * spacing, 0.0, 5.0], d=[0.0, 0.0, -1.0], time=0.0, wavelengths=[]
    )
    assert not dr.any(scene.ray_test(ray))
    si = scene.ray_intersect(ray)
    assert not dr.any(si.is_valid())


def test04_shadow_ray(variants_all_backends_once):
    n_inst = 3
    spacing = 10.0
    to_worlds = translations_along_x(n_inst, spacing)

    scene = mi.load_dict(
        {
            "type": "scene",
            "group": {"type": "shapegroup", "sphere": {"type": "sphere"}},
            "instances": {
                "type": "instancelist",
                "shapegroup": {"type": "ref", "id": "group"},
                "transforms": mi.TensorXf(to_worlds),
            },
        }
    )

    # Ray that passes straight through instance 1's sphere: occluded.
    occluded_ray = mi.Ray3f(
        o=[spacing, 0.0, -5.0], d=[0.0, 0.0, 1.0], time=0.0, wavelengths=[]
    )
    assert dr.all(scene.ray_test(occluded_ray))

    # Ray that passes between instances: not occluded.
    clear_ray = mi.Ray3f(
        o=[0.5 * spacing, 0.0, -5.0], d=[0.0, 0.0, 1.0], time=0.0, wavelengths=[]
    )
    assert not dr.any(scene.ray_test(clear_ray))


def test05_multi_shape_shapegroup(variants_all_backends_once):
    """Instance a shapegroup with more than one child shape, to exercise the
    (instance_index, child_shape_index) packing/unpacking path."""
    n_inst = 3
    spacing = 20.0
    to_worlds = translations_along_x(n_inst, spacing)

    scene = mi.load_dict(
        {
            "type": "scene",
            "group": {
                "type": "shapegroup",
                # Child 0: unit sphere at the group's local origin.
                "sphere": {"type": "sphere"},
                # Child 1: unit cube translated far away within the group,
                # so the two children never overlap.
                "cube": {
                    "type": "cube",
                    "to_world": mi.ScalarTransform4f().translate([0.0, 20.0, 0.0]),
                },
            },
            "instances": {
                "type": "instancelist",
                "shapegroup": {"type": "ref", "id": "group"},
                "transforms": mi.TensorXf(to_worlds),
            },
        }
    )

    for i in range(n_inst):
        x = i * spacing

        # Hits the sphere (child 0) of instance i.
        ray = mi.Ray3f(o=[x, 0.0, 5.0], d=[0.0, 0.0, -1.0], time=0.0, wavelengths=[])
        si = scene.ray_intersect(ray)
        assert dr.all(si.is_valid())
        assert dr.allclose(si.p, mi.Point3f(x, 0.0, 1.0), atol=1e-4)

        # Hits the cube (child 1) of instance i.
        ray = mi.Ray3f(o=[x, 20.0, 5.0], d=[0.0, 0.0, -1.0], time=0.0, wavelengths=[])
        si = scene.ray_intersect(ray)
        assert dr.all(si.is_valid())
        assert dr.allclose(si.p, mi.Point3f(x, 20.0, 1.0), atol=1e-4)
