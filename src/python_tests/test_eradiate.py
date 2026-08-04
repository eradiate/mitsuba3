import mitsuba as mi
import numpy as np
import pytest


SCENE_DICTS = {
    "referenced_bsdf": {
        "type": "scene",
        "bsdf": {"type": "diffuse", "id": "my_bsdf"},
        # Leading underscores ensures that shapes will be traverse first
        "_rectangle_1": {
            "type": "arectangle",
            "bsdf": {"type": "ref", "id": "my_bsdf"},
        },
        "_rectangle_2": {
            "type": "arectangle",
            "bsdf": {"type": "ref", "id": "my_bsdf"},
        },
        "_disk_1": {
            "type": "disk",
            "bsdf": {"type": "ref", "id": "my_bsdf"},
        },
        "_disk_2": {
            "type": "disk",
            "bsdf": {"type": "diffuse"},
        },
    },
    "aliases": {
        "type": "scene",
        "bsdf": {"type": "diffuse", "id": "diffuse_1"},
        "_rectangle_1": {
            "type": "arectangle",
            "bsdf": {"type": "ref", "id": "diffuse_1"},
        },
        "_rectangle_2": {
            "type": "arectangle",
            "bsdf": {"type": "diffuse", "id": "diffuse_2"},
        },
    },
}


@pytest.mark.parametrize(
    "scene_dict, name_id_override, expected",
    [
        (
            "aliases",
            None,
            {
                "allow_thread_reordering",
                "_rectangle_1.silhouette_sampling_weight",
                "_rectangle_1.to_world",
                "_rectangle_2.bsdf.reflectance.value",
                "_rectangle_2.silhouette_sampling_weight",
                "_rectangle_2.to_world",
                "diffuse_1.reflectance.value",
            },
        ),
        (
            "aliases",
            "diffuse_2",
            {
                "allow_thread_reordering",
                "_rectangle_1.silhouette_sampling_weight",
                "_rectangle_1.to_world",
                "_rectangle_2.silhouette_sampling_weight",
                "_rectangle_2.to_world",
                "diffuse_1.reflectance.value",
                "diffuse_2.reflectance.value",
            },
        ),
    ],
    ids=["referenced_bsdf-no_override", "referenced_bsdf-selected_override"],
)
def test_mi_traverse_name_id_override(
    variants_all_backends_once, scene_dict, name_id_override, expected
):
    mi_scene = mi.load_dict(SCENE_DICTS[scene_dict], optimize=False)

    parameters = mi.eradiate.traverse(mi_scene, name_id_override=name_id_override)
    assert isinstance(parameters, mi.SceneParameters)
    assert set(parameters.keys()) == expected


def test_mi_traverse_multi_parent_update(variants_all_backends_once):
    # "my_bsdf" is referenced by three shapes: updating it must notify all of
    # them, not just the first one encountered during traversal.
    mi_scene = mi.load_dict(SCENE_DICTS["referenced_bsdf"], optimize=False)
    parameters = mi.eradiate.traverse(mi_scene)

    out = parameters.update({"my_bsdf.reflectance.value": 0.3})
    updated = {node.id(): keys for node, keys in out if node.id()}

    for shape_id in ("_rectangle_1", "_rectangle_2", "_disk_1"):
        assert "bsdf" in updated[shape_id]

    # "_disk_2" has its own, unrelated BSDF and must not be notified
    assert "_disk_2" not in updated
