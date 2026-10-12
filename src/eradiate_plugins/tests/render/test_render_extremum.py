import argparse
import glob
import sys
from os.path import dirname, exists, join

import mitsuba as mi
import pytest

sys.path.insert(0, join(dirname(__file__), "..", ".."))
from render.tests import test_renders

SCENE_DIR = test_renders.find_resource(
    "resources/data_eradiate/tests/scenes/extremum_media"
)
SCENES = glob.glob(join(SCENE_DIR, "*_volpath.xml"))


def reference_scene(scene_fname):
    """The scene a reference is rendered from: ``*_reference.xml`` if present.

    A reference scene describes the same radiance without the feature under
    test (``repeat``'s tiles as separate media, say), so the reference does
    not depend on it. It shares the ``*_ref_*.exr`` names of its scene and is
    only used to generate references, never by the tests.
    """
    reference = scene_fname.replace("_volpath.xml", "_reference.xml")
    return reference if exists(reference) else scene_fname


def list_all_render_test_configs():
    configs = []
    for variant in mi.variants():
        is_jit = "cuda" in variant or "llvm" in variant
        flag_keys = [
            k for k in test_renders.JIT_FLAG_OPTIONS if (k != "scalar") == is_jit
        ]
        for scene_fname in SCENES:
            for k in flag_keys:
                configs.append((variant, scene_fname, k))
    return configs


@pytest.fixture(autouse=True)
def skip_slow_jit_flags(request):
    # Without symbolic loops a render takes minutes to over 30 min, against
    # seconds for jit_flag_option_3 (the Dr.Jit default)
    if request.config.getoption("--all-jit-flags"):
        return
    key = request.node.callspec.params["jit_flags_key"]
    if key not in ("scalar", "jit_flag_option_3"):
        pytest.skip(f"{key} is opt-in: pass --all-jit-flags")


@pytest.mark.slow
@pytest.mark.parametrize(
    "variant, scene_fname, jit_flags_key", list_all_render_test_configs()
)
def test_render_extremum(variant, scene_fname, jit_flags_key):
    """Render with ``eovolpath`` against ``volpath`` references."""
    test_renders.test_render(variant, scene_fname, "eovolpath", jit_flags_key)


@pytest.mark.slow
@pytest.mark.parametrize(
    "variant, scene_fname, jit_flags_key",
    [
        c
        for c in list_all_render_test_configs()
        if "repeat" in c[1] and "phase" not in c[1]
    ],
)
def test_render_repeat_point_queries(variant, scene_fname, jit_flags_key):
    """Render ``repeat`` with ``volpath``, which samples collisions in the gaps
    between tiles from point queries alone. ``volpath`` evaluates phase
    functions at the unfolded point, hence no position-dependent phase."""
    test_renders.test_render(variant, scene_fname, "volpath", jit_flags_key)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--spp", default=32000, type=int)
    parser.add_argument("--scene", default=None, type=str)
    parser.add_argument("--variant", default=None, type=str)
    args = parser.parse_args()
    test_renders.render_ref_images(list(map(reference_scene, SCENES)), **vars(args))
