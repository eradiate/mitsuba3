import argparse
import glob
import sys
from os.path import dirname, join

import mitsuba as mi
import pytest

sys.path.insert(0, join(dirname(__file__), "..", ".."))
from render.tests import test_renders

SCENES = glob.glob(join(dirname(__file__), "scenes", "extremum_media", "*.xml"))


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


@pytest.mark.slow
@pytest.mark.parametrize(
    "variant, scene_fname, jit_flags_key", list_all_render_test_configs()
)
def test_render_extremum(variant, scene_fname, jit_flags_key):
    """Render with ``eovolpath`` against ``volpath`` references."""
    test_renders.test_render(variant, scene_fname, "eovolpath", jit_flags_key)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--spp", default=32000, type=int)
    parser.add_argument("--scene", default=None, type=str)
    parser.add_argument("--variant", default=None, type=str)
    args = parser.parse_args()
    test_renders.render_ref_images(SCENES, **vars(args))
