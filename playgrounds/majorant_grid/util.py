import os
import mitsuba as mi


def render_cache(fname, overwrite=False, verbose=True):
    """Cache results of long-running rendering functions."""

    def decorator(fn):
        def decorated(*args, **kwargs):
            if (not overwrite) and os.path.exists(fname):
                if verbose:
                    print(f"[↑] {fname}")
                return mi.Bitmap(fname)
            else:
                result = fn(*args, **kwargs)
                mi.Bitmap(result).write(fname)
                if verbose:
                    print(f"[+] {fname}")
                return result

        return decorated

    return decorator
