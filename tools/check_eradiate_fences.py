#!/usr/bin/env python3
"""Verify Eradiate change fences in upstream Mitsuba files.

All modifications to upstream Mitsuba source files (src/ and include/, excluding
src/eradiate_plugins/ and include/eradiate/) must be enclosed within matching
fence markers:

    ERADIATE_CHANGE_BEGIN: <short description>
    ...
    ERADIATE_CHANGE_END

The marker can appear in any comment style (// #ERADIATE_CHANGE_BEGIN,
# ERADIATE_CHANGE_BEGIN, etc.) — only the keyword itself is matched.

Usage
-----
As a pre-commit hook (no arguments):
    Checks staged modifications to upstream files. Newly added files are
    excluded since they are wholly Eradiate-specific contributions.

Manual check on specific files:
    tools/check_eradiate_fences.py <file> [<file> ...]
    Only fence pairing is verified (no git staging context is available).
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

RE_BEGIN = re.compile(r"ERADIATE_CHANGE_BEGIN")
RE_END = re.compile(r"ERADIATE_CHANGE_END")

# Directories containing upstream Mitsuba code (relative to repo root)
UPSTREAM_PREFIXES = ("src/", "include/")

# Eradiate-specific subtrees that live inside the upstream directories above
ERADIATE_PREFIXES = ("src/eradiate_plugins/", "include/eradiate/")


# ---------------------------------------------------------------------------
# File classification
# ---------------------------------------------------------------------------


def is_upstream(path: str) -> bool:
    return any(path.startswith(p) for p in UPSTREAM_PREFIXES) and not any(
        path.startswith(p) for p in ERADIATE_PREFIXES
    )


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------


def staged_modified_files() -> list[str]:
    """Return files that are staged as modified (not newly added)."""
    result = subprocess.run(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=M"],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.splitlines()


def added_line_ranges(filepath: str) -> list[tuple[int, int]]:
    """Return (start, end) line ranges added in the staged diff for filepath."""
    result = subprocess.run(
        ["git", "diff", "--cached", "-U0", "--", filepath],
        capture_output=True,
        text=True,
        check=True,
    )
    ranges: list[tuple[int, int]] = []
    for line in result.stdout.splitlines():
        m = re.match(r"^@@ -\S+ \+(\d+)(?:,(\d+))? @@", line)
        if m:
            start = int(m.group(1))
            count = int(m.group(2)) if m.group(2) is not None else 1
            if count > 0:
                ranges.append((start, start + count - 1))
    return ranges


# ---------------------------------------------------------------------------
# Fence analysis (operates on already-read file lines)
# ---------------------------------------------------------------------------


def fence_ranges(lines: list[str]) -> list[tuple[int, int]]:
    """Return (start, end) 1-based line ranges enclosed by fence markers."""
    ranges: list[tuple[int, int]] = []
    stack: list[int] = []
    for lineno, line in enumerate(lines, 1):
        if RE_BEGIN.search(line):
            stack.append(lineno)
        elif RE_END.search(line) and stack:
            ranges.append((stack.pop(), lineno))
    return ranges


def check_fence_pairing(filepath: str, lines: list[str]) -> list[str]:
    """Return errors for unmatched BEGIN/END markers."""
    errors: list[str] = []
    depth = 0
    for lineno, line in enumerate(lines, 1):
        if RE_BEGIN.search(line):
            depth += 1
        elif RE_END.search(line):
            if depth == 0:
                errors.append(
                    f"{filepath}:{lineno}: ERADIATE_CHANGE_END without matching BEGIN"
                )
            else:
                depth -= 1
    if depth > 0:
        errors.append(f"{filepath}: {depth} unclosed ERADIATE_CHANGE_BEGIN fence(s)")
    return errors


def check_additions_fenced(
    filepath: str,
    lines: list[str],
    fences: list[tuple[int, int]],
) -> list[str]:
    """Return errors for staged additions that fall outside fence markers."""
    hunks = added_line_ranges(filepath)
    if not hunks:
        return []

    errors: list[str] = []
    for hstart, hend in hunks:
        for lineno in range(hstart, hend + 1):
            line = lines[lineno - 1] if lineno <= len(lines) else ""
            # The fence markers themselves are always allowed
            if RE_BEGIN.search(line) or RE_END.search(line):
                continue
            if not any(s <= lineno <= e for s, e in fences):
                errors.append(
                    f"{filepath}:{lineno}: staged change outside ERADIATE_CHANGE fence"
                )
                break  # one error per hunk is sufficient
    return errors


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "files",
        nargs="*",
        metavar="FILE",
        help="Files to check. If omitted, checks staged modifications via git.",
    )
    args = parser.parse_args()

    errors: list[str] = []

    if args.files:
        # Manual mode: check fence pairing only (no git staging context)
        for filepath in args.files:
            lines = Path(filepath).read_text().splitlines()
            errors.extend(check_fence_pairing(filepath, lines))
    else:
        # Pre-commit hook mode: check staged upstream modifications
        upstream_modified = [f for f in staged_modified_files() if is_upstream(f)]
        if not upstream_modified:
            return 0
        for filepath in upstream_modified:
            lines = Path(filepath).read_text().splitlines()
            fences = fence_ranges(lines)
            errors.extend(check_fence_pairing(filepath, lines))
            errors.extend(check_additions_fenced(filepath, lines, fences))

    if errors:
        print("Eradiate fence violations:")
        for err in errors:
            print(f"  {err}")
        print()
        print("All modifications to upstream Mitsuba files must be enclosed in:")
        print("  // #ERADIATE_CHANGE_BEGIN: <description>")
        print("  // #ERADIATE_CHANGE_END")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
