import subprocess
import os

VERSION = (1, 8, 17)

_base = ".".join(map(str, VERSION))


def _git_version():
    """Return a PEP 440 version string augmented with git state when available.

    Clean tag:        1.8.17
    Ahead of tag:     1.8.17.dev91+g8079d44
    Dirty tree:       1.8.17.dev91+g8079d44.dirty
    No git available: 1.8.17
    """
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        desc = subprocess.check_output(
            ["git", "describe", "--tags", "--dirty", "--always", "--long"],
            cwd=here,
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        # desc format: v1.8.17-91-gabcdef[-dirty]
        parts = desc.lstrip("v").split("-")
        if len(parts) < 3:
            # only a bare hash (no tags) — just append it
            return "{}+{}".format(_base, parts[-1])
        _tag, distance, ghash = parts[0], parts[1], parts[2]
        dirty = len(parts) == 4 and parts[3] == "dirty"
        if int(distance) == 0 and not dirty:
            return _base
        suffix = ".dev{}+{}".format(distance, ghash)
        if dirty:
            suffix += ".dirty"
        return _base + suffix
    except Exception:
        return _base


__version__ = _git_version()
