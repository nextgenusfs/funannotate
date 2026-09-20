import subprocess
import os

VERSION = (1, 9, 0)

# Local fix series on top of upstream v1.9.0-beta.12. Bump LOCAL_FIX (f1 -> f2 ...)
# for each further local change; drop it entirely when these land upstream as
# v1.9.0-beta.13. Kept out of VERSION so the numeric tuple stays PEP 440-parseable.
LOCAL_FIX = "beta.12-f1"

_base = ".".join(map(str, VERSION))
if LOCAL_FIX:
    _base = "{}-{}".format(_base, LOCAL_FIX)


def _git_version():
    """Return a PEP 440 version string augmented with git state when available.

    Clean tag:        1.8.17
    Ahead of tag:     1.8.17.dev91+g8079d44
    Dirty tree:       1.8.17.dev91+g8079d44.dirty
    No git available: 1.8.17
    """
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        if not os.path.exists(os.path.join(here, "..", ".git")):
            version_txt = os.path.join(here, "_version.txt")
            if os.path.isfile(version_txt):
                with open(version_txt) as _f:
                    return _f.read().strip()
            return _base
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
