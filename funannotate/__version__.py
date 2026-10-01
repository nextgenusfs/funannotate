import subprocess
import os
import re

VERSION = (1, 9, 0)

# Release label. Hard-coded on purpose: the package is routinely run from a
# bind-mount or an egg with no .git, where the only other source is a baked
# _version.txt that goes stale. Version 1.9.0-beta.12 shipped reporting plain
# "1.9.0" for exactly that reason, so logs could not tell betas apart. Set to ""
# for a final release.
PRERELEASE = "rc.4"

_base = ".".join(map(str, VERSION))
if PRERELEASE:
    _base = "{}-{}".format(_base, PRERELEASE)

# `git describe --long` output: <tag>-<distance>-g<hash>[-dirty]. The tag itself
# contains hyphens ("v1.9.0-beta.12"), so split from the right with a regex; a
# naive str.split("-") mis-assigns the fields and used to fall through to the
# bare base version.
_DESCRIBE = re.compile(r"^v?(?P<tag>.+)-(?P<dist>\d+)-g(?P<hash>[0-9a-f]+)(?P<dirty>-dirty)?$")


def _parse_describe(desc):
    """Turn `git describe --tags --dirty --always --long` output into a version.

    On the release tag, clean:        1.9.0-beta.13
    N commits past a tag, or dirty:   1.9.0-beta.13.dev3+gf03c1ec[.dirty]
    Not describable (bare hash etc.): 1.9.0-beta.13+<hash>

    The label always comes from VERSION/PRERELEASE, never from the tag, so
    development on a branch cut from an older tag still reports the release it
    is heading towards.
    """
    m = _DESCRIBE.match(desc.strip())
    if not m:
        return "{}+{}".format(_base, desc.strip().split("-")[-1])
    dirty = bool(m.group("dirty"))
    if m.group("tag") == _base and int(m.group("dist")) == 0 and not dirty:
        return _base
    suffix = ".dev{}+g{}".format(m.group("dist"), m.group("hash"))
    if dirty:
        suffix += ".dirty"
    return _base + suffix


def _git_version():
    """Return the version string, augmented with git state when available."""
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        if not os.path.exists(os.path.join(here, "..", ".git")):
            version_txt = os.path.join(here, "_version.txt")
            if os.path.isfile(version_txt):
                with open(version_txt) as _f:
                    baked = _f.read().strip()
                # Only trust the baked string if it belongs to THIS release; a
                # stale file from an older build must not override the label.
                if baked == _base or baked.startswith(_base + ".dev") or baked.startswith(_base + "+"):
                    return baked
            return _base
        desc = subprocess.check_output(
            ["git", "describe", "--tags", "--dirty", "--always", "--long"],
            cwd=here,
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        return _parse_describe(desc)
    except Exception:
        return _base


__version__ = _git_version()
