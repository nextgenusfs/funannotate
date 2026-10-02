#!/usr/bin/env python3
"""Release labels and tag checks for funannotate.

funannotate/__version__.py is the single source of the release label:

    VERSION = (1, 9, 0)
    PRERELEASE = "rc.4"     # "" for a final release

The label is "1.9.0-rc.4" (or "1.9.0" when PRERELEASE is empty) and the
release tag is "v" + label: v1.9.0-rc.4 for a release candidate, v1.9.0 for
the final release.

Commands:
  current                 print the label from the working tree
  next  BUMP [--preid ID] print the next label without changing files
  bump  BUMP [--preid ID] [--notes-file F]
                          write the next label to __version__.py,
                          CITATION.cff and CHANGES.md; print the label and
                          write the CHANGES.md section to F
  check-tag TAG [--ref R] fail if version tag TAG (v*) does not match the
                          label in __version__.py at commit R (default: TAG);
                          other tags (e.g. db-mirror) are ignored

BUMP is one of:
  prerelease   next number of the current pre-release (rc.4 -> rc.5); with
               --preid of a later stage, start it (beta.12 --preid rc -> rc.1)
  final        drop the pre-release label (1.9.0-rc.4 -> 1.9.0)
  patch, minor, major            next final release from a final release
  prepatch, preminor, premajor   next version as a pre-release
                                 (--preid, default beta; numbering starts at 1)

Used by .github/workflows/publish_release.yml, the tag check in
.github/workflows/container.yml and the .githooks/pre-push hook.
"""

import argparse
import ast
import datetime
import os
import re
import subprocess
import sys

VERSION_FILE = os.path.join("funannotate", "__version__.py")
STAGES = ["alpha", "beta", "rc"]
BUMPS = ["prerelease", "final", "patch", "minor", "major", "prepatch", "preminor", "premajor"]
VERSION_TAG = re.compile(r"^v\d+\.\d+\.\d+")


class VersionError(Exception):
    pass


def label(version, pre):
    base = ".".join(str(x) for x in version)
    return "{}-{}".format(base, pre) if pre else base


def parse_version_py(text):
    """Return (VERSION tuple, PRERELEASE str) from the text of __version__.py."""
    found = {}
    for node in ast.parse(text).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id in ("VERSION", "PRERELEASE"):
                found[target.id] = ast.literal_eval(node.value)
    if "VERSION" not in found:
        raise VersionError("VERSION is not set in {}".format(VERSION_FILE))
    # files from before the PRERELEASE label (1.8.x) are final releases
    return tuple(found["VERSION"]), found.get("PRERELEASE", "")


def split_pre(pre):
    m = re.fullmatch(r"([a-z]+)\.(\d+)", pre)
    if not m:
        raise VersionError("PRERELEASE {!r} is not of the form <id>.<N>, e.g. rc.4".format(pre))
    return m.group(1), int(m.group(2))


def stage_rank(preid):
    if preid not in STAGES:
        raise VersionError("unknown pre-release id {!r}; use one of {}".format(preid, ", ".join(STAGES)))
    return STAGES.index(preid)


def next_version(version, pre, bump, preid=None):
    """Return (VERSION, PRERELEASE) after BUMP."""
    major, minor, patch = version
    if bump == "prerelease":
        if not pre:
            raise VersionError("{} is a final release; use prepatch, preminor or premajor".format(label(version, pre)))
        cur_id, num = split_pre(pre)
        new_id = preid or cur_id
        if new_id == cur_id:
            return version, "{}.{}".format(cur_id, num + 1)
        if stage_rank(new_id) < stage_rank(cur_id):
            raise VersionError("cannot go from {} back to {}".format(cur_id, new_id))
        return version, "{}.1".format(new_id)
    if bump == "final":
        if not pre:
            raise VersionError("{} is already a final release".format(label(version, pre)))
        return version, ""
    if bump in ("patch", "minor", "major"):
        if pre:
            raise VersionError(
                "{} is a pre-release; use final to release {}".format(label(version, pre), label(version, "")))
    if bump in ("patch", "prepatch"):
        new = (major, minor, patch + 1)
    elif bump in ("minor", "preminor"):
        new = (major, minor + 1, 0)
    elif bump in ("major", "premajor"):
        new = (major + 1, 0, 0)
    else:
        raise VersionError("unknown bump {!r}; use one of {}".format(bump, ", ".join(BUMPS)))
    if bump.startswith("pre"):
        new_id = preid or "beta"
        stage_rank(new_id)
        return new, "{}.1".format(new_id)
    return new, ""


def set_version_py(text, version, pre):
    out, n1 = re.subn(r"(?m)^VERSION = \([^)]*\)$", "VERSION = ({}, {}, {})".format(*version), text)
    out, n2 = re.subn(r'(?m)^PRERELEASE = (["\']).*\1$', 'PRERELEASE = "{}"'.format(pre), out)
    if n1 != 1 or n2 != 1:
        raise VersionError("could not find one VERSION and one PRERELEASE line in {}".format(VERSION_FILE))
    return out


def update_citation(text, new_label, today):
    text = re.sub(r"(?m)^version: .*$", "version: {}".format(new_label), text)
    return re.sub(r"(?m)^date-released: .*$", "date-released: '{}'".format(today), text)


def update_changes(text, new_label, today):
    """Retitle the top "## " section as the release, or add one if the top
    section is already a released version. Return (text, release notes)."""
    heading = "## {} ({})".format(new_label, today)
    m = re.search(r"(?m)^## (.*)$", text)
    if m and not re.match(r"\d+\.\d+\.\d+", m.group(1)):
        text = text[:m.start()] + heading + text[m.end():]
        body = re.search(r"(?ms)^{}\n(.*?)(?=^## |\Z)".format(re.escape(heading)), text)
        notes = body.group(1).strip() if body else ""
        return text, notes or "Release {}".format(new_label)
    insert = heading + "\n\nRelease {}\n\n".format(new_label)
    if m:
        text = text[:m.start()] + insert + text[m.start():]
    else:
        text = text.rstrip("\n") + "\n\n" + insert
    return text, "Release {}".format(new_label)


def read(root, rel):
    with open(os.path.join(root, rel)) as fh:
        return fh.read()


def write(root, rel, text):
    with open(os.path.join(root, rel), "w") as fh:
        fh.write(text)


def bump(root, bump_type, preid, today=None):
    """Apply BUMP to the files under ROOT. Return (new label, release notes)."""
    today = today or datetime.date.today().isoformat()
    version, pre = parse_version_py(read(root, VERSION_FILE))
    new_version, new_pre = next_version(version, pre, bump_type, preid)
    new_label = label(new_version, new_pre)
    write(root, VERSION_FILE, set_version_py(read(root, VERSION_FILE), new_version, new_pre))
    notes = "Release {}".format(new_label)
    if os.path.isfile(os.path.join(root, "CITATION.cff")):
        write(root, "CITATION.cff", update_citation(read(root, "CITATION.cff"), new_label, today))
    if os.path.isfile(os.path.join(root, "CHANGES.md")):
        text, notes = update_changes(read(root, "CHANGES.md"), new_label, today)
        write(root, "CHANGES.md", text)
    return new_label, notes


def tag_problem(tag, version, pre):
    """Return an error message if version tag TAG does not match, else None."""
    if not VERSION_TAG.match(tag):
        return None
    want = "v" + label(version, pre)
    if tag == want:
        return None
    return (
        "tag {} does not match {}: VERSION = {} and PRERELEASE = {!r} give {}. "
        "Bump PRERELEASE (or VERSION) in {} and commit it before tagging, "
        "or tag {} instead.".format(tag, VERSION_FILE, version, pre, want, VERSION_FILE, want)
    )


def check_tag(root, tag, ref=None):
    """Check TAG against __version__.py at commit REF (default: the tag)."""
    if not VERSION_TAG.match(tag):
        return None
    ref = ref or tag
    proc = subprocess.run(
        ["git", "-C", root, "show", "{}:{}".format(ref, VERSION_FILE.replace(os.sep, "/"))],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if proc.returncode != 0:
        return "cannot read {} at {}: {}".format(VERSION_FILE, ref, proc.stderr.strip())
    version, pre = parse_version_py(proc.stdout)
    return tag_problem(tag, version, pre)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=".", help="repository root (default: .)")
    sub = p.add_subparsers(dest="cmd")
    sub.required = True
    sub.add_parser("current")
    for name in ("next", "bump"):
        s = sub.add_parser(name)
        s.add_argument("bump_type", choices=BUMPS)
        s.add_argument("--preid", default=None, choices=STAGES)
        if name == "bump":
            s.add_argument("--notes-file", default=None, help="write the release notes here")
    s = sub.add_parser("check-tag")
    s.add_argument("tag")
    s.add_argument("--ref", default=None, help="commit to check (default: the tag)")
    args = p.parse_args(argv)
    try:
        if args.cmd == "current":
            print(label(*parse_version_py(read(args.root, VERSION_FILE))))
        elif args.cmd == "next":
            version, pre = parse_version_py(read(args.root, VERSION_FILE))
            print(label(*next_version(version, pre, args.bump_type, args.preid)))
        elif args.cmd == "bump":
            new_label, notes = bump(args.root, args.bump_type, args.preid)
            if args.notes_file:
                with open(args.notes_file, "w") as fh:
                    fh.write(notes + "\n")
            print(new_label)
        elif args.cmd == "check-tag":
            problem = check_tag(args.root, args.tag, args.ref)
            if problem:
                print("ERROR: " + problem, file=sys.stderr)
                return 1
    except VersionError as e:
        print("ERROR: {}".format(e), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
