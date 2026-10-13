"""Tests for scripts/release_version.py (release labels and tag checks).

Scheme: funannotate/__version__.py holds VERSION and PRERELEASE. The release
label is "X.Y.Z" when PRERELEASE is empty, else "X.Y.Z-<id>.<N>" (rc.4,
beta.12). A release tag is "v" + label: v1.9.0-rc.4 during release
candidates, v1.9.0 for the final release.
"""
import importlib.util
import os
import subprocess
import tempfile
import unittest


HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "release_version.py")
spec = importlib.util.spec_from_file_location("release_version", SCRIPT)
rv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rv)

VERSION_PY = '''import os

VERSION = (1, 9, 0)

# comment
PRERELEASE = "rc.4"

_base = ".".join(map(str, VERSION))
'''

CITATION = """cff-version: 1.2.0
title: funannotate
version: 1.9.0
date-released: '2026-08-10'
"""

CHANGES = """# Changes

## Branch: release/1.9 — work in progress

- fixed a thing

## 1.9.0-rc.3

- older
"""


class LabelTests(unittest.TestCase):
    def test_label(self):
        self.assertEqual(rv.label((1, 9, 0), "rc.4"), "1.9.0-rc.4")
        self.assertEqual(rv.label((1, 9, 0), ""), "1.9.0")

    def test_parse_version_py(self):
        self.assertEqual(rv.parse_version_py(VERSION_PY), ((1, 9, 0), "rc.4"))


class NextVersionTests(unittest.TestCase):
    def check(self, version, pre, bump, expected, preid=None):
        self.assertEqual(rv.next_version(version, pre, bump, preid), expected)

    def test_prerelease_same_id(self):
        self.check((1, 9, 0), "rc.4", "prerelease", ((1, 9, 0), "rc.5"))
        self.check((1, 9, 0), "beta.12", "prerelease", ((1, 9, 0), "beta.13"))
        self.check((1, 9, 0), "rc.4", "prerelease", ((1, 9, 0), "rc.5"), preid="rc")

    def test_prerelease_moves_to_higher_id(self):
        self.check((1, 9, 0), "beta.12", "prerelease", ((1, 9, 0), "rc.1"), preid="rc")

    def test_prerelease_cannot_go_back(self):
        with self.assertRaises(rv.VersionError):
            rv.next_version((1, 9, 0), "rc.4", "prerelease", "beta")

    def test_prerelease_needs_a_prerelease(self):
        with self.assertRaises(rv.VersionError):
            rv.next_version((1, 9, 0), "", "prerelease", None)

    def test_final(self):
        self.check((1, 9, 0), "rc.4", "final", ((1, 9, 0), ""))

    def test_final_twice_fails(self):
        with self.assertRaises(rv.VersionError):
            rv.next_version((1, 9, 0), "", "final", None)

    def test_patch_minor_major_from_final(self):
        self.check((1, 9, 0), "", "patch", ((1, 9, 1), ""))
        self.check((1, 9, 3), "", "minor", ((1, 10, 0), ""))
        self.check((1, 9, 3), "", "major", ((2, 0, 0), ""))

    def test_patch_from_prerelease_fails(self):
        with self.assertRaises(rv.VersionError):
            rv.next_version((1, 9, 0), "rc.4", "patch", None)

    def test_pre_bumps_start_at_1(self):
        self.check((1, 9, 0), "", "prepatch", ((1, 9, 1), "beta.1"))
        self.check((1, 9, 0), "", "preminor", ((1, 10, 0), "rc.1"), preid="rc")
        self.check((1, 9, 0), "rc.4", "premajor", ((2, 0, 0), "alpha.1"), preid="alpha")


class WriteFilesTests(unittest.TestCase):
    def test_bump_writes_all_files(self):
        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "funannotate"))
            for rel, text in (("funannotate/__version__.py", VERSION_PY),
                              ("CITATION.cff", CITATION), ("CHANGES.md", CHANGES)):
                with open(os.path.join(root, rel), "w") as fh:
                    fh.write(text)
            new, notes = rv.bump(root, "final", None, today="2026-11-01")
            self.assertEqual(new, "1.9.0")
            with open(os.path.join(root, "funannotate/__version__.py")) as fh:
                text = fh.read()
            self.assertIn('PRERELEASE = ""\n', text)
            self.assertIn("VERSION = (1, 9, 0)\n", text)
            self.assertIn("# comment\n", text)
            with open(os.path.join(root, "CITATION.cff")) as fh:
                cff = fh.read()
            self.assertIn("version: 1.9.0\n", cff)
            self.assertIn("date-released: '2026-11-01'\n", cff)
            with open(os.path.join(root, "CHANGES.md")) as fh:
                changes = fh.read()
            self.assertIn("## 1.9.0 (2026-11-01)\n\n- fixed a thing", changes)
            self.assertIn("## 1.9.0-rc.3\n", changes)
            self.assertEqual(notes, "- fixed a thing")

    def test_changes_top_already_released_gets_new_section(self):
        out, notes = rv.update_changes(
            "# Changes\n\n## 1.9.0-rc.4 (2026-10-01)\n\n- old\n", "1.9.0-rc.5", "2026-10-02")
        self.assertTrue(out.startswith("# Changes\n\n## 1.9.0-rc.5 (2026-10-02)\n\n"))
        self.assertIn("## 1.9.0-rc.4 (2026-10-01)\n\n- old\n", out)
        self.assertEqual(notes, "Release 1.9.0-rc.5")


class CheckTagTests(unittest.TestCase):
    def test_matching_tags(self):
        self.assertIsNone(rv.tag_problem("v1.9.0-rc.4", (1, 9, 0), "rc.4"))
        self.assertIsNone(rv.tag_problem("v1.9.0", (1, 9, 0), ""))

    def test_rc_tag_without_prerelease_bump(self):
        msg = rv.tag_problem("v1.9.0-rc.5", (1, 9, 0), "rc.4")
        self.assertIn("1.9.0-rc.4", msg)
        self.assertIn("PRERELEASE", msg)

    def test_final_tag_with_prerelease_still_set(self):
        self.assertIsNotNone(rv.tag_problem("v1.9.0", (1, 9, 0), "rc.4"))

    def test_non_version_tags_are_ignored(self):
        self.assertIsNone(rv.tag_problem("db-mirror", (1, 9, 0), "rc.4"))
        self.assertIsNone(rv.tag_problem("busco-odb9-mirror", (1, 9, 0), "rc.4"))

    def test_file_without_prerelease_is_final(self):
        self.assertEqual(rv.parse_version_py("VERSION = (1, 8, 17)\n"), ((1, 8, 17), ""))
        self.assertIsNone(rv.tag_problem("v1.8.17", (1, 8, 17), ""))

    def test_non_version_tag_skips_reading(self):
        # no git repo at all: a non-version tag must not even be looked up
        with tempfile.TemporaryDirectory() as root:
            self.assertIsNone(rv.check_tag(root, "db-mirror"))

    def test_check_tag_reads_file_at_commit(self):
        with tempfile.TemporaryDirectory() as root:
            git = ["git", "-C", root, "-c", "user.name=t", "-c", "user.email=t@t"]
            subprocess.run(["git", "init", "-q", root], check=True)
            os.makedirs(os.path.join(root, "funannotate"))
            path = os.path.join(root, "funannotate", "__version__.py")
            with open(path, "w") as fh:
                fh.write(VERSION_PY)
            subprocess.run(git + ["add", "-A"], check=True)
            subprocess.run(git + ["commit", "-qm", "rc4"], check=True)
            rc4 = subprocess.run(git + ["rev-parse", "HEAD"], check=True,
                                 capture_output=True, text=True).stdout.strip()
            with open(path, "w") as fh:
                fh.write(VERSION_PY.replace('"rc.4"', '"rc.5"'))
            subprocess.run(git + ["commit", "-qam", "rc5"], check=True)
            self.assertIsNone(rv.check_tag(root, "v1.9.0-rc.5", "HEAD"))
            self.assertIsNotNone(rv.check_tag(root, "v1.9.0-rc.5", rc4))


if __name__ == "__main__":
    unittest.main()
