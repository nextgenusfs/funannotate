import os
import tempfile
import unittest
from unittest import mock

from funannotate import evm_cmds
from funannotate.evm_cmds import build_convert_cmd, build_recombine_cmd, collect_gff3

ALL = lambda tool, flag: True  # noqa: E731  new Rust build
NONE = lambda tool, flag: False  # noqa: E731  Rust build from 1.9.0


class RecombineTests(unittest.TestCase):
    def test_new_rust_gets_threads(self):
        cmd = build_recombine_cmd("rust", "recomb", "/x/partitions_list.out", "evm.out", 4, ALL)
        self.assertEqual(
            cmd, ["recomb", "--partitions", "partitions_list.out", "-O", "evm.out", "--threads", "4"]
        )

    def test_old_rust_gets_no_new_options(self):
        cmd = build_recombine_cmd("rust", "recomb", "/x/p.out", "evm.out", 4, NONE)
        self.assertEqual(cmd, ["recomb", "--partitions", "p.out", "-O", "evm.out"])

    def test_zero_threads_leaves_default(self):
        cmd = build_recombine_cmd("rust", "recomb", "/x/p.out", "evm.out", 0, ALL)
        self.assertNotIn("--threads", cmd)

    def test_perl_unchanged(self):
        cmd = build_recombine_cmd("perl", "r.pl", "/x/p.out", "evm.out", 4, ALL)
        self.assertEqual(
            cmd, ["perl", "r.pl", "--partitions", "p.out", "--output_file_name", "evm.out"]
        )


class ConvertTests(unittest.TestCase):
    def test_new_rust_gets_all_options(self):
        with mock.patch("shutil.which", return_value="/bin/conv"):
            cmd = build_convert_cmd(
                "rust", "conv", "/x/p.out", "evm.out", "g.fa", 4, "evm.concat.gff3", ALL
            )
        self.assertEqual(
            cmd,
            ["conv", "--partitions", "p.out", "-O", "evm.out", "--skip_empty",
             "--threads", "4", "--concat_out", "evm.concat.gff3"],
        )

    def test_old_rust_gets_plain_command(self):
        with mock.patch("shutil.which", return_value="/bin/conv"):
            cmd = build_convert_cmd(
                "rust", "conv", "/x/p.out", "evm.out", "g.fa", 4, "evm.concat.gff3", NONE
            )
        self.assertEqual(cmd, ["conv", "--partitions", "p.out", "-O", "evm.out"])
        self.assertNotIn("--concat_out", cmd)

    def test_each_option_is_probed_separately(self):
        only_skip = lambda tool, flag: flag == "--skip_empty"  # noqa: E731
        with mock.patch("shutil.which", return_value="/bin/conv"):
            cmd = build_convert_cmd(
                "rust", "conv", "/x/p.out", "evm.out", "g.fa", 4, "c.gff3", only_skip
            )
        self.assertIn("--skip_empty", cmd)
        self.assertNotIn("--threads", cmd)
        self.assertNotIn("--concat_out", cmd)

    def test_missing_rust_binary_returns_none(self):
        with mock.patch("shutil.which", return_value=None):
            self.assertIsNone(
                build_convert_cmd("rust", "conv", "/x/p.out", "evm.out", "g.fa", 4, "c", ALL)
            )

    def test_perl_unchanged(self):
        cmd = build_convert_cmd("perl", "c.pl", "/x/p.out", "evm.out", "g.fa", 4, "c", ALL)
        self.assertEqual(cmd[:6], ["perl", "c.pl", "--partitions", "p.out", "--output", "evm.out"])
        self.assertNotIn("--concat_out", cmd)


class ProbeTests(unittest.TestCase):
    def setUp(self):
        evm_cmds.tool_supports.cache_clear()

    def test_flag_found_in_help(self):
        with tempfile.TemporaryDirectory() as d:
            tool = os.path.join(d, "tool")
            with open(tool, "w") as fh:
                fh.write("#!/bin/sh\necho 'Options: --threads <N> --skip_empty'\n")
            os.chmod(tool, 0o755)
            self.assertTrue(evm_cmds.tool_supports(tool, "--threads"))
            self.assertFalse(evm_cmds.tool_supports(tool, "--concat_out"))

    def test_unrunnable_tool_supports_nothing(self):
        self.assertFalse(evm_cmds.tool_supports("/no/such/tool", "--threads"))


class CollectTests(unittest.TestCase):
    def _tree(self, d):
        for name, text in (("c1", "a\n"), ("c2", "b\n")):
            os.makedirs(os.path.join(d, name))
            with open(os.path.join(d, name, "evm.out.gff3"), "w") as fh:
                fh.write(text)

    def test_uses_concat_file_when_present(self):
        with tempfile.TemporaryDirectory() as d:
            self._tree(d)
            concat = os.path.join(d, "evm.concat.gff3")
            with open(concat, "w") as fh:
                fh.write("from-concat\n")
            out = os.path.join(d, "final.gff3")
            self.assertEqual(collect_gff3(d, out, concat), "concat")
            with open(out) as fh:
                self.assertEqual(fh.read(), "from-concat\n")
            self.assertFalse(os.path.exists(concat))

    def test_walks_tree_when_no_concat_file(self):
        with tempfile.TemporaryDirectory() as d:
            self._tree(d)
            out = os.path.join(d, "final.gff3")
            self.assertEqual(collect_gff3(d, out, os.path.join(d, "missing")), "walk")
            with open(out) as fh:
                self.assertEqual(sorted(fh.read().split()), ["a", "b"])

    def test_walk_ignores_other_files(self):
        with tempfile.TemporaryDirectory() as d:
            self._tree(d)
            with open(os.path.join(d, "evm.concat.gff3"), "w") as fh:
                fh.write("stale\n")  # present but not passed: must not be copied
            out = os.path.join(d, "final.gff3")
            collect_gff3(d, out, None)
            with open(out) as fh:
                self.assertNotIn("stale", fh.read())


if __name__ == "__main__":
    unittest.main()
