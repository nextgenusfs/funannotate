"""funannotate 1.10 removed compare and outgroups (moved to funannotate-compare)."""
import importlib
import sys
import unittest
from unittest import mock

from funannotate import funannotate as cli


class RemovedCommandTests(unittest.TestCase):
    def run_cli(self, *argv):
        with mock.patch.object(sys, "argv", ["funannotate", *argv]), \
                mock.patch("sys.stderr") as err:
            with self.assertRaises(SystemExit) as cm:
                cli.main()
        out = "".join(c.args[0] for c in err.write.call_args_list)
        return cm.exception.code, out

    def test_compare_and_outgroups_exit_nonzero_and_name_new_package(self):
        for cmd in ("compare", "outgroups"):
            code, out = self.run_cli(cmd, "-i", "x")
            self.assertEqual(code, 2)
            self.assertIn("funannotate-compare", out)
            self.assertIn("removed in", out)

    def test_outgroups_message_names_the_replacement_commands(self):
        _, out = self.run_cli("outgroups")
        self.assertIn("funannotate-compare setup -i outgroups", out)
        self.assertIn("funannotate-compare outgroups", out)

    def test_not_dispatched(self):
        self.assertNotIn("compare", cli.info)
        self.assertNotIn("outgroups", cli.info)

    def test_modules_are_gone(self):
        for mod in ("compare", "outgroups", "stackedBarGraph"):
            with self.assertRaises(ImportError):
                importlib.import_module("funannotate." + mod)

    def test_setup_busco_outgroups_refused(self):
        from funannotate import setupDB
        self.assertFalse(hasattr(setupDB, "outgroupsDB"))

    def test_database_show_outgroups_refused(self):
        from funannotate import database
        with mock.patch.dict("os.environ", {"FUNANNOTATE_DB": "/x"}), \
                mock.patch("builtins.print") as pr:
            with self.assertRaises(SystemExit) as cm:
                database.main(["--show-outgroups"])
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("funannotate-compare", pr.call_args[0][0])


if __name__ == "__main__":
    unittest.main()
