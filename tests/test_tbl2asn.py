import os
import unittest
from unittest import mock

from funannotate.library import build_tbl2asn_cmd, resolve_tbl2asn_binary


class BuildTbl2asnCmdDialectTests(unittest.TestCase):
    """Tests for build_tbl2asn_cmd dialect-aware flag selection."""

    def _cmd(self, dialect, **kwargs):
        defaults = dict(
            folder="/tmp/tbl2asn_dir",
            template="/tmp/template.sbt",
            organism="Aspergillus fumigatus",
            binary="tbl2asn" if dialect == "tbl2asn" else "table2asn",
            dialect=dialect,
        )
        defaults.update(kwargs)
        return build_tbl2asn_cmd(**defaults)

    # --- indir flag ---

    def test_tbl2asn_uses_dash_p(self):
        cmd = self._cmd("tbl2asn")
        self.assertIn("-p", cmd)
        idx = cmd.index("-p")
        self.assertEqual(cmd[idx + 1], "/tmp/tbl2asn_dir")

    def test_table2asn_uses_indir(self):
        cmd = self._cmd("table2asn")
        self.assertIn("-indir", cmd)
        idx = cmd.index("-indir")
        self.assertEqual(cmd[idx + 1], "/tmp/tbl2asn_dir")

    def test_table2asn_does_not_use_dash_p(self):
        cmd = self._cmd("table2asn")
        self.assertNotIn("-p", cmd)

    def test_tbl2asn_does_not_use_indir(self):
        cmd = self._cmd("tbl2asn")
        self.assertNotIn("-indir", cmd)

    # --- legacy-only flags (-c, -a) ---

    def test_tbl2asn_includes_cleanup_and_assembly_flags(self):
        cmd = self._cmd("tbl2asn")
        self.assertIn("-c", cmd)
        self.assertIn("-a", cmd)

    def test_table2asn_omits_cleanup_and_assembly_flags(self):
        cmd = self._cmd("table2asn")
        self.assertNotIn("-c", cmd)
        self.assertNotIn("-a", cmd)

    # --- organism / modifier string ---

    def test_organism_in_j_flag(self):
        cmd = self._cmd("tbl2asn", organism="Neurospora crassa")
        j_val = cmd[cmd.index("-j") + 1]
        self.assertIn("organism=Neurospora crassa", j_val)

    def test_isolate_in_j_flag(self):
        cmd = self._cmd("tbl2asn", isolate="Af293")
        j_val = cmd[cmd.index("-j") + 1]
        self.assertIn("isolate=Af293", j_val)

    def test_strain_in_j_flag(self):
        cmd = self._cmd("tbl2asn", strain="CEA10")
        j_val = cmd[cmd.index("-j") + 1]
        self.assertIn("strain=CEA10", j_val)

    # --- genetic code injection ---

    def test_gcode_injected_when_not_1(self):
        cmd = self._cmd("tbl2asn", gcode=4)
        j_val = cmd[cmd.index("-j") + 1]
        self.assertIn("gcode=4", j_val)

    def test_gcode_omitted_when_1(self):
        cmd = self._cmd("tbl2asn", gcode=1)
        j_val = cmd[cmd.index("-j") + 1]
        self.assertNotIn("gcode", j_val)

    def test_mgcode_injected_when_not_1(self):
        cmd = self._cmd("tbl2asn", mgcode=4)
        j_val = cmd[cmd.index("-j") + 1]
        self.assertIn("mgcode=4", j_val)

    def test_mgcode_omitted_when_1(self):
        cmd = self._cmd("tbl2asn", mgcode=1)
        j_val = cmd[cmd.index("-j") + 1]
        self.assertNotIn("mgcode", j_val)

    # --- discrepancy ---

    def test_discrepancy_flag_present_when_set(self):
        cmd = self._cmd("tbl2asn", discrepancy="/tmp/discrep.txt")
        self.assertIn("-Z", cmd)
        self.assertEqual(cmd[cmd.index("-Z") + 1], "/tmp/discrep.txt")

    def test_discrepancy_flag_absent_when_none(self):
        cmd = self._cmd("tbl2asn", discrepancy=None)
        self.assertNotIn("-Z", cmd)

    # --- extra parameters ---

    def test_extra_parameters_appended(self):
        cmd = self._cmd("tbl2asn", parameters="-l paired-ends")
        self.assertIn("-l", cmd)
        self.assertIn("paired-ends", cmd)

    # --- missing organism raises ---

    def test_missing_organism_raises(self):
        with self.assertRaises((ValueError, SystemExit)):
            build_tbl2asn_cmd(
                folder="/tmp/x",
                template="/tmp/t.sbt",
                organism=None,
                binary="tbl2asn",
                dialect="tbl2asn",
            )


class ResolveTbl2asnBinaryTests(unittest.TestCase):
    """Tests for resolve_tbl2asn_binary precedence logic."""

    def test_env_var_override_table2asn(self):
        with mock.patch.dict(os.environ, {"FUN_TBL2ASN": "/opt/bin/table2asn"}):
            binary, dialect = resolve_tbl2asn_binary()
        self.assertEqual(binary, "/opt/bin/table2asn")
        self.assertEqual(dialect, "table2asn")

    def test_env_var_override_tbl2asn(self):
        with mock.patch.dict(os.environ, {"FUN_TBL2ASN": "/opt/bin/tbl2asn"}):
            binary, dialect = resolve_tbl2asn_binary()
        self.assertEqual(binary, "/opt/bin/tbl2asn")
        self.assertEqual(dialect, "tbl2asn")

    def test_table2asn_preferred_over_tbl2asn(self):
        env = {k: v for k, v in os.environ.items() if k != "FUN_TBL2ASN"}

        def fake_which(name):
            return "/usr/bin/" + name if name in ("table2asn", "tbl2asn") else None

        with mock.patch.dict(os.environ, env, clear=True):
            with mock.patch("funannotate.library.which", side_effect=fake_which):
                binary, dialect = resolve_tbl2asn_binary()
        self.assertEqual(dialect, "table2asn")

    def test_falls_back_to_tbl2asn(self):
        env = {k: v for k, v in os.environ.items() if k != "FUN_TBL2ASN"}

        def fake_which(name):
            return "/usr/bin/tbl2asn" if name == "tbl2asn" else None

        with mock.patch.dict(os.environ, env, clear=True):
            with mock.patch("funannotate.library.which", side_effect=fake_which):
                binary, dialect = resolve_tbl2asn_binary()
        self.assertEqual(dialect, "tbl2asn")

    def test_raises_when_neither_found(self):
        env = {k: v for k, v in os.environ.items() if k != "FUN_TBL2ASN"}
        with mock.patch.dict(os.environ, env, clear=True):
            with mock.patch("funannotate.library.which", return_value=None):
                with self.assertRaises(RuntimeError):
                    resolve_tbl2asn_binary()


if __name__ == "__main__":
    unittest.main()
