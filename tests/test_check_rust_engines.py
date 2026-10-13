import os
import unittest
from unittest import mock

from funannotate import check


class CheckTrinityRustTests(unittest.TestCase):
    def test_all_tools_present(self):
        # Every Rust utility resolves to a path -> all values truthy.
        with mock.patch.object(
            check.lib, "which_path", side_effect=lambda t: "/env/bin/" + t
        ):
            status = check.check_trinity_rust()
        self.assertEqual(set(status), set(check.TRINITY_RUST_TOOLS))
        self.assertTrue(all(status.values()))
        self.assertEqual(
            status["sam_to_read_coords"], "/env/bin/sam_to_read_coords"
        )

    def test_all_tools_absent(self):
        # which_path returns None for everything -> Rust acceleration inactive.
        with mock.patch.object(check.lib, "which_path", return_value=None):
            status = check.check_trinity_rust()
        self.assertEqual(set(status), set(check.TRINITY_RUST_TOOLS))
        self.assertFalse(any(status.values()))

    def test_partial_install(self):
        # Only one tool on PATH -> mixed truthy/falsey, others None.
        def fake_which(tool):
            return "/env/bin/sam_to_read_coords" \
                if tool == "sam_to_read_coords" else None

        with mock.patch.object(check.lib, "which_path", side_effect=fake_which):
            status = check.check_trinity_rust()
        found = [t for t, p in status.items() if p]
        missing = [t for t, p in status.items() if not p]
        self.assertEqual(found, ["sam_to_read_coords"])
        self.assertEqual(
            set(missing), set(check.TRINITY_RUST_TOOLS) - {"sam_to_read_coords"}
        )

    def test_custom_tool_list(self):
        # Explicit tools argument is honored (and drives the query).
        with mock.patch.object(
            check.lib, "which_path", return_value=None
        ) as m:
            status = check.check_trinity_rust(tools=["only_tool"])
        self.assertEqual(status, {"only_tool": None})
        m.assert_called_once_with("only_tool")


class CheckEvmRustTests(unittest.TestCase):
    def test_present(self):
        with mock.patch.object(
            check.lib, "which_path", return_value="/env/bin/evidence_modeler"
        ) as m:
            self.assertEqual(check.check_evm_rust(), "/env/bin/evidence_modeler")
        m.assert_called_once_with("evidence_modeler")

    def test_absent(self):
        with mock.patch.object(check.lib, "which_path", return_value=None):
            self.assertIsNone(check.check_evm_rust())


class CheckPasaRustTests(unittest.TestCase):
    def test_found_in_pasahome_bin(self):
        # Prefer $PASAHOME/bin/pasa_rust even when it is not on PATH.
        with mock.patch.dict(os.environ, {"PASAHOME": "/opt/pasa/src"}):
            with mock.patch.object(check.os.path, "isfile", return_value=True):
                with mock.patch.object(check.os, "access", return_value=True):
                    with mock.patch.object(
                        check.lib, "which_path", return_value=None
                    ) as m:
                        result = check.check_pasa_rust()
        self.assertEqual(result, os.path.join("/opt/pasa/src", "bin", "pasa_rust"))
        m.assert_not_called()  # resolved via PASAHOME, never touched PATH

    def test_falls_back_to_path(self):
        # No PASAHOME hit -> fall back to PATH lookup.
        with mock.patch.dict(os.environ, {}, clear=True):
            with mock.patch.object(
                check.lib, "which_path", return_value="/env/bin/pasa_rust"
            ) as m:
                result = check.check_pasa_rust()
        self.assertEqual(result, "/env/bin/pasa_rust")
        m.assert_called_once_with("pasa_rust")

    def test_absent_everywhere(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            with mock.patch.object(check.lib, "which_path", return_value=None):
                self.assertIsNone(check.check_pasa_rust())


if __name__ == "__main__":
    unittest.main()
