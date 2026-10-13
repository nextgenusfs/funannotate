"""Regression test: predict crashed with UnboundLocalError on AUGUSTUS_BASE.

predict.py set AUGUSTUS_BASE only when $AUGUSTUS_CONFIG_PATH is a directory
named "config". With a writable copy under another name (e.g. a container user
copying the read-only config to $SCRATCH/augustus_config), BUSCO-trained
Augustus crashed at trainAugustus(AUGUSTUS_BASE, ...).

Found testing the v1.9.0-rc.5 container with `funannotate test -t busco`.
"""
import importlib
import os
import stat
import tempfile
import unittest
from unittest import mock


predict = importlib.import_module("funannotate.predict")


class AugustusBaseTests(unittest.TestCase):
    def test_config_dir_named_config(self):
        self.assertEqual(predict.find_augustus_base("/opt/augustus/config"), "/opt/augustus")
        self.assertEqual(predict.find_augustus_base("/opt/augustus/config/"), "/opt/augustus")

    def test_other_name_uses_augustus_binary(self):
        with tempfile.TemporaryDirectory() as d:
            bindir = os.path.join(d, "env", "bin")
            os.makedirs(bindir)
            exe = os.path.join(bindir, "augustus")
            with open(exe, "w") as fh:
                fh.write("#!/bin/sh\n")
            os.chmod(exe, os.stat(exe).st_mode | stat.S_IEXEC)
            with mock.patch.dict(os.environ, {"PATH": bindir}):
                self.assertEqual(predict.find_augustus_base("/scratch/x/augustus_config"),
                                 os.path.realpath(os.path.join(d, "env")))

    def test_other_name_without_binary_is_still_defined(self):
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.dict(os.environ, {"PATH": d}):
                base = predict.find_augustus_base("/scratch/x/augustus_config")
        self.assertEqual(base, "/scratch/x")


if __name__ == "__main__":
    unittest.main()
