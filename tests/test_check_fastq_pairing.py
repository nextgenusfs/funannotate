"""Regression test for GitHub issue #1212.

CheckFASTQandFix rejected correctly paired reads whose header has a
space-separated description and ends in /1 and /2, e.g. reads from
fasterq-dump with original names:

    @SRR14126314.1 A00953:17:HTJKFDSXX:4:1101:7148:1000/1

A header with a space was only accepted in the Illumina form (description
starting with 1/2); the /1 /2 test was an elif that was never reached.
"""
import importlib
import os
import tempfile
import unittest
from unittest import mock


lib = importlib.import_module("funannotate.library")


def run_check(h1, h2):
    with tempfile.TemporaryDirectory() as tmp:
        r1 = os.path.join(tmp, "R1.fastq")
        r2 = os.path.join(tmp, "R2.fastq")
        for path, header in ((r1, h1), (r2, h2)):
            with open(path, "w") as fh:
                fh.write("@%s\nACGT\n+\nIIII\n" % header)
        with mock.patch.object(lib, "log", mock.MagicMock(), create=True):
            try:
                return lib.CheckFASTQandFix(r1, r2)
            except SystemExit:
                return "exit"


class CheckFASTQPairingTests(unittest.TestCase):
    def test_slash_suffix_with_description(self):           # issue #1212
        self.assertEqual(run_check(
            "SRR14126314.1 A00953:17:HTJKFDSXX:4:1101:7148:1000/1",
            "SRR14126314.1 A00953:17:HTJKFDSXX:4:1101:7148:1000/2"), 0)

    def test_slash_suffix_on_name_before_description(self):
        self.assertEqual(run_check("read1/1 extra", "read1/2 extra"), 0)

    def test_slash_suffix_no_description(self):
        self.assertEqual(run_check("read1/1", "read1/2"), 0)

    def test_illumina_description(self):
        self.assertEqual(run_check(
            "A00953:17:HTJKFDSXX:4:1101:7148:1000 1:N:0:ACGT",
            "A00953:17:HTJKFDSXX:4:1101:7148:1000 2:N:0:ACGT"), 0)

    def test_no_pairing_information_exits(self):
        self.assertEqual(run_check("read1", "read1"), "exit")

    def test_description_without_pairing_exits(self):
        self.assertEqual(run_check("read1 sample=x", "read1 sample=x"), "exit")

    def test_swapped_pairs_exit(self):
        self.assertEqual(run_check("read1/2", "read1/1"), "exit")


if __name__ == "__main__":
    unittest.main()
