"""Regression test for GitHub issue #1196.

tRNAscan-SE 2.x reports pseudogenes in the Note (last) column and keeps the
isotype in the Type column, e.g. "Arg ... pseudo". trnascan2gff3.pl only
dropped rows whose Type was "Pseudo" (tRNAscan-SE 1.x), so pseudo tRNAs went
into the GFF3 as normal tRNA genes.

Rows below are copied from tRNAscan-SE 2.0.12 output for S. cerevisiae R64.
"""
import os
import shutil
import subprocess
import tempfile
import unittest


SCRIPT = os.path.join(os.path.dirname(__file__), "..", "funannotate", "aux_scripts", "trnascan2gff3.pl")

TRNASCAN_OUT = (
    "Sequence\t\ttRNA  \tBounds\ttRNA\tAnti\tIntron Bounds\tInf\t      \n"
    "Name    \ttRNA #\tBegin \tEnd   \tType\tCodon\tBegin\tEnd\tScore\tNote\n"
    "--------\t------\t----- \t------\t----\t-----\t-----\t----\t------\t------\n"
    "I       \t1\t139152\t139254\tPro\tTGG\t139188\t139218\t62.1\t\n"
    "Mito    \t8\t69289 \t69359 \tArg\tACG\t0\t0\t38.9\tpseudo\n"
    "Mito    \t1\t9374  \t9444  \tSup\tTCA\t0\t0\t48.3\t\n"
    "chrX    \t3\t5000  \t5072  \tPseudo\tAAA\t0\t0\t20.1\t\n"
)


@unittest.skipUnless(shutil.which("perl"), "perl not found")
class TrnaScan2Gff3Tests(unittest.TestCase):
    def run_script(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "tRNAscan.out")
            with open(path, "w") as fh:
                fh.write(TRNASCAN_OUT)
            out = subprocess.run(["perl", SCRIPT, "--input", path], check=True,
                                 capture_output=True, text=True).stdout
        return [line.split("\t") for line in out.splitlines() if line and not line.startswith("#")]

    def test_pseudo_note_is_dropped(self):                   # issue #1196
        rows = self.run_script()
        self.assertFalse([r for r in rows if r[0] == "Mito" and r[3] == "69289"])

    def test_normal_trna_is_kept(self):
        rows = self.run_script()
        trna = [r for r in rows if r[2] == "tRNA"]
        self.assertEqual(len(trna), 1)
        self.assertEqual((trna[0][0], trna[0][3], trna[0][4]), ("I", "139152", "139254"))
        self.assertIn("product=tRNA-Pro", trna[0][8])

    def test_v1_pseudo_type_and_suppressor_still_dropped(self):
        rows = self.run_script()
        self.assertFalse([r for r in rows if r[0] == "chrX" or (r[0] == "Mito" and r[3] == "9374")])


if __name__ == "__main__":
    unittest.main()
