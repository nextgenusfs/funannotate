"""`funannotate util gff2prot` and `gff-rename` pass --table to library.gff2dict (whose argument is transl_table).

Both crashed with a TypeError (unexpected keyword 'table') before the fix. Uses the two-gene fixture of
test_genetic_code_tables: g1 reads MLKWG* under table 4 and MLK*G* under table 1.
"""
import contextlib
import io
import os

from test_genetic_code_tables import GeneModelFixture

from funannotate.utilities import gff2prot, gff_reformat


class UtilTableOptionTests(GeneModelFixture):
    def run_gff2prot(self, *extra):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            gff2prot.main(["-g", self.gff3, "-f", self.fasta, *extra])
        return out.getvalue()

    def test_gff2prot_default_table_is_standard_code(self):
        text = self.run_gff2prot()
        self.assertIn("MLK*G", text)
        self.assertNotIn("MLKWG", text)

    def test_gff2prot_table_4_reads_tga_as_tryptophan(self):
        self.assertIn("MLKWG*", self.run_gff2prot("--table", "4"))

    def test_gff2prot_table_12_translates_cug_as_serine(self):
        self.assertIn("MSK", self.run_gff2prot("--table", "12"))

    def test_gff2prot_rejects_an_unassigned_table(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            self.run_gff2prot("--table", "7")

    def test_gff_rename_accepts_table(self):
        out = os.path.join(self.tmp, "renamed.gff3")
        with contextlib.redirect_stdout(io.StringIO()):
            gff_reformat.main(["-g", self.gff3, "-f", self.fasta, "-o", out, "-l", "TST", "--table", "4"])
        with open(out) as f:
            text = f.read()
        self.assertIn("TST_000001", text)
        self.assertIn("TST_000002", text)
