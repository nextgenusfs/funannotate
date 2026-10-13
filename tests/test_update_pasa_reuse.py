"""Regression test for GitHub issue #1213.

`funannotate update` reused training/pasa/alignAssembly.txt without checking
that the PASA database it names still exists. A MySQL database on a server
that only lived for the train job, or a deleted SQLite file, made update fail
late in PASA.

Wanted behaviour:
  - database gone            -> build a new one (check returns False)
  - database present, OK     -> reuse it (check returns True)
  - database present but empty, or its contigs are not in the genome,
    or its state cannot be checked -> stop with an error
"""
import importlib
import os
import subprocess
import tempfile
import unittest
from unittest import mock


update = importlib.import_module("funannotate.update")


class FakeDump:
    """Stand-in for subprocess.run of pasa_asmbl_genes_to_GFF3.dbi."""

    def __init__(self, gff="", returncode=0, stderr=""):
        self.gff, self.returncode, self.stderr = gff, returncode, stderr
        self.calls = []

    def __call__(self, cmd, cwd=None, stdout=None, stderr=None, universal_newlines=None):
        self.calls.append(cmd)
        stdout.write(self.gff)
        return subprocess.CompletedProcess(cmd, self.returncode, None, self.stderr)


GFF_OK = (
    "scaffold_1\tPASA\tgene\t1\t900\t.\t+\t.\tID=g1\n"
    "scaffold_1\tPASA\tmRNA\t1\t900\t.\t+\t.\tID=m1;Parent=g1\n"
    "scaffold_2\tPASA\tgene\t1\t500\t.\t-\t.\tID=g2\n"
)


class CheckExistingPasaDbTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = self.tmp.name
        self.folder = d
        self.genome = os.path.join(d, "genome.fa")
        with open(self.genome, "w") as fh:
            fh.write(">scaffold_1\nACGT\n>scaffold_2\nACGT\n")
        pasahome = os.path.join(d, "pasa")
        os.makedirs(os.path.join(pasahome, "pasa_conf"))
        with open(os.path.join(pasahome, "pasa_conf", "conf.txt"), "w") as fh:
            fh.write("MYSQLSERVER=localhost\nMYSQL_RW_USER=u\nMYSQL_RW_PASSWORD=p\n")
        self.sqlite = os.path.join(d, "Org_pasa")
        self.patches = [
            mock.patch.object(update, "PASA", pasahome, create=True),
            mock.patch.object(update.lib, "log", mock.MagicMock(), create=True),
            mock.patch.dict(os.environ, {}, clear=False),
        ]
        os.environ.pop("PASACONF", None)
        for p in self.patches:
            p.start()

    def tearDown(self):
        for p in self.patches:
            p.stop()
        self.tmp.cleanup()

    def run_check(self, dbname, fake):
        with mock.patch.object(update.subprocess, "run", fake):
            try:
                return update.check_existing_pasa_db(dbname, self.folder, self.genome)
            except SystemExit:
                return "exit"

    def test_sqlite_missing_rebuilds(self):
        fake = FakeDump()
        self.assertIs(self.run_check(self.sqlite, fake), False)
        self.assertEqual(fake.calls, [])

    def test_sqlite_empty_file_rebuilds(self):
        open(self.sqlite, "w").close()
        self.assertIs(self.run_check(self.sqlite, FakeDump()), False)

    def test_sqlite_present_and_matching_reused(self):
        with open(self.sqlite, "w") as fh:
            fh.write("db")
        fake = FakeDump(gff=GFF_OK)
        self.assertIs(self.run_check(self.sqlite, fake), True)
        self.assertEqual(fake.calls[0][2], self.sqlite + ":localhost")

    def test_sqlite_present_wrong_genome_exits(self):
        with open(self.sqlite, "w") as fh:
            fh.write("db")
        gff = GFF_OK + "other_contig\tPASA\tgene\t1\t10\t.\t+\t.\tID=g3\n"
        self.assertEqual(self.run_check(self.sqlite, FakeDump(gff=gff)), "exit")

    def test_sqlite_present_no_gene_models_exits(self):
        with open(self.sqlite, "w") as fh:
            fh.write("db")
        self.assertEqual(self.run_check(self.sqlite, FakeDump(gff="")), "exit")

    def test_mysql_unknown_database_rebuilds(self):              # issue #1213
        fake = FakeDump(returncode=255, stderr=(
            "DBI connect('database=Org_pasa;host=127.0.0.1:4652','u',...) failed: "
            "Unknown database 'Org_pasa'"))
        self.assertIs(self.run_check("Org_pasa", fake), False)
        self.assertEqual(fake.calls[0][2], "Org_pasa:localhost")

    def test_mysql_server_unreachable_exits(self):
        fake = FakeDump(returncode=255, stderr="Can't connect to MySQL server on '127.0.0.1'")
        self.assertEqual(self.run_check("Org_pasa", fake), "exit")

    def test_mysql_present_and_matching_reused(self):
        self.assertIs(self.run_check("Org_pasa", FakeDump(gff=GFF_OK)), True)

    def test_mysql_present_wrong_genome_exits(self):
        gff = "other_contig\tPASA\tgene\t1\t10\t.\t+\t.\tID=g3\n"
        self.assertEqual(self.run_check("Org_pasa", FakeDump(gff=gff)), "exit")


class StopRun(Exception):
    pass


class RunPasaReuseWiringTests(unittest.TestCase):
    """runPASA must build a new database when the reused one is gone."""

    def run_pasa(self, db_exists):
        with tempfile.TemporaryDirectory() as d:
            pasahome = os.path.join(d, "pasa")
            os.makedirs(os.path.join(pasahome, "pasa_conf"))
            with open(os.path.join(pasahome, "pasa_conf", "pasa.alignAssembly.Template.txt"), "w") as fh:
                fh.write("DATABASE=<__DATABASE__>\n")
            config = os.path.join(d, "alignAssembly.txt")
            with open(config, "w") as fh:
                fh.write("DATABASE=/old/train/Org_pasa\n")
            calls = []

            def fake_run(cmd, *_a, **_k):
                calls.append(cmd)
                raise StopRun()

            with mock.patch.object(update, "PASA", pasahome, create=True), \
                    mock.patch.object(update, "LAUNCHPASA", "Launch_PASA_pipeline.pl", create=True), \
                    mock.patch.object(update, "tmpdir", d, create=True), \
                    mock.patch.object(update, "check_existing_pasa_db", lambda *_a: db_exists), \
                    mock.patch.object(update.lib, "log", mock.MagicMock(), create=True), \
                    mock.patch.object(update.lib, "runSubprocess", fake_run), \
                    mock.patch.object(update.lib, "countfasta", lambda *_a: 1), \
                    mock.patch.object(update.lib, "which_path", lambda *_a: "cdbfasta"):
                with self.assertRaises(StopRun):
                    update.runPASA("genome.fa", "t.fa", "t.clean.fa", "align.gff3", "st.gtf",
                                   "no", 3000, 1, "prev.gff3", "Org", "out.gff3", config,
                                   pasa_db="sqlite", aligners=["blat"])
            with open(os.path.join(d, "pasa", "alignAssembly.txt")) as fh:
                align_config = fh.read()
        return calls[0], align_config, d

    def test_missing_db_builds_new_database(self):
        cmd, align_config, d = self.run_pasa(db_exists=False)
        self.assertEqual(cmd[0], "Launch_PASA_pipeline.pl")
        self.assertIn("-C", cmd)
        self.assertEqual(align_config, "DATABASE=%s\n" % os.path.join(d, "pasa", "Org_pasa"))

    def test_present_db_is_reused(self):
        cmd, align_config, _d = self.run_pasa(db_exists=True)
        self.assertEqual(cmd[0], "cdbfasta")
        self.assertEqual(align_config, "DATABASE=/old/train/Org_pasa\n")


if __name__ == "__main__":
    unittest.main()
