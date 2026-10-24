"""Gates for GeneMark and --table in predict.

* check_tables: what predict does with --table / --mtable (errors stop the run, a stop-codon
  reassigning table warns that SNAP/GlimmerHMM/CodingQuarry cannot use it);
* genemark_gcode: whether GeneMark is told the genetic code (only tables 6 and 26, only versions
  that accept --gcode); GeneMark runs inside funannotate only in local (non-container) builds;
* resolve_genemark_path: --GENEMARK_PATH, then $GENEMARK_PATH, then PATH, else skip GeneMark;
* filter_gff_contigs: rows on contigs the genome lacks are dropped and reported. The container
  route passes a GTF from a separate GeneMark/BRAKER run (--genemark_gtf), which can name contigs
  from a different version of the assembly.
"""
import os
import shutil
import tempfile
import unittest

from funannotate import library as lib
from funannotate.genetic_codes import check_tables


class CheckTablesTests(unittest.TestCase):
    def test_standard_table_is_silent(self):
        self.assertEqual(check_tables(1), ([], []))
        self.assertEqual(check_tables(1, None), ([], []))

    def test_cug_table_has_no_warning(self):
        # table 12 changes a sense codon only; every predictor can use it
        self.assertEqual(check_tables(12), ([], []))

    def test_stop_reassigning_table_warns_about_predictors(self):
        errors, warnings = check_tables(4)
        self.assertEqual(errors, [])
        self.assertEqual(len(warnings), 1)
        for word in ("--table 4", "SNAP", "GlimmerHMM", "CodingQuarry", "Augustus", "GeneMark"):
            self.assertIn(word, warnings[0])

    def test_table_6_warns_too(self):
        self.assertEqual(len(check_tables(6)[1]), 1)

    def test_invalid_table_is_an_error(self):
        errors, warnings = check_tables(99)
        self.assertEqual(errors, ["--table 99 is not a valid NCBI translation table id"])
        self.assertEqual(warnings, [])

    def test_invalid_mtable_is_an_error(self):
        errors, _ = check_tables(1, 7)
        self.assertEqual(errors, ["--mtable 7 is not a valid NCBI translation table id"])

    def test_both_errors_are_reported_together(self):
        errors, _ = check_tables(99, 98)
        self.assertEqual(len(errors), 2)

    def test_valid_mtable_is_accepted(self):
        self.assertEqual(check_tables(1, 4), ([], []))

    def test_an_error_suppresses_the_warning(self):
        # table 4 would warn, but the invalid mtable stops the run first
        errors, warnings = check_tables(4, 99)
        self.assertEqual(len(errors), 1)
        self.assertEqual(warnings, [])


class GenemarkGcodeTests(unittest.TestCase):
    @staticmethod
    def probe(result, calls):
        def supports():
            calls.append(1)
            return result
        return supports

    def test_standard_table_adds_nothing_and_does_not_probe(self):
        calls = []
        self.assertEqual(lib.genemark_gcode(1, self.probe(True, calls)), ([], None))
        self.assertEqual(calls, [])

    def test_tables_6_and_26_pass_gcode_when_supported(self):
        for table in (6, 26):
            args, warning = lib.genemark_gcode(table, self.probe(True, []))
            self.assertEqual(args, ["--gcode", str(table)])
            self.assertIsNone(warning)

    def test_tables_6_and_26_warn_when_genemark_is_too_old(self):
        calls = []
        args, warning = lib.genemark_gcode(6, self.probe(False, calls))
        self.assertEqual(args, [])
        self.assertIn("does not support --gcode", warning)
        self.assertEqual(calls, [1])

    def test_other_tables_are_ignored_with_a_warning_and_no_probe(self):
        for table in (4, 12):
            calls = []
            args, warning = lib.genemark_gcode(table, self.probe(True, calls))
            self.assertEqual(args, [])
            self.assertIn("only supports --gcode 6 or 26", warning)
            self.assertIn("ignoring table {}".format(table), warning)
            self.assertEqual(calls, [])

    def test_table_given_as_string(self):
        self.assertEqual(lib.genemark_gcode("26", self.probe(True, []))[0], ["--gcode", "26"])


class ResolveGenemarkPathTests(unittest.TestCase):
    @staticmethod
    def finder(path, calls=None):
        def find():
            if calls is not None:
                calls.append(1)
            return path
        return find

    def test_option_wins_over_everything(self):
        calls = []
        result = lib.resolve_genemark_path(" /opt/gm ", "/env/gm", self.finder("/usr/bin/gmes_petap.pl", calls))
        self.assertEqual(result, ("/opt/gm", "--GENEMARK_PATH"))
        self.assertEqual(calls, [])

    def test_environment_is_used_without_the_option(self):
        calls = []
        result = lib.resolve_genemark_path(None, "/env/gm ", self.finder("/usr/bin/gmes_petap.pl", calls))
        self.assertEqual(result, ("/env/gm", "$GENEMARK_PATH"))
        self.assertEqual(calls, [])

    def test_empty_option_falls_through_to_environment(self):
        self.assertEqual(lib.resolve_genemark_path("", "/env/gm", self.finder(None))[0], "/env/gm")

    def test_empty_but_set_environment_still_beats_path(self):
        # unchanged behavior: a set-but-empty $GENEMARK_PATH is used as is (no PATH search)
        calls = []
        result = lib.resolve_genemark_path(None, "", self.finder("/usr/bin/gmes_petap.pl", calls))
        self.assertEqual(result, ("", "$GENEMARK_PATH"))
        self.assertEqual(calls, [])

    def test_path_lookup_returns_the_directory(self):
        result = lib.resolve_genemark_path(None, None, self.finder("/opt/genemark/gmes_petap.pl"))
        self.assertEqual(result, ("/opt/genemark", "PATH"))

    def test_nothing_found_means_skip_genemark(self):
        self.assertEqual(lib.resolve_genemark_path(None, None, self.finder(None)), (False, None))


class FilterGffContigsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.src = os.path.join(self.tmp, "genemark.gff3")
        self.dst = os.path.join(self.tmp, "filtered.gff3")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def write(self, text):
        with open(self.src, "w") as f:
            f.write(text)

    def read(self):
        with open(self.dst) as f:
            return f.read()

    def test_keeps_known_contigs_and_comments_and_reports_the_rest(self):
        self.write("##gff-version 3\nctg1\tGeneMark\tgene\t1\t9\t.\t+\t.\tID=a\n"
                   "scaffold_7\tGeneMark\tgene\t1\t9\t.\t+\t.\tID=b\n\nctg2\tGeneMark\tgene\t1\t9\t.\t-\t.\tID=c\n")
        missing = lib.filter_gff_contigs(self.src, self.dst, {"ctg1": 100, "ctg2": 50})
        self.assertEqual(missing, {"scaffold_7"})
        out = self.read()
        self.assertIn("ID=a", out)
        self.assertIn("ID=c", out)
        self.assertNotIn("ID=b", out)
        self.assertIn("##gff-version 3", out)

    def test_nothing_missing_copies_everything(self):
        text = "ctg1\tGeneMark\tgene\t1\t9\t.\t+\t.\tID=a\n"
        self.write(text)
        self.assertEqual(lib.filter_gff_contigs(self.src, self.dst, {"ctg1": 10}), set())
        self.assertEqual(self.read(), text)

    def test_gtf_from_a_renamed_assembly_loses_every_row(self):
        # a sidecar GTF made before sort/clean renamed the headers: all rows are dropped
        self.write("contig_1\tGeneMark\tCDS\t1\t9\t.\t+\t0\tgene_id \"1\";\n"
                   "contig_2\tGeneMark\tCDS\t1\t9\t.\t+\t0\tgene_id \"2\";\n")
        missing = lib.filter_gff_contigs(self.src, self.dst, {"scaffold_1": 10, "scaffold_2": 10})
        self.assertEqual(missing, {"contig_1", "contig_2"})
        self.assertEqual(self.read(), "")

    def test_empty_file(self):
        self.write("")
        self.assertEqual(lib.filter_gff_contigs(self.src, self.dst, {"ctg1": 10}), set())
        self.assertEqual(self.read(), "")


if __name__ == "__main__":
    unittest.main()
