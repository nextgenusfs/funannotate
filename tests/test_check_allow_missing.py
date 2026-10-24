"""`funannotate check --strict --allow-missing a,b`: expected-missing tools must not fail strict.

The container image deliberately has no GeneMark (license), so a strict check of the image has to
name gmes_petap.pl and GENEMARK_PATH as allowed to be missing.
"""
import io
import unittest

from funannotate.check import ErrorCounter, parse_allow_missing


class ParseAllowMissingTests(unittest.TestCase):
    def test_separate_value(self):
        self.assertEqual(parse_allow_missing(["--strict", "--allow-missing", "a,b"]), {"a", "b"})

    def test_equals_form(self):
        self.assertEqual(parse_allow_missing(["--allow-missing=a,b"]), {"a", "b"})

    def test_repeated_option_adds_up(self):
        self.assertEqual(parse_allow_missing(["--allow-missing", "a", "--allow-missing=b"]), {"a", "b"})

    def test_spaces_and_empty_items_are_ignored(self):
        self.assertEqual(parse_allow_missing(["--allow-missing", " a , ,b,"]), {"a", "b"})

    def test_absent_or_valueless_option(self):
        self.assertEqual(parse_allow_missing(["--strict"]), set())
        self.assertEqual(parse_allow_missing(["--allow-missing"]), set())


class ErrorCounterAllowedTests(unittest.TestCase):
    def counter(self, allowed=()):
        return ErrorCounter(io.StringIO(), allowed)

    def test_allowed_tool_is_not_counted(self):
        c = self.counter({"gmes_petap.pl"})
        c.write("\tERROR: gmes_petap.pl not installed\n")
        self.assertEqual(c.count, 0)

    def test_other_tools_are_still_counted(self):
        c = self.counter({"gmes_petap.pl"})
        c.write("\tERROR: blat not installed\n")
        c.write("\tERROR: gmes_petap.pl not installed\n")
        self.assertEqual(c.count, 1)

    def test_allowed_environment_variable(self):
        c = self.counter({"GENEMARK_PATH"})
        c.write("\tERROR: GENEMARK_PATH not set. export GENEMARK_PATH=/path/to/dir\n")
        self.assertEqual(c.count, 0)

    def test_found_but_failing_tool_is_counted_unless_named(self):
        text = "\tERROR: tbl2asn found but error running tbl2asn\n"
        c = self.counter()
        c.write(text)
        self.assertEqual(c.count, 1)
        c = self.counter({"tbl2asn"})
        c.write(text)
        self.assertEqual(c.count, 0)

    def test_several_lines_in_one_write(self):
        c = self.counter({"a"})
        c.write("ERROR: a missing\nERROR: b missing\nERROR: c missing\n")
        self.assertEqual(c.count, 2)

    def test_default_counts_everything_and_passes_text_through(self):
        out = io.StringIO()
        c = ErrorCounter(out)
        c.write("ERROR: x missing\n")
        self.assertEqual(c.count, 1)
        self.assertEqual(out.getvalue(), "ERROR: x missing\n")

    def test_allow_matches_the_whole_name(self):
        c = self.counter({"blat"})
        c.write("\tERROR: blatx not installed\n")
        self.assertEqual(c.count, 1)


if __name__ == "__main__":
    unittest.main()
