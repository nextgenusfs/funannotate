"""annotate.parse_version replaces distutils.LooseVersion for emapper and diamond version checks."""
import unittest

from funannotate.annotate import parse_version


class ParseVersionTests(unittest.TestCase):
    def test_plain_versions_order_numerically(self):
        self.assertLess(parse_version("2.0.5"), parse_version("2.1.2"))
        self.assertLess(parse_version("2.1.4"), parse_version("2.1.12"))  # not a string comparison
        self.assertGreaterEqual(parse_version("2.1.12"), parse_version("2.1.4"))

    def test_trailing_git_describe_is_ignored(self):
        self.assertEqual(parse_version("2.1.4-5-gabc1234"), parse_version("2.1.4"))

    def test_leading_v_and_whitespace(self):
        self.assertEqual(parse_version(" v2.0.11\n"), parse_version("2.0.11"))

    def test_unparseable_or_missing_is_zero(self):
        self.assertEqual(parse_version(""), parse_version("0"))
        self.assertEqual(parse_version(None), parse_version("0"))
        self.assertEqual(parse_version("unknown"), parse_version("0"))

    def test_eggnog_window_check(self):
        # annotate refuses emapper versions above 2.0.0 and below 2.0.5
        bad = lambda v: parse_version("2.0.0") < parse_version(v) < parse_version("2.0.5")
        self.assertTrue(bad("2.0.1"))
        self.assertFalse(bad("2.0.5"))
        self.assertFalse(bad("1.0.3"))
        self.assertFalse(bad("2.1.12"))


if __name__ == "__main__":
    unittest.main()
