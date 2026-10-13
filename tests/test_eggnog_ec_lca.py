"""Regression test for GitHub issue #1207.

The eggNOG parser reduced comma-separated EC numbers with os.path.commonprefix,
which compares characters, not EC levels. It could invent an unrelated
4-level EC (1.1.1.10,1.1.1.1 -> 1.1.1.1) or a malformed one
(1.14.18.5,1.14.19.17 -> 1.14.1).
"""
import importlib
import unittest


annotate = importlib.import_module("funannotate.annotate")


class ECCommonAncestorTests(unittest.TestCase):
    def check(self, ecs, expected):
        self.assertEqual(annotate.ec_common_ancestor(ecs), expected)

    def test_does_not_invent_4_level_ec(self):
        self.check("1.1.1.10,1.1.1.1", "1.1.1")

    def test_does_not_split_a_level(self):
        self.check("1.14.18.5,1.14.19.17", "1.14")
        self.check("2.7.11.1,2.7.1.1", "2.7")

    def test_identical(self):
        self.check("1.1.1.1,1.1.1.1", "1.1.1.1")

    def test_single(self):
        self.check("3.2.1.4", "3.2.1.4")

    def test_no_shared_class(self):
        self.check("1.1.1.1,2.1.1.1", "")

    def test_stops_at_dash(self):
        self.check("3.2.1.-,3.2.1.-", "3.2.1")

    def test_whitespace_and_empty_items(self):
        self.check(" 1.2.3.4, 1.2.3.5,", "1.2.3")


if __name__ == "__main__":
    unittest.main()
