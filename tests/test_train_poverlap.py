import unittest

from funannotate.train import pOverlap


class POverlapTests(unittest.TestCase):
    def test_full_overlap(self):
        # two identical intervals -> 100% overlap
        self.assertAlmostEqual(pOverlap([100, 200], [100, 200]), 1.0)

    def test_no_overlap(self):
        # disjoint intervals -> 0% overlap
        self.assertAlmostEqual(pOverlap([100, 200], [300, 400]), 0.0)

    def test_partial_overlap(self):
        # [100,200) vs [150,250): overlap is [150,200) = 50 bp of a 100 bp query
        self.assertAlmostEqual(pOverlap([100, 200], [150, 250]), 0.5)

    def test_asymmetric_overlap(self):
        # overlap fraction is relative to the FIRST interval's length only
        # query [100,200) fully contained in target [50,300) -> 100% of query overlaps
        self.assertAlmostEqual(pOverlap([100, 200], [50, 300]), 1.0)
        # but the reverse direction differs: query [50,300) (250 bp) only
        # overlaps [100,200) (100 bp) of itself -> 100/250 = 0.4
        self.assertAlmostEqual(pOverlap([50, 300], [100, 200]), 0.4)

    def test_touching_intervals_do_not_overlap(self):
        # half-open intervals sharing only a boundary point -> 0% overlap
        self.assertAlmostEqual(pOverlap([100, 200], [200, 300]), 0.0)

    def test_large_interval_performance_and_correctness(self):
        # exercises the O(1) interval-arithmetic path (formerly built a
        # ~1e6-element Python set per call for intervals this size)
        one = [0, 1_000_000]
        two = [500_000, 1_500_000]
        self.assertAlmostEqual(pOverlap(one, two), 0.5)


if __name__ == "__main__":
    unittest.main()
