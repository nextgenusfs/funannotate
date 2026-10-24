"""update.pairwiseAlign: percent identity of a global alignment with match=1, mismatch=0 and no gap cost.

The score is the longest common subsequence length. The expected values here are worked by hand, and a
random check compares against an independent LCS implementation (the former pairwise2.globalxx result).
"""
import random
import unittest

from funannotate.update import pairwiseAlign


def lcs(a, b):
    prev = [0] * (len(b) + 1)
    for x in a:
        cur = [0]
        for j, y in enumerate(b):
            cur.append(prev[j] + 1 if x == y else max(prev[j + 1], cur[j]))
        prev = cur
    return prev[-1]


class PairwiseAlignTests(unittest.TestCase):
    def test_identical_is_100(self):
        self.assertEqual(pairwiseAlign("MKVLA", "MKVLA"), 100.0)

    def test_one_substitution(self):
        self.assertAlmostEqual(pairwiseAlign("MKVLA", "MKALA"), 80.0)

    def test_gap_is_not_penalised(self):
        # the shorter protein is a subsequence of the longer: 4 matches over the longer length 5
        self.assertAlmostEqual(pairwiseAlign("MKVA", "MKVLA"), 80.0)

    def test_nothing_in_common(self):
        self.assertEqual(pairwiseAlign("AAAA", "CCCC"), 0.0)

    def test_matches_longest_common_subsequence(self):
        rng = random.Random(7)
        for _ in range(100):
            a = "".join(rng.choice("ACDEFGHIK") for _ in range(rng.randint(1, 40)))
            b = "".join(rng.choice("ACDEFGHIK") for _ in range(rng.randint(1, 40)))
            expected = lcs(a, b) / float(max(len(a), len(b))) * 100
            self.assertAlmostEqual(pairwiseAlign(a, b), expected, msg=(a, b))


if __name__ == "__main__":
    unittest.main()
