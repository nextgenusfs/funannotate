import os
import random
import shutil
import signal
import tempfile
import unittest
from concurrent.futures.process import BrokenProcessPool

import funannotate.library as lib


def _list2groups(L):
    if len(L) < 1:
        return
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last


def _reference(seq):
    """The pre-2026 per-base implementation, kept here as the oracle."""
    masked, gaps, size = [], [], 0
    for i, c in enumerate(seq):
        if c == "N" or c == "n":
            masked.append(i)
            size += 1
            gaps.append(i)
        elif c.islower():
            masked.append(i)
            size += 1
    return list(_list2groups(masked)), list(_list2groups(gaps)), size


def _new(seq):
    rep = [(m.start(), m.end() - 1) for m in lib._MASKED_RUN.finditer(seq)]
    gaps = [(m.start(), m.end() - 1) for m in lib._GAP_RUN.finditer(seq)]
    return rep, gaps, sum(b - a + 1 for a, b in rep)


def _suicidal_worker(path):
    os.kill(os.getpid(), signal.SIGKILL)


def _random_seq(rng, n):
    out = []
    while sum(map(len, out)) < n:
        kind = rng.choice(["ACGT", "acgt", "N", "n", "ACGT", "acgt"])
        out.append("".join(rng.choice(kind) for _ in range(rng.randint(1, 40))))
    return "".join(out)[:n]


class RunScanEquivalenceTests(unittest.TestCase):
    EDGE_CASES = [
        "",
        "acgt" * 10,
        "N" * 7,
        "n" * 7,
        "ACGT" * 5,
        "ACGnnnNNNacgtACGT",  # gap run touching a lowercase run on both sides
        "acgNNNacg",
        "a",
        "N",
        "ACGrykmswbdhvACG",  # lowercase IUPAC ambiguity codes count as masked
        "ACGRYKMSWBDHVACG",  # uppercase ones do not
        "NNacgtNN",
        "nNnNACGTnN",
        "ACG-*acg",
    ]

    def test_edge_cases_match_reference(self):
        for seq in self.EDGE_CASES:
            self.assertEqual(_new(seq), _reference(seq), repr(seq))

    def test_random_sequences_match_reference(self):
        rng = random.Random(20260906)
        for _ in range(200):
            seq = _random_seq(rng, rng.randint(0, 3000))
            self.assertEqual(_new(seq), _reference(seq), repr(seq[:60]))


class MaskingStats2BedTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _write(self, name, seq):
        path = os.path.join(self.tmp, name + ".fasta")
        with open(path, "w") as fh:
            fh.write(">{} some description\n{}\n".format(name, seq))
        return path

    def test_writes_bed_and_gaps_and_returns_size(self):
        path = self._write("scaf1", "ACGTacgtNNNNACGTnnACGTacgt")
        size = lib.maskingstats2bed(path)
        # acgt(4) + NNNN(4) + nn(2) + acgt(4)
        self.assertEqual(size, 14)
        with open(path.replace(".fasta", ".bed")) as fh:
            bed = [l.rstrip("\n").split("\t") for l in fh]
        # the lowercase run at 4-7 abuts the N run at 8-11, so they form ONE
        # masked run (N bases are masked too) -- same as the historical list2groups
        self.assertEqual(
            bed,
            [
                ["scaf1", "4", "11", "Repeat_"],
                ["scaf1", "16", "17", "Repeat_"],
                ["scaf1", "22", "25", "Repeat_"],
            ],
        )
        with open(path.replace(".fasta", ".gaps")) as fh:
            gaps = [l.rstrip("\n").split("\t") for l in fh]
        self.assertEqual(
            gaps,
            [["scaf1", "8", "11", "assembly-gap_"], ["scaf1", "16", "17", "assembly-gap_"]],
        )

    def test_unmasked_scaffold_writes_nothing(self):
        path = self._write("scaf2", "ACGTACGT")
        self.assertEqual(lib.maskingstats2bed(path), 0)
        self.assertFalse(os.path.exists(path.replace(".fasta", ".bed")))
        self.assertFalse(os.path.exists(path.replace(".fasta", ".gaps")))


class CheckMasklowMemTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genome = os.path.join(self.tmp, "genome.fa")
        with open(self.genome, "w") as fh:
            fh.write(">s10 desc\nACGTacgtNN\n>s2\nacgt\n>s1\nACGT\n")
        self.bed = os.path.join(self.tmp, "repeats.bed")
        self.gaps = os.path.join(self.tmp, "gaps.bed")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_end_to_end_totals_and_numbering(self):
        sizes, length, masked, pct = lib.checkMasklowMem(
            self.genome, self.bed, self.gaps, 2, tmpdir=os.path.join(self.tmp, "mask")
        )
        self.assertEqual(sizes, {"s10": 10, "s2": 4, "s1": 4})
        self.assertEqual(length, 18)
        self.assertEqual(masked, 10)
        self.assertAlmostEqual(pct, 10 / 18.0)
        self.assertFalse(os.path.isdir(os.path.join(self.tmp, "mask")))
        with open(self.bed) as fh:
            names = [l.split("\t")[3].strip() for l in fh]
        # natsorted file order (s1, s2, s10) and sequential numbering across files;
        # s1 is unmasked, s2 is one run, s10's acgt+NN abut into one run
        self.assertEqual(names, ["Repeat_1", "Repeat_2"])
        with open(self.gaps) as fh:
            self.assertEqual([l.split("\t")[3].strip() for l in fh], ["assembly-gap_1"])

    def test_killed_worker_raises_instead_of_hanging(self):
        # A worker taken out by SIGKILL (what the cgroup OOM killer does) must
        # surface as an error, not leave the parent waiting forever.
        original = lib.maskingstats2bed
        lib.maskingstats2bed = _suicidal_worker  # module-level, so it pickles
        signal.signal(signal.SIGALRM, lambda *a: self.fail("checkMasklowMem hung"))
        signal.alarm(60)
        try:
            with self.assertRaises(BrokenProcessPool):
                lib.checkMasklowMem(
                    self.genome, self.bed, self.gaps, 2, tmpdir=os.path.join(self.tmp, "m")
                )
        finally:
            signal.alarm(0)
            lib.maskingstats2bed = original


if __name__ == "__main__":
    unittest.main()
