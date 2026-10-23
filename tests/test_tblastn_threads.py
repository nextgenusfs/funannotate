"""Characterization tests for the tblastn thread gate used in BUSCO training.

funannotate-BUSCO2.py and predict.py both gate on ``tblastn_version > "2.2.31"``.
That is a *string* comparison, so it does not mean "newer than 2.2.31":

    2.2.32 - 2.9.x   -> True  (tblastn forced to 1 thread)
    2.10 and later   -> False (tblastn uses all cpus)

These tests pin that existing behavior so any change is deliberate. If the gate
is ever switched to numeric comparison, update the expectations here (and expect
BUSCO training to become single-threaded on every current BLAST+).

The integration test needs BLAST+ and is skipped when it is not installed.
"""
import importlib
import os
import random
import shutil
import subprocess
import tempfile
import unittest

lib = importlib.import_module("funannotate.library")

AUX = os.path.join(
    os.path.dirname(importlib.import_module("funannotate").__file__),
    "aux_scripts",
    "funannotate-BUSCO2.py",
)
PREDICT = os.path.join(
    os.path.dirname(importlib.import_module("funannotate").__file__), "predict.py"
)


class ThreadGateTests(unittest.TestCase):
    def test_old_blast_is_not_single_threaded(self):
        for v in ["2.2.18", "2.2.28", "2.2.31"]:
            self.assertFalse(lib.tblastn_single_thread(v), v)

    def test_blast_2_2_32_to_2_9_is_single_threaded(self):
        for v in ["2.2.32", "2.3.0", "2.6.0", "2.8.1", "2.9.0"]:
            self.assertTrue(lib.tblastn_single_thread(v), v)

    def test_modern_blast_is_multithreaded(self):
        # lexicographic quirk: "2.10.x" < "2.2.31" as strings
        for v in ["2.10.0", "2.10.1", "2.12.0", "2.14.1", "2.16.0", "2.17.0"]:
            self.assertFalse(lib.tblastn_single_thread(v), v)


class GateSourceSyncTests(unittest.TestCase):
    """BUSCO2.py is not importable in tests (it exits without Augustus tools),
    so check by source that it still uses the same gate as the helper."""

    def test_busco2_uses_same_gate(self):
        with open(AUX) as fh:
            src = fh.read()
        self.assertIn("Analysis._vers_tblastn() > '2.2.31'", src)

    def test_predict_uses_helper(self):
        with open(PREDICT) as fh:
            src = fh.read()
        self.assertIn("lib.tblastn_single_thread(tblastn_version)", src)
        self.assertNotIn('tblastn_version > "2.2.31"', src)


@unittest.skipUnless(
    shutil.which("tblastn") and shutil.which("makeblastdb"), "BLAST+ not installed"
)
class TblastnThreadingTests(unittest.TestCase):
    """Is multi-threaded tblastn output identical to single-threaded?

    Evidence for whether the single-thread workaround is still needed on the
    installed BLAST+. Uses a small synthetic genome with planted protein hits.
    """

    CODONS = {
        "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
        "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
        "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
        "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
    }

    def _make_inputs(self, d, nprot=40, nplanted=40):
        rng = random.Random(42)
        aas = list(self.CODONS)
        prots = {
            "p%d" % i: "M" + "".join(rng.choice(aas) for _ in range(150))
            for i in range(nprot)
        }
        genome = []
        for i in range(nplanted):
            pad1 = "".join(rng.choice("ACGT") for _ in range(rng.randint(500, 2000)))
            pad2 = "".join(rng.choice("ACGT") for _ in range(rng.randint(500, 2000)))
            cds = "".join(self.CODONS[a] for a in prots["p%d" % (i % nprot)])
            genome.append((pad1, cds, pad2))
        with open(os.path.join(d, "genome.fa"), "w") as fh:
            for i, parts in enumerate(genome):
                fh.write(">ctg%d\n%s\n" % (i, "".join(parts)))
        with open(os.path.join(d, "query.fa"), "w") as fh:
            for k, v in prots.items():
                fh.write(">%s\n%s\n" % (k, v))

    def _run(self, d, threads):
        out = os.path.join(d, "hits_%d.tsv" % threads)
        subprocess.run(
            [
                "tblastn", "-evalue", "1e-3", "-num_threads", str(threads),
                "-query", os.path.join(d, "query.fa"),
                "-db", os.path.join(d, "genome"),
                "-out", out, "-outfmt", "6",
            ],
            check=True,
        )
        with open(out) as fh:
            return sorted(line.rstrip("\n") for line in fh)

    def test_threads_give_same_hits(self):
        with tempfile.TemporaryDirectory() as d:
            self._make_inputs(d)
            subprocess.run(
                [
                    "makeblastdb", "-in", os.path.join(d, "genome.fa"),
                    "-dbtype", "nucl", "-out", os.path.join(d, "genome"),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
            )
            one = self._run(d, 1)
            four = self._run(d, 4)
            self.assertGreater(len(one), 0)
            self.assertEqual(one, four)


if __name__ == "__main__":
    unittest.main()
