"""BUSCO training must pass the translation table to tblastn (-db_gencode).

Previously only Augustus received --translation_table, so candidate regions were
found with the standard genetic code and genes using TGA=Trp (tables 4, 3, ...)
were split at in-frame TGA "stops".

funannotate-BUSCO2.py is not importable in tests (it exits without Augustus
tools on PATH), so its plumbing is checked by source. The tblastn integration
test needs BLAST+ and is skipped when it is not installed.
"""
import importlib
import os
import random
import shutil
import subprocess
import tempfile
import unittest

lib = importlib.import_module("funannotate.library")
PKG = os.path.dirname(importlib.import_module("funannotate").__file__)
BUSCO2 = os.path.join(PKG, "aux_scripts", "funannotate-BUSCO2.py")
PREDICT = os.path.join(PKG, "predict.py")


class TranslationArgsTests(unittest.TestCase):
    def test_standard_code_adds_nothing(self):
        self.assertEqual(lib.busco_translation_args(1), [])
        self.assertEqual(lib.busco_translation_args("1"), [])

    def test_alternate_code_sets_augustus_and_tblastn(self):
        self.assertEqual(
            lib.busco_translation_args(4),
            ["-a", "--translation_table=4", "--db_gencode", "4"],
        )
        self.assertEqual(lib.busco_translation_args("12")[-1], "12")


class PlumbingTests(unittest.TestCase):
    def test_predict_uses_helper(self):
        with open(PREDICT) as fh:
            src = fh.read()
        self.assertIn("lib.busco_translation_args(args.table)", src)

    def test_busco2_passes_db_gencode_to_tblastn(self):
        with open(BUSCO2) as fh:
            src = fh.read()
        self.assertIn("'--db_gencode'", src)
        self.assertIn("'-db_gencode'", src)
        self.assertIn("gencode_args", src)
        self.assertIn("params['db_gencode']", src)
        # restart command keeps the option
        self.assertIn("--db_gencode %s", src)


@unittest.skipUnless(
    shutil.which("tblastn") and shutil.which("makeblastdb"), "BLAST+ not installed"
)
class TblastnGencodeTests(unittest.TestCase):
    """A gene with in-frame TGA (Trp in table 4) aligns perfectly only with -db_gencode 4.

    BLAST scores stop codons as mismatches and extends through them, so with the
    standard code the hit is not split; it loses identity and bit score at every
    in-frame TGA instead.
    """

    CODONS = {
        "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
        "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
        "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
        "S": "TCT", "T": "ACT", "V": "GTT", "Y": "TAT",
    }

    def _setup(self, d):
        rng = random.Random(7)
        aas = list(self.CODONS)
        # a protein with a Trp every ~25 residues, encoded with TGA
        prot = ["M"]
        for i in range(300):
            prot.append("W" if i % 25 == 24 else rng.choice(aas))
        prot = "".join(prot)
        cds = "".join("TGA" if a == "W" else self.CODONS[a] for a in prot)
        pad = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        with open(os.path.join(d, "genome.fa"), "w") as fh:
            fh.write(">ctg1\n%s%s%s\n" % (pad(800), cds, pad(800)))
        with open(os.path.join(d, "query.fa"), "w") as fh:
            fh.write(">q\n%s\n" % prot)
        subprocess.run(
            ["makeblastdb", "-in", os.path.join(d, "genome.fa"), "-dbtype", "nucl",
             "-out", os.path.join(d, "genome")],
            check=True, stdout=subprocess.DEVNULL,
        )
        return len(prot)

    def _best_hit(self, d, extra):
        """Return (identical positions, bitscore) of the best hit."""
        out = os.path.join(d, "hits.tsv")
        subprocess.run(
            ["tblastn", "-evalue", "1e-3", "-query", os.path.join(d, "query.fa"),
             "-db", os.path.join(d, "genome"), "-out", out,
             "-outfmt", "6 nident bitscore"] + extra,
            check=True,
        )
        with open(out) as fh:
            hits = [line.split() for line in fh if line.strip()]
        self.assertTrue(hits, "no tblastn hits")
        return max((int(n), float(b)) for n, b in hits)

    def test_db_gencode_recovers_full_identity(self):
        with tempfile.TemporaryDirectory() as d:
            plen = self._setup(d)
            std_ident, std_score = self._best_hit(d, [])
            alt_ident, alt_score = self._best_hit(d, ["-db_gencode", "4"])
            self.assertEqual(alt_ident, plen)
            self.assertLess(std_ident, plen)
            self.assertGreater(alt_score, std_score)


if __name__ == "__main__":
    unittest.main()
