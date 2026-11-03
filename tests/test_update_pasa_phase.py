"""Regression test: funannotate update must not write PASA models with wrong CDS phases.

A PASA build with a Gene_obj.pm get_exons() cache bug (hyphaltip/PASApipeline
f7cede7, 2.6.1-rc.2/rc.3) wrote CDS phases for minus-strand multi-exon genes as
if they were on the plus strand (counted from the 3' end). funannotate took the
5' CDS phase as codon_start, so 6-22% of mRNAs translated with internal stops
(Absidia caerulea RSA 138: 3,535 tbl2asn InternalStop errors). The log said
these models were "dropped", but they were written to the TBL/GBK.

Wanted behaviour:
  - correct PASA models (both strands)       -> TBL written, no internal stops
  - models with the bug's phases             -> update stops (SystemExit)
  - check_internal_stop_rate                 -> stops only above the limit

The gene structure follows FBN33_000004 (6 CDS segments, 630 bp; phases
5'->3' = 0,2,1,0,0,0); the sequence is synthetic.
"""
import importlib
import os
import random
import tempfile
import unittest
from unittest import mock


update = importlib.import_module("funannotate.update")

CDS_LENS = [49, 166, 220, 48, 60, 87]  # 5'->3', 630 bp = 210 codons
GOOD_PHASES = [0, 2, 1, 0, 0, 0]
INTRON = "GT" + "C" * 56 + "AG"
FLANK = "C" * 100
STOPS = {"TAA", "TAG", "TGA"}


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def make_orf(seed):
    rnd = random.Random(seed)
    codons = []
    while len(codons) < 208:
        c = "".join(rnd.choice("ACGT") for _ in range(3))
        if c not in STOPS:
            codons.append(c)
    return "ATG" + "".join(codons) + "TAA"


def build_gene(contig, strand, seed):
    """Return (contig_seq, gene GFF3 lines) for one 6-exon gene."""
    orf = make_orf(seed)
    plus, exons, off = FLANK, [], 0
    for i, n in enumerate(CDS_LENS):
        start = len(plus) + 1
        plus += orf[off:off + n]
        off += n
        exons.append((start, len(plus)))
        if i < len(CDS_LENS) - 1:
            plus += INTRON
    plus += FLANK
    if strand == "-":
        L = len(plus)
        seq = revcomp(plus)
        exons = [(L - e + 1, L - s + 1) for s, e in exons]  # 5'->3' order kept
    else:
        seq = plus
    return seq, exons


def gff_lines(contig, gid, strand, exons, phases):
    """exons and phases in 5'->3' order; written in genomic order as PASA does."""
    lo = min(s for s, _ in exons)
    hi = max(e for _, e in exons)
    out = [
        "{}\tPASA\tgene\t{}\t{}\t.\t{}\t.\tID={};".format(contig, lo, hi, strand, gid),
        "{}\tPASA\tmRNA\t{}\t{}\t.\t{}\t.\tID={}-T1;Parent={};".format(
            contig, lo, hi, strand, gid, gid),
    ]
    for k, ((s, e), ph) in enumerate(zip(exons, phases), 1):
        out.append("{}\tPASA\texon\t{}\t{}\t.\t{}\t.\tID={}-T1.exon{};Parent={}-T1;".format(
            contig, s, e, strand, gid, k, gid))
        out.append("{}\tPASA\tCDS\t{}\t{}\t.\t{}\t{}\tID={}-T1.cds;Parent={}-T1;".format(
            contig, s, e, strand, ph, gid, gid))
    return out


def ascending_phases(exons):
    """Phases as the buggy PASA wrote them: counted from the lowest coordinate."""
    order = sorted(range(len(exons)), key=lambda i: exons[i][0])
    ph, acc = [0] * len(exons), 0
    for i in order:
        ph[i] = (3 - acc % 3) % 3
        acc += exons[i][1] - exons[i][0] + 1
    return ph


class CheckInternalStopRateTests(unittest.TestCase):
    def setUp(self):
        p = mock.patch.object(update.lib, "log", mock.MagicMock(), create=True)
        p.start()
        self.addCleanup(p.stop)

    def test_predict_like_rate_passes(self):
        # Absidia predict output: 6 of 13,344
        self.assertLess(update.check_internal_stop_rate(6, 13344), 0.01)

    def test_zero_mrna_passes(self):
        self.assertEqual(update.check_internal_stop_rate(0, 0), 0.0)

    def test_bug_rate_stops(self):
        # Absidia update with the PASA bug: 3,535 of 15,836
        with self.assertRaises(SystemExit):
            update.check_internal_stop_rate(3535, 15836)

    def test_limit_is_inclusive(self):
        self.assertEqual(update.check_internal_stop_rate(1, 100, limit=0.01), 0.01)
        with self.assertRaises(SystemExit):
            update.check_internal_stop_rate(2, 100, limit=0.01)


class GFF2tblPhaseTests(unittest.TestCase):
    """Run GFF2tblCombinedNEW on a small genome with one minus and one plus gene."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d = self.tmp.name
        self.genome = os.path.join(d, "genome.fa")
        self.trna = os.path.join(d, "trna.gff3")
        self.tbl = os.path.join(d, "genome.tbl")
        open(self.trna, "w").close()
        self.minus_seq, self.minus_exons = build_gene("ctg1", "-", 1)
        self.plus_seq, self.plus_exons = build_gene("ctg2", "+", 2)
        with open(self.genome, "w") as fh:
            fh.write(">ctg1\n{}\n>ctg2\n{}\n".format(self.minus_seq, self.plus_seq))
        p = mock.patch.object(update.lib, "log", mock.MagicMock(), create=True)
        p.start()
        self.addCleanup(p.stop)
        update.ExpressionValues = {}

    def run_gff(self, minus_phases):
        gff = os.path.join(self.tmp.name, "pasa_final.gff3")
        lines = ["##gff-version 3"]
        lines += gff_lines("ctg1", "FUN_000001", "-", self.minus_exons, minus_phases)
        lines += gff_lines("ctg2", "FUN_000002", "+", self.plus_exons, GOOD_PHASES)
        with open(gff, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        update.GFF2tblCombinedNEW(
            gff, self.genome, self.trna, "FUN_", 2, 6, "CFMR", "12345", self.tbl,
            alt_transcripts="0.10", transl_table=1)

    def test_bug_phases_differ_from_correct(self):
        # Guard against a test that cannot fail: the buggy phases for this
        # structure must not equal the correct ones.
        self.assertNotEqual(ascending_phases(self.minus_exons), GOOD_PHASES)
        self.assertEqual(ascending_phases(self.minus_exons)[0], 1)

    def test_correct_phases_write_tbl(self):
        self.run_gff(GOOD_PHASES)
        with open(self.tbl) as fh:
            tbl = fh.read()
        self.assertIn("FUN_000001", tbl)
        self.assertIn("FUN_000002", tbl)
        # complete models: codon_start 1 is the default and is not written
        self.assertNotIn("codon_start\t2", tbl)
        self.assertNotIn("codon_start\t3", tbl)

    def test_bug_phases_stop_update(self):
        with self.assertRaises(SystemExit):
            self.run_gff(ascending_phases(self.minus_exons))


class FallbackToPredictTests(GFF2tblPhaseTests):
    """A PASA locus with an internal stop is replaced by the predict model."""

    UTR5 = 30  # bp of 5' UTR PASA adds to the first exon

    def write(self, name, lines):
        path = os.path.join(self.tmp.name, name)
        with open(path, "w") as fh:
            fh.write("\n".join(["##gff-version 3"] + lines) + "\n")
        return path

    def run_fallback(self, minus_id, minus_phases):
        # predict models (input to PASA): CDS only, correct phases
        predict = gff_lines("ctg1", "FUN_000001", "-", self.minus_exons, GOOD_PHASES)
        predict += gff_lines("ctg2", "FUN_000002", "+", self.plus_exons, GOOD_PHASES)
        fallback = self.write("genome.gff3", predict)
        # PASA model for the minus gene: 5' UTR added (first exon 5' end moved
        # up on the minus strand), CDS unchanged, phases as given
        utr_exons = list(self.minus_exons)
        s, e = utr_exons[0]
        pasa = gff_lines("ctg1", minus_id, "-", self.minus_exons, minus_phases)
        pasa = [
            l.replace("\t{}\t".format(e), "\t{}\t".format(e + self.UTR5), 1)
            if ("\texon\t" in l and "exon1;" in l) or "\tgene\t" in l or "\tmRNA\t" in l
            else l
            for l in pasa
        ]
        pasa += gff_lines("ctg2", "FUN_000002", "+", self.plus_exons, GOOD_PHASES)
        evm = self.write("pasa_final.gff3", pasa)
        with mock.patch.object(update, "check_internal_stop_rate", lambda n, m: 0.0):
            update.GFF2tblCombinedNEW(
                evm, self.genome, self.trna, "FUN_", 2, 6, "CFMR", "12345", self.tbl,
                alt_transcripts="0.10", transl_table=1, fallback_gff=fallback)
        with open(self.tbl) as fh:
            return fh.read(), e

    def test_good_pasa_model_is_kept(self):
        tbl, e = self.run_fallback("FUN_000001", GOOD_PHASES)
        self.assertIn("FUN_000001", tbl)
        self.assertIn(str(e + self.UTR5), tbl)  # PASA 5' UTR end kept

    def test_bad_pasa_model_reverts_to_predict(self):
        tbl, e = self.run_fallback("FUN_000001", ascending_phases(self.minus_exons))
        self.assertIn("FUN_000001", tbl)        # locus kept
        self.assertNotIn(str(e + self.UTR5), tbl)  # PASA UTR gone: predict model
        self.assertNotIn("codon_start\t2", tbl)    # no shifted frame
        self.assertIn("FUN_000002", tbl)        # other locus untouched

    def test_bad_new_pasa_gene_is_dropped(self):
        tbl, _ = self.run_fallback("novel_gene_1", ascending_phases(self.minus_exons))
        # novel gene had no original model: dropped, and it does not get a new
        # locus tag; the predict model it overlaps was not in the PASA output.
        self.assertNotIn("FUN_000003", tbl)
        self.assertIn("FUN_000002", tbl)


if __name__ == "__main__":
    unittest.main()
