"""Tests for minimap2 SAM -> GFF3 / hints conversion (bam2gff3, bam2ExonsHints).

Each SAM record is hand-built so that its CIGAR and cs tag agree, and the
expected genome/transcript coordinates are worked out by hand in the comment.
The original parser advanced the genome position only on cs ':' runs, so
substitutions, deletions and soft clips shifted every later coordinate, and it
judged intron motifs against the SAM flag instead of either orientation.
"""
import os
import tempfile
import unittest
from unittest import mock

import funannotate.library as lib


def sam(qname, flag, pos, cigar, seqlen, cs, nm=0, rname="chr1"):
    seq = "A" * seqlen
    tags = ["NM:i:{}".format(nm), "cs:Z:{}".format(cs)]
    cols = [qname, str(flag), rname, str(pos), "60", cigar, "*", "0", "0", seq, "*"]
    return "\t".join(cols + tags) + "\n"


def gff_rows(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            target = c[8].split("Target=")[1].split()
            rows.append(
                {
                    "chr": c[0],
                    "feature": c[2],
                    "start": int(c[3]),
                    "end": int(c[4]),
                    "strand": c[6],
                    "id": c[8].split("ID=")[1].split(";")[0],
                    "tstart": int(target[1]),
                    "tend": int(target[2]),
                }
            )
    return rows


def run_bam2gff3(lines):
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "out.gff3")
        with mock.patch.object(lib, "execute", return_value=iter(lines)):
            count = lib.bam2gff3("dummy.bam", out)
        return count, gff_rows(out)


def run_hints(lines):
    with tempfile.TemporaryDirectory() as d:
        gff = os.path.join(d, "out.gff3")
        hints = os.path.join(d, "hints.gff")
        with mock.patch.object(lib, "execute", return_value=iter(lines)):
            count = lib.bam2ExonsHints("dummy.bam", gff, hints)
        with open(hints) as fh:
            hint_rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
        return count, gff_rows(gff), hint_rows


def coords(rows):
    return [(r["start"], r["end"], r["tstart"], r["tend"]) for r in rows]


class Bam2Gff3CoordinateTests(unittest.TestCase):
    def test_substitution_before_intron_does_not_shift_exons(self):
        # 10M5N10M at 101: exon1 101-110, intron 111-115, exon2 116-125
        count, rows = run_bam2gff3(
            [sam("t1", 0, 101, "10M5N10M", 20, ":4*ac:5~gt5ag:10", nm=1)]
        )
        self.assertEqual(count, 1)
        self.assertEqual(coords(rows), [(101, 110, 1, 10), (116, 125, 11, 20)])

    def test_deletion_advances_genome_position(self):
        # 5M3D5M5N10M at 101: exon1 101-113 (5+3+5 ref bases), exon2 119-128
        count, rows = run_bam2gff3(
            [sam("t2", 0, 101, "5M3D5M5N10M", 20, ":5-acg:5~gt5ag:10", nm=3)]
        )
        self.assertEqual(coords(rows), [(101, 113, 1, 10), (119, 128, 11, 20)])

    def test_insertion_advances_transcript_position_only(self):
        # 5M2I5M5N10M at 101: exon1 101-110 / query 1-12, exon2 116-125 / 13-22
        count, rows = run_bam2gff3(
            [sam("t3", 0, 101, "5M2I5M5N10M", 22, ":5+tt:5~gt5ag:10", nm=2)]
        )
        self.assertEqual(coords(rows), [(101, 110, 1, 12), (116, 125, 13, 22)])

    def test_soft_clips_are_not_counted_as_aligned(self):
        # 5S10M5N10M5S, read length 30: aligned query is 6-25, not 1-30
        count, rows = run_bam2gff3(
            [sam("t4", 0, 101, "5S10M5N10M5S", 30, ":10~gt5ag:10")]
        )
        self.assertEqual(coords(rows), [(101, 110, 6, 15), (116, 125, 16, 25)])

    def test_reverse_alignment_targets_are_in_transcript_orientation(self):
        # flag 16, 3S10M5N10M, read length 23. Read-orientation query spans
        # exon1 4-13, exon2 14-23; flipped to transcript orientation:
        # exon1 11-20, exon2 1-10. Column 7 is the SAM strand for PASA.
        count, rows = run_bam2gff3(
            [sam("t5", 16, 101, "3S10M5N10M", 23, ":10~ct5ac:10")]
        )
        self.assertEqual(coords(rows), [(101, 110, 11, 20), (116, 125, 1, 10)])
        self.assertEqual({r["strand"] for r in rows}, {"-"})


class Bam2Gff3SpliceMotifTests(unittest.TestCase):
    def test_flag16_with_forward_motif_is_kept(self):
        # antisense-oriented contig of a plus-strand gene: flag 16, gt..ag
        count, rows = run_bam2gff3([sam("t6", 16, 101, "10M5N10M", 20, ":10~gt5ag:10")])
        self.assertEqual(count, 1)

    def test_flag0_with_reverse_motif_is_kept(self):
        # sense-oriented contig of a minus-strand gene: flag 0, ct..ac
        count, rows = run_bam2gff3([sam("t7", 0, 101, "10M5N10M", 20, ":10~ct5ac:10")])
        self.assertEqual(count, 1)

    def test_gc_ag_intron_is_kept(self):
        count, rows = run_bam2gff3([sam("t8", 0, 101, "10M5N10M", 20, ":10~gc5ag:10")])
        self.assertEqual(count, 1)

    def test_any_noncanonical_intron_rejects_alignment(self):
        # first intron non-canonical, last intron canonical: old code kept it
        count, rows = run_bam2gff3(
            [sam("t9", 0, 101, "10M5N10M5N10M", 30, ":10~aa5tt:10~gt5ag:10")]
        )
        self.assertEqual(count, 0)
        self.assertEqual(rows, [])

    def test_mixed_orientation_introns_reject_alignment(self):
        count, rows = run_bam2gff3(
            [sam("t10", 0, 101, "10M5N10M5N10M", 30, ":10~gt5ag:10~ct5ac:10")]
        )
        self.assertEqual(count, 0)


class Bam2Gff3FilterTests(unittest.TestCase):
    def test_secondary_and_supplementary_records_are_skipped(self):
        lines = [
            sam("s1", 256, 101, "10M", 10, ":10"),
            sam("s2", 2048, 101, "10M", 10, ":10"),
            sam("s3", 272, 101, "10M", 10, ":10"),
        ]
        count, rows = run_bam2gff3(lines)
        self.assertEqual(count, 0)

    def test_low_identity_alignment_is_skipped(self):
        # 3 substitutions in 10 aligned bases = 70% identity (< 80%)
        count, rows = run_bam2gff3([sam("s4", 0, 101, "10M", 10, "*ac*ac*ac:7", nm=3)])
        self.assertEqual(count, 0)

    def test_single_exon_alignment(self):
        count, rows = run_bam2gff3([sam("s5", 0, 101, "2S20M", 22, ":20")])
        self.assertEqual(coords(rows), [(101, 120, 3, 22)])
        self.assertEqual(rows[0]["id"], "s5")


class Bam2Gff3ReviewFollowupTests(unittest.TestCase):
    """Cases added after the Fable 5.1 pre-merge review (DECISIONS D16)."""

    def test_default_strand_is_alignment_strand_for_pasa(self):
        count, rows = run_bam2gff3([sam("r1", 16, 101, "10M5N10M", 20, ":10~gt5ag:10")])
        self.assertEqual({r["strand"] for r in rows}, {"-"})

    def test_splice_strand_option_for_evidence_and_hints(self):
        # predict's harmonize_transcripts feeds EVM/Augustus: strand must be the
        # transcribed strand implied by the intron motif (gt..ag -> '+')
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "out.gff3")
            line = sam("r2", 16, 101, "10M5N10M", 20, ":10~gt5ag:10")
            with mock.patch.object(lib, "execute", return_value=iter([line])):
                lib.bam2gff3("dummy.bam", out, strand="splice")
            rows = gff_rows(out)
        self.assertEqual({r["strand"] for r in rows}, {"+"})

    def test_identity_counts_gap_events_not_gap_bases(self):
        # 100 matched bases and one 6-bp deletion: gap-compressed identity is
        # 100/101 = 99.01%, not 100/106 (as with blat-style per-id, a single
        # indel event is one difference)
        count, rows = run_bam2gff3(
            [sam("r3", 0, 101, "50M6D50M", 100, ":50-aaaaaa:50", nm=6)]
        )
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "o.gff3")
            with mock.patch.object(
                lib, "execute",
                return_value=iter([sam("r3", 0, 101, "50M6D50M", 100, ":50-aaaaaa:50", nm=6)]),
            ):
                lib.bam2gff3("dummy.bam", out)
            score = float(open(out).read().splitlines()[1].split("\t")[5])
        self.assertAlmostEqual(score, 99.01, places=2)

    def test_duplicate_and_qc_fail_flags_are_skipped(self):
        lines = [
            sam("f1", 1024, 101, "10M", 10, ":10"),
            sam("f2", 512, 101, "10M", 10, ":10"),
        ]
        count, rows = run_bam2gff3(lines)
        self.assertEqual(count, 0)

    def test_insertion_right_after_intron_keeps_target_contiguous(self):
        # 10M5N2I8M: the 2 inserted bases belong to exon 2, so transcript
        # coordinates stay contiguous (PASA rejects gaps as "Incontiguous
        # alignment"; seen in 137/137 such failures on A. nidulans, D32)
        count, rows = run_bam2gff3(
            [sam("i1", 0, 101, "10M5N2I8M", 20, ":10~gt5ag+tt:8", nm=2)]
        )
        self.assertEqual(coords(rows), [(101, 110, 1, 10), (116, 123, 11, 20)])

    def test_cigar_intron_without_cs_intron_is_rejected(self):
        # CIGAR says spliced, cs does not: inconsistent record, skip it
        count, rows = run_bam2gff3([sam("c1", 0, 101, "10M5N10M", 20, ":20")])
        self.assertEqual(count, 0)


class Bam2ExonsHintsTests(unittest.TestCase):
    def test_intron_hint_coordinates_after_substitution(self):
        count, rows, hints = run_hints(
            [sam("h1", 0, 101, "10M5N10M", 20, ":4*ac:5~gt5ag:10", nm=1)]
        )
        self.assertEqual(count, 1)
        self.assertEqual(coords(rows), [(101, 110, 1, 10), (116, 125, 11, 20)])
        introns = [h for h in hints if h[2] == "intron"]
        self.assertEqual([(int(h[3]), int(h[4])) for h in introns], [(111, 115)])

    def test_hint_strand_follows_intron_motif(self):
        # flag 16 with gt..ag: the gene is on the plus strand
        count, rows, hints = run_hints([sam("h2", 16, 101, "10M5N10M", 20, ":10~gt5ag:10")])
        self.assertEqual(count, 1)
        self.assertEqual({h[6] for h in hints}, {"+"})
        self.assertEqual({r["strand"] for r in rows}, {"+"})

    def test_single_exon_hint_uses_alignment_strand(self):
        count, rows, hints = run_hints([sam("h3", 16, 101, "20M", 20, ":20")])
        self.assertEqual(count, 1)
        self.assertEqual({h[6] for h in hints}, {"-"})

    def test_hint_ids_count_every_sam_record(self):
        # existing behaviour: IDs are minimap2_<record number>, counting skipped records
        lines = [
            sam("x1", 256, 101, "10M", 10, ":10"),
            sam("x2", 0, 201, "20M", 20, ":20"),
        ]
        count, rows, hints = run_hints(lines)
        self.assertEqual({r["id"] for r in rows}, {"minimap2_2"})

    def test_terminal_and_internal_exon_hint_types(self):
        count, rows, hints = run_hints(
            [sam("h4", 0, 101, "10M5N10M5N10M", 30, ":10~gt5ag:10~gt5ag:10")]
        )
        types = [h[2] for h in hints if h[2] != "intron"]
        self.assertEqual(types, ["ep", "exon", "ep"])


if __name__ == "__main__":
    unittest.main()
