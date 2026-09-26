"""Tests for PASA training-model selection (code review R3, R5, F8, F12).

R5/F8: getBestModel must keep ONE model per same-strand locus, ranked by
completeness, then CDS exon count, then CDS length, and only then TPM. The
old code ranked by TPM first and used a one-sided overlap, so a short,
high-TPM fragment inside a full-length model survived next to it.

R3/F12: selectTrainingModels must train only on complete ORFs, and its
overlap removal must be transitive (cluster-based).
"""
import logging
import os
import random
import shutil
import tempfile
import unittest

from funannotate import library as lib
from funannotate import train

# funannotate subcommands assign lib.log at startup; do the same here.
if not hasattr(lib, "log"):
    lib.log = logging.getLogger("funannotate-test")

SENSE = [
    a + b + c
    for a in "ACGT"
    for b in "ACGT"
    for c in "ACGT"
    if a + b + c not in ("TAA", "TAG", "TGA", "ATG")
]


def _rand(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _orf(ncodons, rng):
    return "ATG" + "".join(rng.choice(SENSE) for _ in range(ncodons)) + "TAA"


class IsCompleteCdsTests(unittest.TestCase):
    def test_complete(self):
        self.assertTrue(lib.is_complete_cds("ATGAAACCCTAA"))

    def test_lowercase_is_accepted(self):
        self.assertTrue(lib.is_complete_cds("atgaaaccctag"))

    def test_missing_start(self):
        self.assertFalse(lib.is_complete_cds("CTGAAACCCTAA"))

    def test_missing_stop(self):
        self.assertFalse(lib.is_complete_cds("ATGAAACCCAAA"))

    def test_not_multiple_of_three(self):
        self.assertFalse(lib.is_complete_cds("ATGAAACCTAA"))

    def test_internal_stop(self):
        self.assertFalse(lib.is_complete_cds("ATGTAACCCTAA"))

    def test_empty(self):
        self.assertFalse(lib.is_complete_cds(""))


class ClusterOverlappingTests(unittest.TestCase):
    def _ids(self, clusters):
        return sorted(sorted(c) for c in clusters)

    def test_containment_clusters_regardless_of_order(self):
        # short model inside a long one: 100% of the short one is covered,
        # 10% of the long one. Symmetric overlap (vs shorter) must cluster them.
        models = [("long", "c1", "+", 1, 1000), ("short", "c1", "+", 400, 500)]
        for m in (models, models[::-1]):
            self.assertEqual(
                self._ids(lib.cluster_overlapping(m, min_frac=0.3, strand_aware=True)),
                [["long", "short"]],
            )

    def test_overlap_is_transitive(self):
        # A-B overlap and B-C overlap, A-C do not: one cluster, not two.
        models = [
            ("A", "c1", "+", 1, 100),
            ("B", "c1", "+", 90, 200),
            ("C", "c1", "+", 190, 300),
        ]
        self.assertEqual(
            self._ids(lib.cluster_overlapping(models, min_frac=0.0, strand_aware=False)),
            [["A", "B", "C"]],
        )

    def test_strand_aware_keeps_opposite_strands_apart(self):
        models = [("p", "c1", "+", 1, 100), ("m", "c1", "-", 1, 100)]
        self.assertEqual(
            self._ids(lib.cluster_overlapping(models, min_frac=0.0, strand_aware=True)),
            [["m"], ["p"]],
        )
        self.assertEqual(
            self._ids(lib.cluster_overlapping(models, min_frac=0.0, strand_aware=False)),
            [["m", "p"]],
        )

    def test_below_threshold_stays_apart(self):
        # 10 bp overlap of 100 bp models = 10% < 30%
        models = [("a", "c1", "+", 1, 100), ("b", "c1", "+", 91, 190)]
        self.assertEqual(
            self._ids(lib.cluster_overlapping(models, min_frac=0.3, strand_aware=True)),
            [["a"], ["b"]],
        )

    def test_touching_is_not_overlap(self):
        # 1-based closed intervals: 1-100 and 101-200 share no base
        models = [("a", "c1", "+", 1, 100), ("b", "c1", "+", 101, 200)]
        self.assertEqual(
            self._ids(lib.cluster_overlapping(models, min_frac=0.0, strand_aware=False)),
            [["a"], ["b"]],
        )

    def test_different_contigs_never_cluster(self):
        models = [("a", "c1", "+", 1, 100), ("b", "c2", "+", 1, 100)]
        self.assertEqual(
            self._ids(lib.cluster_overlapping(models, min_frac=0.0, strand_aware=False)),
            [["a"], ["b"]],
        )


class LocusFixture(unittest.TestCase):
    """One contig with:

    gLong : complete two-exon ORF (+), TPM 5
    gFrag : 3' fragment of gLong's second exon (+), no start codon, TPM 50
    gAlt  : complete single-exon ORF spanning only gLong's second exon, i.e. a
            complete but structurally poorer model at the same locus, TPM 80
    gFar  : complete single-exon ORF at a separate locus (+), TPM 1
    gPart : 5'-partial two-exon model alone at its locus (+), TPM 20
    gRev  : complete two-exon ORF on the MINUS strand, alone at its locus,
            TPM 3. gff2dict keeps cds_transcript in genome orientation for
            minus-strand genes, so a sequence-based check wrongly calls these
            incomplete; completeness must not depend on strand.
    """

    def setUp(self):
        rng = random.Random(3)
        self.tmp = tempfile.mkdtemp()
        orf = _orf(200, rng)  # 606 bp
        ex1, ex2 = orf[:180], orf[180:]
        intron = "GT" + _rand(96, rng) + "AG"
        far = _orf(150, rng)
        # gPart: two-exon model alone at its own locus, start codon destroyed
        # -> 5'-partial. The old selection kept such models (no overlap to
        # remove it); R3 must drop it.
        porf = "CTG" + _orf(160, rng)[3:]
        pex1, pex2 = porf[:150], porf[150:]
        pintron = "GT" + _rand(90, rng) + "AG"
        rorf = _orf(170, rng)
        rex1, rex2 = rorf[:210], rorf[210:]
        rintron = "GT" + _rand(84, rng) + "AG"
        rev_block = _revcomp(rex1 + rintron + rex2)
        parts = [_rand(300, rng), ex1, intron, ex2, _rand(400, rng), far, _rand(300, rng),
                 pex1, pintron, pex2, _rand(300, rng), rev_block, _rand(300, rng)]
        seq = "".join(parts)
        # make gAlt a real complete ORF: put an in-frame ATG at the start of
        # its span inside ex2 (codon-aligned within the ORF)
        off = sum(len(p) for p in parts[:3])
        alt_start = off + 1 + 60  # 1-based, codon aligned in ex2 (ex2 offset 180 is codon 60)
        seq = seq[: alt_start - 1] + "ATG" + seq[alt_start + 2 :]
        self.genome = seq
        self.fasta = os.path.join(self.tmp, "genome.fa")
        with open(self.fasta, "w") as f:
            f.write(">ctg1\n" + seq + "\n")

        def span(i):
            s = sum(len(p) for p in parts[:i]) + 1
            return s, s + len(parts[i]) - 1

        e1s, e1e = span(1)
        e2s, e2e = span(3)
        fs, fe = span(5)
        # fragment: last 150 bp of ex2 starting mid-codon-aligned, no ATG
        frag_s = e2e - 149
        rows = []

        def gene(gid, exons, strand="+"):
            lo, hi = exons[0][0], exons[-1][1]
            rows.append(("gene", lo, hi, strand, "ID={}".format(gid)))
            rows.append(("mRNA", lo, hi, strand, "ID={}.t1;Parent={}".format(gid, gid)))
            for n, (s, e) in enumerate(exons, 1):
                rows.append(("exon", s, e, strand, "ID={}.t1.exon{};Parent={}.t1".format(gid, n, gid)))
                rows.append(("CDS", s, e, strand, "ID={}.t1.cds;Parent={}.t1".format(gid, gid)))

        gene("gLong", [(e1s, e1e), (e2s, e2e)])
        gene("gFrag", [(frag_s, e2e)])
        gene("gAlt", [(alt_start, e2e)])
        gene("gFar", [(fs, fe)])
        gene("gPart", [span(7), span(9)])
        rs, re_ = span(11)
        # genome coords of the two minus-strand exons inside rev_block
        gene("gRev", [(rs, rs + len(rex2) - 1), (re_ - len(rex1) + 1, re_)], strand="-")
        self.gff3 = os.path.join(self.tmp, "pasa.gff3")
        with open(self.gff3, "w") as f:
            f.write("##gff-version 3\n")
            for feat, s, e, strand, attr in rows:
                f.write("\t".join(["ctg1", "PASA", feat, str(s), str(e), ".", strand, ".", attr]) + "\n")
        self.tpm = os.path.join(self.tmp, "kallisto.tsv")
        with open(self.tpm, "w") as f:
            f.write("#mRNA-ID\tgene-ID\tLocation\tTPM\n")
            for gid, t in (("gLong", 5), ("gFrag", 50), ("gAlt", 80), ("gFar", 1), ("gPart", 20), ("gRev", 3)):
                f.write("{0}.t1\t{0}\tctg1\t{1}\n".format(gid, t))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _gene_ids(self, gff3):
        ids = set()
        with open(gff3) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) > 8 and cols[2] == "gene":
                    ids.add(cols[8].split("ID=")[1].split(";")[0])
        return ids


class GetBestModelTests(LocusFixture):
    def test_fixture_models_have_expected_completeness(self):
        _, genes = lib.gff2interlap(self.gff3, self.fasta)
        self.assertTrue(lib.is_complete_model(genes["gLong"]))
        self.assertFalse(lib.is_complete_model(genes["gFrag"]))
        self.assertTrue(lib.is_complete_model(genes["gAlt"]))
        self.assertTrue(lib.is_complete_model(genes["gFar"]))
        self.assertFalse(lib.is_complete_model(genes["gPart"]))
        self.assertEqual(len(genes["gPart"]["CDS"][0]), 2)
        self.assertTrue(lib.is_complete_model(genes["gRev"]))
        self.assertEqual(len(genes["gRev"]["CDS"][0]), 2)

    def test_keeps_one_structure_ranked_model_per_locus(self):
        out = os.path.join(self.tmp, "best.gff3")
        train.getBestModel(self.gff3, self.fasta, self.tpm, out, pasa_alignment_overlap=30)
        # gLong beats the higher-TPM fragment and the higher-TPM single-exon
        # model at its locus. gFar and gPart are alone at their loci and are
        # kept: getBestModel ranks within a locus, R3 filtering happens later.
        self.assertEqual(self._gene_ids(out), {"gLong", "gFar", "gPart", "gRev"})


class GuardedRankingTests(unittest.TestCase):
    """Review finding 1 (DECISIONS D41/D55/D68/D69). complete_min_frac=0 (default)
    is complete-first; complete_min_frac=0.8 is the guarded rule (a complete model
    wins only if its CDS is >= 80% of the locus's longest CDS). Guarded won on
    RefSeq evidence metrics but was prediction-neutral, so it is opt-in."""

    def test_default_is_complete_first(self):
        # D69: the default reproduces complete-first (what every predict arm used);
        # guarded ranking was prediction-neutral on 3 RefSeq genomes.
        feats = {"long_partial": (False, 3, 1500), "short_complete": (True, 1, 300)}
        self.assertEqual(train.pick_locus_model(feats, {}), "short_complete")

    def test_guarded_long_partial_beats_short_complete(self):
        feats = {"long_partial": (False, 3, 1500), "short_complete": (True, 1, 300)}
        self.assertEqual(train.pick_locus_model(feats, {}, complete_min_frac=0.8), "long_partial")

    def test_guarded_near_full_length_complete_beats_partial(self):
        feats = {"partial": (False, 3, 1500), "complete": (True, 3, 1290)}   # 86%
        self.assertEqual(train.pick_locus_model(feats, {}, complete_min_frac=0.8), "complete")

    def test_tpm_breaks_ties_last(self):
        feats = {"a": (True, 2, 900), "b": (True, 2, 900)}
        self.assertEqual(train.pick_locus_model(feats, {"a": 1.0, "b": 5.0}), "b")


class SelectTrainingModelsTests(LocusFixture):
    def test_trains_only_on_complete_non_overlapping_models(self):
        out = os.path.join(self.tmp, "final_training_models.gff3")
        empty_gtf = os.path.join(self.tmp, "keepers.gtf")
        open(empty_gtf, "w").close()
        n = lib.selectTrainingModels(
            self.gff3, self.fasta, empty_gtf, out, self.tmp, min_models=200
        )
        _, kept = lib.gff2interlap(out, self.fasta)
        spans = sorted((v["location"][0], v["location"][1]) for v in kept.values())
        # gFrag and gPart are incomplete -> dropped (gPart has no overlap, so
        # only a completeness filter removes it); gLong and gAlt overlap ->
        # only the two-exon gLong survives; gFar and the minus-strand gRev are
        # separate complete models -> kept.
        self.assertEqual(n, 3)
        self.assertEqual(len(spans), 3)
        for v in kept.values():
            self.assertTrue(lib.is_complete_model(v))
        self.assertEqual(sorted(v["strand"] for v in kept.values()), ["+", "+", "-"])


def _write_align(path, rows):
    """rows: (contig, start, end, strand, match_id, target, tstart, tend)"""
    with open(path, "w") as f:
        for c, s, e, st, mid, t, ts, te in rows:
            f.write("\t".join([c, "exonerate", "nucleotide_to_protein_match", str(s), str(e),
                               "100.00", st, ".", "ID={};Target={} {} {}".format(mid, t, ts, te)]) + "\n")


class ProteinSingleExonShareTests(unittest.TestCase):
    """Single-exon share estimated from exonerate protein2genome alignments (R6 option b)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.fa = os.path.join(self.tmp, "prot.fa")
        with open(self.fa, "w") as f:
            for name, n in (("P1", 100), ("P2", 200), ("P3", 150), ("P4", 300)):
                f.write(">{}\n{}\n".format(name, "M" * n))
        self.gff = os.path.join(self.tmp, "aln.gff3")
        _write_align(self.gff, [
            ("c1", 1000, 1300, "+", "m1", "P1", 1, 100),        # single segment, full length
            ("c1", 5000, 5600, "-", "m2", "P2", 1, 200),        # single segment, full length
            ("c1", 9000, 9200, "+", "m3", "P3", 1, 70),         # 2 exons, intron 100 bp
            ("c1", 9301, 9540, "+", "m3", "P3", 71, 150),
            ("c1", 20000, 20150, "+", "m4", "P4", 1, 50),       # partial (17% of P4) -> ignored
        ])

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_share_counts_only_near_full_length_alignments(self):
        share, n = lib.protein_single_exon_share(self.gff, self.fa, min_cov=0.9, min_intron=20)
        self.assertEqual(n, 3)
        self.assertAlmostEqual(share, 2.0 / 3.0)

    def test_small_gap_is_not_an_intron(self):
        # a 5-bp gap between segments (frameshift/indel), not an intron
        _write_align(self.gff, [
            ("c1", 1000, 1150, "+", "m1", "P1", 1, 50),
            ("c1", 1156, 1306, "+", "m1", "P1", 51, 100),
        ])
        share, n = lib.protein_single_exon_share(self.gff, self.fa, min_cov=0.9, min_intron=20)
        self.assertEqual((share, n), (1.0, 1))

    def test_missing_file_returns_none(self):
        self.assertEqual(lib.protein_single_exon_share(os.path.join(self.tmp, "x"), self.fa), (None, 0))


class SelectionDecisionLogTests(LocusFixture):
    def test_selection_records_each_filter_step(self):
        path = os.path.join(self.tmp, "training_decisions.tsv")
        lib.set_training_decision_log(path, command="predict")
        try:
            gtf = os.path.join(self.tmp, "keepers.gtf")
            open(gtf, "w").close()
            n = lib.selectTrainingModels(self.gff3, self.fasta, gtf,
                                         os.path.join(self.tmp, "t.gff3"), self.tmp, min_models=200)
        finally:
            lib.set_training_decision_log(None)
        with open(path) as f:
            stages = [line.split("\t")[1] for line in f.readlines()[1:]]
        for st in ("select_complete_orf", "select_keeper_filter", "select_multi_cds",
                   "select_redundancy", "select_overlap", "select_final"):
            self.assertIn(st, stages)
        with open(path) as f:
            final = [l for l in f if "\tselect_final\t" in l][0].split("\t")
        self.assertEqual(final[3], str(n))  # columns: command, stage, decision, value, ...


class AllPartialInputTests(LocusFixture):
    def test_returns_zero_instead_of_exiting_when_no_complete_models(self):
        # keep only the incomplete models gFrag and gPart -> R3 leaves nothing.
        # Must return 0 so predict can fall back to BUSCO (review finding 2:
        # an empty proteins FASTA made diamond makedb fail -> sys.exit(1)).
        keep = ("ID=gFrag", "Parent=gFrag", "ID=gPart", "Parent=gPart")
        part = os.path.join(self.tmp, "partial.gff3")
        with open(self.gff3) as f, open(part, "w") as o:
            for line in f:
                if line.startswith("#") or any(k in line for k in keep):
                    o.write(line)
        gtf = os.path.join(self.tmp, "keepers.gtf")
        open(gtf, "w").close()
        out = os.path.join(self.tmp, "train.gff3")
        self.assertEqual(lib.selectTrainingModels(part, self.fasta, gtf, out, self.tmp), 0)


class ModelSingleExonShareTests(unittest.TestCase):
    """Single-exon share from ab-initio gene models (GeneMark-ES), R6 option (b).
    On the 3 RefSeq test genomes GeneMark-ES was within +1.3..+6.5 points of the
    RefSeq single-CDS share, the protein-alignment estimate -2.5..-13.6 (D44)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_gtf_counts_cds_per_gene_id(self):
        p = os.path.join(self.tmp, "gm.gtf")
        with open(p, "w") as f:
            for gid, n in (("g1", 1), ("g2", 3), ("g3", 1), ("g4", 2)):
                for i in range(n):
                    f.write("c1\tGeneMark.hmm\tCDS\t{0}\t{1}\t.\t+\t0\tgene_id \"{2}\"; transcript_id \"{2}.t\";\n".format(100 * i + 1, 100 * i + 50, gid))
        self.assertEqual(lib.model_single_exon_share(p), (0.5, 4))

    def test_gff3_counts_cds_per_parent(self):
        p = os.path.join(self.tmp, "gm.gff3")
        with open(p, "w") as f:
            f.write("c1\tGeneMark\tCDS\t1\t90\t.\t+\t0\tID=a.cds;Parent=a.t1\n")
            f.write("c1\tGeneMark\tCDS\t200\t290\t.\t+\t0\tID=b.cds;Parent=b.t1\n")
            f.write("c1\tGeneMark\tCDS\t400\t490\t.\t+\t0\tID=b.cds;Parent=b.t1\n")
        self.assertEqual(lib.model_single_exon_share(p), (0.5, 2))

    def test_missing_file(self):
        self.assertEqual(lib.model_single_exon_share(os.path.join(self.tmp, "x")), (None, 0))


class SingleExonTrainingTests(LocusFixture):
    """selectTrainingModels with R6 option (b). min_models=1 makes the multi-CDS
    requirement active (gLong and gRev are complete multi-exon models), which is
    the production situation that excludes every single-exon gene today."""

    def _run(self, **kw):
        out = os.path.join(self.tmp, "train.gff3")
        gtf = os.path.join(self.tmp, "keepers.gtf")
        open(gtf, "w").close()
        n = lib.selectTrainingModels(self.gff3, self.fasta, gtf, out, self.tmp, min_models=1, **kw)
        _, kept = lib.gff2interlap(out, self.fasta)
        return n, kept

    def _far_align(self, strand="+", shift=0):
        _, genes = lib.gff2interlap(self.gff3, self.fasta)
        s, e = genes["gFar"]["location"]
        p = os.path.join(self.tmp, "aln.gff3")
        _write_align(p, [("ctg1", s + shift, e + shift, strand, "m1", "PX", 1, 150)])
        return p

    def test_default_excludes_single_exon_genes(self):
        n, kept = self._run()
        self.assertEqual(n, 2)
        self.assertTrue(all(len(v["CDS"][0]) > 1 for v in kept.values()))

    def test_supported_single_exon_gene_is_admitted_within_cap(self):
        n, kept = self._run(single_exon_align=self._far_align(), single_exon_share=0.5)
        self.assertEqual(n, 3)
        self.assertEqual(sum(len(v["CDS"][0]) == 1 for v in kept.values()), 1)

    def test_cap_of_zero_admits_none(self):
        # share 0.2 -> cap floor(0.25 * 2 multi-exon models) = 0
        n, _ = self._run(single_exon_align=self._far_align(), single_exon_share=0.2)
        self.assertEqual(n, 2)

    def test_unsupported_single_exon_gene_is_not_admitted(self):
        n, _ = self._run(single_exon_align=self._far_align(shift=100000), single_exon_share=0.5)
        self.assertEqual(n, 2)

    def test_opposite_strand_alignment_is_not_support(self):
        n, _ = self._run(single_exon_align=self._far_align(strand="-"), single_exon_share=0.5)
        self.assertEqual(n, 2)


if __name__ == "__main__":
    unittest.main()
