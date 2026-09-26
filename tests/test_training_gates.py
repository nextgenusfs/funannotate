"""Tests for the evidence-quality gates used by train and predict.

Two gates guard against RNA-seq evidence that does not belong to the genome:

* train: sample reads and map them to the genome; too few mapped reads means
  the RNA-seq is from another organism (host tissue, wrong species, stale
  file) and PASA should not run.
* predict: count complete-ORF PASA training models; too few means Augustus
  and SNAP would train on fragments, so training falls back to BUSCO.
"""
import gzip
import os
import random
import shutil
import tempfile
import unittest

from funannotate import library as lib

RNG = random.Random(42)


def _randseq(n, rng=RNG):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _orf(ncodons, rng=RNG):
    # ATG + sense codons (no stops) + TAA
    sense = [
        a + b + c
        for a in "ACGT"
        for b in "ACGT"
        for c in "ACGT"
        if a + b + c not in ("TAA", "TAG", "TGA")
    ]
    return "ATG" + "".join(rng.choice(sense) for _ in range(ncodons)) + "TAA"


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class TrainingModelFixture(unittest.TestCase):
    """Genome with four PASA-style models on one contig:

    g1: complete single-exon ORF, + strand
    g2: complete two-exon ORF with a GT..AG intron, + strand
    g3: complete single-exon ORF, - strand
    g4: partial (no start codon), + strand
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        orf1 = _orf(100)
        orf2 = _orf(120)
        ex2a, ex2b = orf2[:150], orf2[150:]
        intron = "GT" + _randseq(80) + "AG"
        orf3 = _orf(90)
        orf4 = "CCC" + _orf(80)[3:]  # start codon replaced -> partial
        parts = [
            _randseq(200),
            orf1,
            _randseq(200),
            ex2a,
            intron,
            ex2b,
            _randseq(200),
            _revcomp(orf3),
            _randseq(200),
            orf4,
            _randseq(200),
        ]
        genome = "".join(parts)
        self.fasta = os.path.join(self.tmp, "genome.fa")
        with open(self.fasta, "w") as f:
            f.write(">ctg1\n" + genome + "\n")

        def span(i):
            start = sum(len(p) for p in parts[:i]) + 1
            return start, start + len(parts[i]) - 1

        s1, e1 = span(1)
        a_s, a_e = span(3)
        b_s, b_e = span(5)
        s3, e3 = span(7)
        s4, e4 = span(9)
        rows = []

        def gene(gid, strand, exons):
            lo, hi = exons[0][0], exons[-1][1]
            rows.append(("ctg1", "PASA", "gene", lo, hi, strand, "ID={}".format(gid)))
            rows.append(
                ("ctg1", "PASA", "mRNA", lo, hi, strand,
                 "ID={}.t1;Parent={}".format(gid, gid))
            )
            for n, (s, e) in enumerate(exons, 1):
                rows.append(
                    ("ctg1", "PASA", "exon", s, e, strand,
                     "ID={}.t1.exon{};Parent={}.t1".format(gid, n, gid))
                )
                rows.append(
                    ("ctg1", "PASA", "CDS", s, e, strand,
                     "ID={}.t1.cds;Parent={}.t1".format(gid, gid))
                )

        gene("g1", "+", [(s1, e1)])
        gene("g2", "+", [(a_s, a_e), (b_s, b_e)])
        gene("g3", "-", [(s3, e3)])
        gene("g4", "+", [(s4, e4)])
        self.gff3 = os.path.join(self.tmp, "pasa.gff3")
        with open(self.gff3, "w") as f:
            f.write("##gff-version 3\n")
            for r in rows:
                f.write("\t".join(map(str, r[:5])) + "\t.\t" + r[5] + "\t.\t" + r[6] + "\n")

    def tearDown(self):
        shutil.rmtree(self.tmp)


class CountCompleteOrfModelsTests(TrainingModelFixture):
    def test_counts_complete_and_partial_models(self):
        counts = lib.count_complete_orf_models(self.gff3, self.fasta)
        self.assertEqual(counts["total"], 4)
        self.assertEqual(counts["complete"], 3)
        self.assertEqual(counts["no_start"], 1)
        self.assertEqual(counts["no_stop"], 0)
        self.assertEqual(counts["not_mult3"], 0)

    def test_missing_file_counts_zero(self):
        counts = lib.count_complete_orf_models(
            os.path.join(self.tmp, "absent.gff3"), self.fasta
        )
        self.assertEqual(counts["total"], 0)
        self.assertEqual(counts["complete"], 0)


class RunPasaTrainingGateTests(TrainingModelFixture):
    def test_writes_report_and_returns_verdict(self):
        report = os.path.join(self.tmp, "predict_training_gate.tsv")
        passed = lib.run_pasa_training_gate(self.gff3, self.fasta, 3, report)
        self.assertTrue(passed)
        with open(report) as f:
            header, row = [line.rstrip("\n").split("\t") for line in f]
        rec = dict(zip(header, row))
        self.assertEqual(rec["gate"], "pasa_complete_orf")
        self.assertEqual(rec["total"], "4")
        self.assertEqual(rec["complete"], "3")
        self.assertEqual(rec["min_complete"], "3")
        self.assertEqual(rec["passed"], "True")

    def test_fails_when_below_minimum(self):
        report = os.path.join(self.tmp, "predict_training_gate.tsv")
        self.assertFalse(lib.run_pasa_training_gate(self.gff3, self.fasta, 4, report))
        with open(report) as f:
            self.assertIn("False", f.read())


class PasaTrainingGateTests(unittest.TestCase):
    def test_passes_when_complete_models_meet_minimum(self):
        ok, msg = lib.pasa_training_gate(
            {"total": 900, "complete": 600, "no_start": 200, "no_stop": 100, "not_mult3": 0},
            min_complete=500,
        )
        self.assertTrue(ok)
        self.assertIn("600", msg)

    def test_fails_below_minimum_and_names_busco_fallback(self):
        ok, msg = lib.pasa_training_gate(
            {"total": 844, "complete": 167, "no_start": 603, "no_stop": 432, "not_mult3": 0},
            min_complete=500,
        )
        self.assertFalse(ok)
        self.assertIn("167", msg)
        self.assertIn("844", msg)
        self.assertIn("--min_pasa_complete_models", msg)
        self.assertIn("BUSCO", msg)

    def test_zero_minimum_disables_gate(self):
        ok, _ = lib.pasa_training_gate(
            {"total": 10, "complete": 0, "no_start": 10, "no_stop": 10, "not_mult3": 0},
            min_complete=0,
        )
        self.assertTrue(ok)


class RnaseqConcordanceGateTests(unittest.TestCase):
    def test_passes_at_or_above_minimum_rate(self):
        ok, msg = lib.rnaseq_concordance_gate(200000, 194600, min_rate=10.0)
        self.assertTrue(ok)
        self.assertIn("97.3", msg)

    def test_fails_below_minimum_rate_with_actionable_message(self):
        ok, msg = lib.rnaseq_concordance_gate(200000, 1000, min_rate=10.0)
        self.assertFalse(ok)
        self.assertIn("0.5", msg)
        self.assertIn("--min_rnaseq_map_rate", msg)
        self.assertIn("host", msg)

    def test_no_reads_sampled_fails(self):
        ok, msg = lib.rnaseq_concordance_gate(0, 0, min_rate=10.0)
        self.assertFalse(ok)
        self.assertIn("no reads", msg.lower())

    def test_zero_minimum_disables_gate(self):
        ok, _ = lib.rnaseq_concordance_gate(200000, 0, min_rate=0)
        self.assertTrue(ok)

    def test_gate_exit_code_is_distinct(self):
        # 1 = generic failure, 75 = OOM/tempfail; the gate needs its own code.
        self.assertEqual(lib.RNASEQ_GATE_EXIT, 3)


@unittest.skipUnless(lib.which("minimap2") and lib.which("samtools"), "needs minimap2 + samtools")
class SampleReadMapRateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        rng = random.Random(7)
        self.genome_seq = _randseq(60000, rng)
        self.genome = os.path.join(self.tmp, "g.fa")
        with open(self.genome, "w") as f:
            f.write(">chr\n" + self.genome_seq + "\n")
        self.reads = os.path.join(self.tmp, "r1.fq.gz")
        with gzip.open(self.reads, "wt") as f:
            for i in range(400):
                if i % 2 == 0:  # on-target read
                    p = rng.randint(0, len(self.genome_seq) - 151)
                    s = self.genome_seq[p:p + 150]
                else:  # off-target read
                    s = _randseq(150, rng)
                f.write("@r{}\n{}\n+\n{}\n".format(i, s, "I" * 150))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_counts_sampled_and_mapped_reads(self):
        sampled, mapped = lib.sample_read_map_rate(
            self.reads, self.genome, n_reads=400, cpus=1, tmpdir=self.tmp
        )
        self.assertEqual(sampled, 400)
        self.assertGreaterEqual(mapped, 190)
        self.assertLessEqual(mapped, 210)

    def test_run_gate_writes_report_and_passes(self):
        report = os.path.join(self.tmp, "train_rnaseq_gate.tsv")
        passed = lib.run_rnaseq_concordance_gate(
            self.reads, self.genome, min_rate=10.0, n_reads=400, cpus=1,
            tmpdir=self.tmp, report=report,
        )
        self.assertTrue(passed)
        with open(report) as f:
            header, row = [line.rstrip("\n").split("\t") for line in f]
        rec = dict(zip(header, row))
        self.assertEqual(rec["gate"], "rnaseq_read_concordance")
        self.assertEqual(rec["sampled"], "400")
        self.assertEqual(rec["min_rate_pct"], "10.0")
        self.assertEqual(rec["passed"], "True")
        self.assertTrue(rec["reads"].endswith("r1.fq.gz"))

    def test_run_gate_fails_above_observed_rate(self):
        report = os.path.join(self.tmp, "train_rnaseq_gate.tsv")
        self.assertFalse(
            lib.run_rnaseq_concordance_gate(
                self.reads, self.genome, min_rate=90.0, n_reads=400, cpus=1,
                tmpdir=self.tmp, report=report,
            )
        )

    def test_respects_sample_size(self):
        sampled, _ = lib.sample_read_map_rate(
            self.reads, self.genome, n_reads=100, cpus=1, tmpdir=self.tmp
        )
        self.assertEqual(sampled, 100)


if __name__ == "__main__":
    unittest.main()


class ReviewFixTests(unittest.TestCase):
    """Fixes from the Fable 5.1 review of gates + R3/R5 (DECISIONS D41)."""

    def test_is_complete_model_rejects_cds_not_multiple_of_three(self):
        # 13-bp CDS "ATGAAACTTTAAG": translation drops the trailing base -> "MKL*"
        gene = {"protein": ["MKL*"], "codon_start": [1], "CDS": [[(1, 13)]], "strand": "+"}
        self.assertFalse(lib.is_complete_model(gene))
        gene["CDS"] = [[(1, 12)]]
        self.assertTrue(lib.is_complete_model(gene))

    def test_gate_applies_when_any_predictor_trains_from_pasa(self):
        modes = {"augustus": "pretrained", "snap": "pasa", "glimmerhmm": "pasa"}
        self.assertTrue(lib.pasa_gate_applies(modes, run_busco=False, augustus_done=False))

    def test_gate_skipped_when_nothing_trains_from_pasa(self):
        modes = {"augustus": "pretrained", "snap": "pretrained"}
        self.assertFalse(lib.pasa_gate_applies(modes, run_busco=False, augustus_done=False))

    def test_gate_skipped_when_busco_training_already_runs(self):
        modes = {"augustus": "busco", "snap": "pasa"}
        self.assertFalse(lib.pasa_gate_applies(modes, run_busco=True, augustus_done=False))

    def test_gate_skipped_on_resume_or_supplied_augustus(self):
        # augustus.gff3 checkpoint exists, or --augustus_gff given: flipping
        # snap to BUSCO now would mix PASA- and BUSCO-trained predictors
        modes = {"augustus": "pasa", "snap": "pasa"}
        self.assertFalse(lib.pasa_gate_applies(modes, run_busco=False, augustus_done=True))


class PasaFeatureFlagTests(unittest.TestCase):
    """R2/F4 opt-in PASA flags (PASApipeline v2.6.1-rc.2, DECISIONS D18/D62) are
    passed only when the installed Launch_PASA_pipeline.pl supports them."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.new = os.path.join(self.tmp, "Launch_new.pl")
        with open(self.new, "w") as f:
            f.write('GetOptions("UNSPLICED_JOIN_SPLICED" => \\$x, "ONE_ALIGNMENT_PER_CDNA" => \\$y);\n')
        self.old = os.path.join(self.tmp, "Launch_old.pl")
        with open(self.old, "w") as f:
            f.write('GetOptions("stringent_alignment_overlap=f" => \\$z);\n')

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_supported_flags_are_passed(self):
        from funannotate import train
        self.assertEqual(
            train.pasa_feature_flags(self.new, unspliced_join=True, one_alignment=True),
            ["--UNSPLICED_JOIN_SPLICED", "--ONE_ALIGNMENT_PER_CDNA"],
        )

    def test_nothing_requested_passes_nothing(self):
        from funannotate import train
        self.assertEqual(train.pasa_feature_flags(self.new, False, False), [])

    def test_unsupported_flags_are_dropped(self):
        from funannotate import train
        self.assertEqual(train.pasa_feature_flags(self.old, True, True), [])

    def test_missing_launcher_passes_nothing(self):
        from funannotate import train
        self.assertEqual(train.pasa_feature_flags(os.path.join(self.tmp, "nope.pl"), True, True), [])


class TrainingDecisionLogTests(unittest.TestCase):
    """Auditable record of every training-data decision (logfiles/training_decisions.tsv)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, "logfiles", "training_decisions.tsv")
        lib.set_training_decision_log(self.path, command="predict")

    def tearDown(self):
        lib.set_training_decision_log(None)
        shutil.rmtree(self.tmp)

    def _rows(self):
        with open(self.path) as f:
            header = f.readline().rstrip("\n").split("\t")
            return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]

    def test_records_rows_with_threshold_and_reason(self):
        lib.record_training_decision("pasa_gate", "complete PASA models", 167, ">=500",
                                     "BUSCO training", "fewer complete models than threshold")
        rows = self._rows()
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual(r["command"], "predict")
        self.assertEqual(r["stage"], "pasa_gate")
        self.assertEqual(r["value"], "167")
        self.assertEqual(r["threshold"], ">=500")
        self.assertEqual(r["outcome"], "BUSCO training")
        self.assertIn("timestamp", r)

    def test_appends_across_commands(self):
        lib.record_training_decision("a", "x", 1, "", "ok", "")
        lib.set_training_decision_log(self.path, command="train")
        lib.record_training_decision("b", "y", 2, "", "ok", "")
        self.assertEqual([r["command"] for r in self._rows()], ["predict", "train"])

    def test_tabs_and_newlines_in_text_are_sanitized(self):
        lib.record_training_decision("s", "d\twith tab", 3, "", "ok", "line1\nline2")
        r = self._rows()[0]
        self.assertEqual(r["decision"], "d with tab")
        self.assertEqual(r["reason"], "line1 line2")

    def test_summary_lists_recorded_decisions(self):
        lib.record_training_decision("pasa_gate", "complete PASA models", 2146, ">=500", "PASA training", "")
        text = lib.training_decision_summary()
        self.assertIn("pasa_gate", text)
        self.assertIn("2146", text)
        self.assertIn("PASA training", text)

    def test_no_log_set_does_not_fail(self):
        lib.set_training_decision_log(None)
        lib.record_training_decision("s", "d", 1, "", "ok", "")  # must not raise
