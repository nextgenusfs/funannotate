"""Tests for choosing the training source: PASA/transcript models vs BUSCO.

predict.py makes two decisions, both reduced here to pure functions:

1. resolve_pasa_gate: the PASA complete-ORF gate and the RNA-seq identity gate
   say whether PASA models are fit to train from; if not, switch to BUSCO,
   unless BUSCO cannot run (then keep PASA rather than exit).
2. resolve_min_training_models: even a passing training set must reach
   --min_training_models; if not, fall back to BUSCO or exit.

The chain tests feed real PASA-style GFF3 files through the gate and into
those decisions, including the poor-evidence cases (no output, no CDS, all
partial models) that must hand training to BUSCO alone.
"""
import os
import unittest

from funannotate import library as lib
from test_training_gates import TrainingModelFixture

MIN_COMPLETE = 3  # the fixture has 3 complete + 1 partial model
MIN_MODELS = 200  # --min_training_models default


def source_for(gff3, fasta, report, min_complete=MIN_COMPLETE, identity_ok=True,
               busco_fallback=True, augustus_functional=True):
    """Run the PASA gate on a GFF3 and return the resolved training source."""
    passed = lib.run_pasa_training_gate(gff3, fasta, min_complete, report)
    return lib.resolve_pasa_gate(passed, identity_ok, busco_fallback, augustus_functional)


class ResolvePasaGateTests(unittest.TestCase):
    def test_both_gates_pass_keeps_pasa(self):
        self.assertEqual(lib.resolve_pasa_gate(True, True, True, True), "pasa")

    def test_pasa_gate_fails_switches_to_busco(self):
        self.assertEqual(lib.resolve_pasa_gate(False, True, True, True), "busco")

    def test_identity_gate_fails_switches_to_busco(self):
        self.assertEqual(lib.resolve_pasa_gate(True, False, True, True), "busco")

    def test_both_fail_switches_to_busco(self):
        self.assertEqual(lib.resolve_pasa_gate(False, False, True, True), "busco")

    def test_no_busco_fallback_keeps_pasa_despite_failure(self):
        self.assertEqual(lib.resolve_pasa_gate(False, True, False, True), "pasa_override")

    def test_broken_augustus_keeps_pasa_despite_failure(self):
        # BUSCO training needs augustus --proteinprofile
        self.assertEqual(lib.resolve_pasa_gate(False, True, True, False), "pasa_override")

    def test_passing_gates_ignore_fallback_settings(self):
        self.assertEqual(lib.resolve_pasa_gate(True, True, False, False), "pasa")


class ResolveMinTrainingModelsTests(unittest.TestCase):
    def resolve(self, total, source="pasa", fallback=True, functional=True, minimum=MIN_MODELS):
        return lib.resolve_min_training_models(total, minimum, source, fallback, functional)

    def test_enough_models(self):
        self.assertEqual(self.resolve(500), "enough")

    def test_exactly_at_minimum_is_enough(self):
        self.assertEqual(self.resolve(MIN_MODELS), "enough")

    def test_one_below_minimum_falls_back_to_busco(self):
        self.assertEqual(self.resolve(MIN_MODELS - 1), "busco_fallback")

    def test_zero_models_falls_back_to_busco(self):
        self.assertEqual(self.resolve(0), "busco_fallback")

    def test_pretrained_needs_no_training_set(self):
        self.assertEqual(self.resolve(0, source="pretrained"), "skip")

    def test_too_few_busco_models_exit_instead_of_looping(self):
        # BUSCO is already the source; falling back to it again makes no sense
        self.assertEqual(self.resolve(10, source="busco"), "exit")

    def test_no_busco_fallback_exits(self):
        self.assertEqual(self.resolve(10, fallback=False), "exit")

    def test_fallback_needs_working_augustus(self):
        self.assertEqual(self.resolve(10, functional=False), "exit_no_augustus")

    def test_enough_models_do_not_need_augustus_check(self):
        self.assertEqual(self.resolve(MIN_MODELS, functional=False), "enough")

    def test_minimum_zero_never_triggers(self):
        self.assertEqual(self.resolve(0, minimum=0), "enough")


class PoorEvidenceChainTests(TrainingModelFixture):
    """PASA output that is empty, unusable or fragmentary trains from BUSCO."""

    def setUp(self):
        super().setUp()
        self.report = os.path.join(self.tmp, "gate.tsv")

    def write(self, name, text):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as f:
            f.write(text)
        return path

    def test_good_evidence_trains_from_pasa(self):
        self.assertEqual(source_for(self.gff3, self.fasta, self.report), "pasa")

    def test_pasa_produced_no_file(self):
        absent = os.path.join(self.tmp, "no_pasa_output.gff3")
        self.assertEqual(source_for(absent, self.fasta, self.report), "busco")

    def test_pasa_output_is_empty_file(self):
        # PASA found no meaningful alignments: zero-byte GFF3
        empty = self.write("empty.gff3", "")
        self.assertEqual(source_for(empty, self.fasta, self.report), "busco")

    def test_pasa_output_has_header_only(self):
        header = self.write("header.gff3", "##gff-version 3\n")
        self.assertEqual(source_for(header, self.fasta, self.report), "busco")

    def test_models_without_cds_count_as_nothing(self):
        # genes/mRNAs/exons but no CDS (no ORF could be called)
        noncoding = self.write(
            "noncoding.gff3",
            "##gff-version 3\n"
            "ctg1\tPASA\tgene\t1\t900\t.\t+\t.\tID=g1\n"
            "ctg1\tPASA\tmRNA\t1\t900\t.\t+\t.\tID=g1.t1;Parent=g1\n"
            "ctg1\tPASA\texon\t1\t900\t.\t+\t.\tID=g1.t1.exon1;Parent=g1.t1\n",
        )
        counts = lib.count_complete_orf_models(noncoding, self.fasta)
        self.assertEqual(counts["total"], 0)
        self.assertEqual(source_for(noncoding, self.fasta, self.report), "busco")

    def test_all_models_partial_trains_from_busco(self):
        # keep only g4 (no start codon): plenty of models in principle, none usable
        with open(self.gff3) as f:
            lines = [ln for ln in f if ln.startswith("#") or "g4" in ln]
        partial = self.write("partial.gff3", "".join(lines))
        counts = lib.count_complete_orf_models(partial, self.fasta)
        self.assertEqual((counts["total"], counts["complete"]), (1, 0))
        self.assertEqual(source_for(partial, self.fasta, self.report, min_complete=1), "busco")

    def test_boundary_exactly_min_complete_passes(self):
        self.assertEqual(source_for(self.gff3, self.fasta, self.report, min_complete=3), "pasa")

    def test_boundary_one_above_min_complete_fails(self):
        self.assertEqual(source_for(self.gff3, self.fasta, self.report, min_complete=4), "busco")

    def test_low_identity_overrides_good_models(self):
        # good models from reads of another strain/species: identity gate decides
        self.assertEqual(
            source_for(self.gff3, self.fasta, self.report, identity_ok=False), "busco")

    def test_poor_evidence_without_busco_keeps_pasa(self):
        empty = self.write("empty.gff3", "")
        self.assertEqual(
            source_for(empty, self.fasta, self.report, busco_fallback=False), "pasa_override")

    def test_disabled_gate_passes_empty_pasa_but_min_models_still_catches_it(self):
        # --min_pasa_complete_models 0 turns the first gate off ...
        empty = self.write("empty.gff3", "")
        self.assertEqual(source_for(empty, self.fasta, self.report, min_complete=0), "pasa")
        # ... and the second decision then routes zero models to BUSCO
        total = lib.countGFFgenes(empty)
        self.assertEqual(total, 0)
        self.assertEqual(
            lib.resolve_min_training_models(total, MIN_MODELS, "pasa", True, True),
            "busco_fallback")

    def test_gate_passes_but_set_too_small_for_augustus(self):
        # 3 complete models clear the gate yet are far below --min_training_models
        self.assertEqual(source_for(self.gff3, self.fasta, self.report), "pasa")
        total = lib.countGFFgenes(self.gff3)
        self.assertEqual(
            lib.resolve_min_training_models(total, MIN_MODELS, "pasa", True, True),
            "busco_fallback")


class GateAppliesChainTests(unittest.TestCase):
    """The PASA gate is only consulted when it can change something."""

    def test_pretrained_augustus_skips_gate_and_training_size_check(self):
        modes = {"augustus": "pretrained", "snap": "pretrained", "glimmerhmm": "pretrained"}
        self.assertFalse(lib.pasa_gate_applies(modes, False, False))
        self.assertEqual(
            lib.resolve_min_training_models(0, MIN_MODELS, modes["augustus"], True, True),
            "skip")

    def test_busco_source_is_final_even_when_too_small(self):
        modes = {"augustus": "busco", "snap": "busco", "glimmerhmm": "busco"}
        self.assertFalse(lib.pasa_gate_applies(modes, True, False))
        self.assertEqual(
            lib.resolve_min_training_models(50, MIN_MODELS, modes["augustus"], True, True),
            "exit")


if __name__ == "__main__":
    unittest.main()
