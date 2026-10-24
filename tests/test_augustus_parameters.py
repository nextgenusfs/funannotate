"""Tests for the Augustus parameter decisions in predict/train.

* where the Augustus parameters come from: pre-trained paths from a
  parameters.json (absolute, relative to the file, relative to the working
  directory, or missing -> train from scratch);
* how many training models are set aside to test Augustus;
* how the held-out accuracy printed by `augustus` is read, labelled and judged.

The accuracy fixture copies the layout of real `augustus --species=X test.gb`
output (Augustus 3.x): sensitivity comes before specificity on every line.
"""
import math
import os
import shutil
import tempfile
import unittest

from funannotate import library as lib
from funannotate.predict import resolveTrainingPaths

EVAL_TEMPLATE = """\

---------------------------------------------\\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       {n_sn} |       {n_sp} |
---------------------------------------------/

----------------------------------------------------------------------------------------------------------\\
           |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
           | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
           | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
----------------------------------------------------------------------------------------------------------|
           |        |        |      |                  0 |                  2 |             |             |
exon level |    100 |    120 |   80 | ------------------ | ------------------ |       {e_sn} |       {e_sp} |
           |    100 |    120 |      |    0 |    0 |    0 |    0 |    0 |    2 |             |             |
----------------------------------------------------------------------------------------------------------/

----------------------------------------------------------------------------\\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |    50 |    60 |   30 |   20 |   30 |       {g_sn} |       {g_sp} |
----------------------------------------------------------------------------/
# total time: 0.0358
"""


def write_eval(path, n=(0.91, 0.85), e=(0.72, 0.66), g=(0.41, 0.35)):
    with open(path, "w") as f:
        f.write(
            EVAL_TEMPLATE.format(
                n_sn=n[0], n_sp=n[1], e_sn=e[0], e_sp=e[1], g_sn=g[0], g_sp=g[1]
            )
        )


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cwd = os.getcwd()

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.tmp)

    def mkdir(self, *parts):
        path = os.path.join(self.tmp, *parts)
        os.makedirs(path, exist_ok=True)
        return path


class ResolveTrainingPathsTests(TempDirCase):
    def data(self, **paths):
        return {name: [{"path": path}] for name, path in paths.items()}

    def test_absolute_path_is_kept(self):
        aug = self.mkdir("params", "augustus")
        out, missing = resolveTrainingPaths(self.data(augustus=aug), self.tmp)
        self.assertEqual(out["augustus"][0]["path"], os.path.abspath(aug))
        self.assertEqual(missing, [])

    def test_relative_path_resolves_against_parameters_dir(self):
        self.mkdir("store", "augustus")
        out, missing = resolveTrainingPaths(
            self.data(augustus="augustus"), os.path.join(self.tmp, "store")
        )
        self.assertEqual(
            out["augustus"][0]["path"], os.path.join(self.tmp, "store", "augustus")
        )
        self.assertEqual(missing, [])

    def test_relative_path_falls_back_to_working_directory(self):
        self.mkdir("elsewhere")
        self.mkdir("work", "snap")
        os.chdir(os.path.join(self.tmp, "work"))
        out, missing = resolveTrainingPaths(
            self.data(snap="snap"), os.path.join(self.tmp, "elsewhere")
        )
        self.assertEqual(
            os.path.realpath(out["snap"][0]["path"]),
            os.path.realpath(os.path.join(self.tmp, "work", "snap")),
        )
        self.assertEqual(missing, [])

    def test_anchor_wins_over_working_directory(self):
        self.mkdir("store", "augustus")
        self.mkdir("work", "augustus")
        os.chdir(os.path.join(self.tmp, "work"))
        out, _ = resolveTrainingPaths(
            self.data(augustus="augustus"), os.path.join(self.tmp, "store")
        )
        self.assertEqual(
            out["augustus"][0]["path"], os.path.join(self.tmp, "store", "augustus")
        )

    def test_missing_path_is_reported_for_retraining(self):
        out, missing = resolveTrainingPaths(
            self.data(augustus=os.path.join(self.tmp, "gone")), self.tmp
        )
        self.assertEqual([p for p, _ in missing], ["augustus"])

    def test_only_the_missing_predictor_is_reported(self):
        good = self.mkdir("good")
        out, missing = resolveTrainingPaths(
            self.data(augustus=good, snap=os.path.join(self.tmp, "gone")), self.tmp
        )
        self.assertEqual([p for p, _ in missing], ["snap"])
        self.assertEqual(out["augustus"][0]["path"], os.path.abspath(good))

    def test_untrained_predictors_are_ignored(self):
        data = {"augustus": [{}], "genemark": [], "snap": [{"path": ""}]}
        out, missing = resolveTrainingPaths(data, self.tmp)
        self.assertEqual(missing, [])

    def test_unknown_predictor_names_are_left_alone(self):
        data = {"mystery": [{"path": "/does/not/exist"}]}
        out, missing = resolveTrainingPaths(data, self.tmp)
        self.assertEqual(missing, [])
        self.assertEqual(out["mystery"][0]["path"], "/does/not/exist")


class AugustusTestSetSizeTests(unittest.TestCase):
    def test_small_sets_hold_out_a_fifth(self):
        self.assertEqual(lib.augustus_test_set_size(200), 40)
        self.assertEqual(lib.augustus_test_set_size(499), 100)

    def test_middle_range_holds_out_one_hundred(self):
        self.assertEqual(lib.augustus_test_set_size(500), 100)
        self.assertEqual(lib.augustus_test_set_size(1000), 100)

    def test_large_sets_hold_out_a_tenth(self):
        self.assertEqual(lib.augustus_test_set_size(1001), 100)
        self.assertEqual(lib.augustus_test_set_size(5000), 500)

    def test_default_minimum_leaves_training_models(self):
        # --min_training_models default is 200: 40 held out, 160 train
        total = 200
        self.assertLess(lib.augustus_test_set_size(total), total)


class TrainResultsParsingTests(TempDirCase):
    def test_reads_sensitivity_then_specificity(self):
        path = os.path.join(self.tmp, "aug.txt")
        write_eval(path)
        self.assertEqual(
            lib.getTrainResults(path), (0.91, 0.85, 0.72, 0.66, 0.41, 0.35)
        )

    def test_nan_specificity_is_parsed(self):
        # Augustus prints nan when it predicts nothing
        path = os.path.join(self.tmp, "aug.txt")
        write_eval(path, n=(0, "nan"), e=(0, "nan"), g=(0, "nan"))
        res = lib.getTrainResults(path)
        self.assertEqual(res[4], 0.0)
        self.assertTrue(math.isnan(res[5]))


class ResultsTableTests(unittest.TestCase):
    def test_headers_match_the_order_of_the_values(self):
        table = lib.augustus_results_table((0.91, 0.85, 0.72, 0.66, 0.41, 0.35))
        self.assertEqual(table[0], ["Feature", "Sensitivity", "Specificity"])
        self.assertEqual(table[1], ["nucleotides", "91.0%", "85.0%"])
        self.assertEqual(table[2], ["exons", "72.0%", "66.0%"])
        self.assertEqual(table[3], ["genes", "41.0%", "35.0%"])


class TrainingAssessmentTests(unittest.TestCase):
    def results(self, gene_sn):
        return (0.9, 0.8, 0.7, 0.6, gene_sn, 0.5)

    def test_good_accuracy_is_ok(self):
        self.assertEqual(lib.augustus_training_assessment(self.results(0.62)), "ok")

    def test_low_gene_sensitivity_is_low(self):
        self.assertEqual(lib.augustus_training_assessment(self.results(0.31)), "low")

    def test_boundary_at_threshold_is_ok(self):
        self.assertEqual(
            lib.augustus_training_assessment(self.results(lib.AUGUSTUS_LOW_GENE_SN)),
            "ok",
        )

    def test_uses_gene_sensitivity_not_specificity(self):
        res = (0.9, 0.8, 0.7, 0.6, 0.80, 0.10)  # high sensitivity, low specificity
        self.assertEqual(lib.augustus_training_assessment(res), "ok")

    def test_nothing_predicted_is_low(self):
        res = (0.0, float("nan"), 0.0, float("nan"), 0.0, float("nan"))
        self.assertEqual(lib.augustus_training_assessment(res), "low")

    def test_nan_gene_sensitivity_is_unmeasured(self):
        res = (0.9, 0.8, 0.7, 0.6, float("nan"), 0.5)
        self.assertEqual(lib.augustus_training_assessment(res), "unmeasured")

    def test_missing_results_are_unmeasured(self):
        self.assertEqual(lib.augustus_training_assessment(None), "unmeasured")

    def test_custom_threshold(self):
        self.assertEqual(
            lib.augustus_training_assessment(self.results(0.62), min_gene_sn=0.7), "low"
        )


if __name__ == "__main__":
    unittest.main()
