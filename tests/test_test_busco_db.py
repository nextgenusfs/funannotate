"""funannotate test --busco_db: run the busco and annotate tests with another lineage."""
from funannotate import test as ftest


def test_default_lineage_adds_no_arguments():
    assert ftest.busco_test_args("dikarya") == []
    assert ftest.busco_test_args(None) == []
    assert ftest.busco_test_args("fungi_odb10") == ["--busco_db", "fungi_odb10"]


def test_predict_arguments_lower_the_training_minimum_for_other_lineages():
    assert ftest.busco_test_predict_args("dikarya") == []
    assert ftest.busco_test_predict_args("fungi_odb10") == [
        "--busco_db", "fungi_odb10", "--min_training_models", "100"]


def test_gene_count_range():
    assert ftest.busco_test_gene_range("dikarya") == (1500, 1800)
    assert ftest.busco_test_gene_range(None) == (1500, 1800)
    low, high = ftest.busco_test_gene_range("fungi_odb10")
    assert (low, high) == (1400, 1750)
    low, high = ftest.busco_test_gene_range("ascomycota_odb10")  # not calibrated: broad window
    assert low < 1400 and high > 1800


def test_option_is_parsed(monkeypatch, capsys):
    import argparse
    seen = {}
    real = argparse.ArgumentParser.parse_args

    def spy(self, args=None, namespace=None):
        ns = real(self, args, namespace)
        seen["busco_db"] = ns.busco_db
        raise SystemExit(0)

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", spy)
    try:
        ftest.main(["-t", "annotate", "--busco_db", "fungi_odb10"])
    except SystemExit:
        pass
    assert seen["busco_db"] == "fungi_odb10"
