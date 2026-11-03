"""A BUSCO that became Complete only after the retrained second pass has no GFF; it must be skipped,
not crash the BUSCO training step (found with fungi_odb10: 38 of 199 models)."""
import os

from funannotate import library as lib


def make(tmp_path, ids):
    d = tmp_path / "busco" / "run_aug" / "augustus_output" / "gffs"
    d.mkdir(parents=True)
    for i in ids:
        (d / (i + ".gff")).write_text("x\n")
    return str(tmp_path / "busco")


def test_split_into_found_and_missing(tmp_path):
    loc = make(tmp_path, ["771at4751", "EOG0926004Z"])
    found, missing = lib.busco_gff_paths(loc, "aug", ["771at4751", "119833at4751", "EOG0926004Z"])
    assert [os.path.basename(f) for f in found] == ["771at4751.gff", "EOG0926004Z.gff"]
    assert missing == ["119833at4751"]


def test_no_gff_folder_at_all(tmp_path):
    found, missing = lib.busco_gff_paths(str(tmp_path / "nothing"), "aug", ["a", "b"])
    assert found == [] and missing == ["a", "b"]


def test_predict_uses_the_helper_and_does_not_open_missing_files():
    src = open(os.path.join(os.path.dirname(lib.__file__), "predict.py")).read()
    assert "lib.busco_gff_paths(" in src
    assert "for file in busco_gffs:" in src
