"""#1069: getTrainResults raised UnboundLocalError when Augustus printed no accuracy report."""
from funannotate import library as lib

REPORT = (
    "nucleotide level | 0.9 | 0.8 |\n"
    "exon level | 1 | 2 | 3 | 4 | 5 | 0.5 | 0.6 |\n"
    "gene level | 1 | 2 | 3 | 4 | 5 | 0.3 | 0.4 |\n"
)


def test_full_report(tmp_path):
    f = tmp_path / "r.txt"
    f.write_text(REPORT)
    assert lib.getTrainResults(str(f)) == (0.9, 0.8, 0.5, 0.6, 0.3, 0.4)


def test_failed_augustus_run_gives_none(tmp_path):
    f = tmp_path / "r.txt"
    f.write_text("augustus: ERROR\n\tNo genbank sequences found.\n")
    assert lib.getTrainResults(str(f)) is None


def test_partial_report_gives_none(tmp_path):
    f = tmp_path / "r.txt"
    f.write_text(REPORT.splitlines()[0] + "\n")
    assert lib.getTrainResults(str(f)) is None
