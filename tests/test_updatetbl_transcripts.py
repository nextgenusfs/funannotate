"""updateTBL: a second transcript without annotations must not inherit the first transcript's."""
import logging

import pytest

from funannotate import library as lib


@pytest.fixture(autouse=True)
def _log(monkeypatch):
    monkeypatch.setattr(lib, "log", logging.getLogger("test"), raising=False)

TBL = (
    ">Feature ctg1\n"
    "1\t900\tgene\n"
    "\t\t\tlocus_tag\tFUN_000001\n"
    "1\t300\tmRNA\n"
    "\t\t\tproduct\thypothetical protein\n"
    "\t\t\tprotein_id\tgnl|ncbi|FUN_000001-T1\n"
    "\t\t\ttranscript_id\tgnl|ncbi|FUN_000001-T1_mrna\n"
    "1\t300\tCDS\n"
    "\t\t\tproduct\thypothetical protein\n"
    "\t\t\tprotein_id\tgnl|ncbi|FUN_000001-T1\n"
    "\t\t\tcodon_start\t1\n"
    "1\t900\tmRNA\n"
    "\t\t\tproduct\thypothetical protein\n"
    "\t\t\tprotein_id\tgnl|ncbi|FUN_000001-T2\n"
    "\t\t\ttranscript_id\tgnl|ncbi|FUN_000001-T2_mrna\n"
    "1\t900\tCDS\n"
    "\t\t\tproduct\thypothetical protein\n"
    "\t\t\tprotein_id\tgnl|ncbi|FUN_000001-T2\n"
    "\t\t\tcodon_start\t1\n"
)


def run(tmp_path, annot):
    src, dst = tmp_path / "in.tbl", tmp_path / "out.tbl"
    src.write_text(TBL)
    lib.updateTBL(str(src), annot, str(dst))
    return dst.read_text()


def test_only_the_annotated_transcript_gets_notes(tmp_path):
    out = run(tmp_path, {"FUN_000001-T1": {"note": ["PFAM:PF00001"]}})
    assert out.count("PFAM:PF00001") == 1


def test_both_annotated_each_get_their_own(tmp_path):
    out = run(tmp_path, {"FUN_000001-T1": {"note": ["PFAM:PF00001"]},
                         "FUN_000001-T2": {"note": ["PFAM:PF00002"]}})
    assert out.count("PFAM:PF00001") == 1 and out.count("PFAM:PF00002") == 1
