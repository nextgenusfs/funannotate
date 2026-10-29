"""#1110: SignalP 6 IDs keep the FASTA header description; #1120: which() on unreadable PATH entries."""
import os
import stat

from funannotate import library as lib


def test_signalp6_id_is_first_word(tmp_path):
    sp = tmp_path / "sp.txt"
    sp.write_text(
        "# SignalP-6.0\tOrganism: Eukarya\n"
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        "gene1-T1 gene1\tSP\t0.01\t0.99\tCS pos: 20-21. Pr: 0.9\n"
        "gene2-T1\tOTHER\t0.99\t0.01\t\n"
    )
    out = tmp_path / "sec.txt"
    lib.parseSignalP(str(sp), str(out))
    lines = out.read_text().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("gene1-T1\tnote\tSECRETED:SignalP(1-20")


def test_signalp5_id_is_first_word(tmp_path):
    sp = tmp_path / "sp.txt"
    sp.write_text(
        "# SignalP-5.0\tOrganism: Eukarya\n"
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        "gene1-T1 gene1\tSP(Sec/SPI)\t0.01\t0.99\tCS pos: 20-21. AB-CD. Pr: 0.9\n"
    )
    out = tmp_path / "sec.txt"
    lib.parseSignalP(str(sp), str(out))
    assert out.read_text().startswith("gene1-T1\tnote\tSECRETED:SignalP(")


def test_which_is_false_when_path_entry_unreadable(tmp_path, monkeypatch):
    blocked = tmp_path / "blocked"
    blocked.mkdir()
    blocked.chmod(0)
    try:
        monkeypatch.setenv("PATH", str(blocked))
        assert lib.which("definitely-not-a-program-xyz") is False
    finally:
        blocked.chmod(stat.S_IRWXU)
