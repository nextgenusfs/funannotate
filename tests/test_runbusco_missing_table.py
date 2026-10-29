"""#1108: funannotate-BUSCO2.py exits 0 when it fails; annotate then died with FileNotFoundError."""
import logging

import pytest

from funannotate import library as lib


def test_missing_table_is_a_clear_error(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(lib, "log", logging.getLogger("t"), raising=False)
    monkeypatch.setattr(lib, "runSubprocess", lambda *a, **k: None)  # "succeeds", writes nothing
    with caplog.at_level(logging.ERROR):
        with pytest.raises(SystemExit) as e:
            lib.runBUSCO("p.fa", "/db/dikarya", 2, str(tmp_path), str(tmp_path / "out.txt"))
    assert e.value.code == 1
    assert "full_table_busco.tsv" in caplog.text and "/db/dikarya" in caplog.text


def test_table_is_parsed(tmp_path, monkeypatch):
    monkeypatch.setattr(lib, "log", logging.getLogger("t"), raising=False)
    run = tmp_path / "run_busco"
    run.mkdir()
    (run / "full_table_busco.tsv").write_text(
        "# header\nEOG1\tComplete\tp1\t1\t1\nEOG2\tMissing\nEOG3\tDuplicated\tp3\t1\t1\nEOG4\tFragmented\tp4\t1\t1\n"
    )
    monkeypatch.setattr(lib, "runSubprocess", lambda *a, **k: None)
    out = tmp_path / "out.txt"
    lib.runBUSCO("p.fa", "/db/dikarya", 2, str(tmp_path), str(out))
    assert out.read_text() == "p1\tnote\tBUSCO:EOG1\np3\tnote\tBUSCO:EOG3\n"
