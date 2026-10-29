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


def test_script_runs_with_this_interpreter(tmp_path, monkeypatch):
    """Not through the `#!/usr/bin/env python` shebang, which may find no python, or another one."""
    import sys

    monkeypatch.setattr(lib, "log", logging.getLogger("t"), raising=False)
    calls = []
    monkeypatch.setattr(lib, "runSubprocess", lambda cmd, *a, **k: calls.append(cmd))
    run = tmp_path / "run_busco"
    run.mkdir()
    (run / "full_table_busco.tsv").write_text("EOG1\tComplete\tp1\t1\t1\n")
    lib.runBUSCO("p.fa", "/db/dikarya", 2, str(tmp_path), str(tmp_path / "out.txt"))
    assert calls[0][0] == sys.executable
    assert calls[0][1].endswith("funannotate-BUSCO2.py")


def test_no_complete_buscos_is_a_warning(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(lib, "log", logging.getLogger("t"), raising=False)
    run = tmp_path / "run_busco"
    run.mkdir()
    (run / "full_table_busco.tsv").write_text(
        "# header\nEOG2\tMissing\nEOG4\tFragmented\tp4\t1\t1\n"
    )
    monkeypatch.setattr(lib, "runSubprocess", lambda *a, **k: None)
    out = tmp_path / "out.txt"
    with caplog.at_level(logging.WARNING):
        lib.runBUSCO("p.fa", "/db/dikarya", 2, str(tmp_path), str(out))
    assert out.read_text() == ""
    assert "no Complete or Duplicated dikarya" in caplog.text
    assert [r.levelno for r in caplog.records] == [logging.WARNING]


def test_no_warning_when_buscos_found(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(lib, "log", logging.getLogger("t"), raising=False)
    run = tmp_path / "run_busco"
    run.mkdir()
    (run / "full_table_busco.tsv").write_text("EOG1\tComplete\tp1\t1\t1\n")
    monkeypatch.setattr(lib, "runSubprocess", lambda *a, **k: None)
    with caplog.at_level(logging.WARNING):
        lib.runBUSCO("p.fa", "/db/dikarya", 2, str(tmp_path), str(tmp_path / "out.txt"))
    assert caplog.text == ""
