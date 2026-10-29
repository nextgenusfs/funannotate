"""#1102: a failing trainGlimmerHMM must skip GlimmerHMM, not end `funannotate predict`."""
import logging
import subprocess

from funannotate import library as lib


def test_training_failure_returns_none_and_writes_no_predictions(tmp_path, monkeypatch):
    monkeypatch.setattr(lib, "log", logging.getLogger("t"), raising=False)
    monkeypatch.setattr(lib, "gff2dict", lambda *a, **k: {})
    monkeypatch.setattr(lib, "dict2glimmer", lambda genes, out: open(out, "w").close())
    calls = []

    def fake_run(cmd, directory, logfile, **kw):
        calls.append(cmd[0])
        assert kw.get("raise_not_exit") is True
        raise subprocess.CalledProcessError(255, cmd)

    monkeypatch.setattr(lib, "runSubprocess", fake_run)
    out = tmp_path / "glimmer.gff3"
    assert lib.runGlimmerHMM("g.fa", "m.gff3", str(tmp_path), str(out)) is None
    assert calls == ["trainGlimmerHMM"]  # prediction was not attempted
    assert not out.exists()
