"""`funannotate check` looks for the programs aux_scripts/funannotate-BUSCO2.py needs.

The script exits with its own message when etraining or one of the Augustus perl scripts
(new_species.pl, gff2gbSmallDNA.pl, optimize_augustus.pl) is missing; tblastn, makeblastdb,
hmmsearch and augustus were already checked.
"""
import os
import stat

import funannotate
from funannotate import check

AUX = os.path.join(os.path.dirname(funannotate.__file__), "aux_scripts", "funannotate-BUSCO2.py")
CHECK = os.path.join(os.path.dirname(funannotate.__file__), "check.py")


def _exe(path):
    path.write_text("#!/bin/sh\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def test_every_busco2_program_is_checked():
    with open(AUX) as fh:
        script = fh.read()
    with open(CHECK) as fh:
        checker = fh.read()
    for prog in ["etraining", "new_species.pl", "gff2gbSmallDNA.pl", "optimize_augustus.pl"]:
        assert prog in script
        assert prog in check.PRESENCE_ONLY
    for prog in ["tblastn", "makeblastdb", "hmmsearch", "augustus"]:
        assert prog in script
        assert "'%s'" % prog in checker  # in programs1..programs6


def test_found_on_path(tmp_path, monkeypatch):
    exe = _exe(tmp_path / "etraining")
    monkeypatch.setenv("PATH", str(tmp_path))
    assert check.find_presence_only("etraining") == str(exe)


def test_augustus_script_in_config_sibling_scripts_dir(tmp_path, monkeypatch):
    (tmp_path / "config").mkdir()
    (tmp_path / "scripts").mkdir()
    exe = _exe(tmp_path / "scripts" / "new_species.pl")
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", str(tmp_path / "config") + os.sep)
    assert check.find_presence_only("new_species.pl") == str(exe)
    # the script uses that fallback only when the config dir is named config
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", str(tmp_path / "scripts"))
    assert check.find_presence_only("new_species.pl") is False


def test_etraining_has_no_scripts_fallback(tmp_path, monkeypatch):
    (tmp_path / "config").mkdir()
    (tmp_path / "scripts").mkdir()
    _exe(tmp_path / "scripts" / "etraining")
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", str(tmp_path / "config"))
    assert check.find_presence_only("etraining") is False


def test_missing(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.delenv("AUGUSTUS_CONFIG_PATH", raising=False)
    assert check.find_presence_only("optimize_augustus.pl") is False
