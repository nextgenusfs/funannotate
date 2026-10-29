"""1.10 BUSCO hygiene: the Python 2 copy and the functions without callers are gone."""
import os

from funannotate import library, setupDB

AUX = os.path.join(os.path.dirname(library.__file__), "aux_scripts")


def test_py2_copy_removed_and_py3_script_present():
    assert not os.path.exists(os.path.join(AUX, "funannotate-BUSCO2-py2.py"))
    assert os.path.isfile(os.path.join(AUX, "funannotate-BUSCO2.py"))


def test_functions_without_callers_removed():
    for name in ("runBUSCOv5", "genemark2busco", "evidence2busco"):
        assert not hasattr(library, name)
    assert not hasattr(setupDB, "download_busco_v5")
    assert hasattr(library, "runBUSCO")


def test_training_source_label_is_spelled_busco():
    path = os.path.join(os.path.dirname(library.__file__), "predict.py")
    with open(path) as handle:
        text = handle.read()
    assert "BUCSCO" not in text
    assert 'sourceText = "BUSCO {:}"' in text
