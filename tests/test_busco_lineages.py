"""BUSCO lineage names (downloads.json, setup, resources) and the provenance recorded for a run."""
import json
import logging
import os
import re

import pytest

from funannotate import library as lib
from funannotate import resources, setupDB

PKG = os.path.dirname(lib.__file__)

# user-facing lineage name -> dataset folder inside the archive, where it is not NAME_odb9 / NAME_ensembl
RENAMED = {"sordariomycetes": "sordariomyceta_odb9", "saccharomycetes": "saccharomyceta_odb9"}
ENSEMBL = {"protists", "alveolata_stramenophiles"}


@pytest.fixture(scope="module")
def busco_links():
    with open(os.path.join(PKG, "downloads.json")) as handle:
        return json.load(handle)["busco"]


def test_lineage_lists_agree(busco_links):
    tree = set(re.findall(r"(\w+) \(\d+\)", resources.buscoTree))
    assert set(setupDB.BUSCO_LINEAGES) == set(busco_links) == tree
    assert len(setupDB.BUSCO_LINEAGES) == len(set(setupDB.BUSCO_LINEAGES))


def test_dataset_folder_names_follow_the_documented_pattern(busco_links):
    for name, (url, folder) in busco_links.items():
        expected = RENAMED.get(name) or (name + ("_ensembl" if name in ENSEMBL else "_odb9"))
        assert folder == expected, name
        assert os.path.basename(url).startswith(folder + ".tar")
    # the one archive with an odd asset name on the db-mirror release
    assert busco_links["sordariomycetes"][0].endswith("sordariomyceta_odb9.tar-2.gz")


def test_dataset_cfg_provenance(tmp_path):
    d = tmp_path / "dikarya"
    d.mkdir()
    (d / "dataset.cfg").write_text(
        "name=dikarya_odb9\nspecies=aspergillus_nidulans\ndomain=eukaryota\n"
        "creation_date=2016-02-13\nnumber_of_BUSCOs=1312\nnumber_of_species=75\n")
    info = lib.busco_lineage_info(str(d))
    assert info["name"] == "dikarya_odb9" and info["number_of_BUSCOs"] == 1312
    assert lib.busco_lineage_label(str(d)) == "dikarya_odb9 (created 2016-02-13, 1,312 BUSCOs)"


def test_label_without_dataset_cfg(tmp_path):
    d = tmp_path / "mystery"
    d.mkdir()
    assert lib.busco_lineage_info(str(d)) == {}
    assert lib.busco_lineage_label(str(d)) == "mystery (no dataset.cfg)"


def make_dataset(root, folder, cfg_name):
    d = root / folder
    d.mkdir()
    if cfg_name is not None:
        (d / "dataset.cfg").write_text("name={}\ncreation_date=2016-02-13\nnumber_of_BUSCOs=1759\n".format(cfg_name))


def test_finish_install_renames_to_the_lineage_name(tmp_path, caplog):
    make_dataset(tmp_path, "saccharomyceta_odb9", "saccharomyceta_odb9")
    with caplog.at_level(logging.INFO):
        info = setupDB.finish_busco_install(str(tmp_path), "saccharomycetes", "saccharomyceta_odb9",
                                            logging.getLogger("t"))
    assert (tmp_path / "saccharomycetes" / "dataset.cfg").is_file()
    assert not (tmp_path / "saccharomyceta_odb9").exists()
    assert info["number_of_BUSCOs"] == 1759
    assert "saccharomycetes: dataset saccharomyceta_odb9 (created 2016-02-13, 1,759 BUSCOs)" in caplog.text


def test_finish_install_same_name_is_untouched(tmp_path):
    make_dataset(tmp_path, "dikarya_odb9", "dikarya_odb9")
    setupDB.finish_busco_install(str(tmp_path), "dikarya_odb9", "dikarya_odb9", logging.getLogger("t"))
    assert (tmp_path / "dikarya_odb9").is_dir()


def test_finish_install_missing_folder_is_a_clear_error(tmp_path):
    with pytest.raises(RuntimeError, match='did not contain the folder "fungi_odb9"'):
        setupDB.finish_busco_install(str(tmp_path), "fungi", "fungi_odb9", logging.getLogger("t"))


def test_finish_install_warns_on_unexpected_dataset(tmp_path, caplog):
    make_dataset(tmp_path, "fungi_odb9", "something_else")
    with caplog.at_level(logging.WARNING):
        setupDB.finish_busco_install(str(tmp_path), "fungi", "fungi_odb9", logging.getLogger("t"))
    assert 'dataset.cfg names "something_else", downloads.json expects "fungi_odb9"' in caplog.text
    make_dataset(tmp_path, "microsporidia_odb9", None)
    with caplog.at_level(logging.WARNING):
        setupDB.finish_busco_install(str(tmp_path), "microsporidia", "microsporidia_odb9", logging.getLogger("t"))
    assert "no dataset.cfg" in caplog.text


def test_predict_and_annotate_record_the_lineage():
    with open(os.path.join(PKG, "predict.py")) as handle:
        predict = handle.read()
    assert '"busco_training", "BUSCO lineage"' in predict
    assert "lib.busco_lineage_label(os.path.join(FUNDB, args.busco_db))" in predict
    with open(os.path.join(PKG, "annotate.py")) as handle:
        annotate = handle.read()
    assert annotate.count("lib.busco_lineage_label(buscoDB)") == 2


def test_setup_choices_come_from_the_lineage_list(monkeypatch):
    import argparse
    seen = {}
    real = argparse.ArgumentParser.add_argument

    def spy(self, *a, **k):
        if "--busco_db" in a:
            seen["choices"] = k.get("choices")
        return real(self, *a, **k)

    monkeypatch.setattr(argparse.ArgumentParser, "add_argument", spy)
    with pytest.raises(SystemExit):
        setupDB.main(["-h"])
    assert seen["choices"] == ["all"] + setupDB.BUSCO_LINEAGES
