"""`funannotate setup -i busco -b fungi_odb10`: OrthoDB 10 lineages from the BUSCO download server."""
import hashlib
import io
import logging
import os
import tarfile

import pytest

from funannotate import resources, setupDB

INDEX = (
    "fungi_odb10\t2024-01-08\td389cb627ded374b47db05630191494e\tEukaryota\tlineages\n"
    "microsporidia_odb10\t2024-01-08\tc6c96c9fddc451b403ba83630506101e\tEukaryota\tlineages\n"
    "bacteria_odb10\t2024-01-08\t00000000000000000000000000000000\tProkaryota\tlineages\n"
    "fungi_odb12\t2026-05-22\t11111111111111111111111111111111\tEukaryota\tlineages\n"
    "fungi_odb10\t2024-01-08\tbroken line\n"
    "eukaryota_odb10\t2024-01-08\t93e2a6ec05708dbf0f23db078b88c5a8\tEukaryota\tlineages\n"
)


def test_name_validation():
    for ok in ("all", "dikarya", "sordariomycetes", "fungi_odb10", "alveolata_odb10", "mucoromycota_odb10"):
        assert setupDB.busco_name_problem(ok) is None, ok
    odb12 = setupDB.busco_name_problem("fungi_odb12")
    assert "lengths_cutoff" in odb12 and "_odb10" in odb12
    assert "odb12" in setupDB.busco_name_problem("fungi_odb12.2")
    for bad in ("fungi_odb11", "Fungi_odb10", "fungi odb10", "../fungi_odb10", "fungi", "nonsense"):
        msg = setupDB.busco_name_problem(bad)
        if bad == "fungi":  # the odb9 name is valid
            assert msg is None
        else:
            assert msg and "fungi_odb10" in msg, bad


def test_parse_file_versions_skips_malformed_lines():
    index = setupDB.parse_file_versions(INDEX)
    # the malformed second fungi_odb10 line (3 columns) is skipped, the valid one is kept
    assert index["fungi_odb10"] == ("2024-01-08", "d389cb627ded374b47db05630191494e", "Eukaryota", "lineages")
    assert "bacteria_odb10" in index and "fungi_odb12" in index
    assert index["microsporidia_odb10"][1] == "c6c96c9fddc451b403ba83630506101e"


def test_odb10_lineage_url_and_refusals():
    index = setupDB.parse_file_versions(INDEX)
    url, md5 = setupDB.odb10_lineage("microsporidia_odb10", index)
    assert url == "https://busco-data.ezlab.org/v5/data/lineages/microsporidia_odb10.2024-01-08.tar.gz"
    assert md5 == "c6c96c9fddc451b403ba83630506101e"
    with pytest.raises(ValueError, match="Prokaryota"):
        setupDB.odb10_lineage("bacteria_odb10", index)
    with pytest.raises(ValueError, match="not in the BUSCO dataset list"):
        setupDB.odb10_lineage("nothing_odb10", index)


def make_archive(path, folder="microsporidia_odb10", extra=None):
    cfg = "name={0}\ncreation_date=2024-01-08\nnumber_of_BUSCOs=600\n".format(folder).encode()
    with tarfile.open(path, "w:gz") as tar:
        info = tarfile.TarInfo(folder + "/dataset.cfg")
        info.size = len(cfg)
        tar.addfile(info, io.BytesIO(cfg))
        for name, data in (extra or {}).items():
            member = tarfile.TarInfo(name)
            member.size = len(data)
            tar.addfile(member, io.BytesIO(data))


def md5_of(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def test_install_verifies_extracts_and_cleans_up(tmp_path):
    src = tmp_path / "server.tar.gz"
    make_archive(str(src))
    index = {"microsporidia_odb10": ("2024-01-08", md5_of(str(src)), "Eukaryota", "lineages")}
    fundb = tmp_path / "db"
    fundb.mkdir()
    urls = []

    def fetch(url, path):
        urls.append(url)
        with open(src, "rb") as a, open(path, "wb") as b:
            b.write(a.read())

    md5 = setupDB.install_odb10_lineage("microsporidia_odb10", str(fundb), index, fetch, logging.getLogger("t"))
    assert md5 == md5_of(str(src))
    assert urls == ["https://busco-data.ezlab.org/v5/data/lineages/microsporidia_odb10.2024-01-08.tar.gz"]
    assert (fundb / "microsporidia_odb10" / "dataset.cfg").is_file()
    assert not (fundb / "microsporidia_odb10.tar.gz").exists()
    row = setupDB.busco_info_row(str(fundb / "microsporidia_odb10"), md5)
    assert row == ("busco", str(fundb / "microsporidia_odb10"), "microsporidia_odb10", "2024-01-08", 600, md5)


def test_install_rejects_a_bad_checksum_and_keeps_nothing(tmp_path):
    src = tmp_path / "server.tar.gz"
    make_archive(str(src))
    index = {"microsporidia_odb10": ("2024-01-08", "0" * 32, "Eukaryota", "lineages")}
    fundb = tmp_path / "db"
    fundb.mkdir()

    def fetch(url, path):
        with open(src, "rb") as a, open(path, "wb") as b:
            b.write(a.read())

    with pytest.raises(RuntimeError, match="md5"):
        setupDB.install_odb10_lineage("microsporidia_odb10", str(fundb), index, fetch, logging.getLogger("t"))
    assert os.listdir(str(fundb)) == []


def test_install_rejects_a_path_traversal_member(tmp_path):
    src = tmp_path / "evil.tar.gz"
    make_archive(str(src), extra={"../escape.txt": b"x"})
    index = {"microsporidia_odb10": ("2024-01-08", md5_of(str(src)), "Eukaryota", "lineages")}
    fundb = tmp_path / "db"
    fundb.mkdir()

    def fetch(url, path):
        with open(src, "rb") as a, open(path, "wb") as b:
            b.write(a.read())

    with pytest.raises(RuntimeError, match="unsafe path"):
        setupDB.install_odb10_lineage("microsporidia_odb10", str(fundb), index, fetch, logging.getLogger("t"))
    assert not (tmp_path / "escape.txt").exists()


def test_install_archive_with_the_wrong_folder(tmp_path):
    src = tmp_path / "server.tar.gz"
    make_archive(str(src), folder="something_else")
    index = {"microsporidia_odb10": ("2024-01-08", md5_of(str(src)), "Eukaryota", "lineages")}
    fundb = tmp_path / "db"
    fundb.mkdir()

    def fetch(url, path):
        with open(src, "rb") as a, open(path, "wb") as b:
            b.write(a.read())

    with pytest.raises(RuntimeError, match='did not contain the folder "microsporidia_odb10"'):
        setupDB.install_odb10_lineage("microsporidia_odb10", str(fundb), index, fetch, logging.getLogger("t"))


def test_calcmd5_matches_hashlib_for_a_large_file(tmp_path):
    f = tmp_path / "big.bin"
    f.write_bytes(os.urandom(3 * (1 << 20) + 17))
    assert setupDB.calcmd5(str(f)) == md5_of(str(f))


def test_setup_rejects_an_odb12_name_before_doing_anything(capsys):
    with pytest.raises(SystemExit) as e:
        setupDB.main(["-i", "busco", "-b", "fungi_odb12", "-d", "/nonexistent/should/not/be/created"])
    assert e.value.code == 2
    assert "lengths_cutoff" in capsys.readouterr().err
    assert not os.path.exists("/nonexistent")


def test_listed_odb10_lineages_are_valid_names():
    for name, count in resources.buscoOdb10:
        assert setupDB.busco_name_problem(name) is None
        assert count > 0
