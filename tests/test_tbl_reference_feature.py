"""The REFERENCE feature is written for tbl2asn but not for table2asn.

NCBI's table2asn rejects the sequencing-center REFERENCE feature and exits with
code 2, which hid real failures behind a routine one. The legacy tbl2asn keeps its
historic feature table.
"""
import pytest

import funannotate.library as lib


def write_table(tmp_path):
    out = tmp_path / "genome.tbl"
    lib.dicts2tbl({}, {"scaffold_1": []}, {"scaffold_1": 1000}, "CFMR", "12345", [], str(out))
    return out.read_text()


def use_dialect(monkeypatch, dialect):
    monkeypatch.setattr(lib, "resolve_tbl2asn_binary", lambda: (dialect, dialect))


def test_table2asn_gets_no_reference_feature(tmp_path, monkeypatch):
    use_dialect(monkeypatch, "table2asn")
    text = write_table(tmp_path)
    assert ">Feature scaffold_1" in text
    assert "REFERENCE" not in text
    assert "CFMR" not in text


def test_tbl2asn_keeps_the_reference_feature(tmp_path, monkeypatch):
    use_dialect(monkeypatch, "tbl2asn")
    text = write_table(tmp_path)
    assert "1\t1000\tREFERENCE\n" in text
    assert "\t\t\tCFMR\t12345\n" in text


def test_no_tool_found_keeps_the_historic_output(tmp_path, monkeypatch):
    def missing():
        raise RuntimeError("no tbl2asn or table2asn")

    monkeypatch.setattr(lib, "resolve_tbl2asn_binary", missing)
    assert lib.tbl_writes_reference_feature() is True
    assert "REFERENCE" in write_table(tmp_path)
