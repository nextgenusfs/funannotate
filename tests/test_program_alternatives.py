"""Dependency checks accept any installed name of a program.

bioconda's fasta3 installs only fasta36, and macOS arm64 has no tbl2asn, only
NCBI's table2asn. A requirement for "fasta" or "tbl2asn" must be met by any of
the alternative names on PATH, and a missing program must still stop the run.
"""
import os
import stat

import pytest

import funannotate.library as lib


def make_program(directory, name, output="version: 36.3.8\n"):
    path = os.path.join(str(directory), name)
    with open(path, "w") as fh:
        fh.write("#!/bin/sh\necho '%s'\n" % output.strip())
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


@pytest.fixture
def only_path(tmp_path, monkeypatch):
    """PATH holds only the temporary directory."""
    monkeypatch.setenv("PATH", str(tmp_path))
    return tmp_path


@pytest.mark.parametrize("installed", ["fasta", "fasta36", "fasta3"])
def test_fasta_requirement_met_by_any_alternative(only_path, installed):
    make_program(only_path, installed)
    name, path = lib.resolve_program("fasta")
    assert name == installed
    assert path.endswith(installed)


@pytest.mark.parametrize("installed", ["table2asn", "tbl2asn"])
def test_tbl2asn_requirement_met_by_either_tool(only_path, installed):
    make_program(only_path, installed, "%s 1.0 arguments:" % installed)
    name, _ = lib.resolve_program("tbl2asn")
    assert name == installed


def test_table2asn_is_preferred_when_both_are_installed(only_path):
    make_program(only_path, "tbl2asn")
    make_program(only_path, "table2asn")
    assert lib.resolve_program("tbl2asn")[0] == "table2asn"


def test_fasta_prefers_fasta_then_fasta36(only_path):
    make_program(only_path, "fasta3")
    make_program(only_path, "fasta36")
    assert lib.resolve_program("fasta")[0] == "fasta36"
    make_program(only_path, "fasta")
    assert lib.resolve_program("fasta")[0] == "fasta"


def test_program_without_alternatives_is_resolved_by_its_own_name(only_path):
    make_program(only_path, "diamond")
    assert lib.resolve_program("diamond")[0] == "diamond"
    assert lib.resolve_program("hisat2") == (None, None)


def test_missing_program_with_alternatives_returns_none(only_path):
    assert lib.resolve_program("fasta") == (None, None)
    assert lib.resolve_program("tbl2asn") == (None, None)


def test_check_dependencies_passes_with_only_the_alternatives(only_path):
    make_program(only_path, "fasta36")
    make_program(only_path, "table2asn", "table2asn 1.0 arguments:")
    # No SystemExit: both requirements are met.
    lib.CheckDependencies(["fasta", "tbl2asn"])


def test_check_dependencies_still_stops_when_nothing_is_installed(only_path, monkeypatch):
    messages = []

    class Recorder(object):
        def error(self, message):
            messages.append(message)

    # CheckDependencies logs through lib.log when it exists, and to stderr
    # otherwise; record the message so the test does not depend on which.
    monkeypatch.setattr(lib, "log", Recorder(), raising=False)
    with pytest.raises(SystemExit):
        lib.CheckDependencies(["fasta", "tbl2asn"])
    text = " ".join(messages)
    assert "fasta/fasta36/fasta3" in text
    assert "table2asn/tbl2asn" in text
