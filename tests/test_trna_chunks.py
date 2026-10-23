import os
import shutil
import stat
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from funannotate import trna_chunks as tc  # noqa: E402

HEAD = (
    "Sequence        \t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\t      \n"
    "Name            \ttRNA #\tBegin\tEnd \tType\tCodon\tBegin\tEnd\tScore\tNote\n"
    "--------        \t------\t-----\t------\t----\t-----\t-----\t----\t------\t----\n"
)


def write_fasta(path, n, seqlen=60):
    with open(path, "w") as fh:
        for i in range(n):
            fh.write(">ctg{:04d} desc\n{}\n{}\n".format(i, "ACGT" * (seqlen // 4), "ACGT" * 5))
    return path


def contigs(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def test_fasta_stats(tmp_path):
    f = write_fasta(str(tmp_path / "g.fa"), 5)
    assert tc.fasta_stats(f) == (5, 5 * 80)


@pytest.mark.parametrize("n,k", [(10, 4), (3, 8), (1, 4), (7, 1), (8, 4)])
def test_split_every_contig_once_no_empty(tmp_path, n, k):
    f = write_fasta(str(tmp_path / "g.fa"), n)
    paths = tc.split_fasta(f, k, str(tmp_path))
    assert len(paths) == min(k, n)
    seen = []
    sizes = []
    for p in paths:
        c = contigs(p)
        assert c, "empty chunk"
        seen += c
        sizes.append(len(c))
    assert sorted(seen) == contigs(f)
    assert max(sizes) - min(sizes) <= 1


def test_split_keeps_sequence_with_header(tmp_path):
    f = write_fasta(str(tmp_path / "g.fa"), 4)
    p = tc.split_fasta(f, 2, str(tmp_path))[0]
    lines = open(p).read().splitlines()
    assert lines[0].startswith(">ctg0000") and not lines[1].startswith(">")
    assert lines[3].startswith(">ctg0002")


def test_should_chunk(tmp_path):
    many = write_fasta(str(tmp_path / "m.fa"), 50)
    assert tc.should_chunk(many, 4, min_contigs=50)
    assert not tc.should_chunk(many, 1, min_contigs=50)
    assert not tc.should_chunk(many, 4, min_contigs=51)
    assert not tc.should_chunk(many, 4, min_contigs=50, max_mean=10)


def test_merge_one_header_skips_empty(tmp_path):
    a = tmp_path / "a.out"
    b = tmp_path / "b.out"
    c = tmp_path / "c.out"
    a.write_text(HEAD + "ctg1   \t1\t5 \t80 \tMet\tCAT\t0\t0\t70.1\t\n")
    b.write_text("")
    c.write_text(HEAD + "ctg2                \t1\t7 \t90 \tGly\tGCC\t0\t0\t60.0\t\n")
    out = tmp_path / "m.out"
    assert tc.merge_outputs([str(a), str(b), str(c), str(tmp_path / "none")], str(out)) == 2
    lines = out.read_text().splitlines()
    assert len(lines) == 5
    assert lines[0].startswith("Sequence") and sum(l.startswith("Name") for l in lines) == 1
    assert lines[3].startswith("ctg1") and lines[4].startswith("ctg2")


def test_merge_no_rows_gives_zero_byte_file(tmp_path):
    a = tmp_path / "a.out"
    a.write_text("")
    out = tmp_path / "m.out"
    assert tc.merge_outputs([str(a)], str(out)) == 0
    assert out.stat().st_size == 0


FAKE = """#!/usr/bin/env python3
import os, sys
a = sys.argv[1:]
out = a[a.index("-o") + 1]
assert "-Q" in a and a[a.index("--thread") + 1] == "1"
tmp = os.environ["TMPDIR"]
assert os.path.isdir(tmp), tmp
rows = []
for l in open(a[-1]):
    if l.startswith(">"):
        name = l[1:].split()[0]
        if int(name[3:]) % 3 == 0:
            rows.append(name + "\\t1\\t5\\t80\\tMet\\tCAT\\t0\\t0\\t70.1\\t\\n")
if rows:
    with open(out, "w") as o:
        o.write(@HEAD@)
        o.writelines(rows)
else:
    open(out, "w").close()
"""


def fake_exe(path):
    with open(path, "w") as fh:
        fh.write(FAKE.replace("@HEAD@", repr(HEAD)))
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
    return path


def parsed(path):
    rows = []
    for l in open(path):
        if l.startswith(tc.HEADER_STARTS):
            continue
        c = [x.strip() for x in l.split("\t")]
        rows.append((c[0], int(c[2]), int(c[3])))
    return sorted(rows)


def test_run_chunked_with_fake_exe(tmp_path):
    f = write_fasta(str(tmp_path / "g.fa"), 20)
    exe = fake_exe(str(tmp_path / "fake-trnascan"))
    out = str(tmp_path / "t.out")
    assert tc.run_chunked(f, out, 4, str(tmp_path), exe=exe)
    assert [r[0] for r in parsed(out)] == ["ctg{:04d}".format(i) for i in range(0, 20, 3)]
    assert not [d for d in os.listdir(str(tmp_path)) if d.startswith("trnascan_chunks_")]


def test_run_chunked_failure_returns_false(tmp_path):
    f = write_fasta(str(tmp_path / "g.fa"), 6)
    exe = str(tmp_path / "bad")
    with open(exe, "w") as fh:
        fh.write("#!/bin/sh\nexit 3\n")
    os.chmod(exe, 0o755)
    out = str(tmp_path / "t.out")
    assert not tc.run_chunked(f, out, 2, str(tmp_path), exe=exe)
    assert not os.path.exists(out)


@pytest.mark.skipif(not shutil.which("tRNAscan-SE"), reason="tRNAscan-SE not installed")
def test_chunked_equals_whole_real(tmp_path):
    f = str(tmp_path / "g.fa")
    sys.path.insert(0, os.path.dirname(__file__))
    # a tRNA-Met-like sequence between filler contigs
    trna = ("GCCTCGTTAGCGCAGTAGGTAGCGCGTCAGTCTCATAATCTGAAGGTCGTGAGTTCGATCCTCACACGGGGCA")
    with open(f, "w") as fh:
        for i in range(9):
            fh.write(">ctg{:04d}\n{}\n".format(i, ("ACGTTGCA" * 20) + (trna if i % 2 == 0 else "") + "TTGACA" * 10))
    whole = str(tmp_path / "whole.out")
    subprocess.check_call(["tRNAscan-SE", "-q", "-Q", "-o", whole, "--thread", "1", f])
    chunked = str(tmp_path / "chunked.out")
    assert tc.run_chunked(f, chunked, 3, str(tmp_path))
    assert parsed(whole) == parsed(chunked)
