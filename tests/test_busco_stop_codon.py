"""#1078/#1071/#1069: BUSCO Augustus GTF stop_codon rows must not reach the EVM converter."""
from funannotate import library as lib

GTF = (
    "ctg1\tAUGUSTUS\tgene\t283\t2340\t.\t+\t.\tg1\n"
    "ctg1\tAUGUSTUS\ttranscript\t283\t2340\t.\t+\t.\tg1.t1\n"
    "ctg1\tAUGUSTUS\tCDS\t283\t2340\t.\t+\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n"
    "ctg1\tAUGUSTUS\tstop_codon\t2338\t2340\t.\t+\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n"
    "ctg1\tAUGUSTUS\tstart_codon\t283\t285\t.\t+\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n"
    "# not a feature line\n"
)


def test_drops_only_stop_codon_rows(tmp_path):
    src, dst = tmp_path / "in.gff", tmp_path / "out.gff"
    src.write_text(GTF)
    assert lib.augustus_gtf_drop_stop_codon(str(src), str(dst)) == 1
    out = dst.read_text().splitlines()
    assert len(out) == 5
    assert not any("\tstop_codon\t" in line for line in out)
    assert any("\tCDS\t283\t2340\t" in line for line in out)
    assert any("\tstart_codon\t" in line for line in out)
    assert out[-1] == "# not a feature line"


def test_no_stop_codon_rows_is_identity(tmp_path):
    src, dst = tmp_path / "in.gff", tmp_path / "out.gff"
    src.write_text("a\tb\tCDS\t1\t3\t.\t+\t0\tx\n")
    assert lib.augustus_gtf_drop_stop_codon(str(src), str(dst)) == 0
    assert dst.read_text() == src.read_text()
