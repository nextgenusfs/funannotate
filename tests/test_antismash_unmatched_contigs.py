"""#1147: antiSMASH records renamed away from the genome's contig names matched no gene, silently."""
import logging

from funannotate import library as lib


def test_unmatched_contig_names_are_warned(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(lib, "log", logging.getLogger("t"), raising=False)
    bed = tmp_path / "clusters.bed"
    bed.write_text("renamed_contig_1\t100\t900\tCluster_1\n")
    monkeypatch.setattr(lib, "gff2dict", lambda *a, **k: {
        "g1": {"type": "mRNA", "location": (150, 500), "contig": "contig_1", "ids": ["g1-T1"]}})
    out = tmp_path / "clusters.txt"
    with caplog.at_level(logging.WARNING):
        result = lib.GetClusterGenes(str(bed), "x.gff3", "x.fa", str(out))
    assert result == {}
    assert "renamed_contig_1" in caplog.text and "--allow-long-headers" in caplog.text


def test_matching_names_give_clusters_and_no_warning(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(lib, "log", logging.getLogger("t"), raising=False)
    bed = tmp_path / "clusters.bed"
    bed.write_text("contig_1\t100\t900\tCluster_1\n")
    monkeypatch.setattr(lib, "gff2dict", lambda *a, **k: {
        "g1": {"type": "mRNA", "location": (150, 500), "contig": "contig_1", "ids": ["g1-T1"]}})
    out = tmp_path / "clusters.txt"
    with caplog.at_level(logging.WARNING):
        result = lib.GetClusterGenes(str(bed), "x.gff3", "x.fa", str(out))
    assert result == {"Cluster_1": ["g1-T1"]}
    assert "not found" not in caplog.text
