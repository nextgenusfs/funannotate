"""#1089: a locus tag prefix that contains '-T' (for example --name Al-Tr-26) was cut at the first '-T',
so every gene collapsed into one called `Al`. Only the trailing transcript suffix -T<n> is stripped."""
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from funannotate import library as lib


def make_record(tag, n):
    seq = "ATG" + "GCT" * 20 + "TAA"
    rec = SeqRecord(Seq(seq * n), id="ctg1")
    for i in range(n):
        off = i * len(seq)
        loc = FeatureLocation(off, off + len(seq), strand=1)
        locus = "%s_%06d" % (tag, i + 1)
        rec.features.append(SeqFeature(loc, type="gene", qualifiers={"locus_tag": [locus]}))
        rec.features.append(SeqFeature(loc, type="mRNA", qualifiers={
            "locus_tag": [locus], "transcript_id": ["ncbi:%s-T1_mrna" % locus]}))
        rec.features.append(SeqFeature(loc, type="CDS", qualifiers={
            "locus_tag": [locus], "protein_id": ["ncbi:%s-T1" % locus], "codon_start": ["1"],
            "translation": ["M" + "A" * 20]}))
    return rec


def collect(tag, n=3):
    rec = make_record(tag, n)
    genes = {}
    for f in rec.features:
        lib.gb_feature_add2dict(f, rec, genes)
    return genes


def test_plain_locus_tag():
    assert sorted(collect("FUN")) == ["FUN_000001", "FUN_000002", "FUN_000003"]


def test_locus_tag_prefix_containing_dash_T():
    assert sorted(collect("Al-Tr-26")) == ["Al-Tr-26_000001", "Al-Tr-26_000002", "Al-Tr-26_000003"]
