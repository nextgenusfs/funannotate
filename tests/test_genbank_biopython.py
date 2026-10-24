"""The GenBank feature parser must work on current Biopython (nofuzzy_start/nofuzzy_end were removed in 1.87)."""
import io

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from funannotate import library as lib


def make_record():
    # gene 1-30 (+), two exons 1-12 and 19-30, CDS ATG AAA GTT | TAA GCC CTG ... in frame over both exons
    seq = Seq("ATGAAAGTTTAAGGGATGCTGCCCTGA" + "A" * 20)
    rec = SeqRecord(seq, id="ctg1", name="ctg1", description="t", annotations={"molecule_type": "DNA"})
    q = {"locus_tag": ["GN_000001"]}
    rec.features.append(SeqFeature(FeatureLocation(0, 27, strand=1), type="gene", qualifiers=dict(q)))
    loc = CompoundLocation([FeatureLocation(0, 9, strand=1), FeatureLocation(12, 27, strand=1)])
    rec.features.append(SeqFeature(loc, type="mRNA", qualifiers=dict(q, product=["hypothetical protein"], transcript_id=["GN_000001-T1"])))
    rec.features.append(SeqFeature(loc, type="CDS", qualifiers=dict(q, product=["hypothetical protein"], protein_id=["GN_000001-T1"], codon_start=["1"])))
    return rec


def test_gb_feature_add2dict_reads_exons_with_current_biopython():
    handle = io.StringIO()
    SeqIO.write(make_record(), handle, "genbank")
    handle.seek(0)
    genes = {}
    for rec in SeqIO.parse(handle, "genbank"):
        for f in rec.features:
            lib.gb_feature_add2dict(f, rec, genes)
    gene = genes["GN_000001"]
    assert gene["type"] == "mRNA"
    assert gene["location"] == (1, 27)
    assert gene["mRNA"][0] == [(1, 9), (13, 27)]
    assert gene["strand"] == "+"
