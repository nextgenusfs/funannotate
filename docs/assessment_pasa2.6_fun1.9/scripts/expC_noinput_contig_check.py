#!/usr/bin/python3.12
"""Final-missed RefSeq genes with no EVM input overlap: contig length and contig-end distance."""
import collections, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from error_breakdown import Index, gene_overlap, genes_from
from masked_score import mask_spans, overlaps
sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/refseq_benchmark")
from predict_scorer import chrom_lengths, coding_cds
ref_gff, run, genome = sys.argv[1:4]; masks = sys.argv[4:]
L = {}
for l in open(f"{run}/out/predict_misc/genome.softmasked.fa.fai" if os.path.exists(f"{run}/out/predict_misc/genome.softmasked.fa.fai") else "/dev/null"):
    c = l.split("\t"); L[c[0]] = int(c[1])
if not L:
    name = None
    for l in open(f"{run}/out/predict_misc/genome.softmasked.fa"):
        if l[0] == ">": name = l[1:].split()[0]; L[name] = 0
        else: L[name] += len(l.strip())
chroms = set(chrom_lengths(ref_gff)); mask = mask_spans(masks)
ref = genes_from({t: v for t, v in coding_cds(ref_gff, True).items() if v[0] in chroms})
inp = genes_from({t: v for t, v in coding_cds(f"{run}/out/predict_misc/gene_predictions.gff3", False).items() if v[0] in chroms})
gm = {k: v for k, v in inp.items()}
ii = Index(inp)
sc = [g for g in ref.values() if not overlaps(mask, g["chrom"], g["start"], g["end"])]
def noinp(g): return not any(gene_overlap(g, inp[p]) > 0 for p in ii.hits(g["chrom"], g["strand"], g["start"], g["end"]))
# any-strand input overlap too
def anyinp(g):
    return any(gene_overlap(g, inp[p]) > 0 for st in "+-" for p in ii.hits(g["chrom"], st, g["start"], g["end"]))
bins = [(0, 50000), (50000, 500000), (500000, 10**10)]
tot = collections.Counter(); miss = collections.Counter(); anti = 0; n = 0
for g in sc:
    b = next(i for i, (a, z) in enumerate(bins) if a <= L.get(g["chrom"], 0) < z)
    tot[b] += 1
    if noinp(g):
        miss[b] += 1; n += 1
        if anyinp(g): anti += 1
print(genome, f"noinput={n} (opposite-strand input overlap {anti})",
      " ".join(f"contig{bins[b][0]//1000}-{bins[b][1]//1000 if bins[b][1]<10**10 else 'inf'}kb: {miss[b]}/{tot[b]}" for b in range(3)), sep="\t")
