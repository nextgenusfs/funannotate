#!/usr/bin/python3.12
"""For final-missed RefSeq genes: repeat (repeatmasker.bed) coverage of CDS, split by whether
any EVM input overlaps; and for genes in evm.round1 but not final, whether the gene's EVM model
is listed in bad_models.gff / repeat.gene.models.txt. Prints one line per run."""
import collections, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from error_breakdown import Index, gene_overlap, genes_from
from masked_score import mask_spans, overlaps
sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/refseq_benchmark")
from predict_scorer import chrom_lengths, coding_cds
ref_gff, run, genome, variant = sys.argv[1:5]; masks = sys.argv[5:]
chroms = set(chrom_lengths(ref_gff)); mask = mask_spans(masks); misc = f"{run}/out/predict_misc"
ref = genes_from({t: v for t, v in coding_cds(ref_gff, True).items() if v[0] in chroms})
final = genes_from({t: v for t, v in coding_cds(glob.glob(f"{run}/out/predict_results/*.gff3")[0], False).items() if v[0] in chroms})
inp = genes_from({t: v for t, v in coding_cds(f"{misc}/gene_predictions.gff3", False).items() if v[0] in chroms})
evm = genes_from({t: v for t, v in coding_cds(f"{misc}/evm.round1.gff3", False).items() if v[0] in chroms})
rep = collections.defaultdict(list)
for l in open(f"{misc}/repeatmasker.bed"):
    c = l.split("\t")
    if c[0] in chroms: rep[c[0]].append((int(c[1]) + 1, int(c[2])))
for k in rep: rep[k].sort()
def rep_frac(g):
    cds = max(g["tx"].values(), key=len); tot = sum(e - s + 1 for s, e in cds); cov = 0
    for s, e in cds:
        for a, b in rep.get(g["chrom"], ()):
            if a > e: break
            if b >= s: cov += min(b, e) - max(a, s) + 1
    return min(cov / tot, 1.0)
bad = set()
for f in ("bad_models.gff", "repeat.gene.models.txt"):
    p = f"{misc}/{f}"
    if os.path.exists(p):
        for l in open(p):
            for tok in l.replace(";", "\t").replace("=", "\t").split():
                bad.add(tok)
fi, ii, ei = Index(final), Index(inp), Index(evm)
def hit(ix, genes, g): return [p for p in ix.hits(g["chrom"], g["strand"], g["start"], g["end"]) if gene_overlap(g, genes[p]) > 0]
c = collections.Counter()
for g in ref.values():
    if overlaps(mask, g["chrom"], g["start"], g["end"]) or hit(fi, final, g): continue
    k = "withinput" if hit(ii, inp, g) else "noinput"
    c[k] += 1; c[k + "_rep50"] += rep_frac(g) >= 0.5
    e = hit(ei, evm, g)
    if e:
        c["evm_removed"] += 1
        c["evm_removed_listed_bad"] += any(p in bad or any(t in bad for t in evm[p]["tx"]) for p in e)
        c["evm_removed_rep50"] += rep_frac(g) >= 0.5
print(genome, variant, " ".join(f"{k}={v}" for k, v in sorted(c.items())), sep="\t")
