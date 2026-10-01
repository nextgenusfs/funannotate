#!/usr/bin/python3.12
"""GeneMark-alone vs final models, intron-chain level (ends ignored for multi-exon genes;
single-exon genes need the exact span). Same mask and gene set as error_breakdown.py."""
import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from error_breakdown import genes_from, introns
from masked_score import mask_spans, overlaps
from source_accuracy import split_by_source
sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/refseq_benchmark")
from predict_scorer import chrom_lengths, coding_cds
ref_gff, run, genome, variant = sys.argv[1:5]; masks = sys.argv[5:]
chroms = set(chrom_lengths(ref_gff)); mask = mask_spans(masks)
def key(g, ch): return (g["chrom"], g["strand"], introns(ch) if len(ch) > 1 else ch)
def sc(g): return not overlaps(mask, g["chrom"], g["start"], g["end"])
ref = {k: v for k, v in genes_from({t: v for t, v in coding_cds(ref_gff, True).items() if v[0] in chroms}).items() if sc(v)}
rk = {key(g, ch) for g in ref.values() for ch in g["tx"].values()}
sets = split_by_source(f"{run}/out/predict_misc/gene_predictions.gff3", chroms)
sets["final"] = genes_from({t: v for t, v in coding_cds(glob.glob(f"{run}/out/predict_results/*.gff3")[0], False).items() if v[0] in chroms})
out = [genome, variant]
for s in ("GeneMark", "final"):
    gs = {k: v for k, v in sets[s].items() if sc(v)}
    ks = {key(g, ch) for g in gs.values() for ch in g["tx"].values()}
    sn = sum(any(key(g, ch) in ks for ch in g["tx"].values()) for g in ref.values()) / len(ref)
    pr = sum(any(key(g, ch) in rk for ch in g["tx"].values()) for g in gs.values()) / len(gs)
    out.append(f"{s} {100*sn:.1f}/{100*pr:.1f} F1 {200*sn*pr/(sn+pr):.1f}")
print("\t".join(out))
