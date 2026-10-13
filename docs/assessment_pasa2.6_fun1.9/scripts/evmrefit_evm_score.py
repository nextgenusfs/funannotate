#!/usr/bin/python3.12
"""Score EVM outputs on experiment B holdout chromosomes (EVM refit, 2026-09-30).

Gene-level scores against RefSeq protein-coding genes on the holdout chromosomes:
  ic : intron-chain match (multi-exon: same introns, ends free, as gffcompare; single-exon:
       exact span)
  ex : exact CDS chain (start and stop must also match)
Sn = RefSeq genes with a matching chain / RefSeq genes; Pr = predicted genes with a matching
chain / predicted genes. Loads the RefSeq file once and scores many predictions.
Usage: evm_score.py --ref-gff3 REF --holdout CHROMS --genome G --out TSV  LABEL=PRED.gff3 ...
"""
import argparse
import os
import sys

sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_C")
from error_breakdown import genes_from, introns  # noqa: E402

sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/refseq_benchmark")
from predict_scorer import coding_cds  # noqa: E402


def keys(g, exact):
    for ch in g["tx"].values():
        yield (g["chrom"], g["strand"], ch if exact or len(ch) == 1 else introns(ch))


def score(ref, pred, exact):
    rk = {k for g in ref.values() for k in keys(g, exact)}
    pk = {k for g in pred.values() for k in keys(g, exact)}
    sn = sum(any(k in pk for k in keys(g, exact)) for g in ref.values()) / len(ref)
    pr = sum(any(k in rk for k in keys(g, exact)) for g in pred.values()) / max(len(pred), 1)
    f1 = 0 if sn + pr == 0 else 2 * sn * pr / (sn + pr)
    return 100 * sn, 100 * pr, 100 * f1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff3", required=True)
    ap.add_argument("--holdout", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("preds", nargs="+")
    a = ap.parse_args()
    hold = {l.strip() for l in open(a.holdout) if l.strip()}
    ref = genes_from({t: v for t, v in coding_cds(a.ref_gff3, True).items() if v[0] in hold})
    new = not os.path.exists(a.out) or os.path.getsize(a.out) == 0
    with open(a.out, "a") as o:
        if new:
            o.write("genome\tlabel\tref_genes\tpred_genes\tic_sn\tic_pr\tic_f1\tex_sn\tex_pr\tex_f1\n")
        for lp in a.preds:
            label, path = lp.split("=", 1)
            pred = genes_from({t: v for t, v in coding_cds(path, False).items() if v[0] in hold})
            ic = score(ref, pred, False)
            ex = score(ref, pred, True)
            o.write(f"{a.genome}\t{label}\t{len(ref)}\t{len(pred)}\t"
                    + "\t".join(f"{x:.2f}" for x in ic + ex) + "\n")


if __name__ == "__main__":
    main()
