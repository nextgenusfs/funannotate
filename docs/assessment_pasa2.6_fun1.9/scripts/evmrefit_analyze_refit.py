#!/usr/bin/python3.12
"""Summarize EVM refit scores: each weight set minus the rerun baseline ('base'), per genome/arm,
for intron-chain F1 (ic_f1) and exact-CDS F1 (ex_f1). Also shows each input source alone
(GeneMark / Augustus+HiQ / pasa) from the original gene_predictions.gff3 for reference.
Usage: analyze_refit.py [GENOME ...]   (default: every genome under runs/)"""
import csv
import glob
import os
import sys

E = os.path.dirname(os.path.abspath(__file__))
W = {r["set"]: r for r in csv.DictReader(open(f"{E}/weight_sets.tsv"), delimiter="\t")}
rows = []
for f in sorted(glob.glob(f"{E}/runs/*/*/scores.tsv")):
    arm = os.path.basename(os.path.dirname(f))
    g = os.path.basename(os.path.dirname(os.path.dirname(f)))
    if len(sys.argv) > 1 and g not in sys.argv[1:]:
        continue
    for r in csv.DictReader(open(f), delimiter="\t"):
        r["arm"] = arm
        r["g"] = g
        rows.append(r)
runs = sorted({(r["g"], r["arm"]) for r in rows})
by = {(r["g"], r["arm"], r["label"]): r for r in rows}
sets = [s for s in W if all((g, a, s) in by for g, a in runs)]
print("delta vs base (rerun with current weights); columns = genome/arm; ic = intron-chain F1, ex = exact-CDS F1")
hdr = ["set", "weights A/H/G/S/P/Pr/T"] + [f"{g.split('_')[0][:4]}.{g.split('_')[1][:4]}.{a.split('.')[0][:5]}" for g, a in runs]
print("\t".join(hdr + ["mean_ic", "mean_ex", "min_ic"]))
out = []
for s in sets:
    w = W[s]
    wt = "/".join(w[k] for k in ("Augustus", "HiQ", "GeneMark", "snap", "pasa", "proteins", "transcripts"))
    dic = [float(by[(g, a, s)]["ic_f1"]) - float(by[(g, a, "base")]["ic_f1"]) for g, a in runs]
    dex = [float(by[(g, a, s)]["ex_f1"]) - float(by[(g, a, "base")]["ex_f1"]) for g, a in runs]
    cells = [f"{i:+.2f}/{x:+.2f}" for i, x in zip(dic, dex)]
    out.append((sum(dic) / len(dic), s, wt, cells, sum(dex) / len(dex), min(dic)))
for m, s, wt, cells, mx, mn in sorted(out, reverse=True):
    print("\t".join([s, wt] + cells + [f"{m:+.2f}", f"{mx:+.2f}", f"{mn:+.2f}"]))
print("\nbase absolute (ic_f1 / ex_f1):", "  ".join(
    f"{g.split('_')[0][:4]}.{a.split('.')[0]}: {by[(g, a, 'base')]['ic_f1']}/{by[(g, a, 'base')]['ex_f1']}" for g, a in runs))
