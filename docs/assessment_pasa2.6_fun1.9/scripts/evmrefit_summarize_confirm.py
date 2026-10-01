#!/usr/bin/python3.12
"""Full-predict confirmation (D127): base vs new EVM weights, final models and evm.round1,
holdout gene-level intron-chain (ic) and exact-CDS (ex) F1, plus gffcompare locus F1."""
import csv, glob, os, statistics as st
C = os.path.dirname(os.path.abspath(__file__))
out = ["genome\tweights\tn_reps\tpred_genes\tfinal_ic\tfinal_ex\tevm_ic\tgffcmp_locus_f1"]
res = {}
for d in sorted(glob.glob(f"{C}/*/*_r*.B")):
    g = os.path.basename(os.path.dirname(d)); w = os.path.basename(d).split("_r")[0]
    gs = f"{d}/gene_score.tsv"; sc = f"{d}/score.tsv"
    if not os.path.exists(gs):
        continue
    r = {x["label"]: x for x in csv.DictReader(open(gs), delimiter="\t")}
    s = next(csv.DictReader(open(sc), delimiter="\t"))
    lf = 2 * float(s["locus_sn"]) * float(s["locus_pr"]) / (float(s["locus_sn"]) + float(s["locus_pr"]))
    res.setdefault((g, w), []).append((int(r["final"]["pred_genes"]), float(r["final"]["ic_f1"]),
                                       float(r["final"]["ex_f1"]), float(r["evm_round1"]["ic_f1"]), lf))
for (g, w), v in sorted(res.items()):
    m = [st.mean(x[i] for x in v) for i in range(5)]
    out.append(f"{g}\t{w}\t{len(v)}\t{m[0]:.0f}\t{m[1]:.2f}\t{m[2]:.2f}\t{m[3]:.2f}\t{m[4]:.2f}")
out.append("\ndelta new - base (means)")
for g in sorted({g for g, _ in res}):
    if (g, "base") in res and (g, "new") in res:
        b = [st.mean(x[i] for x in res[(g, "base")]) for i in range(5)]
        n = [st.mean(x[i] for x in res[(g, "new")]) for i in range(5)]
        rng = lambda k, i: max(x[i] for x in res[(g, k)]) - min(x[i] for x in res[(g, k)])
        out.append(f"{g}\tfinal_ic {n[1]-b[1]:+.2f}\tfinal_ex {n[2]-b[2]:+.2f}\tlocus {n[4]-b[4]:+.2f}\tgenes {n[0]-b[0]:+.0f}\trep range ic base {rng('base',1):.2f} new {rng('new',1):.2f}")
txt = "\n".join(out); open(f"{C}/confirm_summary.txt", "w").write(txt + "\n"); print(txt)
