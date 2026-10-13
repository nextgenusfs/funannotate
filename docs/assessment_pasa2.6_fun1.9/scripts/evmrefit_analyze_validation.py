#!/usr/bin/python3.12
"""EVM refit validation on held-out experiment B genomes (not used to choose the weights).
For each finalist weight set: per-genome delta vs 'base' (rerun with current weights), for
intron-chain F1 (ic) and exact-CDS F1 (ex), per arm. Reports the mean with a 95% bootstrap CI
over genomes, the median, the number of genomes improved / worsened by > 0.1 (EVM rerun noise
is about +/-0.04), and the same split by yeast (Saccharomycotina) vs other fungi.
Writes validation_deltas.tsv and validation_summary.txt."""
import csv
import os
import random
import statistics as st

E = os.path.dirname(os.path.abspath(__file__))
genomes = [l.strip() for l in open(f"{E}/validation_genomes.txt") if l.strip()]
sets = [l.strip() for l in open(f"{E}/finalists.txt") if l.strip() and l.strip() != "base"]
YEAST = {"Williopsis", "Metschnikowia", "Brettanomyces", "Hyphopichia", "Clavispora", "Naumovozyma",
         "Kluyveromyces", "Kuraishia", "Henningerozyma", "Saccharomyces", "Scheffersomyces", "Debaryomyces"}
rows = []
for g in genomes:
    for arm in ("pasa.B", "busco_r1.B"):
        f = f"{E}/runs/{g}/{arm}/scores.tsv"
        if not os.path.exists(f):
            continue
        sc = {r["label"]: r for r in csv.DictReader(open(f), delimiter="\t")}
        if "base" not in sc:
            continue
        for s in sets:
            if s in sc:
                rows.append(dict(genome=g, arm=arm, set=s, yeast=g.split("_")[0] in YEAST,
                                 base_ic=float(sc["base"]["ic_f1"]), base_ex=float(sc["base"]["ex_f1"]),
                                 d_ic=float(sc[s]["ic_f1"]) - float(sc["base"]["ic_f1"]),
                                 d_ex=float(sc[s]["ex_f1"]) - float(sc["base"]["ex_f1"])))
with open(f"{E}/validation_deltas.tsv", "w", newline="") as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)


def boot(x, n=5000, seed=1):
    r = random.Random(seed)
    m = sorted(st.mean(r.choices(x, k=len(x))) for _ in range(n))
    return m[int(0.025 * n)], m[int(0.975 * n)]


out = []
p = out.append
p("Held-out validation: delta vs current weights (A1/H2/G1/S1/P6/Pr1/T1); n = genomes")
for arm in ("pasa.B", "busco_r1.B"):
    for grp, sel in (("all", lambda r: True), ("yeast", lambda r: r["yeast"]), ("non-yeast", lambda r: not r["yeast"])):
        p(f"\n== arm {arm}, {grp} ==")
        p("set\tn\tmetric\tmean [95% CI]\tmedian\tbetter>0.1\tworse>0.1\tmin\tmax")
        for s in sets:
            rr = [r for r in rows if r["arm"] == arm and r["set"] == s and sel(r)]
            if len(rr) < 2:
                continue
            for k in ("d_ic", "d_ex"):
                x = [r[k] for r in rr]
                lo, hi = boot(x)
                p(f"{s}\t{len(x)}\t{k[2:]}\t{st.mean(x):+.2f} [{lo:+.2f}, {hi:+.2f}]\t{st.median(x):+.2f}\t"
                  f"{sum(v > 0.1 for v in x)}\t{sum(v < -0.1 for v in x)}\t{min(x):+.2f}\t{max(x):+.2f}")
txt = "\n".join(out)
open(f"{E}/validation_summary.txt", "w").write(txt + "\n")
print(txt)
