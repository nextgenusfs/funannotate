"""Experiment B (D81): complete-model threshold for PASA vs BUSCO training, across 40 genomes.

Per genome: PASA - BUSCO holdout F1 (locus, exon, intron chain), BUSCO = mean of 3 repeats.
Predictors: complete-ORF PASA models on the train chromosomes (= what the gate counts in step A),
and final PASA training models after selection.
For each threshold T, policy "PASA if complete >= T else BUSCO":
  gain vs always-PASA and vs always-BUSCO (mean over genomes, bootstrap 95% CI over genomes);
  mean PASA - BUSCO among genomes >= T (bootstrap CI); conservative T* = smallest T where that
  lower bound >= 0 for T and every larger T with >= 5 genomes.
"""
import csv, math, os, random, statistics

E = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B"
OUT = os.path.join(E, "read_identity")
LEVELS = [("locus", "locus_sn", "locus_pr"), ("exon", "exon_sn", "exon_pr"),
          ("intron_chain", "intron_chain_sn", "intron_chain_pr")]
TS = [100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000, 2500, 3000]
B = 2000


def scores(path):
    r = next(csv.DictReader(open(path), delimiter="\t"))
    out = {}
    for lev, a, b in LEVELS:
        sn, pr = float(r[a]), float(r[b])
        out[lev] = 2 * sn * pr / (sn + pr) if sn + pr else 0.0
    return out


def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rank(x), rank(y); mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den


def ci(fn, rows, seed=1):
    rng = random.Random(seed)
    vals = sorted(fn([rng.choice(rows) for _ in rows]) for _ in range(B))
    vals = [v for v in vals if v is not None]
    return vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1]


cm = {r["name"]: r for r in csv.DictReader(open(os.path.join(OUT, "complete_models.tsv")), delimiter="\t")}
rows = []
for g in csv.DictReader(open(os.path.join(E, "genomes.tsv")), delimiter="\t"):
    n = g["name"]
    sp = os.path.join(E, n, "pasa.B", "score.tsv")
    bs = [os.path.join(E, n, "busco_r%d.B" % i, "score.tsv") for i in (1, 2, 3)]
    if n not in cm or not os.path.isfile(sp) or not all(os.path.isfile(b) for b in bs):
        continue
    p = scores(sp); bl = [scores(b) for b in bs]
    row = dict(name=n, complete=int(cm[n]["complete_train"]), complete_genome=int(cm[n]["complete_genome"]),
               final=int(cm[n]["final_training_models"]))
    for lev, _, _ in LEVELS:
        row["pasa_" + lev] = p[lev]; row["busco_" + lev] = statistics.mean(b[lev] for b in bl)
        row["diff_" + lev] = p[lev] - row["busco_" + lev]
    rows.append(row)
rows.sort(key=lambda r: r["complete"])

with open(os.path.join(OUT, "complete_threshold_by_genome.tsv"), "w") as f:
    cols = ["name", "complete", "complete_genome", "final"] + [k + "_" + l for l, _, _ in LEVELS for k in ("pasa", "busco", "diff")]
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(("%.2f" % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

print("genomes:", len(rows))
for lev, _, _ in LEVELS:
    d = [r["diff_" + lev] for r in rows]
    print("\n=== %s F1 ===" % lev)
    print("Spearman(log complete_train, diff) %.3f; (log final, diff) %.3f; (log complete_genome, diff) %.3f" % (
        spearman([math.log(r["complete"]) for r in rows], d), spearman([math.log(r["final"]) for r in rows], d),
        spearman([math.log(r["complete_genome"]) for r in rows], d)))
    mp = statistics.mean(r["pasa_" + lev] for r in rows); mb = statistics.mean(r["busco_" + lev] for r in rows)
    print("mean F1 always-PASA %.2f, always-BUSCO %.2f" % (mp, mb))
    print("T\tn>=T\tpos>=T\tmean diff >=T [95%CI]\tn<T\tmean diff <T\tgain vs PASA [CI]\tgain vs BUSCO [CI]")
    lows = {}
    for t in TS:
        ab = [r for r in rows if r["complete"] >= t]; be = [r for r in rows if r["complete"] < t]
        def mean_ab(s, t=t):
            x = [r["diff_" + lev] for r in s if r["complete"] >= t]
            return statistics.mean(x) if x else None
        def gain_p(s, t=t):
            return statistics.mean((0.0 if r["complete"] >= t else -r["diff_" + lev]) for r in s)
        def gain_b(s, t=t):
            return statistics.mean((r["diff_" + lev] if r["complete"] >= t else 0.0) for r in s)
        m = mean_ab(rows); lo, hi = ci(mean_ab, rows) if len(ab) >= 2 else (float("nan"), float("nan"))
        lows[t] = (lo, len(ab))
        gp = gain_p(rows); gpl, gph = ci(gain_p, rows); gb = gain_b(rows); gbl, gbh = ci(gain_b, rows)
        mbe = statistics.mean(r["diff_" + lev] for r in be) if be else float("nan")
        print("%d\t%d\t%d\t%+.2f [%+.2f, %+.2f]\t%d\t%+.2f\t%+.2f [%+.2f, %+.2f]\t%+.2f [%+.2f, %+.2f]" % (
            t, len(ab), sum(r["diff_" + lev] > 0 for r in ab), m, lo, hi, len(be), mbe, gp, gpl, gph, gb, gbl, gbh))
    star = None
    for t in TS:
        if all(lows[u][0] >= 0 for u in TS if u >= t and lows[u][1] >= 5):
            star = t; break
    print("conservative T* (lower CI of mean diff among >=T is >= 0 here and above, n>=5):", star)
