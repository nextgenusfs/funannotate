"""Experiment B: does median read identity predict PASA vs BUSCO training outcome?

Joins read_identity.tsv with holdout scores (step B) per genome:
  pasa  = pasa.B locus F1
  busco = mean locus F1 of busco_r1..r3.B (sd reported)
  norna = norna.B locus F1 (BUSCO training, no RNA-seq evidence; 15 genomes)
  complete = complete-ORF PASA models on the train chromosomes (pasa.A pasa_gate record)
Writes per-genome table and a threshold sweep. Genome-level bootstrap (2000) for mean diffs.
"""
import csv, math, os, random, statistics, sys

E = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B"
OUT = os.path.join(E, "read_identity")


def f1(path):
    with open(path) as f:
        r = next(csv.DictReader(f, delimiter="\t"))
    sn, pr = float(r["locus_sn"]), float(r["locus_pr"])
    return 2 * sn * pr / (sn + pr) if sn + pr else 0.0


def complete_models(name):
    p = os.path.join(E, name, "pasa.A", "out", "logfiles", "training_decisions.tsv")
    if not os.path.isfile(p):
        return None
    for line in open(p):
        c = line.rstrip("\n").split("\t")
        if len(c) > 3 and c[1] == "pasa_gate":
            try:
                return int(c[3])
            except ValueError:
                pass
    return None


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
    rx, ry = rank(x), rank(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def boot_ci(vals, n=2000, seed=1):
    if len(vals) < 2:
        return (float("nan"), float("nan"))
    rng = random.Random(seed)
    m = sorted(statistics.mean(rng.choice(vals) for _ in vals) for _ in range(n))
    return m[int(0.025 * n)], m[int(0.975 * n) - 1]


rows = []
for g in csv.DictReader(open(os.path.join(OUT, "read_identity.tsv")), delimiter="\t"):
    if g["status"] != "ok" or not g["median_identity_pct"]:
        continue
    n = g["name"]
    sp = os.path.join(E, n, "pasa.B", "score.tsv")
    bs = [os.path.join(E, n, "busco_r%d.B" % i, "score.tsv") for i in (1, 2, 3)]
    if not os.path.isfile(sp) or not all(os.path.isfile(b) for b in bs):
        continue
    bf = [f1(b) for b in bs]
    npath = os.path.join(E, n, "norna.B", "score.tsv")
    rows.append(dict(
        name=n, stratum=g["stratum"], map_rate=float(g["map_rate_pct"]),
        identity=float(g["median_identity_pct"]), p10=float(g["p10_identity_pct"]),
        complete=complete_models(n), pasa=f1(sp), busco=statistics.mean(bf),
        busco_sd=statistics.stdev(bf), norna=f1(npath) if os.path.isfile(npath) else None))
for r in rows:
    r["diff"] = r["pasa"] - r["busco"]
    r["ev_gain"] = None if r["norna"] is None else r["busco"] - r["norna"]
rows.sort(key=lambda r: r["identity"])

with open(os.path.join(OUT, "identity_vs_training.tsv"), "w") as out:
    cols = ["name", "stratum", "map_rate", "identity", "p10", "complete", "pasa", "busco",
            "busco_sd", "diff", "norna", "ev_gain"]
    out.write("\t".join(cols) + "\n")
    for r in rows:
        out.write("\t".join("" if r[c] is None else (("%.2f" % r[c]) if isinstance(r[c], float) else str(r[c])) for c in cols) + "\n")

print("genomes analysed:", len(rows))
print("\nname\tidentity\tmap%\tcomplete\tPASA-BUSCO locusF1\tbusco_sd\tBUSCO-norna")
for r in rows:
    print("%s\t%.2f\t%.1f\t%s\t%+.2f\t%.2f\t%s" % (r["name"], r["identity"], r["map_rate"], r["complete"],
          r["diff"], r["busco_sd"], "" if r["ev_gain"] is None else "%+.2f" % r["ev_gain"]))

x = [r["identity"] for r in rows]; y = [r["diff"] for r in rows]
print("\nSpearman(identity, PASA-BUSCO), all: %.3f (n=%d)" % (spearman(x, y), len(rows)))
hi = [r for r in rows if r["complete"] is not None and r["complete"] >= 500]
print("Spearman, complete>=500 only: %.3f (n=%d)" % (spearman([r["identity"] for r in hi], [r["diff"] for r in hi]), len(hi)))
cc = [r for r in rows if r["complete"]]
print("Spearman(log complete, PASA-BUSCO): %.3f (n=%d)" % (spearman([math.log(r["complete"]) for r in cc], [r["diff"] for r in cc]), len(cc)))
print("Spearman(identity, log complete): %.3f" % spearman([r["identity"] for r in cc], [math.log(r["complete"]) for r in cc]))

print("\nThreshold sweep (identity gate; genomes with complete>=500, i.e. not already sent to BUSCO by the PASA gate)")
print("T\tn_below\tBUSCO_better_below\tmean(PASA-BUSCO)_below [95%CI]\tn_above\tmean_above [95%CI]\tnet_gain_if_gated")
for t in [x / 2.0 for x in range(180, 200)]:
    below = [r["diff"] for r in hi if r["identity"] < t]
    above = [r["diff"] for r in hi if r["identity"] >= t]
    if not below:
        continue
    lo, up = boot_ci(below); alo, aup = boot_ci(above)
    print("%.1f\t%d\t%d\t%+.2f [%+.2f, %+.2f]\t%d\t%s\t%+.2f" % (
        t, len(below), sum(d < 0 for d in below), statistics.mean(below), lo, up, len(above),
        ("%+.2f [%+.2f, %+.2f]" % (statistics.mean(above), alo, aup)) if above else "-",
        -sum(below) / len(hi)))

ev = [r for r in rows if r["ev_gain"] is not None]
if ev:
    g = [r["ev_gain"] for r in ev]; lo, up = boot_ci(g)
    print("\nRNA-seq as evidence with BUSCO training (busco - norna locus F1): mean %+.2f [%+.2f, %+.2f], n=%d, positive in %d"
          % (statistics.mean(g), lo, up, len(g), sum(v > 0 for v in g)))
    lowid = [r["ev_gain"] for r in ev if r["identity"] < 99]
    if lowid:
        print("  identity < 99%%: mean %+.2f, n=%d, positive in %d" % (statistics.mean(lowid), len(lowid), sum(v > 0 for v in lowid)))
