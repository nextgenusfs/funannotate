#!/usr/bin/env python3
"""Experiment A analysis (DECISIONS D81): calibrate the complete-model gate.

Input: titration_scores.tsv, one row per genome x N x draw, with columns
  genome, N, draw, n_complete_used, n_keepers,
  locus_sn, locus_pr, exon_sn, exon_pr, intron_chain_sn, intron_chain_pr, busco_c
and one comparator row per genome with N == "busco" (busco-forced training,
same evidence).

For each genome and N:
  F1 = 2*Sn*Pr/(Sn+Pr) at locus level (also exon and intron-chain);
  delta = F1(PASA-trained, N) - F1(busco-forced);
  bootstrap 95% CI of mean delta by resampling draws (with replacement).
  With one draw (N > 500) the CI collapses to the point value.
Conservative N* (user decision D81): the smallest N such that the lower 95%
bound of delta is >= 0 at that N and at every larger N tested.
Also reports whether n_keepers or n_complete_used tracks delta better
(Spearman correlation across all PASA rows of a genome).

Limitation: the bootstrap covers subsample (draw) variation only. A bootstrap
over holdout genes needs per-gene match tables from gffcompare, which
predict_scorer.py does not keep; add them if the draw-level CIs are too narrow
to be credible.

Usage: titration_analysis.py titration_scores.tsv [--level locus|exon|intron_chain]
           [--boot 2000] [--seed 1] [--out summary.tsv]
"""
import argparse
import collections
import csv
import random
import sys


def f1(sn, pr):
    sn, pr = float(sn), float(pr)
    return 0.0 if sn + pr == 0 else 2 * sn * pr / (sn + pr)


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0 + 1
        i = j + 1
    return r


def spearman(x, y):
    if len(x) < 3:
        return float("nan")
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scores")
    ap.add_argument("--level", default="locus", choices=["locus", "exon", "intron_chain"])
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out")
    a = ap.parse_args()
    rnd = random.Random(a.seed)
    sn_col, pr_col = a.level + "_sn", a.level + "_pr"

    busco = collections.defaultdict(list)
    pasa = collections.defaultdict(lambda: collections.defaultdict(list))
    rows_by_genome = collections.defaultdict(list)
    with open(a.scores) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r[sn_col] in ("", "NA") or r[pr_col] in ("", "NA"):
                continue
            if r.get("status", "ok") not in ("ok", ""):
                continue
            # skip failed runs: a crashed step must not count as low accuracy
            if any(r.get(k, "0") not in ("0", "", None) for k in ("exit_stepA", "exit_stepB")):
                sys.stderr.write("skipping failed run %s N=%s draw=%s\n" % (r["genome"], r["N"], r["draw"]))
                continue
            val = f1(r[sn_col], r[pr_col])
            if r["N"] == "busco":
                busco[r["genome"]].append(val)
            elif not r["N"].isdigit():
                continue  # reference rows such as N="busco_code_new" are not comparators
            else:
                pasa[r["genome"]][int(r["N"])].append(val)
                rows_by_genome[r["genome"]].append(r)

    out = []
    header = ["genome", "N", "draws_pasa/busco", "f1_mean", "busco_f1", "delta_mean", "delta_lo95", "delta_hi95"]
    nstar = {}
    for g in sorted(pasa):
        if g not in busco:
            sys.stderr.write("no busco comparator for %s; skipped\n" % g)
            continue
        bvals = busco[g]
        b = sum(bvals) / len(bvals)
        per_n = []
        for n in sorted(pasa[g]):
            vals = pasa[g][n]
            mean = sum(vals) / len(vals) - b
            if len(vals) > 1 or len(bvals) > 1:
                # resample PASA draws and BUSCO repeats independently, so the
                # comparator's run-to-run noise is part of the interval
                boots = sorted(
                    sum(rnd.choice(vals) for _ in vals) / len(vals)
                    - sum(rnd.choice(bvals) for _ in bvals) / len(bvals)
                    for _ in range(a.boot))
                lo, hi = boots[int(0.025 * a.boot)], boots[int(0.975 * a.boot) - 1]
            else:
                lo = hi = mean
            d = vals
            per_n.append((n, lo))
            out.append([g, n, "%d/%d" % (len(d), len(bvals)), "%.3f" % (sum(vals) / len(vals)), "%.3f" % b,
                        "%.3f" % mean, "%.3f" % lo, "%.3f" % hi])
        # conservative N*: lower bound >= 0 at this N and all larger N
        cand = None
        for i, (n, lo) in enumerate(per_n):
            if all(l >= 0 for _, l in per_n[i:]):
                cand = n
                break
        nstar[g] = cand
        # which count tracks accuracy better
        rs = rows_by_genome[g]
        y = [f1(r[sn_col], r[pr_col]) for r in rs]
        rho_c = spearman([float(r["n_complete_used"]) for r in rs], y)
        rho_k = spearman([float(r["n_keepers"]) for r in rs if r["n_keepers"] not in ("", "NA")],
                         [f1(r[sn_col], r[pr_col]) for r in rs if r["n_keepers"] not in ("", "NA")])
        out.append([g, "N*", "", "", "", "", "", "conservative_N*=%s rho_complete=%.2f rho_keepers=%.2f"
                    % (cand if cand is not None else "none", rho_c, rho_k)])

    text = "\t".join(header) + "\n" + "\n".join("\t".join(str(x) for x in r) for r in out) + "\n"
    if a.out:
        open(a.out, "w").write(text)
    sys.stdout.write(text)
    vals = [v for v in nstar.values() if v is not None]
    sys.stdout.write("\nper-genome conservative N*: %s\n" % nstar)
    if vals:
        sys.stdout.write("proposed threshold (max over genomes, most conservative): %d\n" % max(vals))


if __name__ == "__main__":
    main()
