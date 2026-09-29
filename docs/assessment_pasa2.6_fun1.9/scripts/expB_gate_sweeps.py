"""Experiment B: final PASA training models as a gate, alone and combined with complete>=500.

Reads complete_threshold_by_genome.tsv (from analyze_complete_threshold.py).
Writes final_models_sweep.txt and combined_gate_sweep.txt next to it.
"""
import csv, os, random, statistics

OUT = os.path.dirname(os.path.abspath(__file__))
R = list(csv.DictReader(open(os.path.join(OUT, "complete_threshold_by_genome.tsv")), delimiter="\t"))
for r in R:
    r["final"] = int(r["final"]); r["complete"] = int(r["complete"]); r["d"] = float(r["diff_locus"])


def ci(fn, rows, B=2000, seed=1):
    rng = random.Random(seed); v = sorted(fn([rng.choice(rows) for _ in rows]) for _ in range(B)); return v[50], v[1949]


with open(os.path.join(OUT, "final_models_sweep.txt"), "w") as f:
    f.write("Final PASA training models (after selection) as gate variable, locus F1, n=%d\n" % len(R))
    f.write("T\tn>=T\tpos>=T\tmean diff >=T [CI]\tgain vs always-PASA [CI]\n")
    for t in [150, 200, 250, 275, 300, 350, 400, 450, 500, 600, 700, 800, 1000]:
        ab = [r for r in R if r["final"] >= t]
        def m(s, t=t):
            x = [r["d"] for r in s if r["final"] >= t]; return statistics.mean(x) if x else 0.0
        def g(s, t=t):
            return statistics.mean(0.0 if r["final"] >= t else -r["d"] for r in s)
        lo, hi = ci(m, R); gl, gh = ci(g, R)
        f.write("%d\t%d\t%d\t%+.2f [%+.2f, %+.2f]\t%+.2f [%+.2f, %+.2f]\n" % (
            t, len(ab), sum(r["d"] > 0 for r in ab), m(R), lo, hi, g(R), gl, gh))
    f.write("\nGenomes with final < 300:\n")
    for r in sorted(R, key=lambda r: r["final"]):
        if r["final"] < 300:
            f.write("%s %d %d %+.2f\n" % (r["name"], r["complete"], r["final"], r["d"]))

with open(os.path.join(OUT, "combined_gate_sweep.txt"), "w") as f:
    f.write("Combined policy: PASA if complete>=500 AND final>=F, else BUSCO. Gain in mean F1 over "
            "always-PASA; and over the complete>=500 gate alone. n=%d\n" % len(R))
    for lev in ("locus", "exon", "intron_chain"):
        k = "diff_" + lev
        base = lambda s, k=k: statistics.mean(0.0 if r["complete"] >= 500 else -float(r[k]) for r in s)
        f.write("\n%s: complete>=500 alone vs always-PASA %+.2f\n" % (lev, base(R)))
        f.write("F\tn_switched_extra\tgain vs always-PASA [CI]\tgain vs complete-gate [CI]\n")
        for F in (0, 150, 200, 250, 275, 300, 350, 400):
            pol = lambda s, F=F, k=k: statistics.mean(
                0.0 if (r["complete"] >= 500 and r["final"] >= F) else -float(r[k]) for r in s)
            extra = lambda s, F=F: pol(s, F) - base(s)
            n_extra = sum(1 for r in R if r["complete"] >= 500 and r["final"] < F)
            a = ci(pol, R); b = ci(extra, R)
            f.write("%d\t%d\t%+.2f [%+.2f, %+.2f]\t%+.2f [%+.2f, %+.2f]\n" % (
                F, n_extra, pol(R), a[0], a[1], extra(R), b[0], b[1]))
    f.write("\nExtra genomes switched at F=300 (complete>=500, final<300): name complete final PASA-BUSCO locus\n")
    for r in R:
        if r["complete"] >= 500 and r["final"] < 300:
            f.write("%s %d %d %s\n" % (r["name"], r["complete"], r["final"], r["diff_locus"]))

# Same policy sweep with the whole-genome complete count as the gate variable (production
# counts on the whole genome; experiment B trained on the train chromosomes). Added 2026-09-28
# after an independent review; see D117.
with open(os.path.join(OUT, "whole_genome_gate_sweep.txt"), "w") as f:
    rat = sorted(int(r["complete_genome"]) / r["complete"] for r in R)
    f.write("whole-genome / train-chromosome complete-model ratio: min %.2f, median %.2f, max %.2f (n=%d)\n"
            % (rat[0], statistics.median(rat), rat[-1], len(R)))
    for var in ("complete", "complete_genome"):
        f.write("\ngate variable: %s\nT\tgenomes_below\tgain vs always-PASA locus [CI]\n" % var)
        for t in (300, 500, 700, 1000, 1500, 2000):
            g = lambda s, t=t, var=var: statistics.mean(0.0 if int(r[var]) >= t else -r["d"] for r in s)
            lo, hi = ci(g, R)
            f.write("%d\t%d\t%+.2f [%+.2f, %+.2f]\n" % (t, sum(int(r[var]) < t for r in R), g(R), lo, hi))
    f.write("\nGenomes with < 500 complete models on the train chromosomes: name train whole_genome PASA-BUSCO locus\n")
    for r in sorted(R, key=lambda r: r["complete"]):
        if r["complete"] < 500:
            f.write("%s\t%d\t%s\t%+.2f\n" % (r["name"], r["complete"], r["complete_genome"], r["d"]))
