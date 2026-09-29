"""Experiment C analysis: PASA vs BUSCO training on the whole genome, masked scoring (D118).

Per genome: mean and SD over 3 repeats per arm for locus, exon, intron-chain F1; PASA - BUSCO with a
rough 95% interval from the repeat SDs (t, df=4, two groups of 3); complete-model count and final
training models from pasa_r1's decision log (whole genome = the production gate variable);
experiment B (training chromosomes) PASA - BUSCO for comparison. Writes expC_summary.tsv.
"""
import csv, math, os, statistics
C = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(C, "..", "experiment_B", "read_identity", "complete_threshold_by_genome.tsv")
LEV = [("locus", "locus_sn", "locus_pr"), ("exon", "exon_sn", "exon_pr"), ("intron_chain", "intron_chain_sn", "intron_chain_pr")]
f1 = lambda sn, pr: 2 * sn * pr / (sn + pr) if sn + pr else 0.0
S = list(csv.DictReader(open(os.path.join(C, "scores.tsv")), delimiter="\t"))
expb = {r["name"]: r for r in csv.DictReader(open(B), delimiter="\t")}


def dec(g, stage):
    p = os.path.join(C, g, "pasa_r1", "out", "logfiles", "training_decisions.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["stage"] == stage:
            return r["value"]
    return ""


out = []
for g in sorted({r["genome"] for r in S}):
    row = {"genome": g, "complete_genome": dec(g, "pasa_gate"), "final_pasa_models": dec(g, "select_final"),
           "complete_train_expB": expb[g]["complete"], "ref_mrnas_scored": "", "ref_mrnas_total": ""}
    for lev, a, b in LEV:
        v = {arm: [f1(float(r[a]), float(r[b])) for r in S if r["genome"] == g and r["variant"].startswith(arm)] for arm in ("pasa", "busco")}
        mp, mb = statistics.mean(v["pasa"]), statistics.mean(v["busco"])
        sp, sb = statistics.stdev(v["pasa"]), statistics.stdev(v["busco"])
        half = 2.776 * math.sqrt(sp ** 2 / 3 + sb ** 2 / 3)
        row.update({lev + "_pasa": mp, lev + "_busco": mb, lev + "_pasa_sd": sp, lev + "_busco_sd": sb,
                    lev + "_diff": mp - mb, lev + "_ci_lo": mp - mb - half, lev + "_ci_hi": mp - mb + half,
                    lev + "_diff_expB": float(expb[g]["diff_" + lev])})
    r0 = [r for r in S if r["genome"] == g][0]
    row["ref_mrnas_scored"], row["ref_mrnas_total"] = r0["ref_mrnas_scored"], r0["ref_mrnas_total"]
    out.append(row)
out.sort(key=lambda r: int(r["complete_genome"] or 0))
cols = list(out[0].keys())
with open(os.path.join(C, "expC_summary.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in out:
        f.write("\t".join(("%.2f" % r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
print("genome\tcomplete(whole)\tcomplete(expB train)\tfinal\tscored/total mRNA\tlocus PASA(sd)\tlocus BUSCO(sd)\tPASA-BUSCO [95%]\texpB PASA-BUSCO")
for r in out:
    print("%s\t%s\t%s\t%s\t%s/%s\t%.2f (%.2f)\t%.2f (%.2f)\t%+.2f [%+.2f, %+.2f]\t%+.2f" % (
        r["genome"], r["complete_genome"], r["complete_train_expB"], r["final_pasa_models"], r["ref_mrnas_scored"], r["ref_mrnas_total"],
        r["locus_pasa"], r["locus_pasa_sd"], r["locus_busco"], r["locus_busco_sd"], r["locus_diff"], r["locus_ci_lo"], r["locus_ci_hi"], r["locus_diff_expB"]))
print("\nexon and intron chain PASA-BUSCO (expC | expB):")
for r in out:
    print("%s\texon %+.2f [%+.2f, %+.2f] | %+.2f\tintron_chain %+.2f [%+.2f, %+.2f] | %+.2f" % (
        r["genome"], r["exon_diff"], r["exon_ci_lo"], r["exon_ci_hi"], r["exon_diff_expB"],
        r["intron_chain_diff"], r["intron_chain_ci_lo"], r["intron_chain_ci_hi"], r["intron_chain_diff_expB"]))
