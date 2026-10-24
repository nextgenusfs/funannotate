"""Experiment D: does --optimize_augustus help, and when?

Reads data/expD_optimize_scores.tsv (scripts/expD_collect_optimize.py). Per arm (PASA-trained, BUSCO-trained):
  1. held-out gene sensitivity before and after optimization (self-consistency; same-arm comparison only);
  2. paired holdout locus F1, optimized minus plain, mean over genomes with a genome-level bootstrap 95% CI;
  3. whether the gain depends on the starting held-out gene sensitivity (Spearman) and on the low group;
  4. run-to-run spread of optimized runs (_r2 repeats) against the size of the gain;
  5. cost: step A wall time, plain vs optimized.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

B, SEED = 2000, 1
HERE = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(HERE, "..", "data", "expD_optimize_scores.tsv")
d = pd.read_csv(path, sep="\t")
grp = pd.read_csv(os.path.join(HERE, "..", "data", "expD_genomes.tsv"), sep="\t").set_index("genome").group
d["f1"] = 2 * d.locus_sn * d.locus_pr / (d.locus_sn + d.locus_pr)
rng = np.random.default_rng(SEED)


def boot(v):
    v = np.asarray(v, float)
    m = v[rng.integers(0, len(v), size=(B, len(v)))].mean(axis=1)
    return "{:+.2f} [{:+.2f}, {:+.2f}]".format(v.mean(), *np.percentile(m, [2.5, 97.5]))


for arm in ("pasa", "busco"):
    a = d[d.arm == arm]
    plain = a[a.variant == "plain"].set_index("genome")
    opt = a[a.variant == "opt"].set_index("genome")
    j = plain.join(opt, lsuffix="_p", rsuffix="_o", how="inner").dropna(subset=["f1_p", "f1_o"])
    print("\n== {} arm: {} genomes with both plain and optimized scores".format(arm, len(j)))
    if len(j) < 3:
        continue
    j["dF1"] = j.f1_o - j.f1_p
    print("held-out gene Sn  plain {:.3f} -> optimized {:.3f} (mean over genomes)".format(
        j.init_gene_sn_p.mean(), j.final_gene_sn_o.mean()))
    print("holdout locus F1 gain, optimized - plain: {}   (genomes better/worse by >0.5: {}/{})".format(
        boot(j.dF1), int((j.dF1 > 0.5).sum()), int((j.dF1 < -0.5).sum())))
    for lab, sel in (("low group", grp.reindex(j.index) == "low"), ("spread group", grp.reindex(j.index) == "spread")):
        if sel.sum() >= 2:
            print("  {}: n={} gain {}".format(lab, int(sel.sum()), boot(j.dF1[sel])))
    r, p = spearmanr(j.init_gene_sn_p, j.dF1)
    print("gain vs starting held-out gene Sn: rho={:+.2f} p={:.3g}".format(r, p))
    r, p = spearmanr(j.final_gene_sn_o - j.init_gene_sn_p, j.dF1)
    print("gain vs change in held-out gene Sn: rho={:+.2f} p={:.3g}".format(r, p))
    if j.busco_c_p.notna().all() and j.busco_c_o.notna().all():
        print("proteome BUSCO C change: {}".format(boot(j.busco_c_o - j.busco_c_p)))
    print("step A wall time: plain {:.0f} min, optimized {:.0f} min (median)".format(
        j.wall_A_s_p.median() / 60, j.wall_A_s_o.median() / 60))
    rep = a[a.variant == "opt_r2"].set_index("genome").join(opt, lsuffix="_2", rsuffix="_1", how="inner").dropna(subset=["f1_1", "f1_2"])
    if len(rep):
        print("optimized run-to-run |difference| in locus F1 (n={}): {}".format(
            len(rep), ", ".join("{:.2f}".format(x) for x in (rep.f1_2 - rep.f1_1).abs())))
