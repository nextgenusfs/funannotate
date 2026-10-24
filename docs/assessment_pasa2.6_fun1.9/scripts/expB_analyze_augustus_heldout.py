"""Experiment B: does held-out Augustus accuracy calibrate a training-source gate?

Inputs (data/):
  expB_augustus_heldout.tsv              held-out accuracy printed by `augustus` on the models set aside by
                                         randomSplit.pl, per genome and arm (scripts/expB_collect_augustus_heldout.py)
  expB_complete_threshold_by_genome.tsv  complete-ORF PASA models on the train chromosomes and the holdout
                                         locus F1 of the PASA- and BUSCO-trained predictors (40 genomes)
Questions:
  1. Within an arm, does held-out gene sensitivity track that arm's own holdout locus F1?
  2. Gate "train from PASA if PASA held-out gene Sn >= T, else BUSCO": gain over always-PASA / always-BUSCO
     (mean over genomes, bootstrap 95% CI over genomes, 2,000 resamples).
  3. Is the held-out number comparable across arms (PASA arm vs BUSCO arm)?
  4. Compared with the complete-model count gate, and the two combined.
Caveat: the held-out set of an arm comes from that arm's own training models, so values are not comparable
across arms by construction; question 3 measures how far that matters in practice.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
B, SEED = 2000, 1

aug = pd.read_csv(os.path.join(DATA, "expB_augustus_heldout.tsv"), sep="\t")
gen = pd.read_csv(os.path.join(DATA, "expB_complete_threshold_by_genome.tsv"), sep="\t").rename(columns={"name": "genome"})
pasa = aug[aug.arm == "pasa.A"].set_index("genome")
busco = aug[aug.arm.str.startswith("busco_r")].groupby("genome").mean(numeric_only=True)
d = gen.set_index("genome").join(pasa[["gene_sn", "gene_sp", "exon_sn", "nuc_sn"]].add_prefix("pasa_"), how="inner") \
       .join(busco[["gene_sn", "gene_sp", "exon_sn", "nuc_sn"]].add_prefix("busco_"), how="inner")
n = len(d)
rng = np.random.default_rng(SEED)
idx = rng.integers(0, n, size=(B, n))


def boot(v):
    v = np.asarray(v, float)
    m = v[idx].mean(axis=1)
    return v.mean(), np.percentile(m, 2.5), np.percentile(m, 97.5)


def fmt(t):
    return "{:+.2f} [{:+.2f}, {:+.2f}]".format(*t)


print("genomes:", n)
print("mean holdout locus F1: always-PASA {:.2f}, always-BUSCO {:.2f}".format(d.pasa_locus.mean(), d.busco_locus.mean()))

print("\n1. held-out gene Sn vs the same arm's holdout locus F1 (Spearman)")
for arm in ("pasa", "busco"):
    r, p = spearmanr(d[arm + "_gene_sn"], d[arm + "_locus"])
    print("  {:5s} rho={:+.2f} p={:.3g}".format(arm, r, p))
print("  gene Sp vs locus F1: pasa rho={:+.2f}, busco rho={:+.2f}".format(
    spearmanr(d.pasa_gene_sp, d.pasa_locus)[0], spearmanr(d.busco_gene_sp, d.busco_locus)[0]))

print("\n3. across arms")
print("  held-out gene Sn: PASA mean {:.3f}, BUSCO mean {:.3f}".format(d.pasa_gene_sn.mean(), d.busco_gene_sn.mean()))
r, p = spearmanr(d.pasa_gene_sn - d.busco_gene_sn, d.diff_locus)
print("  (PASA - BUSCO held-out gene Sn) vs (PASA - BUSCO holdout locus F1): rho={:+.2f} p={:.3g}".format(r, p))
r, p = spearmanr(d.pasa_gene_sn, d.diff_locus)
print("  PASA held-out gene Sn vs PASA - BUSCO locus F1: rho={:+.2f} p={:.3g}".format(r, p))
r, p = spearmanr(d.complete, d.diff_locus)
print("  complete-ORF count vs PASA - BUSCO locus F1: rho={:+.2f} p={:.3g}".format(r, p))
r, p = spearmanr(d.complete, d.pasa_gene_sn)
print("  complete-ORF count vs PASA held-out gene Sn: rho={:+.2f} p={:.3g}".format(r, p))


def sweep(title, pick_pasa, thresholds, label):
    print("\n" + title)
    print("  {:>8} {:>8} {:>8} {:>22} {:>22}".format(label, "n PASA", "n BUSCO", "gain vs always-PASA", "gain vs always-BUSCO"))
    for t in thresholds:
        use = pick_pasa(t)
        f1 = np.where(use, d.pasa_locus, d.busco_locus)
        print("  {:>8} {:>8d} {:>8d} {:>22} {:>22}".format(
            t, int(use.sum()), int((~use).sum()),
            fmt(boot(f1 - d.pasa_locus.values)), fmt(boot(f1 - d.busco_locus.values))))


sweep("2. gate on PASA held-out gene sensitivity (PASA if >= T)",
      lambda t: (d.pasa_gene_sn >= t).values, [0.0, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75], "T")
sweep("4a. complete-ORF count gate (PASA if >= T; train-chromosome counts, as in section 5 of the README)",
      lambda t: (d.complete >= t).values, [300, 400, 500, 700, 1000], "T")
sweep("4b. combined: complete >= 500 and PASA held-out gene Sn >= T",
      lambda t: ((d.complete >= 500) & (d.pasa_gene_sn >= t)).values, [0.4, 0.5, 0.6, 0.7], "T")

print("\ngenomes where PASA held-out gene Sn < 0.5:")
low = d[d.pasa_gene_sn < 0.5][["complete", "pasa_gene_sn", "busco_gene_sn", "pasa_locus", "busco_locus", "diff_locus"]]
print(low.round(2).to_string() if len(low) else "  none")
