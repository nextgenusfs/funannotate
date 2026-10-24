"""Experiment D (--optimize_augustus): choose genomes from experiment B.

All genomes whose PASA held-out gene sensitivity is below 0.5 (the low-accuracy group; 6 of 40), plus an even
spread of 10 of the remaining 34 ordered by that sensitivity, so a benefit can be related to starting accuracy.
Three genomes (the low genome with the largest PASA - BUSCO locus F1, the median and the top of the spread)
get a second optimized run of each arm to measure run-to-run spread.
Output: data/expD_genomes.tsv (genome, group, pasa_gene_sn, repeat)
"""
import os
import numpy as np
import pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
a = pd.read_csv(os.path.join(D, "expB_augustus_heldout.tsv"), sep="\t")
g = pd.read_csv(os.path.join(D, "expB_complete_threshold_by_genome.tsv"), sep="\t").rename(columns={"name": "genome"}).set_index("genome")
p = a[a.arm == "pasa.A"].set_index("genome").gene_sn
d = g.join(p.rename("pasa_gene_sn"), how="inner")
low = d[d.pasa_gene_sn < 0.5].sort_values("pasa_gene_sn")
rest = d[d.pasa_gene_sn >= 0.5].sort_values("pasa_gene_sn")
pick = rest.iloc[np.unique(np.linspace(0, len(rest) - 1, 10).round().astype(int))]
rows = [(n, "low", r.pasa_gene_sn) for n, r in low.iterrows()] + [(n, "spread", r.pasa_gene_sn) for n, r in pick.iterrows()]
rep = {low.diff_locus.idxmax(), pick.index[len(pick) // 2], pick.index[-1]}
out = pd.DataFrame(rows, columns=["genome", "group", "pasa_gene_sn"])
out["repeat"] = out.genome.isin(rep).astype(int)
out.to_csv(os.path.join(D, "expD_genomes.tsv"), sep="\t", index=False, float_format="%.3f")
print(out.to_string(index=False))
print(len(out), "genomes;", int(out["repeat"].sum()), "with a repeat")
