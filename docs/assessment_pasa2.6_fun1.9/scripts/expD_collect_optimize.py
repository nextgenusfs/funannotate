"""Experiment D: collect plain and --optimize_augustus runs (run on the cluster).

Usage: python3 expD_collect_optimize.py .../experiment_B > expD_optimize_scores.tsv
One row per genome x arm (pasa, busco) x variant (plain, opt, opt_r2): held-out Augustus accuracy before
optimization (augustus.initial.training.txt) and after (augustus.final.training.txt, opt runs only), step A wall
seconds, and the step B holdout scores and proteome BUSCO. Plain baselines are experiment B's pasa and busco_r1.
"""
import csv, os, re, sys

root = sys.argv[1]
nan = float("nan")
genomes = [r["genome"] for r in csv.DictReader(open(os.path.join(root, "expD_genomes.tsv")), delimiter="\t")]
VARIANTS = [("pasa", "plain", "pasa"), ("busco", "plain", "busco_r1"), ("pasa", "opt", "pasa_opt"),
            ("busco", "opt", "busco_opt"), ("pasa", "opt_r2", "pasa_opt_r2"), ("busco", "opt_r2", "busco_opt_r2")]
COLS = ["genome", "arm", "variant", "step_A_ok", "step_B_ok", "wall_A_s",
        "init_gene_sn", "init_gene_sp", "final_gene_sn", "final_gene_sp",
        "locus_sn", "locus_pr", "exon_sn", "exon_pr", "intron_chain_sn", "intron_chain_pr", "pred_genes", "busco_c"]


def num(x):
    try:
        return float(x)
    except ValueError:
        return nan


def gene_level(path):
    if not os.path.isfile(path):
        return nan, nan
    for line in open(path, errors="replace"):
        c = [x.strip() for x in line.rstrip("\n").split("|")]
        if line.startswith("gene level") and len(c) > 7:
            return num(c[6]), num(c[7])
    return nan, nan


def status(path):
    d = {}
    if os.path.isfile(path):
        for line in open(path):
            k, _, v = line.rstrip("\n").partition("\t")
            d[k] = v
    return d


print("\t".join(COLS))
for g in genomes:
    for arm, variant, dirname in VARIANTS:
        A, B = os.path.join(root, g, dirname + ".A"), os.path.join(root, g, dirname + ".B")
        if not os.path.isdir(A):
            continue
        sa, sb = status(os.path.join(A, "run_status.tsv")), status(os.path.join(B, "run_status.tsv"))
        mis = os.path.join(A, "out", "predict_misc")
        isn, isp = gene_level(os.path.join(mis, "augustus.initial.training.txt"))
        fsn, fsp = gene_level(os.path.join(mis, "augustus.final.training.txt"))
        score = {}
        sp = os.path.join(B, "score.tsv")
        if os.path.isfile(sp):
            score = next(csv.DictReader(open(sp), delimiter="\t"), {})
        busco = nan
        bp = os.path.join(B, "busco_short_summary.txt")
        if os.path.isfile(bp):
            m = re.search(r"C:([0-9.]+)%", open(bp, errors="replace").read())
            busco = float(m.group(1)) if m else nan
        row = [g, arm, variant, sa.get("exit_status", ""), sb.get("exit_status", ""), sa.get("wall_seconds", ""),
               isn, isp, fsn, fsp] + [score.get(k, "") for k in
               ("locus_sn", "locus_pr", "exon_sn", "exon_pr", "intron_chain_sn", "intron_chain_pr", "pred_genes")] + [busco]
        print("\t".join(str(x) for x in row))
