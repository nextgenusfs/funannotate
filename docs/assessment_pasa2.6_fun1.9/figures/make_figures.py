#!/usr/bin/env python3
"""Figures for docs/assessment_pasa2.6_fun1.9 (PNG for the README, PDF for the paper).

Every figure reads only files in ../data/. Run from anywhere:
    /usr/bin/python3.12 docs/assessment_pasa2.6_fun1.9/figures/make_figures.py

Palette: the dataviz reference categorical slots, validated with
validate_palette.js (light mode, surface #fcfcfb): all checks pass; aqua and
yellow are below 3:1 contrast, so every series also has a direct label and its
own marker. A genome keeps the same color and marker in every figure.
"""
import csv
import gzip
import math
import os
import statistics
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

SURFACE = "#fcfcfb"
INK = "#1f1f1e"
INK2 = "#5c5b56"
GRID = "#e4e3dc"
ZERO = "#8a8981"
GENOMES = {  # name -> (label, color, marker)
    "Aspergillus_nidulans_FGSC_A4": ("A. nidulans", "#2a78d6", "o"),
    "Botrytis_cinerea_B05.10": ("B. cinerea", "#eb6834", "s"),
    "Neurospora_crassa_OR74A": ("N. crassa", "#1baf7a", "^"),
    "Cryptococcus_neoformans_H99": ("C. neoformans", "#eda100", "D"),
    "Schizophyllum_commune_H4-8": ("S. commune", "#e87ba4", "v"),
}
BLUE_DARK, BLUE_LIGHT = "#1c5cab", "#86b6ef"  # one-hue ordinal pair (Sn, Pr)

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 9.5,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "lines.linewidth": 2,
})


def read_tsv(name):
    path = os.path.join(DATA, name)
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f1(sn, pr):
    sn, pr = float(sn), float(pr)
    return 2 * sn * pr / (sn + pr) if sn + pr else 0.0


def save(fig, stem, source):
    fig.text(0.01, -0.035, "Source: data/" + source, fontsize=7.5, color=INK2, ha="left", va="top")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, stem + "." + ext), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem)


def fig1_titration():
    rows = read_tsv("titration_analysis_locus.tsv")
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.axhline(0, color=ZERO, lw=1)
    ax.axvline(500, color=INK2, lw=1, ls=(0, (4, 3)))
    by = defaultdict(list)
    for r in rows:
        if r["N"].isdigit():
            by[r["genome"]].append((int(r["N"]), float(r["delta_mean"]),
                                    float(r["delta_lo95"]), float(r["delta_hi95"])))
    offsets = {g: f for g, f in zip(GENOMES, (0.94, 1.0, 1.06, 1.12, 1.18))}
    for g, pts in by.items():
        lab, col, mk = GENOMES[g]
        pts.sort()
        xs = [n * offsets[g] for n, *_ in pts]
        ys = [d for _, d, _, _ in pts]
        lo = [d - l for _, d, l, _ in pts]
        hi = [h - d for _, d, _, h in pts]
        ax.plot(xs, ys, color=col, lw=2, zorder=2)
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt="none", ecolor=col, elinewidth=1.2, capsize=0, zorder=2)
        ax.scatter(xs, ys, s=38, marker=mk, color=col, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        dy = -11 if g == "Neurospora_crassa_OR74A" else 0
        ax.annotate(lab, (xs[-1], ys[-1]), xytext=(8, dy), textcoords="offset points",
                    color=INK, fontsize=8.5, va="center")
    ax.text(520, -3.55, "current gate: 500", color=INK2, fontsize=8, va="center")
    ax.set_xscale("log")
    ax.set_xticks([50, 100, 200, 300, 500, 750, 1000, 2000])
    ax.set_xticklabels(["50", "100", "200", "300", "500", "750", "1000", "2000"])
    ax.set_xlim(40, 3300)
    ax.set_xlabel("Complete PASA training models (N, log scale)")
    ax.set_ylabel("Holdout locus F1: PASA − BUSCO (points)")
    ax.set_title("BUSCO training wins below ~300 PASA models; at 500+ PASA is within about 1 point or ahead")
    ax.text(0.01, 0.98, "Below 0: BUSCO-trained predictors are more accurate.\nWhiskers: 95% bootstrap over subsample draws.\n"
            "BUSCO comparator = mean of 3 repeat runs.\nN. crassa pool is 1,999, so no N = 2000.", transform=ax.transAxes, fontsize=7.5, color=INK2, va="top")
    save(fig, "fig1_titration_pasa_vs_busco", "titration_analysis_locus.tsv")


def fig2_training_source():
    rows = read_tsv("scorecard.tsv")
    val = {(r["genome"], r["arm"], r["evidence"]): f1(r["locus_sn"], r["locus_pr"]) for r in rows}
    order = ["Neurospora_crassa_OR74A", "Aspergillus_nidulans_FGSC_A4", "Botrytis_cinerea_B05.10",
             "Cryptococcus_neoformans_H99", "Schizophyllum_commune_H4-8"]
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    for i, g in enumerate(order):
        lab, col, mk = GENOMES[g]
        p, b = val.get((g, "tx", "fixed")), val.get((g, "busco", "fixed"))
        if p is None or b is None:
            continue
        y = len(order) - 1 - i
        ax.plot([p, b], [y, y], color=GRID, lw=3, zorder=1, solid_capstyle="round")
        ax.scatter([p], [y], s=60, marker=mk, color=col, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.scatter([b], [y], s=60, marker=mk, facecolor=SURFACE, edgecolor=col, linewidth=2, zorder=3)
        ax.annotate("PASA %.1f · BUSCO %.1f" % (p, b), (max(p, b), y), xytext=(10, 0),
                    textcoords="offset points", va="center", fontsize=8, color=INK)
        se = val.get((g, "seR1R2", "fixed"))
        if se is not None:
            ax.scatter([se], [y - 0.28], s=40, marker=mk, color=col, alpha=0.55, edgecolor=SURFACE, zorder=3)
            ax.annotate("PASA + R1/R2 + single-exon %.1f" % se, (se, y - 0.28), xytext=(8, 0),
                        textcoords="offset points", va="center", fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([GENOMES[g][0] for g in reversed(order)])
    ax.set_ylim(-0.7, len(order) - 0.4)
    ax.set_xlim(25, 100)
    ax.set_xlabel("Holdout locus F1 (%), same fixed EVM evidence")
    ax.set_title("Training source: filled = PASA-trained, open = BUSCO-trained")
    ax.grid(axis="y", visible=False)
    save(fig, "fig2_training_source_by_genome", "scorecard.tsv (arms tx, busco, seR1R2; evidence fixed)")


def fig3_evidence_vs_training():
    rows = read_tsv("scorecard.tsv")
    v = {(r["genome"], r["arm"], r["evidence"]): (float(r["locus_sn"]), float(r["locus_pr"])) for r in rows}
    cases = [("Neurospora_crassa_OR74A", "txR1R2", "N. crassa (divergent): R1 + R2"),
             ("Botrytis_cinerea_B05.10", "txR1", "B. cinerea (same strain): R1")]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
    for ax, (g, arm, title) in zip(axes, cases):
        base_f, base_o = v[(g, "tx", "fixed")], v[(g, "tx", "own")]
        new_f, new_o = v[(g, arm, "fixed")], v[(g, arm, "own")]
        eff = {"Training\nalone": [new_f[k] - base_f[k] for k in (0, 1)],
               "Evidence\nalone": [new_o[k] - new_f[k] for k in (0, 1)],
               "Both": [new_o[k] - base_o[k] for k in (0, 1)]}
        xs = range(len(eff))
        w = 0.36
        for k, (name, col) in enumerate((("Locus Sn", BLUE_DARK), ("Locus Pr", BLUE_LIGHT))):
            vals = [e[k] for e in eff.values()]
            bars = ax.bar([x + (k - 0.5) * w for x in xs], vals, width=w - 0.04, color=col, label=name, zorder=2)
            for b_, val in zip(bars, vals):
                ax.annotate("%+.1f" % val, (b_.get_x() + b_.get_width() / 2, val),
                            xytext=(0, 3 if val >= 0 else -10), textcoords="offset points",
                            ha="center", fontsize=7.5, color=INK2)
        ax.axhline(0, color=ZERO, lw=1)
        ax.set_xticks(list(xs))
        ax.set_xticklabels(list(eff.keys()))
        ax.set_title(title, fontsize=9.5)
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("Change in holdout locus score (points)")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("RNA-seq fixes help mostly through the evidence EVM uses, not the training set",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    save(fig, "fig3_evidence_vs_training_effect", "scorecard.tsv (tx vs txR1R2 / txR1, fixed and own evidence)")


def fig4_r1_validation():
    rows = read_tsv("r1_validation.tsv")
    d = defaultdict(dict)
    for r in rows:
        d[r["genome"]][r["stage"]] = 100.0 * int(r["spliced_valid"]) / int(r["spliced_total"])
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for g, st in d.items():
        lab, col, mk = GENOMES[g]
        ys = [st["before_R1"], st["after_R1"]]
        ax.plot([0, 1], ys, color=col, lw=2, zorder=2)
        ax.scatter([0, 1], ys, s=50, marker=mk, color=col, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.annotate("%.1f%%" % ys[0], (0, ys[0]), xytext=(-8, 0), textcoords="offset points", ha="right", va="center", fontsize=8, color=INK2)
        ax.annotate("%.0f%%  %s" % (ys[1], lab), (1, ys[1]), xytext=(8, 0), textcoords="offset points", va="center", fontsize=8.5, color=INK)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Before fix (2018 parser)", "After fix (CIGAR-based)"])
    ax.set_xlim(-0.45, 1.6)
    ax.set_ylim(-3, 103)
    ax.set_ylabel("Spliced minimap2 alignments\npassing PASA validation (%)")
    ax.set_title("The minimap2 conversion fix restores spliced evidence")
    ax.grid(axis="x", visible=False)
    save(fig, "fig4_minimap2_fix_validation", "r1_validation.tsv")


def fig5_intron_accuracy():
    rows = read_tsv("intron_accuracy.tsv")
    bins = [">=99", "97-99", "95-97", "<95"]
    aligners = [("minimap2_fixed", "minimap2 (fixed)", "#2a78d6"), ("gmap", "gmap", "#eb6834"), ("blat", "blat", "#1baf7a")]
    v = {(r["aligner"], r["identity_bin"]): float(r["pct_exact_refseq_intron"]) for r in rows}
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    w = 0.26
    for k, (key, lab, col) in enumerate(aligners):
        xs = [i + (k - 1) * w for i in range(len(bins))]
        ys = [v[(key, b)] for b in bins]
        ax.bar(xs, ys, width=w - 0.04, color=col, label="%s (all: %.1f%%)" % (lab, v[(key, "all")]), zorder=2)
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels(["≥ 99%", "97-99%", "95-97%", "< 95%"])
    ax.set_xlabel("Alignment identity to the genome (N. crassa, divergent RNA-seq)")
    ax.set_ylabel("Introns exactly matching RefSeq (%)")
    ax.set_ylim(0, 100)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=8)
    ax.set_title("minimap2 keeps splice sites exact as read identity falls", pad=26)
    ax.grid(axis="x", visible=False)
    save(fig, "fig5_intron_accuracy_by_aligner", "intron_accuracy.tsv")


def fig6_f1_production():
    rows = read_tsv("production_f1_scan.tsv.gz")
    by = defaultdict(list)
    for r in rows:
        if r["frac_not_mod3"] not in ("NA", "") and r["pasa_gff3_date"] >= "2026-06-01":
            by[r["pasa_gff3_date"]].append(float(r["frac_not_mod3"]))
    dates = sorted(by)
    import datetime as dt
    xs = [dt.date.fromisoformat(d) for d in dates]
    med = [statistics.median(by[d]) for d in dates]
    q1 = [statistics.quantiles(by[d], n=4)[0] if len(by[d]) >= 4 else min(by[d]) for d in dates]
    q3 = [statistics.quantiles(by[d], n=4)[2] if len(by[d]) >= 4 else max(by[d]) for d in dates]
    n = [len(by[d]) for d in dates]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.axvspan(dt.date(2026, 6, 29), dt.date(2026, 9, 24), color="#eceae3", zorder=0)
    ax.text(dt.date(2026, 7, 1), 63, "F1 bug in deployed PASA (2026-06-29 to 09-24)", fontsize=7.5, color=INK2, va="top")
    ax.vlines(xs, [a * 100 for a in q1], [b * 100 for b in q3], color="#86b6ef", lw=1.2, zorder=2)
    ax.scatter(xs, [m * 100 for m in med], s=[8 + 22 * math.log10(1 + k) for k in n], color="#2a78d6",
               edgecolor=SURFACE, linewidth=1, zorder=3)
    ax.set_ylabel("PASA models with CDS length\nnot divisible by 3 (%)")
    ax.set_ylim(-2, 65)
    ax.set_title("Frame-broken PASA training models in 8,007 production genomes, by date")
    fig.autofmt_xdate()
    save(fig, "fig6_f1_production_timeline", "production_f1_scan.tsv.gz. Dot = median per file date (size ~ log number of genomes); whisker = interquartile range.")


def fig7_identity():
    rows = read_tsv("production_identity.tsv.gz")
    vals = [float(r["median_id"]) for r in rows if r["median_id"] not in ("NA", "")]
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    edges = [85 + 0.5 * i for i in range(31)]
    clipped = [min(max(v, 85.0), 99.999) for v in vals]
    ax.hist(clipped, bins=edges, color="#2a78d6", edgecolor=SURFACE, linewidth=0.8, zorder=2)
    ymax = ax.get_ylim()[1]
    for x, lab, yf, ha in ((90, " < 90%: train from BUSCO", 0.45, "left"),
                           (99, "same strain: ≥ 99% ", 0.93, "right")):
        ax.axvline(x, color=INK2, lw=1, ls=(0, (4, 3)))
        ax.text(x, ymax * yf, lab, fontsize=7.5, color=INK2, va="top", ha=ha)
    ax.set_xticks([85, 87.5, 90, 92.5, 95, 97.5, 100])
    ax.set_xticklabels(["≤ 85", "87.5", "90", "92.5", "95", "97.5", "100"])
    n = len(vals)
    ge99 = sum(v >= 99 for v in vals)
    lt90 = sum(v < 90 for v in vals)
    ax.text(0.01, 0.62, "n = %d genomes\n≥ 99%%: %d (%.0f%%)\n90-99%%: %d (%.0f%%)\n< 90%%: %d (%.1f%%)"
            % (n, ge99, 100 * ge99 / n, n - ge99 - lt90, 100 * (n - ge99 - lt90) / n, lt90, 100 * lt90 / n),
            transform=ax.transAxes, fontsize=8, color=INK2)
    ax.set_xlabel("Median transcript-to-genome identity per genome (%); reads 0.4-0.9 points higher")
    ax.set_ylabel("Genomes")
    ax.set_title("Production RNA-seq identity against the gate categories")
    ax.grid(axis="x", visible=False)
    save(fig, "fig7_production_rnaseq_identity", "production_identity.tsv.gz")


def fig8_single_exon():
    rows = read_tsv("single_exon_scores.tsv")
    v = {(r["genome"], r["arm"]): r for r in rows}
    metrics = [("single_sn", "Single-exon Sn"), ("single_pr", "Single-exon Pr"),
               ("multi_sn", "Multi-exon Sn"), ("multi_pr", "Multi-exon Pr")]
    genomes = [g for g in GENOMES if (g, "se") in v and (g, "tx2") in v]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    w = 0.8 / len(genomes)
    for k, g in enumerate(genomes):
        lab, col, mk = GENOMES[g]
        xs = [i + (k - (len(genomes) - 1) / 2) * w for i in range(len(metrics))]
        ys = [float(v[(g, "se")][m]) - float(v[(g, "tx2")][m]) for m, _ in metrics]
        ax.bar(xs, ys, width=w - 0.03, color=col, label=lab, zorder=2)
    ax.axhline(0, color=ZERO, lw=1)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([m[1] for m in metrics])
    ax.set_ylabel("Change with single-exon training (points)")
    ax.set_title("Single-exon training genes: more single-exon genes found, multi-exon precision up")
    ax.legend(ncol=3, fontsize=8, loc="upper right")
    ax.grid(axis="x", visible=False)
    save(fig, "fig8_single_exon_training_effect", "single_exon_scores.tsv (se − tx2, exact CDS-chain match)")


def fig9_identity_gate():
    rows = read_tsv("expB_identity_vs_training.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6), sharey=True)
    hi_col, lo_col = "#2a78d6", "#eb6834"
    for ax, key, xlab in ((axes[0], "identity", "Median read identity to the genome (%)"),
                          (axes[1], "complete", "Complete PASA models, train chromosomes (log)")):
        ax.axhline(0, color=ZERO, lw=1)
        for few, col, mk, lab in ((False, hi_col, "o", "≥ 500 complete"),
                                  (True, lo_col, "s", "< 500: PASA gate → BUSCO")):
            pts = [r for r in rows if (int(r["complete"]) < 500) == few]
            ax.scatter([float(r[key]) for r in pts], [float(r["diff"]) for r in pts], s=26,
                       color=col, marker=mk, edgecolor=SURFACE, linewidth=0.6, zorder=3, label=lab)
        ax.set_xlabel(xlab)
    axes[0].axvline(95, color=INK2, lw=1, ls=(0, (4, 3)))
    axes[0].text(95.1, 5.4, "95%", fontsize=7.5, color=INK2, va="top")
    for r in rows:
        if r["name"].startswith(("Penicillium_antarcticum", "Metschnikowia")):
            axes[0].annotate(r["name"].split("_")[0][0] + ". " + r["name"].split("_")[1],
                             (float(r["identity"]), float(r["diff"])), xytext=(4, 4),
                             textcoords="offset points", fontsize=7, color=INK2)
    axes[1].set_xscale("log")
    ticks = [150, 300, 500, 1000, 2000, 4000]
    axes[1].set_xticks(ticks)
    axes[1].set_xticklabels([format(t, ",") for t in ticks])
    axes[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axes[1].text(520, 5.4, "500", fontsize=7.5, color=INK2, va="top")
    axes[1].axvline(500, color=INK2, lw=1, ls=(0, (4, 3)))
    axes[0].set_ylabel("PASA − BUSCO training, holdout locus F1")
    axes[1].legend(fontsize=7.5, loc="lower right")
    axes[0].set_title("Experiment B: identity does not predict the training outcome", loc="left")
    save(fig, "fig9_identity_vs_training_outcome",
         "expB_identity_vs_training.tsv (40 RefSeq genomes; BUSCO = mean of 3 repeats)")


if __name__ == "__main__":
    fig1_titration()
    fig2_training_source()
    fig3_evidence_vs_training()
    fig4_r1_validation()
    fig5_intron_accuracy()
    fig6_f1_production()
    fig7_identity()
    fig8_single_exon()
    fig9_identity_gate()
