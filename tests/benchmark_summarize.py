#!/usr/bin/env python3
"""Summarize the clean Rust-vs-baseline benchmark.

Reads <bench_root>/<stage>/results/<engine>_rep<N>.{metrics,env} for the stages
trinity, pasa and evm, and prints one table per stage: wall time, peak memory,
output checks, and engine ratios with the spread over replicates.

It refuses to report a comparison when the legs are not comparable:
  - a run with status other than OK (FAILED, INVALID, missing) is listed and
    left out of the statistics;
  - legs that ran on different CPU models are reported as not comparable.

Usage:
    tests/benchmark_summarize.py /path/to/BENCH_ROOT [--markdown]
"""
import argparse
import collections
import glob
import os
import re
import statistics
import sys

STAGES = [
    ("trinity", "Trinity (train --stop_after_trinity)", "wall"),
    ("pasa", "PASA stage of train", "pasa_stage_seconds"),
    ("evm", "EVM (funannotate predict)", "wall"),
]


def parse_kv(path):
    data = {}
    if not os.path.isfile(path):
        return data
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = re.match(r"^\s*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)$", line)
            if m:
                data["wall"] = parse_clock(m.group(1))
                continue
            m = re.match(r"^\s*Maximum resident set size \(kbytes\): (\d+)$", line)
            if m:
                data["max_rss_gb"] = int(m.group(1)) / 1048576.0
                continue
            m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*)$", line)
            if m:
                data[m.group(1)] = m.group(2)
    return data


def parse_clock(text):
    parts = [float(x) for x in text.strip().split(":")]
    seconds = 0.0
    for p in parts:
        seconds = seconds * 60 + p
    return seconds


def load_stage(root, stage):
    runs = []
    for mpath in sorted(glob.glob(os.path.join(root, stage, "results", "*_rep*.metrics"))):
        name = os.path.basename(mpath)[: -len(".metrics")]
        m = re.match(r"^(.+)_rep(\d+)$", name)
        if not m:
            continue
        run = {"engine": m.group(1), "rep": int(m.group(2)), "name": name}
        run.update(parse_kv(mpath))
        env = parse_kv(mpath[: -len(".metrics")] + ".env")
        for k, v in env.items():
            run.setdefault(k, v)
        run.setdefault("status", "UNKNOWN")
        runs.append(run)
    return runs


def num(run, key):
    try:
        return float(run[key])
    except (KeyError, ValueError):
        return None


def fmt_time(sec):
    if sec is None:
        return "NA"
    h, rem = divmod(int(round(sec)), 3600)
    m, s = divmod(rem, 60)
    return "%d:%02d:%02d" % (h, m, s)


def spread(values):
    if not values:
        return "NA"
    med = statistics.median(values)
    return "%s (min %s, max %s, n=%d)" % (fmt_time(med), fmt_time(min(values)), fmt_time(max(values)), len(values))


GFF_GLOBS = {
    "pasa": "work_pasa_{name}/output/training/funannotate_train.pasa.gff3",
    "evm": "predict_work_{name}/predict_results/*.gff3",
}


def structure_lines(path):
    """Feature lines without IDs: sequence, type, start, end, strand."""
    counts = collections.Counter()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7:
                counts[(f[0], f[2], f[3], f[4], f[6])] += 1
    return counts


def divergence_report(root, stage, good):
    """Compare how far apart runs are: within one engine and between engines.

    PASA and EVM number IDs in a run-dependent order, and train is not deterministic, so two
    runs of the same engine already differ. A difference between engines only matters if it is
    larger than that.
    """
    pattern = GFF_GLOBS.get(stage)
    if not pattern:
        return []
    data = {}
    for r in good:
        hits = sorted(glob.glob(os.path.join(root, stage, pattern.format(name=r["name"]))))
        if hits:
            data[r["name"]] = (r["engine"], structure_lines(hits[0]))
    names = sorted(data)
    if len(names) < 2:
        return []
    within, between = [], []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ca, cb = data[a][1], data[b][1]
            diff = sum(((ca - cb) + (cb - ca)).values())
            (within if data[a][0] == data[b][0] else between).append(diff)
    total = statistics.median([sum(c.values()) for _, c in data.values()])
    out = ["", "Output differences (feature lines that differ between two runs, IDs ignored; a run has about %d lines):" % total]
    if within:
        out.append("- same engine, different replicate: median %d (range %d-%d, %d pairs)" % (
            statistics.median(within), min(within), max(within), len(within)))
    if between:
        out.append("- different engines: median %d (range %d-%d, %d pairs)" % (
            statistics.median(between), min(between), max(between), len(between)))
    if within and between:
        w, b = statistics.median(within), statistics.median(between)
        if w > 0:
            out.append("  Between-engine differences are %.2fx the within-engine differences." % (b / w))
    return out


def stage_report(root, stage, title, metric):
    runs = load_stage(root, stage)
    out = ["", "## %s" % title, ""]
    if not runs:
        out.append("No results found in %s/%s/results." % (root, stage))
        return out

    out.append("| run | status | host | CPU | %s | wall | max RSS (GB) | output check |" % ("PASA stage" if metric != "wall" else "stage time"))
    out.append("|---|---|---|---|---|---|---|---|")
    for r in runs:
        check = r.get("transcripts") or r.get("pasa_gff3_genes") or r.get("genes") or ""
        if r.get("transcripts"):
            check = "%s transcripts" % r["transcripts"]
        elif r.get("pasa_gff3_genes"):
            check = "%s genes" % r["pasa_gff3_genes"]
        elif r.get("genes"):
            check = "%s genes" % r["genes"]
        out.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (
            r["name"], r["status"], r.get("host", "NA"), r.get("cpu_model", "NA"),
            fmt_time(num(r, metric)), fmt_time(num(r, "wall")),
            "%.1f" % r["max_rss_gb"] if "max_rss_gb" in r else "NA", check))

    good = [r for r in runs if r["status"] == "OK"]
    bad = [r for r in runs if r["status"] != "OK"]
    out.append("")
    if bad:
        out.append("Left out of the statistics (not OK): " + ", ".join("%s=%s" % (r["name"], r["status"]) for r in bad))
        out.append("")

    cpus = sorted({r.get("cpu_model", "NA") for r in good})
    comparable = len(cpus) == 1
    if not comparable:
        out.append("NOT COMPARABLE: the legs ran on different CPU models: " + "; ".join(cpus))
        return out
    out.append("All compared legs ran on: %s" % cpus[0])
    out.append("")

    by_engine = {}
    for r in good:
        v = num(r, metric)
        if v is not None:
            by_engine.setdefault(r["engine"], []).append(v)
    out.append("| engine | %s, median (range) |" % metric)
    out.append("|---|---|")
    for eng in sorted(by_engine):
        out.append("| %s | %s |" % (eng, spread(by_engine[eng])))

    base = by_engine.get("perl")
    if base:
        base_med = statistics.median(base)
        out.append("")
        for eng in sorted(by_engine):
            if eng == "perl":
                continue
            med = statistics.median(by_engine[eng])
            if med <= 0 or min(by_engine[eng]) <= 0 or min(base) <= 0:
                out.append("%s vs perl: a time is zero; no ratio reported" % eng)
                continue
            # Ratio of medians, plus the worst and best case from the extremes.
            worst = min(base) / max(by_engine[eng])
            best = max(base) / min(by_engine[eng])
            out.append("%s vs perl: median speedup %.2fx (extremes %.2fx to %.2fx)" % (eng, base_med / med, worst, best))
            if worst <= 1.0 <= best:
                out.append("  The ranges overlap. This does not show a difference at this replicate count.")

    # Per-tool comparisons: each isolates one Rust tool by comparing two engines that
    # differ only in that tool.
    pairs = [
        ("rust", "nocdb", "cdbyank_rust + faidx_rust", "nocdb has the C++ cdbyank"),
        ("rust", "noslclust", "slclust_rust", "noslclust has the C++ slclust"),
        ("rustasm", "rust", "the Rust assembler (pasa_rust)", "rust has the C++ assembler"),
    ]
    pair_lines = []
    for fast, slow, label, note in pairs:
        if fast in by_engine and slow in by_engine:
            f, sl = by_engine[fast], by_engine[slow]
            if min(f) <= 0 or min(sl) <= 0:
                continue
            ratio = statistics.median(sl) / statistics.median(f)
            overlap = not (max(f) < min(sl) or max(sl) < min(f))
            pair_lines.append("- %s: %s median %s vs %s median %s: %.3fx (%s)%s" % (
                label, fast, fmt_time(statistics.median(f)), slow, fmt_time(statistics.median(sl)),
                ratio, note, "; ranges overlap, so no difference is shown" if overlap else "; ranges do not overlap"))
    if pair_lines:
        out.append("")
        out.append("Single-tool comparisons (replicates: %s):" % ", ".join("%s n=%d" % (e, len(by_engine[e])) for e in sorted(by_engine)))
        out.extend(pair_lines)

    out.extend(divergence_report(root, stage, good))

    # Output identity checks. The ID-dependent GFF3 checksums are not used: see divergence_report.
    for key, label in (("transcripts", "transcript counts"),):
        vals = {}
        for r in good:
            if key in r:
                vals.setdefault(r["engine"], set()).add(r[key])
        if vals:
            allvals = set().union(*vals.values())
            out.append("")
            out.append("%s: %s" % (label, "identical in all legs" if len(allvals) == 1 else "DIFFER between runs: " +
                                    "; ".join("%s=%s" % (e, sorted(v)) for e, v in sorted(vals.items()))))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bench_root")
    args = ap.parse_args(argv)
    lines = ["# Clean Rust benchmark summary", "", "Root: %s" % args.bench_root]
    for stage, title, metric in STAGES:
        lines.extend(stage_report(args.bench_root, stage, title, metric))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
