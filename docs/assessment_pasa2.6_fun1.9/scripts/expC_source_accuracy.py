#!/usr/bin/python3.12
"""Per-source accuracy of EVM inputs, and where missed RefSeq genes are lost (2026-09-30).

For one predict run and the experiment C mask:
  - per EVM input source (gene_predictions.gff3 column 2) and for evm.round1 and the final
    models: gene count, exact-chain sensitivity (scored RefSeq genes with an exact chain in the
    source) and exact-chain precision (scored source genes with an exact RefSeq chain).
  - for RefSeq genes that the final models miss (no same-strand CDS overlap): whether any EVM
    input overlaps them, whether evm.round1 overlaps them (lost after EVM = filter step), and
    whether transcript or protein evidence overlaps them.
Usage: source_accuracy.py --ref-gff3 REF --run-dir RUN --mask M... --genome G --variant V --out TSV
"""
import argparse
import collections
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from error_breakdown import Index, gene_overlap, genes_from  # noqa: E402
from masked_score import mask_spans, overlaps  # noqa: E402

sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/refseq_benchmark")
from predict_scorer import chrom_lengths, coding_cds  # noqa: E402

SOURCES = ["Augustus", "HiQ", "GeneMark", "snap", "pasa"]


def split_by_source(path, chroms):
    src = {}
    with open(path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) > 8 and c[2] == "mRNA":
                src[c[8].split("ID=")[1].split(";")[0].strip()] = c[1]
    models = coding_cds(path, False)
    out = collections.defaultdict(dict)
    for t, v in models.items():
        if v[0] in chroms:
            out[src.get(t, "?")][t] = v
    return {s: genes_from(m) for s, m in out.items()}


def match_spans(path, chroms):
    """Evidence alignments (match/cDNA_match/nucleotide_to_protein_match) as pseudo-genes."""
    spans = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[0] in chroms:
                spans[(c[0], c[6])].append((int(c[3]), int(c[4])))
    for k in spans:
        spans[k].sort()
    return spans


def span_hit(spans, ch, st, s, e):
    for a, b in spans.get((ch, st), ()):
        if a > e:
            return False
        if b >= s:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    for k in ("--ref-gff3", "--run-dir", "--genome", "--variant", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--mask", nargs="*", default=[])
    a = ap.parse_args()
    chroms = set(chrom_lengths(a.ref_gff3))
    mask = mask_spans(a.mask)
    misc = os.path.join(a.run_dir, "out", "predict_misc")
    ref = genes_from({t: v for t, v in coding_cds(a.ref_gff3, True).items() if v[0] in chroms})
    sets = split_by_source(os.path.join(misc, "gene_predictions.gff3"), chroms)
    sets["evm_round1"] = genes_from({t: v for t, v in coding_cds(os.path.join(misc, "evm.round1.gff3"), False).items() if v[0] in chroms})
    final = [p for p in glob.glob(os.path.join(a.run_dir, "out", "predict_results", "*.gff3"))][0]
    sets["final"] = genes_from({t: v for t, v in coding_cds(final, False).items() if v[0] in chroms})

    def scored(g):
        return not overlaps(mask, g["chrom"], g["start"], g["end"])

    ref_s = {k: v for k, v in ref.items() if scored(v)}
    ref_chain = {(g["chrom"], g["strand"], ch) for g in ref_s.values() for ch in g["tx"].values()}
    rows = []
    idx = {}
    for name, genes in sets.items():
        gs = {k: v for k, v in genes.items() if scored(v)}
        idx[name] = (Index(genes), genes)
        chains = {(g["chrom"], g["strand"], ch) for g in genes.values() for ch in g["tx"].values()}
        sn = sum(1 for g in ref_s.values() if any((g["chrom"], g["strand"], ch) in chains for ch in g["tx"].values()))
        pr = sum(1 for g in gs.values() if any((g["chrom"], g["strand"], ch) in ref_chain for ch in g["tx"].values()))
        rows.append((name, len(gs), sn, pr))

    tr = match_spans(os.path.join(misc, "transcript_alignments.gff3"), chroms)
    pa = os.path.join(misc, "pasa_predictions.gff3")
    pr_al = match_spans(os.path.join(misc, "protein_alignments.gff3"), chroms)

    def any_hit(name, g):
        ix, genes = idx[name]
        return any(gene_overlap(g, genes[p]) > 0 for p in ix.hits(g["chrom"], g["strand"], g["start"], g["end"]))

    miss = collections.Counter()
    for g in ref_s.values():
        if any_hit("final", g):
            continue
        miss["missed"] += 1
        in_any = [s for s in SOURCES if s in idx and any_hit(s, g)]
        if in_any:
            miss["input_overlap"] += 1
            if in_any == ["pasa"]:
                miss["input_only_pasa"] += 1
            if set(in_any) <= {"GeneMark", "snap"}:
                miss["input_only_GeneMark_snap"] += 1
        else:
            miss["no_input"] += 1
        if any_hit("evm_round1", g):
            miss["in_evm_round1"] += 1
        ev_t = span_hit(tr, g["chrom"], g["strand"], g["start"], g["end"])
        ev_p = span_hit(pr_al, g["chrom"], g["strand"], g["start"], g["end"])
        miss["transcript_evidence"] += ev_t
        miss["protein_evidence"] += ev_p
        miss["no_evidence"] += not (ev_t or ev_p)
        L = sum(e - s + 1 for s, e in max(g["tx"].values(), key=len))
        miss["short_lt300bp"] += L < 300
        miss["single_exon"] += all(len(ch) == 1 for ch in g["tx"].values())

    new = not os.path.exists(a.out) or os.path.getsize(a.out) == 0
    with open(a.out, "a") as o:
        if new:
            o.write("genome\tvariant\tset\tgenes_scored\tref_scored\texact_sn_n\texact_pr_n\tsn\tpr\n")
        for name, n, sn, pr in rows:
            o.write(f"{a.genome}\t{a.variant}\t{name}\t{n}\t{len(ref_s)}\t{sn}\t{pr}\t"
                    f"{100 * sn / len(ref_s):.1f}\t{100 * pr / max(n, 1):.1f}\n")
        o.write(f"{a.genome}\t{a.variant}\tMISSED\t" + ";".join(f"{k}={v}" for k, v in sorted(miss.items())) + "\n")
    del pa


if __name__ == "__main__":
    main()
