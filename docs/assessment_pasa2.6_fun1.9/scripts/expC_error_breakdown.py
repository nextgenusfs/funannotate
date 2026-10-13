#!/usr/bin/python3.12
"""Error breakdown of final gene models against RefSeq (recommendation #2, 2026-09-30).

For one predict run, classify every scored RefSeq protein-coding gene (same mask and same
sequences as masked_score.py) by how the final models represent it, and check whether any EVM
input (predict_misc/gene_predictions.gff3: Augustus, HiQ, GeneMark, snap, pasa) already had an
exact CDS chain for it.

RefSeq gene classes (first match wins; overlap = same-strand CDS bp overlap >= 1):
  exact          a final model has the exact CDS chain of a RefSeq isoform
  missed         no final model overlaps the gene on the same strand
                 (missed_masked: the only overlapping final model was removed by the mask)
  merged         an overlapping final model also overlaps another RefSeq gene (fusion)
  split          >= 2 final models overlap the gene (and none is a fusion)
  ends_wrong     1:1; same intron chain, different start or stop
  exons_fewer    1:1; final model has fewer CDS segments than the best RefSeq isoform
  exons_more     1:1; more CDS segments
  splice_diff    1:1; same number of CDS segments, an internal boundary differs
Final-model classes: fp_unmatched = scored final model that overlaps no RefSeq gene on either
strand (may be a real gene missing from RefSeq; not an error by itself).

The candidate check ("input_exact") asks: for a non-exact RefSeq gene, did any EVM input model
have the exact CDS chain? Yes -> the combination/selection step lost it; no -> no input had it.

Usage: error_breakdown.py --ref-gff3 REF --run-dir RUN --mask M1 [M2 ...] --genome G --variant V --out TSV
"""
import argparse
import collections
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from masked_score import mask_spans, overlaps  # noqa: E402

sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/refseq_benchmark")
from predict_scorer import chrom_lengths, coding_cds  # noqa: E402

CLASSES = ["exact", "ends_wrong", "splice_diff", "exons_fewer", "exons_more",
           "split", "merged", "missed", "missed_masked"]


def introns(cds):
    return tuple((cds[i][1], cds[i + 1][0]) for i in range(len(cds) - 1))


class Index:
    """Per (chrom, strand) sorted gene spans for overlap lookups."""

    def __init__(self, genes):
        self.by = collections.defaultdict(list)
        for g, v in genes.items():
            self.by[(v["chrom"], v["strand"])].append((v["start"], v["end"], g))
        for k in self.by:
            self.by[k].sort()

    def hits(self, chrom, strand, s, e):
        out = []
        for gs, ge, g in self.by.get((chrom, strand), ()):
            if gs > e:
                break
            if ge >= s:
                out.append(g)
        return out


def genes_from(models):
    genes = {}
    for t, (ch, strand, cds, gene) in models.items():
        g = genes.setdefault(gene, {"chrom": ch, "strand": strand, "start": cds[0][0],
                                    "end": cds[-1][1], "tx": {}})
        g["start"] = min(g["start"], cds[0][0])
        g["end"] = max(g["end"], cds[-1][1])
        g["tx"][t] = tuple(cds)
    return genes


def cds_overlap(a, b):
    """bp overlap between two CDS chains."""
    n, i, j = 0, 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s <= e:
            n += e - s + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def gene_overlap(g1, g2):
    return max(cds_overlap(x, y) for x in g1["tx"].values() for y in g2["tx"].values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff3", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--mask", nargs="*", default=[])
    ap.add_argument("--genome", required=True)
    ap.add_argument("--variant", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-gene", help="optional per-gene TSV")
    a = ap.parse_args()

    chroms = set(chrom_lengths(a.ref_gff3))
    merged_mask = mask_spans(a.mask)
    pred_gff = [p for p in glob.glob(os.path.join(a.run_dir, "out", "predict_results", "*.gff3"))][0]
    ref = genes_from({t: v for t, v in coding_cds(a.ref_gff3, True).items() if v[0] in chroms})
    pred = genes_from({t: v for t, v in coding_cds(pred_gff, False).items() if v[0] in chroms})
    inputs = genes_from({t: v for t, v in coding_cds(
        os.path.join(a.run_dir, "out", "predict_misc", "gene_predictions.gff3"), False).items()
        if v[0] in chroms})
    input_chains = collections.defaultdict(set)  # (chrom, strand, chain) -> sources
    src_of = {}
    with open(os.path.join(a.run_dir, "out", "predict_misc", "gene_predictions.gff3")) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) > 8 and c[2] == "mRNA":
                tid = c[8].split("ID=")[1].split(";")[0].strip()
                src_of[tid] = c[1]
    for g in inputs.values():
        for t, chain in g["tx"].items():
            input_chains[(g["chrom"], g["strand"], chain)].add(src_of.get(t, "?"))

    def scored(g):
        return not overlaps(merged_mask, g["chrom"], g["start"], g["end"])

    ref_idx, pred_idx = Index(ref), Index(pred)
    counts = collections.Counter()
    input_exact = collections.Counter()
    input_src = collections.Counter()
    per_gene = []
    for rid, r in ref.items():
        if not scored(r):
            continue
        cand = [p for p in pred_idx.hits(r["chrom"], r["strand"], r["start"], r["end"])
                if gene_overlap(r, pred[p]) > 0]
        ref_chains = set(r["tx"].values())
        exact = any(ch in ref_chains for p in cand for ch in pred[p]["tx"].values())
        if exact:
            cls = "exact"
        elif not cand:
            cls = "missed"
        else:
            live = [p for p in cand if scored(pred[p])]
            fusion = any(
                any(o != rid and gene_overlap(ref[o], pred[p]) > 0
                    for o in ref_idx.hits(pred[p]["chrom"], pred[p]["strand"],
                                          pred[p]["start"], pred[p]["end"]))
                for p in cand)
            if not live:
                cls = "missed_masked"
            elif fusion:
                cls = "merged"
            elif len(cand) >= 2:
                cls = "split"
            else:
                pch = max(pred[cand[0]]["tx"].values(), key=len)
                rch = max(ref_chains, key=lambda x: cds_overlap(x, pch))
                if introns(pch) == introns(rch):
                    cls = "ends_wrong"
                elif len(pch) < len(rch):
                    cls = "exons_fewer"
                elif len(pch) > len(rch):
                    cls = "exons_more"
                else:
                    cls = "splice_diff"
        counts[cls] += 1
        had = set()
        if cls != "exact":
            for ch in ref_chains:
                had |= input_chains.get((r["chrom"], r["strand"], ch), set())
            if had:
                input_exact[cls] += 1
                for s in had:
                    input_src[s] += 1
        per_gene.append((rid, cls, ",".join(sorted(had))))

    fp = sum(1 for pid, p in pred.items() if scored(p) and not any(
        gene_overlap(p, ref[o]) > 0
        for st in ("+", "-")
        for o in ref_idx.hits(p["chrom"], st, p["start"], p["end"])))
    n_scored_pred = sum(1 for p in pred.values() if scored(p))

    total = sum(counts.values())
    row = collections.OrderedDict(genome=a.genome, variant=a.variant, ref_genes_scored=total,
                                  pred_genes_scored=n_scored_pred, fp_unmatched=fp)
    for c in CLASSES:
        row[c] = counts[c]
    for c in CLASSES[1:]:
        row[f"{c}_input_exact"] = input_exact[c]
    for s in ("Augustus", "HiQ", "GeneMark", "snap", "pasa"):
        row[f"input_exact_src_{s}"] = input_src[s]
    new = not os.path.exists(a.out) or os.path.getsize(a.out) == 0
    with open(a.out, "a") as o:
        if new:
            o.write("\t".join(row) + "\n")
        o.write("\t".join(str(v) for v in row.values()) + "\n")
    if a.per_gene:
        with open(a.per_gene, "w") as o:
            o.write("ref_gene\tclass\tinput_sources_with_exact_chain\n")
            for x in per_gene:
                o.write("\t".join(x) + "\n")
    print("\t".join(f"{k}={v}" for k, v in row.items()))


if __name__ == "__main__":
    main()
