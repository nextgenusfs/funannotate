#!/usr/bin/env python3
"""Score one or more gene-model training sets against a reference annotation.

For each labelled GFF3 (label=path), report size, accuracy and diversity:
  models, single_cds_pct, complete_pct (ATG..stop, no internal stop, len%3==0,
  from the genome sequence), exact_chain (CDS start, stop and all introns equal
  to a reference mRNA), exact_chain_pct (precision), distinct_ref_genes_exact
  (distinct reference genes with an exact match), distinct_ref_genes_hit
  (distinct reference genes overlapped in frame), redundant_models (models
  hitting a reference gene already hit), median_cds_len, pct_gt3_introns,
  exon_count_ks (Kolmogorov-Smirnov distance between the CDS-segment-count
  distributions of the set and the reference, 0 = identical).
The reference line gives the same statistics for the reference itself.

Usage: training_set_vs_refseq.py --ref-gff3 REF[.gz] --genome FA[.gz]
           LABEL=train1.gff3 [LABEL2=train2.gff3 ...] [--out table.tsv]
"""
import argparse
import bisect
import collections
import gzip
import re
import sys

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def xopen(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_models(path):
    """transcript -> (chrom, strand, segs in transcription order, gene)."""
    cds = collections.defaultdict(list)
    gene_of = {}
    with xopen(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or line.startswith("#"):
                continue
            if c[2] in ("mRNA", "transcript"):
                i = re.search(r"ID=([^;]+)", c[8])
                g = re.search(r"Parent=([^;]+)", c[8])
                if i:
                    gene_of[i.group(1)] = g.group(1) if g else i.group(1)
            elif c[2] == "CDS":
                p = re.search(r"Parent=([^;]+)", c[8])
                if p:
                    for par in p.group(1).split(","):
                        cds[par].append((c[0], int(c[3]), int(c[4]), c[6]))
    out = {}
    for t, v in cds.items():
        st = v[0][3]
        segs = sorted(((s, e) for _, s, e, _ in v), reverse=(st == "-"))
        out[t] = (v[0][0], st, segs, gene_of.get(t, t))
    return out


def read_fasta(p, wanted):
    seqs, name, buf = {}, None, []
    with xopen(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in wanted:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in wanted:
                buf.append(line.strip())
    if name in wanted:
        seqs[name] = "".join(buf).upper()
    return seqs


def cds_seq(genome, chrom, strand, segs):
    g = genome.get(chrom)
    if g is None:
        return None
    s = "".join(g[a - 1:b] for a, b in sorted(segs))
    return s.translate(COMP)[::-1] if strand == "-" else s


def is_complete(seq):
    if not seq or len(seq) % 3 or len(seq) < 6:
        return False
    if seq[:3] != "ATG" or seq[-3:] not in STOPS:
        return False
    return not any(seq[i:i + 3] in STOPS for i in range(0, len(seq) - 3, 3))


def chain_key(chrom, strand, segs):
    return (chrom, strand, tuple(sorted(segs)))


def ks(a, b):
    if not a or not b:
        return float("nan")
    ca, cb = collections.Counter(a), collections.Counter(b)
    fa = fb = 0.0
    d = 0.0
    for k in sorted(set(ca) | set(cb)):
        fa += ca[k] / len(a)
        fb += cb[k] / len(b)
        d = max(d, abs(fa - fb))
    return d


def frame_map(strand, segs):
    idx = {}
    p = 0
    for s, e in segs:
        rng = range(s, e + 1) if strand == "+" else range(e, s - 1, -1)
        for g in rng:
            idx[g] = p % 3
            p += 1
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff3", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out")
    ap.add_argument("sets", nargs="+", help="LABEL=path.gff3")
    a = ap.parse_args()
    ref = read_models(a.ref_gff3)
    sets = [(s.split("=", 1)[0], read_models(s.split("=", 1)[1])) for s in a.sets]
    chroms = {m[0] for m in ref.values()}
    for _, ms in sets:
        chroms |= {m[0] for m in ms.values()}
    genome = read_fasta(a.genome, chroms)

    ref_chain = {}
    by_chrom = collections.defaultdict(list)
    for t, (ch, st, segs, g) in ref.items():
        ref_chain.setdefault(chain_key(ch, st, segs), g)
        by_chrom[ch].append((min(x for x, _ in segs), max(y for _, y in segs), t))
    starts = {}
    for ch in by_chrom:
        by_chrom[ch].sort()
        starts[ch] = [x[0] for x in by_chrom[ch]]
    maxlen = max(hi - lo for L in by_chrom.values() for lo, hi, _ in L)
    ref_counts = [len(m[2]) for m in ref.values()]

    cols = ["set", "models", "single_cds_pct", "complete_pct", "exact_chain",
            "exact_chain_pct", "distinct_ref_genes_exact", "distinct_ref_genes_hit",
            "redundant_models", "median_cds_len", "pct_gt3_introns", "exon_count_ks"]
    rows = []

    def summarize(label, ms, compare):
        n = len(ms)
        counts = [len(m[2]) for m in ms.values()]
        lens = sorted(sum(b - a + 1 for a, b in m[2]) for m in ms.values())
        comp = sum(1 for m in ms.values() if is_complete(cds_seq(genome, m[0], m[1], m[2])))
        exact = 0
        genes_exact, genes_hit = set(), set()
        redundant = 0
        if compare:
            for ch, st, segs, _ in ms.values():
                g = ref_chain.get(chain_key(ch, st, segs))
                if g:
                    exact += 1
                    genes_exact.add(g)
                # in-frame overlap with any reference mRNA
                mf = frame_map(st, segs)
                lo, hi = min(x for x, _ in segs), max(y for _, y in segs)
                hit = None
                L = by_chrom.get(ch, [])
                i = bisect.bisect_left(starts.get(ch, []), lo - maxlen)
                while i < len(L) and L[i][0] <= hi and hit is None:
                    rlo, rhi, rt = L[i]
                    i += 1
                    if rhi < lo or ref[rt][1] != st:
                        continue
                    rf = frame_map(st, ref[rt][2])
                    if any(p in rf and rf[p] == f for p, f in mf.items()):
                        hit = ref[rt][3]
                if hit is not None:
                    if hit in genes_hit:
                        redundant += 1
                    genes_hit.add(hit)
        rows.append([
            label, n,
            "%.1f" % (100.0 * sum(1 for c in counts if c == 1) / n) if n else "NA",
            "%.1f" % (100.0 * comp / n) if n else "NA",
            exact if compare else "NA",
            "%.1f" % (100.0 * exact / n) if (compare and n) else "NA",
            len(genes_exact) if compare else "NA",
            len(genes_hit) if compare else "NA",
            redundant if compare else "NA",
            lens[len(lens) // 2] if lens else "NA",
            "%.1f" % (100.0 * sum(1 for c in counts if c > 4) / n) if n else "NA",
            "%.3f" % ks(counts, ref_counts),
        ])

    summarize("reference", ref, False)
    for label, ms in sets:
        summarize(label, ms, True)
    text = "\t".join(cols) + "\n" + "\n".join("\t".join(str(x) for x in r) for r in rows) + "\n"
    if a.out:
        open(a.out, "w").write(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
