#!/usr/bin/env python3
"""R13c: compare PASA training models with a reference annotation.

For each training model, find the same-strand reference mRNA whose CDS
shares the most bases with it in the same reading frame. Then:
  complete models: does the start, stop and intron chain match the reference?
  incomplete models: how many codons of reference CDS lie beyond each open
      end, and does that missing stretch cross a reference intron? A
      missing end within N codons and in the same exon is what a simple
      N-codon genome extension could rescue.

Usage: r13_reference_compare.py --train-gff3 F --td-genome-gff3 F --ref-gff3 F[.gz] --out PREFIX
"""
import argparse
import bisect
import collections
import gzip
import re
import sys
from urllib.parse import unquote


def xopen(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_cds(path, only_parent_re=None):
    cds = collections.defaultdict(list)
    with xopen(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "CDS":
                continue
            m = re.search(r"Parent=([^;]+)", c[8])
            if not m:
                continue
            cds[m.group(1)].append((c[0], int(c[3]), int(c[4]), c[6]))
    out = {}
    for k, v in cds.items():
        strand = v[0][3]
        segs = sorted(((s, e) for _, s, e, _ in v), reverse=(strand == "-"))
        out[k] = (v[0][0], strand, segs)
    return out


def orf_types(td):
    t = {}
    with xopen(td) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) > 8 and c[2] == "mRNA":
                m = re.search(r"ORF type:(\w+)", unquote(c[8]))
                if m:
                    t[re.search(r"ID=([^;]+)", c[8]).group(1)] = m.group(1)
    return t


def cds_index(strand, segs):
    """genome pos -> 0-based CDS position (transcription order)."""
    idx = {}
    p = 0
    for s, e in segs:
        rng = range(s, e + 1) if strand == "+" else range(e, s - 1, -1)
        for g in rng:
            idx[g] = p
            p += 1
    return idx, p


def introns(segs):
    lo = sorted(segs)
    return tuple((a[1] + 1, b[0] - 1) for a, b in zip(lo, lo[1:]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-gff3", required=True)
    ap.add_argument("--td-genome-gff3", required=True)
    ap.add_argument("--ref-gff3", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    types = orf_types(a.td_genome_gff3)
    models = read_cds(a.train_gff3)
    ref = read_cds(a.ref_gff3)

    # reference single-CDS baseline (one mRNA per gene is not enforced; all mRNAs)
    ref_single = sum(1 for _, _, s in ref.values() if len(s) == 1)

    # bin reference mRNAs by chrom for overlap search
    by_chrom = collections.defaultdict(list)
    for rid, (ch, st, segs) in ref.items():
        lo = min(s for s, _ in segs)
        hi = max(e for _, e in segs)
        by_chrom[ch].append((lo, hi, rid))
    starts = {}
    for ch in by_chrom:
        by_chrom[ch].sort()
        starts[ch] = [x[0] for x in by_chrom[ch]]
    maxlen = max(hi - lo for L in by_chrom.values() for lo, hi, _ in L)

    tally = collections.Counter()
    miss = collections.defaultdict(list)
    rows = []
    for mid, (ch, st, segs) in models.items():
        t = types.get(mid, "NA")
        mlo = min(s for s, _ in segs)
        mhi = max(e for _, e in segs)
        midx, mlen = cds_index(st, segs)
        best = None
        L = by_chrom.get(ch, [])
        i = bisect.bisect_left(starts.get(ch, []), mlo - maxlen)
        while i < len(L) and L[i][0] <= mhi:
            lo, hi, rid = L[i]
            i += 1
            if hi < mlo or ref[rid][1] != st:
                continue
            ridx, rlen = cds_index(st, ref[rid][2])
            inframe = sum(1 for g, p in midx.items() if g in ridx and ridx[g] % 3 == p % 3)
            if inframe and (best is None or inframe > best[0]):
                best = (inframe, rid, ridx, rlen)
        if best is None:
            cat = "no_inframe_reference"
            tally[(t, cat)] += 1
            rows.append((mid, t, "NA", cat, "NA", "NA"))
            continue
        _, rid, ridx, rlen = best
        rsegs = ref[rid][2]
        first_g = segs[0][0] if st == "+" else segs[0][1]
        last_g = segs[-1][1] if st == "+" else segs[-1][0]
        # CDS position of model ends in the reference
        p5 = ridx.get(first_g)
        p3 = ridx.get(last_g)
        if t == "complete":
            same_start = p5 == 0
            same_stop = p3 == rlen - 1
            same_chain = introns(segs) == introns(rsegs)
            cat = "exact" if (same_start and same_stop and same_chain) else (
                "ends_match_chain_differs" if (same_start and same_stop) else "ends_differ")
            tally[(t, cat)] += 1
            rows.append((mid, t, rid, cat, "NA", "NA"))
            continue
        rec5 = rec3 = "NA"
        if t in ("5prime_partial", "internal"):
            if p5 is None:
                rec5 = "end_outside_ref_cds"
            else:
                codons = p5 // 3
                # same exon: all missing reference CDS bases are contiguous with first_g
                seg = next((s, e) for s, e in rsegs if s <= first_g <= e)
                same_exon = (codons * 3) <= ((first_g - seg[0]) if st == "+" else (seg[1] - first_g))
                rec5 = "{}:{}".format(codons, "same_exon" if same_exon else "crosses_intron")
                miss[(t, "start", "same_exon" if same_exon else "crosses_intron")].append(codons)
        if t in ("3prime_partial", "internal"):
            if p3 is None:
                rec3 = "end_outside_ref_cds"
            else:
                codons = (rlen - 1 - p3) // 3
                seg = next((s, e) for s, e in rsegs if s <= last_g <= e)
                same_exon = (codons * 3) <= ((seg[1] - last_g) if st == "+" else (last_g - seg[0]))
                rec3 = "{}:{}".format(codons, "same_exon" if same_exon else "crosses_intron")
                miss[(t, "stop", "same_exon" if same_exon else "crosses_intron")].append(codons)
        tally[(t, "has_inframe_reference")] += 1
        rows.append((mid, t, rid, "incomplete", rec5, rec3))

    with gzip.open(a.out + ".per_model.tsv.gz", "wt") as out:
        out.write("model\torf_type\tref_mrna\tcategory\tmissing_5p\tmissing_3p\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    def q(v, f):
        v = sorted(v)
        return v[int(f * (len(v) - 1))] if v else "NA"

    with open(a.out + ".summary.tsv", "w") as s:
        s.write("reference_mrnas\t{}\nreference_single_cds\t{} ({:.1f}%)\n".format(
            len(ref), ref_single, 100.0 * ref_single / len(ref)))
        tm = collections.Counter(types.get(m, "NA") for m in models)
        ms = sum(1 for _, _, sg in models.values() if len(sg) == 1)
        s.write("training_models\t{}\ntraining_single_cds\t{} ({:.1f}%)\n".format(
            len(models), ms, 100.0 * ms / len(models)))
        s.write("\norf_type\tcategory\tcount\n")
        for (t, c), n in sorted(tally.items()):
            s.write("{}\t{}\t{}\n".format(t, c, n))
        s.write("\norf_type\tmissing_end\tlocation\tn\tmedian_codons\tp25\tp75\tle6\tle9\n")
        for key in sorted(miss):
            v = miss[key]
            s.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                key[0], key[1], key[2], len(v), q(v, .5), q(v, .25), q(v, .75),
                sum(1 for x in v if x <= 6), sum(1 for x in v if x <= 9)))
    sys.stdout.write(open(a.out + ".summary.tsv").read())


if __name__ == "__main__":
    main()
