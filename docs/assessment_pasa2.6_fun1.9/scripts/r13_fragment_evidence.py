#!/usr/bin/env python3
"""R13d: for incomplete training models, is the missing part of the reference
CDS covered by other evidence (union of PASA assemblies, union of transcript
alignments, the single best transcript alignment)?

Coverage is genomic base overlap with the reference CDS, same strand for
PASA assemblies (their strand is PASA's), either strand for transcript
alignments (unstranded data). Note: transcript alignment coordinates from
funannotate's pre-fix bam2gff3 can be shifted by a few bases; that does not
matter at the scale measured here.

Usage: r13_fragment_evidence.py --refcmp PREFIX.per_model.tsv.gz --train-gff3 F
         --ref-gff3 F --pasa-gff3 F --tx-gff3 F --out PREFIX
"""
import argparse, collections, gzip, re, sys

def xopen(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def feats(path, ftype=None, key=r"Parent=([^;]+)", need_type=True):
    d = collections.defaultdict(list)
    with xopen(path) as fh:
        for l in fh:
            c = l.rstrip("\n").split("\t")
            if len(c) < 9 or l.startswith("#"): continue
            if ftype and c[2] != ftype: continue
            m = re.search(key, c[8])
            if m: d[m.group(1)].append((c[0], int(c[3]), int(c[4]), c[6]))
    return d

def bases(segs): 
    s = set()
    for _, a, b, _ in segs: s.update(range(a, b + 1))
    return s

def main():
    ap = argparse.ArgumentParser()
    for x in ("refcmp", "train-gff3", "ref-gff3", "pasa-gff3", "tx-gff3", "out"):
        ap.add_argument("--" + x, required=True)
    a = ap.parse_args()
    ref = feats(a.ref_gff3, "CDS")
    mod = feats(a.train_gff3, "CDS")
    pasa = feats(a.pasa_gff3, None, r"Target=(\S+)")
    tx = feats(a.tx_gff3, None, r"ID=([^;]+)")
    # index evidence by chrom
    def index(d):
        idx = collections.defaultdict(list)
        for k, segs in d.items():
            idx[segs[0][0]].append((min(s[1] for s in segs), max(s[2] for s in segs), k))
        return idx
    pidx, tidx = index(pasa), index(tx)
    cats = collections.Counter(); rows = []
    with xopen(a.refcmp) as fh:
        next(fh)
        for l in fh:
            m, t, rid, cat, m5, m3 = l.rstrip("\n").split("\t")
            if cat != "incomplete" or rid not in ref: continue
            rsegs = ref[rid]; ch, st = rsegs[0][0], rsegs[0][3]
            rb = bases(rsegs); lo, hi = min(rb), max(rb)
            fm = len(rb & bases(mod[m])) / len(rb)
            pu = set(); tu = set(); tbest = 0.0
            for s, e, k in pidx.get(ch, []):
                if e >= lo and s <= hi and pasa[k][0][3] == st: pu |= rb & bases(pasa[k])
            for s, e, k in tidx.get(ch, []):
                if e >= lo and s <= hi:
                    b = rb & bases(tx[k]); tu |= b; tbest = max(tbest, len(b) / len(rb))
            fp, ft = len(pu) / len(rb), len(tu) / len(rb)
            if ft < 0.9: c = "no_evidence_for_missing_part"
            elif tbest >= 0.9: c = "one_transcript_covers_gene"
            else: c = "pieces_cover_gene_no_single_transcript"
            cats[c] += 1
            rows.append((m, t, rid, "%.2f" % fm, "%.2f" % fp, "%.2f" % ft, "%.2f" % tbest, c))
    with gzip.open(a.out + ".per_model.tsv.gz", "wt") as o:
        o.write("model\torf_type\tref\tmodel_cov\tpasa_union_cov\ttx_union_cov\tbest_tx_cov\tclass\n")
        for r in rows: o.write("\t".join(r) + "\n")
    def med(i): 
        v = sorted(float(r[i]) for r in rows); return v[len(v)//2] if v else float("nan")
    with open(a.out + ".summary.tsv", "w") as s:
        s.write("incomplete_models_with_ref\t%d\n" % len(rows))
        s.write("median_ref_cds_cov_model\t%.2f\nmedian_ref_cds_cov_pasa_union\t%.2f\n"
                "median_ref_cds_cov_tx_union\t%.2f\nmedian_ref_cds_cov_best_tx\t%.2f\n" % (med(3), med(4), med(5), med(6)))
        for k, v in sorted(cats.items()): s.write("%s\t%d\t%.1f%%\n" % (k, v, 100.0 * v / len(rows)))
    sys.stdout.write(open(a.out + ".summary.tsv").read())

if __name__ == "__main__": main()
