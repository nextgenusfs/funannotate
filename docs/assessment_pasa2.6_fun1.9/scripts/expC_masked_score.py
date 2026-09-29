#!/usr/bin/env python3
"""Score whole-genome predictions against RefSeq, leaving out every training gene (experiment C).

Experiment B scored held-out chromosomes. Experiment C trains on the whole genome, so there is
no held-out chromosome. Instead, the gene spans of all training models of all arms of a genome
(PASA final_training_models.gff3 and BUSCO busco.final.gff3) form a mask. RefSeq mRNAs and
predicted transcripts whose CDS span overlaps the mask (either strand) are removed; the rest are
scored with gffcompare exactly as predict_scorer.py does. Every arm of a genome is scored on the
same gene set, and no training gene is scored.

Usage:
  masked_score.py --ref-gff3 REF --pred-gff3 PRED --mask T1.gff3 [T2.gff3 ...]
                  --genome NAME --variant NAME [--chroms list.txt] [--out table.tsv]
--chroms restricts scoring to listed sequences (default: all non-mitochondrial RefSeq sequences).
"""
import argparse, bisect, collections, os, subprocess, sys, tempfile

sys.path.insert(0, "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/refseq_benchmark")
from predict_scorer import LEVELS, attrs, chrom_lengths, coding_cds, parse_stats, write_gtf, xopen  # noqa: E402


def mask_spans(paths):
    spans = collections.defaultdict(list)
    for p in paths:
        with xopen(p) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) > 8 and c[2] == "gene":
                    spans[c[0]].append((int(c[3]), int(c[4])))
    merged = {}
    for ch, iv in spans.items():
        iv.sort(); out = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= out[-1][1] + 1:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[ch] = ([s for s, _ in out], [e for _, e in out])
    return merged


def overlaps(merged, ch, s, e):
    if ch not in merged:
        return False
    starts, ends = merged[ch]
    i = bisect.bisect_right(starts, e) - 1
    return i >= 0 and ends[i] >= s


def keep(models, merged, chroms):
    return {t: v for t, v in models.items()
            if v[0] in chroms and not overlaps(merged, v[0], v[2][0][0], v[2][-1][1])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff3", required=True); ap.add_argument("--pred-gff3", required=True)
    ap.add_argument("--mask", nargs="*", default=[]); ap.add_argument("--genome", required=True)
    ap.add_argument("--variant", required=True); ap.add_argument("--chroms"); ap.add_argument("--out")
    ap.add_argument("--gffcompare", default="gffcompare")
    a = ap.parse_args()
    chroms = ({l.strip() for l in open(a.chroms) if l.strip()} if a.chroms
              else set(chrom_lengths(a.ref_gff3)))
    merged = mask_spans(a.mask)
    ref_all = {t: v for t, v in coding_cds(a.ref_gff3, True).items() if v[0] in chroms}
    ref = keep(ref_all, merged, chroms)
    pred = keep(coding_cds(a.pred_gff3, False), merged, chroms)
    with tempfile.TemporaryDirectory() as d:
        rg, pg = os.path.join(d, "ref.gtf"), os.path.join(d, "pred.gtf")
        write_gtf(ref, chroms, rg); write_gtf(pred, chroms, pg)
        subprocess.run([a.gffcompare, "-r", rg, "-o", os.path.join(d, "cmp"), pg], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        sp = os.path.join(d, "cmp")
        stats = parse_stats(sp + ".stats" if os.path.exists(sp + ".stats") else sp)
    row = collections.OrderedDict(genome=a.genome, variant=a.variant,
                                  ref_mrnas_total=len(ref_all), ref_mrnas_scored=len(ref),
                                  pred_genes_scored=len({v[3] for v in pred.values()}),
                                  mask_bp=sum(e - s + 1 for st, en in merged.values() for s, e in zip(st, en)))
    for lv in LEVELS:
        k = lv.lower().replace(" ", "_")
        row[k + "_sn"] = stats.get(k + "_sn", "NA"); row[k + "_pr"] = stats.get(k + "_pr", "NA")
    header = "\t".join(row) + "\n"; line = "\t".join(str(v) for v in row.values()) + "\n"
    if a.out:
        new = not os.path.exists(a.out) or os.path.getsize(a.out) == 0
        with open(a.out, "a") as o:
            if new:
                o.write(header)
            o.write(line)
    sys.stdout.write(header + line)


if __name__ == "__main__":
    main()
