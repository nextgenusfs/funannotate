#!/usr/bin/env python3
"""Scan production PASA training models for the F1 (duplicate GFF3 row) signature.

For genome_annotation_training/<out>/training/funannotate_train.pasa.gff3,
report the file date, number of models, and the fraction of models whose
total CDS length is not a multiple of 3 (a complete or partial TransDecoder
ORF mapped correctly to the genome always has length % 3 == 0; rc.1 runs
with the F1 fix show 0). Needs no genome sequence.
Usage: production_f1_scan.py STAGEA_TSV TRAINING_ROOT OUT_TSV [NPROC]
"""
import sys, os, csv, time, collections
from multiprocessing import Pool
stagea, root, out = sys.argv[1:4]
nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 16
def one(name):
    f = os.path.join(root, name, "training", "funannotate_train.pasa.gff3")
    if not os.path.exists(f):
        return [name, "NA", 0, "NA", "NA"]
    ln = collections.Counter()
    with open(f) as fh:
        for l in fh:
            c = l.split("\t")
            if len(c) > 8 and c[2] == "CDS":
                p = c[8].split("Parent=", 1)[1].split(";", 1)[0].strip()
                ln[p] += int(c[4]) - int(c[3]) + 1
    n = len(ln)
    bad = sum(1 for v in ln.values() if v % 3)
    date = time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(f)))
    return [name, date, n, bad, f"{bad/n:.3f}" if n else "NA"]
names = [r["out"] for r in csv.DictReader(open(stagea), delimiter="\t") if r.get("aug_training") == "PASA"]
with Pool(nproc) as p:
    rows = p.map(one, names, chunksize=20)
with open(out, "w") as o:
    o.write("out\tpasa_gff3_date\tmodels\tcds_len_not_mod3\tfrac_not_mod3\n")
    for r in rows:
        o.write("\t".join(str(x) for x in r) + "\n")
print(f"wrote {len(rows)} rows")
