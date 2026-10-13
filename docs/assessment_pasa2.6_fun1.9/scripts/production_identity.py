#!/usr/bin/env python3
"""Per-genome transcript-to-genome identity for production PASA-trained genomes.

Reads column 6 (percent identity) of cDNA_match rows in
genome_annotation_training/<out>/training/transcript.alignments.gff3, one value
per transcript alignment (rows of one alignment share an ID; the first row's
value is used). Caveat: these files come from funannotate's pre-fix bam2gff3,
which dropped alignments below 80% identity and many spliced alignments whose
intron motif disagreed with the SAM flag; the distribution is truncated at 80%
and n_alignments is reported next to n_transcripts (Trinity FASTA headers).

Usage: production_identity.py STAGEA_TSV TRAINING_ROOT OUT_TSV [NPROC]
"""
import sys, os, csv, statistics
from multiprocessing import Pool
stagea, root, out = sys.argv[1:4]
nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 16

def one(name):
    d = os.path.join(root, name, "training")
    gff = os.path.join(d, "transcript.alignments.gff3")
    if not os.path.exists(gff):
        return [name, "NA"] + ["NA"] * 8
    seen, ids = set(), []
    try:
        with open(gff) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.split("\t", 9)
                if len(c) < 9 or c[2] != "cDNA_match":
                    continue
                aid = c[8].split("ID=", 1)[1].split(";", 1)[0] if "ID=" in c[8] else line
                if aid in seen:
                    continue
                seen.add(aid)
                ids.append(float(c[5]))
    except (OSError, ValueError):
        return [name, "ERR"] + ["NA"] * 8
    ntx = "NA"
    fa = os.path.join(d, "trinity.fasta")
    if os.path.exists(fa):
        with open(fa) as fh:
            ntx = sum(1 for l in fh if l.startswith(">"))
    if not ids:
        return [name, 0, ntx] + ["NA"] * 7
    ids.sort()
    n = len(ids)
    q = lambda f: ids[min(n - 1, int(f * n))]
    ge99 = sum(1 for x in ids if x >= 99) / n
    lt90 = sum(1 for x in ids if x < 90) / n
    med = statistics.median(ids)
    cat = "same_strain" if med >= 99 else ("divergent" if med >= 90 else "other_species")
    return [name, n, ntx, f"{med:.2f}", f"{q(0.10):.2f}", f"{q(0.90):.2f}",
            f"{ge99:.3f}", f"{1-ge99-lt90:.3f}", f"{lt90:.3f}", cat]

names = []
with open(stagea) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("aug_training") == "PASA":
            names.append(r["out"])
with Pool(nproc) as p:
    rows = p.map(one, names, chunksize=20)
with open(out, "w") as o:
    o.write("out\tn_alignments\tn_transcripts\tmedian_id\tp10_id\tp90_id\tfrac_ge99\tfrac_90_99\tfrac_lt90\tgate_category\n")
    for r in rows:
        o.write("\t".join(str(x) for x in r) + "\n")
print(f"wrote {len(rows)} rows to {out}")
