#!/usr/bin/env python3
"""Score gene predictions against RefSeq on held-out chromosomes with gffcompare.

Owner: the PASA-side review session (CODE_REVIEW_20260925). Shared with the
predict arms run by the selection session.

Both RefSeq and the predictions are reduced to protein-coding CDS, written
as GTF with each CDS segment as an 'exon' (CDS-level evaluation: UTRs are
ignored). Only the held-out chromosomes are kept, so genes used to train
Augustus/SNAP are not scored.

Subcommands:
  split  --ref-gff3 REF --out-prefix P [--min-len 1000000]
         Sort chromosomes (>= min-len, excluding mitochondrion) by length and
         alternate them into P.train_chroms.txt and P.holdout_chroms.txt.
         Filter each training set to the train chromosomes before training.
  score  --ref-gff3 REF --pred-gff3 PRED --holdout P.holdout_chroms.txt
         --genome NAME --variant NAME [--out table.tsv] [--gffcompare gffcompare]
         Appends one row: genome, variant, genes, single_exon_pct and
         gffcompare sensitivity/precision at base, exon, intron, intron-chain,
         transcript and locus level.
"""
import argparse
import collections
import gzip
import os
import re
import subprocess
import sys
import tempfile

LEVELS = ["Base", "Exon", "Intron", "Intron chain", "Transcript", "Locus"]


def xopen(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def attrs(col9):
    return dict(x.split("=", 1) for x in col9.strip().strip(";").split(";") if "=" in x)


def coding_cds(path, refseq):
    """transcript -> (chrom, strand, [(s, e)]). For RefSeq keep only mRNAs of
    protein_coding genes; for predictions keep every CDS parent."""
    coding_genes, mrna_gene = set(), {}
    cds = collections.defaultdict(list)
    with xopen(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            a = attrs(c[8])
            if c[2] == "gene" and a.get("gene_biotype") == "protein_coding":
                coding_genes.add(a.get("ID"))
            elif c[2] in ("mRNA", "transcript") and "ID" in a:
                mrna_gene[a["ID"]] = a.get("Parent", a["ID"])
            elif c[2] == "CDS" and "Parent" in a:
                for p in a["Parent"].split(","):
                    cds[p].append((c[0], c[6], int(c[3]), int(c[4])))
    out = {}
    for t, parts in cds.items():
        if refseq and mrna_gene.get(t) not in coding_genes:
            continue
        out[t] = (parts[0][0], parts[0][1], sorted((s, e) for _, _, s, e in parts),
                  mrna_gene.get(t, t))
    return out


def write_gtf(models, chroms, path):
    with open(path, "w") as o:
        for t, (ch, st, segs, g) in sorted(models.items()):
            if ch not in chroms:
                continue
            for s, e in segs:
                o.write('{}\tcds\texon\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}";\n'
                        .format(ch, s, e, st, g, t))


def chrom_lengths(ref_gff3):
    lens = {}
    with xopen(ref_gff3) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                _, ch, s, e = line.split()[:4]
                lens[ch] = int(e)
            elif not line.startswith("#"):
                c = line.split("\t")
                if len(c) > 8 and c[2] == "region":
                    a = attrs(c[8])
                    if a.get("genome") == "mitochondrion":
                        lens.pop(c[0], None)
                        lens[c[0] + "#mito"] = 0
    return {k: v for k, v in lens.items() if not k.endswith("#mito")}


def do_split(a):
    lens = chrom_lengths(a.ref_gff3)
    big = sorted((ch for ch, l in lens.items() if l >= a.min_len), key=lambda c: -lens[c])
    train, hold = big[0::2], big[1::2]
    open(a.out_prefix + ".train_chroms.txt", "w").write("\n".join(train) + "\n")
    open(a.out_prefix + ".holdout_chroms.txt", "w").write("\n".join(hold) + "\n")
    print("train:", ",".join(train))
    print("holdout:", ",".join(hold))


def parse_stats(path):
    res = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*(Base|Exon|Intron chain|Intron|Transcript|Locus) level:\s+([\d.]+)\s+\|\s+([\d.]+)", line)
            if m:
                k = m.group(1).lower().replace(" ", "_")
                res[k + "_sn"] = float(m.group(2))
                res[k + "_pr"] = float(m.group(3))
    return res


def do_score(a):
    hold = {l.strip() for l in open(a.holdout) if l.strip()}
    ref = coding_cds(a.ref_gff3, True)
    pred = coding_cds(a.pred_gff3, False)
    pred_h = {t: v for t, v in pred.items() if v[0] in hold}
    with tempfile.TemporaryDirectory() as d:
        rg, pg = os.path.join(d, "ref.gtf"), os.path.join(d, "pred.gtf")
        write_gtf(ref, hold, rg)
        write_gtf(pred, hold, pg)
        subprocess.run([a.gffcompare, "-r", rg, "-o", os.path.join(d, "cmp"), pg],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        sp = os.path.join(d, "cmp")
        stats = parse_stats(sp + ".stats" if os.path.exists(sp + ".stats") else sp)
    genes = {v[3] for v in pred_h.values()}
    single = sum(1 for v in pred_h.values() if len(v[2]) == 1)
    row = collections.OrderedDict()
    row["genome"] = a.genome
    row["variant"] = a.variant
    row["holdout_ref_mrnas"] = sum(1 for v in ref.values() if v[0] in hold)
    row["pred_genes"] = len(genes)
    row["pred_single_exon_pct"] = round(100.0 * single / len(pred_h), 1) if pred_h else 0.0
    for lv in LEVELS:
        k = lv.lower().replace(" ", "_")
        row[k + "_sn"] = stats.get(k + "_sn", "NA")
        row[k + "_pr"] = stats.get(k + "_pr", "NA")
    header = "\t".join(row) + "\n"
    line = "\t".join(str(v) for v in row.values()) + "\n"
    if a.out:
        new = not os.path.exists(a.out) or os.path.getsize(a.out) == 0
        with open(a.out, "a") as o:
            if new:
                o.write(header)
            o.write(line)
    sys.stdout.write(header + line)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split")
    s.add_argument("--ref-gff3", required=True)
    s.add_argument("--out-prefix", required=True)
    s.add_argument("--min-len", type=int, default=1000000)
    c = sub.add_parser("score")
    c.add_argument("--ref-gff3", required=True)
    c.add_argument("--pred-gff3", required=True)
    c.add_argument("--holdout", required=True)
    c.add_argument("--genome", required=True)
    c.add_argument("--variant", required=True)
    c.add_argument("--out")
    c.add_argument("--gffcompare", default="gffcompare")
    a = ap.parse_args()
    do_split(a) if a.cmd == "split" else do_score(a)


if __name__ == "__main__":
    main()
