#!/usr/bin/env python3
"""R13b: could incomplete training models be completed by extending their
open end(s) a few codons on the genome?

For each model in the training GFF3 whose TransDecoder ORF type is not
'complete', read the genome in frame beyond the open end:
  missing start (5prime_partial, internal): walk upstream codon by codon,
      success = ATG reached before any in-frame stop.
  missing stop  (3prime_partial, internal): walk downstream,
      success = first in-frame stop codon.
The walk assumes the gene continues colinearly on the genome (no intron in
the extension), which is only a fair assumption for short extensions.

Null baseline: the same walk at the same site but in the two other reading
frames (offset +1 and +2 nt). A real rescue should clearly beat the null.

Usage: r13_genome_extension.py --train-gff3 F --td-genome-gff3 F --genome FA[.gz]
           [--per-assembly r13.per_assembly.tsv.gz] [--max-codons 30] --out PREFIX
With --per-assembly (from r13_truncation_analysis.py), results are also split
by R13 class (member_has_complete / no_member_complete / members_unscored).
"""
import argparse
import collections
import gzip
import re
import sys
from urllib.parse import unquote

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def xopen(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


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


def rc(s):
    return s.translate(COMP)[::-1]


def orf_types(td_gff3):
    t = {}
    with xopen(td_gff3) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[2] == "mRNA":
                i = re.search(r"ID=([^;]+)", c[8]).group(1)
                m = re.search(r"ORF type:(\w+)", unquote(c[8]))
                if m:
                    t[i] = m.group(1)
    return t


def models(train_gff3):
    cds = collections.defaultdict(list)
    with xopen(train_gff3) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) > 8 and c[2] == "CDS":
                p = re.search(r"Parent=([^;]+)", c[8]).group(1)
                cds[p].append((c[0], int(c[3]), int(c[4]), c[6], c[7]))
    return cds


def codon_at(seq, strand, pos):
    """Codon whose transcript-5' base is at genome 1-based pos (+ strand: pos..pos+2;
    - strand: pos-2..pos reverse-complemented)."""
    if strand == "+":
        if pos < 1 or pos + 2 > len(seq):
            return None
        return seq[pos - 1:pos + 2]
    if pos - 2 < 1 or pos > len(seq):
        return None
    return rc(seq[pos - 3:pos])


def walk(seq, strand, first_pos, step, want, maxc):
    """Walk codons starting at first_pos, moving `step` nt per codon along the
    transcript direction. Returns ('hit', k) / ('blocked', k) / ('none', maxc)."""
    pos = first_pos
    for k in range(1, maxc + 1):
        cod = codon_at(seq, strand, pos)
        if cod is None or "N" in cod:
            return ("edge", k)
        if want == "start":
            if cod == "ATG":
                return ("hit", k)
            if cod in STOPS:
                return ("blocked", k)
        else:
            if cod in STOPS:
                return ("hit", k)
        pos += step
    return ("none", maxc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-gff3", required=True)
    ap.add_argument("--td-genome-gff3", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--per-assembly")
    ap.add_argument("--max-codons", type=int, default=30)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    types = orf_types(a.td_genome_gff3)
    cds = models(a.train_gff3)
    klass = {}
    if a.per_assembly:
        with xopen(a.per_assembly) as fh:
            next(fh)
            for line in fh:
                c = line.split("\t")
                klass[c[0]] = c[3]
    genome = read_fasta(a.genome, {v[0][0] for v in cds.values()})
    K = (6, 9, a.max_codons)
    # tallies[(group, end, frame)] -> Counter of 'hit<=k'
    tal = collections.defaultdict(collections.Counter)
    rows = []
    for tid, parts in cds.items():
        t = types.get(tid)
        if t is None or t == "complete":
            continue
        chrom, strand = parts[0][0], parts[0][3]
        seq = genome.get(chrom)
        if seq is None:
            continue
        lo = min(p[1] for p in parts)
        hi = max(p[2] for p in parts)
        asmbl = tid.split(".")[0]
        groups = ["all:" + t]
        if klass:
            groups.append("{}:{}".format(klass.get(asmbl, "NA"), t))
        res = {}
        for end in ("start", "stop"):
            if end == "start" and t not in ("5prime_partial", "internal"):
                continue
            if end == "stop" and t not in ("3prime_partial", "internal"):
                continue
            for off in (0, 1, 2):
                if end == "start":
                    # first codon upstream of the CDS 5' end, in frame (+off)
                    if strand == "+":
                        first, step = lo - 3 + off, -3
                    else:
                        first, step = hi + 3 - off, 3
                else:
                    if strand == "+":
                        first, step = hi + 1 + off, 3
                    else:
                        first, step = lo - 1 - off, -3
                st, k = walk(seq, strand, first, step, end, a.max_codons)
                if off == 0:
                    res[end] = (st, k)
                for g in groups:
                    key = (g, end, "inframe" if off == 0 else "null")
                    tal[key]["n"] += 1
                    for kk in K:
                        if st == "hit" and k <= kk:
                            tal[key]["hit<={}".format(kk)] += 1
                    if st == "blocked":
                        tal[key]["blocked"] += 1
        rows.append((tid, t, klass.get(asmbl, "NA"),
                     "{}:{}".format(*res["start"]) if "start" in res else "NA",
                     "{}:{}".format(*res["stop"]) if "stop" in res else "NA"))
    with gzip.open(a.out + ".per_model.tsv.gz", "wt") as out:
        out.write("model\torf_type\tr13_class\tstart_walk\tstop_walk\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    with open(a.out + ".summary.tsv", "w") as s:
        s.write("group\tend\tframe\tn\t" + "\t".join("hit<={}".format(k) for k in K)
                + "\t" + "\t".join("pct<={}".format(k) for k in K) + "\tblocked\n")
        for key in sorted(tal):
            c = tal[key]
            n = c["n"]
            hits = [c["hit<={}".format(k)] for k in K]
            s.write("\t".join([key[0], key[1], key[2], str(n)] + [str(h) for h in hits]
                              + ["{:.1f}".format(100.0 * h / n) if n else "NA" for h in hits]
                              + [str(c["blocked"])]) + "\n")
    sys.stdout.write(open(a.out + ".summary.tsv").read())


if __name__ == "__main__":
    main()
