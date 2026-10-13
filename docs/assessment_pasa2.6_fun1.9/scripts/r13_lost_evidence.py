#!/usr/bin/env python3
"""R13e: where is transcript evidence lost inside PASA?

1. Validation outcome per aligner (prog), split by single- vs multi-segment
   alignments, plus the top failure reasons (from alignment.validations.output).
2. Cross-tab of per-transcript validity across aligners.
3. For genes flagged by r13_fragment_evidence.py as "transcripts cover the
   gene but PASA assemblies do not" (tx_union_cov >= 0.9, pasa_union_cov < 0.9),
   classify the transcripts that overlap the gene:
     all_failed_validation  every overlapping transcript failed in every aligner
     some_valid_not_assembled  at least one valid alignment exists, but the
                               PASA assemblies still do not cover the gene
                               (clustering / assembly / subclustering loss)

Usage: r13_lost_evidence.py --validations alignment.validations.output[.gz]
         --fragev fragev.per_model.tsv.gz --ref-gff3 REF --tx-gff3 trinity.alignments.gff3
         --out PREFIX
"""
import argparse
import collections
import gzip
import re
import sys


def xopen(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def main():
    ap = argparse.ArgumentParser()
    for x in ("validations", "fragev", "ref-gff3", "tx-gff3", "out"):
        ap.add_argument("--" + x, required=True)
    a = ap.parse_args()

    # 1+2: validation outcomes
    outcome = collections.Counter()
    reasons = collections.defaultdict(collections.Counter)
    status = collections.defaultdict(dict)  # acc -> prog -> set of OK/ERROR
    with xopen(a.validations) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) != 15 or c[0].startswith("#"):
                continue
            prog, acc, nseg, valid, comment = c[0], c[1], int(c[5]), c[8], c[14]
            seg = "single" if nseg == 1 else "multi"
            outcome[(prog, seg, valid)] += 1
            if valid != "OK":
                reasons[prog][re.sub(r"\d+(\.\d+)?", "N", comment)[:80]] += 1
            status[acc].setdefault(prog, set()).add(valid)
    progs = sorted({k[0] for k in outcome})

    def acc_state(acc, prog):
        s = status.get(acc, {}).get(prog)
        if not s:
            return "NA"
        return "OK" if "OK" in s else "ERROR"

    xt = collections.Counter(tuple(acc_state(acc, p) for p in progs) for acc in status)

    # 3: lost genes
    ref = collections.defaultdict(list)
    with xopen(a.ref_gff3) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) > 8 and c[2] == "CDS":
                m = re.search(r"Parent=([^;]+)", c[8])
                if m:
                    ref[m.group(1)].append((c[0], int(c[3]), int(c[4])))
    tx = collections.defaultdict(list)
    with xopen(a.tx_gff3) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) > 8 and not line.startswith("#"):
                m = re.search(r"ID=([^;]+)", c[8])
                if m:
                    tx[c[0]].append((int(c[3]), int(c[4]), m.group(1)))
    lost = []
    with xopen(a.fragev) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if float(c[5]) >= 0.9 and float(c[4]) < 0.9:
                lost.append(c[2])
    gene_class = collections.Counter()
    tx_states = collections.Counter()
    for rid in set(lost):
        segs = ref.get(rid)
        if not segs:
            continue
        ch = segs[0][0]
        lo, hi = min(s for _, s, _ in segs), max(e for _, _, e in segs)
        accs = {acc for s, e, acc in tx.get(ch, []) if e >= lo and s <= hi}
        states = [tuple(acc_state(x, p) for p in progs) for x in accs]
        for st in states:
            tx_states[st] += 1
        if any("OK" in st for st in states):
            gene_class["some_valid_not_assembled"] += 1
        else:
            gene_class["all_failed_validation"] += 1

    with open(a.out + ".summary.tsv", "w") as s:
        s.write("## validation outcome: prog\tsegments\tresult\tcount\n")
        for k in sorted(outcome):
            s.write("{}\t{}\t{}\t{}\n".format(*k, outcome[k]))
        for p in progs:
            s.write("\n## top failure reasons: {}\n".format(p))
            for r, n in reasons[p].most_common(5):
                s.write("{}\t{}\n".format(n, r))
        s.write("\n## per-transcript status across aligners ({})\n".format(" | ".join(progs)))
        for k, n in xt.most_common(12):
            s.write("{}\t{}\n".format(n, " | ".join(k)))
        s.write("\n## lost genes (tx_union>=0.9, pasa_union<0.9): {}\n".format(len(set(lost))))
        for k, n in gene_class.most_common():
            s.write("{}\t{}\n".format(k, n))
        s.write("\n## transcripts over lost genes, status ({})\n".format(" | ".join(progs)))
        for k, n in tx_states.most_common(10):
            s.write("{}\t{}\n".format(n, " | ".join(k)))
    sys.stdout.write(open(a.out + ".summary.tsv").read())


if __name__ == "__main__":
    main()
