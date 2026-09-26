#!/usr/bin/env python3
"""R13: why are PASA assemblies missing start/stop codons?

For every PASA assembly whose best TransDecoder ORF is not complete, ask
whether any member transcript (from pasa_assemblies_described.txt) has a
complete ORF in the transcript-level TransDecoder run that PASA does with
--TRANSDECODER.

  member_has_complete  -> the transcript had the whole ORF; the assembly lost
                          it (alignment trimming / clipping / validation /
                          assembly), i.e. an alignment-side cause.
  no_member_complete   -> no member transcript is complete either; the
                          transcripts themselves are fragments.
  members_unscored     -> no member transcript has any TransDecoder ORF.

Optionally restricts to the models funannotate chose for training, and
measures soft clipping per transcript from the minimap2 BAM.

Inputs are read as plain text or .gz. Usage:
  r13_truncation_analysis.py --pasa-dir DIR [--train-gff3 funannotate_train.pasa.gff3]
                             [--bam trinity.alignments.bam] [--samtools samtools]
                             --out PREFIX
Writes PREFIX.summary.tsv and PREFIX.per_assembly.tsv.gz.
"""
import argparse
import collections
import glob
import gzip
import os
import re
import subprocess
import sys
from urllib.parse import unquote

RANK = {"complete": 3, "5prime_partial": 2, "3prime_partial": 2, "internal": 1}


def xopen(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def one(pattern, required=True):
    hits = sorted(glob.glob(pattern)) + sorted(glob.glob(pattern + ".gz"))
    if not hits:
        if required:
            sys.exit("R13: no file matches {}".format(pattern))
        return None
    return hits[0]


def best_orfs(td_gff3):
    """seqid -> (type, orf_len_aa) of its best ORF (rank, then length)."""
    best = {}
    with xopen(td_gff3) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "mRNA":
                continue
            name = unquote(c[8])
            m = re.search(r"ORF type:(\w+)", name)
            if not m:
                continue
            t = m.group(1)
            ln = re.search(r"len:(\d+)", name)
            ln = int(ln.group(1)) if ln else (int(c[4]) - int(c[3]) + 1) // 3
            key = (RANK.get(t, 0), ln)
            if c[0] not in best or key > best[c[0]][2]:
                best[c[0]] = (t, ln, key)
    return {k: (v[0], v[1]) for k, v in best.items()}


def members(described):
    out = {}
    dup_assemblies = 0
    with xopen(described) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            accs = c[3].split(",")
            if len(set(accs)) < len(accs):
                dup_assemblies += 1
            out[c[2]] = sorted(set(accs))
    return out, dup_assemblies


def chosen_assemblies(train_gff3):
    keep = set()
    with xopen(train_gff3) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) > 8 and c[2] == "mRNA":
                m = re.search(r"ID=(asmbl_\d+)", c[8])
                if m:
                    keep.add(m.group(1))
    return keep


def soft_clips(bam, samtools):
    """transcript -> fraction of its length soft-clipped (primary records)."""
    clip = {}
    p = subprocess.Popen([samtools, "view", "-F", "0x904", bam],
                         stdout=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout:
        c = line.split("\t", 6)
        ops = re.findall(r"(\d+)([MIDNSHP=X])", c[5])
        qlen = sum(int(n) for n, o in ops if o in "MIS=XH")
        s = sum(int(n) for n, o in ops if o in "SH")
        if qlen:
            clip[c[0]] = s / qlen
    if p.wait():
        sys.exit("R13: samtools view failed on {}".format(bam))
    return clip


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pasa-dir", required=True)
    ap.add_argument("--train-gff3")
    ap.add_argument("--bam")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = a.pasa_dir
    asm_td = one(os.path.join(d, "*.assemblies.fasta.transdecoder.gff3"))
    described = one(os.path.join(d, "*.pasa_assemblies_described.txt"))
    tx_td = [f for f in glob.glob(os.path.join(d, "*.transdecoder.gff3*"))
             if ".assemblies.fasta." not in f]
    if not tx_td:
        sys.exit("R13: no transcript-level *.transdecoder.gff3 in {} "
                 "(needs PASA --TRANSDECODER)".format(d))
    asm_orf = best_orfs(asm_td)
    tx_orf = {}
    for f in tx_td:
        tx_orf.update(best_orfs(f))
    mem, dup_asm = members(described)
    chosen = chosen_assemblies(a.train_gff3) if a.train_gff3 else None
    clip = soft_clips(a.bam, a.samtools) if a.bam else None

    counts = collections.Counter()
    gain = []  # member complete ORF length - assembly best ORF length
    clip_by_class = collections.defaultdict(list)
    with gzip.open(a.out + ".per_assembly.tsv.gz", "wt") as out:
        out.write("asmbl\tasmbl_orf\tasmbl_len_aa\tclass\tn_members\tbest_member\t"
                  "member_orf\tmember_len_aa\tmax_member_clip_frac\n")
        for asmbl, (t, ln) in sorted(asm_orf.items()):
            if chosen is not None and asmbl not in chosen:
                continue
            counts["assemblies"] += 1
            if t == "complete":
                counts["complete"] += 1
                continue
            accs = mem.get(asmbl, [])
            scored = [(tx_orf[x], x) for x in accs if x in tx_orf]
            complete = [s for s in scored if s[0][0] == "complete"]
            if complete:
                cls = "member_has_complete"
                (mt, ml), mx = max(complete, key=lambda s: s[0][1])
                gain.append(ml - ln)
            elif scored:
                cls = "no_member_complete"
                (mt, ml), mx = max(scored, key=lambda s: (RANK.get(s[0][0], 0), s[0][1]))
            else:
                cls = "members_unscored"
                mt, ml, mx = "NA", 0, "NA"
            counts[cls] += 1
            mc = "NA"
            if clip is not None:
                vals = [clip[x] for x in accs if x in clip]
                if vals:
                    mc = max(vals)
                    clip_by_class[cls].append(mc)
                    mc = "{:.3f}".format(mc)
            out.write("\t".join(str(v) for v in
                                (asmbl, t, ln, cls, len(accs), mx, mt, ml, mc)) + "\n")

    def med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else float("nan")

    with open(a.out + ".summary.tsv", "w") as s:
        s.write("metric\tvalue\n")
        s.write("scope\t{}\n".format("training_models" if chosen is not None else "all_assemblies"))
        for k in ("assemblies", "complete", "member_has_complete",
                  "no_member_complete", "members_unscored"):
            s.write("{}\t{}\n".format(k, counts[k]))
        inc = counts["assemblies"] - counts["complete"]
        if inc:
            s.write("frac_incomplete_with_complete_member\t{:.3f}\n".format(
                counts["member_has_complete"] / inc))
        s.write("median_orf_gain_aa_if_member_used\t{}\n".format(med(gain)))
        s.write("assemblies_listing_same_transcript_twice\t{}\n".format(dup_asm))
        if clip is not None:
            allc = list(clip.values())
            s.write("transcripts_in_bam\t{}\n".format(len(allc)))
            s.write("frac_transcripts_clip_ge_10pct\t{:.3f}\n".format(
                sum(1 for x in allc if x >= 0.10) / len(allc) if allc else float("nan")))
            for cls, v in sorted(clip_by_class.items()):
                s.write("median_max_member_clip_{}\t{:.3f}\n".format(cls, med(v)))
                s.write("frac_clip_ge_10pct_{}\t{:.3f}\n".format(
                    cls, sum(1 for x in v if x >= 0.10) / len(v)))
    sys.stdout.write(open(a.out + ".summary.tsv").read())


if __name__ == "__main__":
    main()
