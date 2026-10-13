#!/usr/bin/env python3
"""Classify spliced-alignment introns that do not exactly match a reference intron.

Input: PASA alignment.validations.output[.gz] (introns from the 'alignment'
column of multi-segment rows for one prog), a RefSeq GFF3 (introns from exon
features) and the genome FASTA (for splice motifs).
Classes for non-matching introns:
  shifted_same_len   both ends shifted by the same 1-10 bp (equivalent placement in a repeat)
  shifted_1_10bp     a reference intron within 10 bp at each end
  noncanonical       donor/acceptor not GT-AG/GC-AG/AT-AC (either orientation)
  in_ref_gene        inside a reference mRNA span, different intron (alt isoform / retained splicing / error)
  outside_ref_gene   no reference mRNA overlaps
Usage: intron_discordance.py VALIDATIONS PROG REF_GFF3 GENOME_FA[.gz]
"""
import sys, gzip, re, bisect, collections
val, prog, refgff, genome = sys.argv[1:5]
op = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)
ex = collections.defaultdict(list); span = collections.defaultdict(list)
for l in op(refgff):
    c = l.split("\t")
    if len(c) > 8 and c[2] == "exon":
        m = re.search(r"Parent=([^;]+)", c[8])
        if m: ex[m.group(1)].append((c[0], int(c[3]), int(c[4])))
ref = collections.defaultdict(set)
for k, v in ex.items():
    v.sort(key=lambda x: x[1]); span[v[0][0]].append((v[0][1], v[-1][2]))
    for a, b in zip(v, v[1:]): ref[a[0]].add((a[2] + 1, b[1] - 1))
refl = {ch: sorted(s) for ch, s in ref.items()}
for ch in span: span[ch].sort()
seq, name, buf = {}, None, []
for l in op(genome):
    if l.startswith(">"):
        if name: seq[name] = "".join(buf).upper()
        name, buf = l[1:].split()[0], []
    else: buf.append(l.strip())
seq[name] = "".join(buf).upper()
CAN = {("GT","AG"),("GC","AG"),("AT","AC"),("CT","AC"),("CT","GC"),("GT","AT")}
cls = collections.Counter(); by_valid = collections.defaultdict(collections.Counter); n = exact = 0
for l in op(val):
    c = l.rstrip("\n").split("\t")
    if len(c) != 15 or c[0] != prog or int(c[5]) < 2: continue
    ok = c[8] == "OK"; ch = c[4]
    segs = sorted((min(a, b), max(a, b)) for a, b in ((int(x), int(y)) for x, y in re.findall(r"(\d+)\(\d+\)-(\d+)\(\d+\)", c[13])))
    for a, b in zip(segs, segs[1:]):
        s, e = a[1] + 1, b[0] - 1; n += 1
        if (s, e) in ref.get(ch, ()): exact += 1; by_valid[ok]["exact"] += 1; continue
        L = refl.get(ch, []); i = bisect.bisect_left(L, (s - 10, 0)); k = None
        while i < len(L) and L[i][0] <= s + 10:
            if abs(L[i][1] - e) <= 10: k = L[i]; break
            i += 1
        if k:
            k2 = "shifted_same_len" if (k[0] - s) == (k[1] - e) else "shifted_1_10bp"
        elif (seq[ch][s-1:s+1], seq[ch][e-2:e]) not in CAN: k2 = "noncanonical"
        elif any(x <= s and e <= y for x, y in span.get(ch, [])[:bisect.bisect_right(span.get(ch, []), (s, 10**12))]): k2 = "in_ref_gene"
        else: k2 = "outside_ref_gene"
        cls[k2] += 1; by_valid[ok][k2] += 1
print(f"introns={n} exact={exact} ({100*exact/n:.1f}%) non-matching={n-exact}")
for k, v in cls.most_common(): print(f"  {k:18s} {v:6d}  {100*v/(n-exact):5.1f}% of non-matching")
for ok in (True, False):
    t = sum(by_valid[ok].values()); print(f"  [{'valid' if ok else 'failed'} alignments] introns={t}: " + ", ".join(f"{k}={v}" for k, v in by_valid[ok].most_common()))
