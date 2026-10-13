# -*- coding: utf-8 -*-
"""Reduce redundant transcripts before PASA alignment assembly.

Both filters act only on the transcripts PASA assembles (the cleaned and raw
transcript FASTA and the imported minimap2 alignment GFF3). The alignment BAM
and GFF3 used later as EVM transcript evidence are not changed.

- contained fragments: a transcript is dropped when another transcript of the
  same Trinity gene, on the same sequence and strand, contains it: the dropped
  transcript's intron chain is a contiguous part of the other's, with identical
  intron coordinates, and its first and last exons lie inside the matching exons
  of the other. A single-exon transcript is dropped when it lies inside one exon
  of another transcript of the same gene. It adds no junction PASA does not
  already see.
- top-N isoforms: keep at most N isoforms per Trinity gene, ranked by abundance
  (TPM, then length).

Trinity gene = transcript ID without its trailing "_i<N>". Names without that
suffix form their own gene, so neither filter removes them. A transcript with
more than one alignment locus is never removed.
"""
import re
from collections import defaultdict

_ISOFORM_RE = re.compile(r"_i\d+$")


def trinity_gene(tid):
    return _ISOFORM_RE.sub("", tid)


def read_alignment_exons(gff3):
    """transcript ID -> list of (chrom, strand, sorted [(start, end), ...]) loci.

    Reads cDNA_match lines (one per exon) as written by lib.bam2gff3. Exons of
    one transcript on a different sequence or strand form a separate locus.
    """
    parts = defaultdict(lambda: defaultdict(list))
    with open(gff3) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "cDNA_match":
                continue
            m = re.search(r"ID=([^;]+)", c[8])
            if not m:
                continue
            parts[m.group(1)][(c[0], c[6])].append((int(c[3]), int(c[4])))
    out = {}
    for tid, loci in parts.items():
        out[tid] = [(ch, st, sorted(ex)) for (ch, st), ex in loci.items()]
    return out


def _introns(exons):
    return [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)]


def _contains(big, small, mode="strict"):
    """True if exon chain `small` is contained in exon chain `big`.

    strict : as in the module doc (introns identical, ends inside matching exons).
    introns: a spliced `small` only needs its intron chain to be a contiguous part
             of big's; its ends are free, so end diversity (UTR length, retained
             terminal intron parts) is given up. Single-exon transcripts as strict.
    """
    if len(small) == 1:
        s, e = small[0]
        return any(bs <= s and e <= be for bs, be in big)
    if mode == "strict" and (small[0][0] < big[0][0] or small[-1][1] > big[-1][1]):
        return False
    bi, si = _introns(big), _introns(small)
    n = len(si)
    for j in range(len(bi) - n + 1):
        if bi[j:j + n] == si:
            if mode == "introns":
                return True
            # first exon of small inside exon j of big, last inside exon j + n
            return small[0][0] >= big[j][0] and small[-1][1] <= big[j + n][1]
    return False


def contained_transcripts(alignments, mode="strict"):
    """Return the set of transcript IDs to drop as contained fragments.

    Among transcripts with identical chains, the one with the smallest ID is kept.
    mode: 'strict' or 'introns' (see _contains).
    """
    if mode not in ("strict", "introns"):
        raise ValueError("unknown contained-fragment mode: %s" % mode)
    groups = defaultdict(list)
    for tid, loci in alignments.items():
        if len(loci) != 1:
            continue
        ch, st, ex = loci[0]
        groups[(trinity_gene(tid), ch, st)].append((tid, ex))
    drop = set()
    for members in groups.values():
        if len(members) < 2:
            continue
        # sort so a transcript is only compared with ones that could contain it,
        # and ties keep one copy: strict -> larger span first; introns -> more
        # exons first (a longer chain can contain a chain with a longer span)
        span = lambda x: x[1][-1][1] - x[1][0][0]
        if mode == "strict":
            members.sort(key=lambda x: (-span(x), -len(x[1]), x[0]))
        else:
            members.sort(key=lambda x: (-len(x[1]), -span(x), x[0]))
        for i, (tid, ex) in enumerate(members):
            for otid, oex in members[:i]:
                if otid in drop:
                    continue
                if _contains(oex, ex, mode):
                    drop.add(tid)
                    break
    return drop


def top_isoforms(ids, abundance, lengths, max_per_gene):
    """Return the set of IDs to drop so each Trinity gene keeps <= max_per_gene.

    abundance: ID -> TPM (missing = 0). lengths: ID -> transcript length.
    max_per_gene <= 0 keeps everything.
    """
    if max_per_gene <= 0:
        return set()
    by_gene = defaultdict(list)
    for tid in ids:
        by_gene[trinity_gene(tid)].append(tid)
    drop = set()
    for members in by_gene.values():
        if len(members) <= max_per_gene:
            continue
        members.sort(key=lambda t: (-abundance.get(t, 0.0), -lengths.get(t, 0), t))
        drop.update(members[max_per_gene:])
    return drop


def fasta_lengths(fasta):
    lengths, name = {}, None
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def write_fasta_without(fasta, drop, output):
    keep, kept = True, 0
    with open(fasta) as fh, open(output, "w") as out:
        for line in fh:
            if line.startswith(">"):
                keep = line[1:].split()[0] not in drop
                kept += keep
            if keep:
                out.write(line)
    return kept


def write_table_without(table, drop, output):
    """Copy a table whose first whitespace field is a transcript ID (e.g. the
    seqclean .cln report), leaving out rows for IDs in drop."""
    with open(table) as fh, open(output, "w") as out:
        for line in fh:
            f = line.split(None, 1)
            if f and f[0] in drop:
                continue
            out.write(line)


def write_gff3_without(gff3, drop, output):
    with open(gff3) as fh, open(output, "w") as out:
        for line in fh:
            if not line.startswith("#"):
                m = re.search(r"ID=([^;\s]+)", line)
                if m and m.group(1) in drop:
                    continue
            out.write(line)
