"""Training-set diversity metrics against RefSeq, for benchmark.py.

Owner: the PASA-side review session (CODE_REVIEW_20260925). benchmark.py
calls diversity(genes, ref), where genes is a funannotate gff2dict dict and
ref is the tuple returned by benchmark.load_refseq(). Only the first
transcript of each gene (v["CDS"][0]) is scored, as in benchmark.score().

Columns returned:
  refseq_genes_hit      distinct RefSeq coding genes overlapped by a model
                        on the same strand and in the same reading frame
  recall_hit_pct        refseq_genes_hit / RefSeq coding genes
  redundant_models      models whose in-frame RefSeq gene was already hit
                        by another model (duplicates of one locus)
  gt3_introns_pct       models with more than 3 introns (> 4 CDS segments)
  ref_gt3_introns_pct   the same share for RefSeq (one chain per gene)
  exon_count_ks         Kolmogorov-Smirnov distance between CDS-segment
                        counts of models and RefSeq (0 = same distribution)
  median_cds_len        median CDS length of models (nt)
  ref_median_cds_len    the same for RefSeq
  cds_len_ks            KS distance between CDS-length distributions
The RefSeq side uses the longest CDS chain per gene.
"""
import bisect
import collections

_CACHE = {}


def _frames(strand, ex):
    """ex sorted ascending -> list of (s, e, offset_of_s_in_CDS) for '+',
    or (s, e, offset_of_e_in_CDS) for '-'."""
    out = []
    p = 0
    order = ex if strand == "+" else list(reversed(ex))
    for s, e in order:
        out.append((s, e, p))
        p += e - s + 1
    return out


def _frame_at(strand, fr, g):
    for s, e, off in fr:
        if s <= g <= e:
            return (off + (g - s if strand == "+" else e - g)) % 3
    return None


def _ks(a, b):
    if not a or not b:
        return float("nan")
    a, b = sorted(a), sorted(b)
    i = j = 0
    d = 0.0
    while i < len(a) and j < len(b):
        x = min(a[i], b[j])
        while i < len(a) and a[i] <= x:
            i += 1
        while j < len(b) and b[j] <= x:
            j += 1
        d = max(d, abs(i / len(a) - j / len(b)))
    return round(d, 3)


def _prepare(ref):
    key = id(ref)
    if key in _CACHE:
        return _CACHE[key]
    chains = ref[0]
    per_gene = {}
    for (ctg, strand, ex), gid in chains.items():
        ln = sum(e - s + 1 for s, e in ex)
        if gid not in per_gene or ln > per_gene[gid][0]:
            per_gene[gid] = (ln, len(ex))
    idx = collections.defaultdict(list)
    for (ctg, strand, ex), gid in chains.items():
        idx[(ctg, strand)].append((ex[0][0], ex[-1][1], gid, _frames(strand, ex)))
    starts = {}
    maxspan = 0
    for k in idx:
        idx[k].sort(key=lambda x: x[0])
        starts[k] = [x[0] for x in idx[k]]
        maxspan = max([maxspan] + [x[1] - x[0] for x in idx[k]])
    prep = (per_gene, idx, starts, maxspan)
    _CACHE[key] = prep
    return prep


def diversity(genes, ref):
    per_gene, idx, starts, maxspan = _prepare(ref)
    coding = ref[3]
    counts, lens = [], []
    hit = set()
    redundant = 0
    for v in genes.values():
        ex = sorted(v["CDS"][0])
        if not ex:
            continue
        ctg, strand = v["contig"], v["strand"]
        counts.append(len(ex))
        lens.append(sum(e - s + 1 for s, e in ex))
        mfr = _frames(strand, ex)
        lo, hi = ex[0][0], ex[-1][1]
        L = idx.get((ctg, strand), [])
        i = bisect.bisect_left(starts.get((ctg, strand), []), lo - maxspan)
        gene_hit = None
        while i < len(L) and L[i][0] <= hi and gene_hit is None:
            rlo, rhi, gid, rfr = L[i]
            i += 1
            if rhi < lo:
                continue
            for s, e, _ in mfr:
                for rs, re_, _ in rfr:
                    g = max(s, rs)
                    if g <= min(e, re_) and _frame_at(strand, mfr, g) == _frame_at(strand, rfr, g):
                        gene_hit = gid
                        break
                if gene_hit:
                    break
        if gene_hit is not None:
            if gene_hit in hit:
                redundant += 1
            hit.add(gene_hit)
    ref_counts = [c for _, c in per_gene.values()]
    ref_lens = [l for l, _ in per_gene.values()]
    pct = lambda a, b: round(100.0 * a / b, 1) if b else 0.0
    med = lambda x: sorted(x)[len(x) // 2] if x else 0
    return {
        "refseq_genes_hit": len(hit),
        "recall_hit_pct": pct(len(hit), len(coding)),
        "redundant_models": redundant,
        "gt3_introns_pct": pct(sum(1 for c in counts if c > 4), len(counts)),
        "ref_gt3_introns_pct": pct(sum(1 for c in ref_counts if c > 4), len(ref_counts)),
        "exon_count_ks": _ks(counts, ref_counts),
        "median_cds_len": med(lens),
        "ref_median_cds_len": med(ref_lens),
        "cds_len_ks": _ks(lens, ref_lens),
    }
