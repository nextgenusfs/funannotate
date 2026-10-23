"""Chunked, parallel tRNAscan-SE for genomes with very many small contigs.

tRNAscan-SE spends most of its time per contig, not per base: for each contig it
writes a temp file and starts one cmsearch process per covariance model. The
``--thread`` option only sets ``cmsearch --cpu`` and does nothing for short contigs.
Splitting the contigs over several tRNAscan-SE processes avoids this.

Contigs are never split, so the tRNA calls are the same as a whole-genome run.
Chunks are balanced by contig count (round robin), because the cost is per contig.
"""
import os
import shutil
import subprocess
import tempfile

# use the chunked path from this many contigs
MIN_CONTIGS = 2000
# below this mean contig length, --thread does not help and chunking does
MAX_MEAN_CONTIG = 100000

HEADER_STARTS = ("Sequence", "Name", "--------")


def fasta_stats(fasta):
    """Return (number of contigs, total bases) of a FASTA file."""
    n = 0
    total = 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
            else:
                total += len(line.strip())
    return n, total


def should_chunk(fasta, cpus, min_contigs=MIN_CONTIGS, max_mean=MAX_MEAN_CONTIG):
    """True when the genome has many short contigs and more than one CPU is available."""
    if cpus < 2:
        return False
    n, total = fasta_stats(fasta)
    return n >= min_contigs and total / max(n, 1) <= max_mean


def split_fasta(fasta, k, outdir):
    """Write k chunk FASTA files with whole contigs dealt round robin.

    Returns the list of chunk paths. Never returns an empty chunk, so k is reduced
    to the number of contigs when needed (tRNAscan-SE exits 25 on an empty FASTA).
    """
    n, _ = fasta_stats(fasta)
    k = max(1, min(k, n))
    paths = [os.path.join(outdir, "chunk{:03d}.fa".format(i)) for i in range(k)]
    handles = [open(p, "w") for p in paths]
    try:
        idx = -1
        out = None
        with open(fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    idx += 1
                    out = handles[idx % k]
                if out is not None:
                    out.write(line)
    finally:
        for h in handles:
            h.close()
    return paths


def merge_outputs(outputs, merged):
    """Merge tRNAscan-SE tabular outputs: one header, then all rows in chunk order.

    A chunk without tRNAs gives a 0-byte file and is skipped. If no chunk has
    rows, a 0-byte file is written, which is what a whole-genome run gives.
    """
    header = None
    nrows = 0
    with open(merged, "w") as out:
        for path in outputs:
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                continue
            head = []
            with open(path) as fh:
                for line in fh:
                    if line.startswith(HEADER_STARTS):
                        head.append(line)
                        continue
                    if header is None:
                        header = head
                        out.writelines(header)
                    out.write(line)
                    nrows += 1
    return nrows


def run_chunked(fasta, output, cpus, workdir, log=None, exe="tRNAscan-SE"):
    """Run tRNAscan-SE in parallel on chunks of fasta and merge into output.

    workdir is the parent for the temporary chunk directory. Each process gets its
    own TMPDIR, because tRNAscan-SE names its temp files from the process id only.
    Returns True when every chunk exited 0.
    """
    tmp = tempfile.mkdtemp(prefix="trnascan_chunks_", dir=workdir)
    ok = True
    try:
        chunks = split_fasta(fasta, cpus, tmp)
        procs = []
        for i, chunk in enumerate(chunks):
            tdir = os.path.join(tmp, "tmp{:03d}".format(i))
            os.makedirs(tdir)
            out = os.path.join(tmp, "chunk{:03d}.out".format(i))
            # -Q: never prompt on overwrite (would hang with no stdin)
            cmd = [exe, "-q", "-Q", "-o", out, "--thread", "1", chunk]
            env = dict(os.environ, TMPDIR=tdir)
            if log:
                log.debug(" ".join(cmd))
            # stdout and stderr go to files: pipes could fill while we wait on another chunk
            fo = open(out + ".stdout", "w")
            fe = open(out + ".stderr", "w")
            p = subprocess.Popen(cmd, stdout=fo, stderr=fe, env=env)
            procs.append((cmd, out, p, fo, fe))
        outputs = []
        for cmd, out, p, fo, fe in procs:
            p.wait()
            fo.close()
            fe.close()
            outputs.append(out)
            if p.returncode != 0:
                ok = False
                if log:
                    log.error("CMD ERROR: {}".format(" ".join(cmd)))
            if log:
                for ext in (".stdout", ".stderr"):
                    with open(out + ext) as fh:
                        text = fh.read()
                    if text:
                        log.debug(text)
        if ok:
            merge_outputs(outputs, output)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return ok
