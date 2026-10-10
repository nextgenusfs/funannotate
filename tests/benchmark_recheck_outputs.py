#!/usr/bin/env python3
"""Add structure-only checksums to benchmark results made before the scripts had them.

PASA and EVM number their IDs in a run-dependent order, so two runs of the same input give
different files and different md5sums, with the same gene models. The structure checksum
ignores IDs and uses only sequence, type, start, end and strand of every feature.

For each results/<engine>_rep<N>.metrics in the pasa and evm stages, read the GFF3 that the run
copied back and append pasa_gff3_struct_md5= or gff3_struct_md5= when it is missing.

Usage: tests/benchmark_recheck_outputs.py /path/to/BENCH_ROOT
"""
import glob
import hashlib
import os
import sys


def structure_md5(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7:
                rows.append("\t".join([f[0], f[2], f[3], f[4], f[6]]))
    rows.sort()
    return hashlib.md5(("\n".join(rows) + "\n").encode()).hexdigest()


def main(root):
    jobs = (
        ("pasa", "work_pasa_%s", "output/training/funannotate_train.pasa.gff3", "pasa_gff3_struct_md5"),
        ("evm", "predict_work_%s", "predict_results/*.gff3", "gff3_struct_md5"),
    )
    for stage, work_fmt, gff_glob, key in jobs:
        for mpath in sorted(glob.glob(os.path.join(root, stage, "results", "*_rep*.metrics"))):
            name = os.path.basename(mpath)[: -len(".metrics")]
            text = open(mpath).read()
            if key + "=" in text or "status=OK" not in text:
                continue
            hits = sorted(glob.glob(os.path.join(root, stage, work_fmt % name, gff_glob)))
            if not hits:
                print("%s/%s: no GFF3 found" % (stage, name))
                continue
            digest = structure_md5(hits[0])
            with open(mpath, "a") as fh:
                fh.write("%s=%s\n" % (key, digest))
            print("%s/%s: %s" % (stage, name, digest))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
