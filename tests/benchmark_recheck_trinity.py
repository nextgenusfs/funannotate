#!/usr/bin/env python3
"""Re-check Trinity benchmark runs that were marked FAILED by a harness bug.

Runs started before the fix looked for training/funannotate_train.trinity-GG.fasta,
which --stop_after_trinity never writes, so valid runs were marked FAILED. This
script redoes the checks from the files that the run copied back:

  - funannotate must have exited 0 (exit_status=0 in the .metrics file);
  - training/trinity.fasta must exist and be non-empty;
  - the engine proof must hold: the rust leg's logs contain calls to
    rust_bio_utils binaries, the perl leg's logs contain none.

If all hold, it sets status=OK and records status_repaired=recheck so the
change is visible. Otherwise it sets status=INVALID or leaves FAILED.

Usage:  tests/benchmark_recheck_trinity.py /path/to/BENCH_ROOT/trinity
"""
import glob
import os
import re
import sys


def read_kv(path):
    kv = {}
    for line in open(path):
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*)$", line.rstrip("\n"))
        if m:
            kv[m.group(1)] = m.group(2)
    return kv


def count_rust_calls(out):
    n = 0
    needle = "rust_bio_utils/target/release"
    for sub in ("logfiles", "training"):
        for dirpath, _, files in os.walk(os.path.join(out, sub)):
            for name in files:
                if not name.endswith(".log"):
                    continue
                try:
                    with open(os.path.join(dirpath, name), errors="replace") as fh:
                        n += sum(1 for line in fh if needle in line)
                except OSError:
                    pass
    return n


def recheck(stage_dir):
    results = os.path.join(stage_dir, "results")
    for mpath in sorted(glob.glob(os.path.join(results, "*_rep*.metrics"))):
        name = os.path.basename(mpath)[: -len(".metrics")]
        kv = read_kv(mpath)
        if kv.get("status") != "FAILED" or kv.get("exit_status") != "0":
            continue
        if kv.get("transcripts") is not None:
            continue  # the fixed script already measured it
        engine = name.split("_rep")[0]
        out = os.path.join(stage_dir, "work_" + name, "output")
        fasta = os.path.join(out, "training", "trinity.fasta")
        if not os.path.isfile(fasta) or os.path.getsize(fasta) == 0:
            print("%s: no trinity.fasta; left FAILED" % name)
            continue
        ntx = 0
        nbases = 0
        for line in open(fasta):
            if line.startswith(">"):
                ntx += 1
            else:
                nbases += len(line.strip())
        calls = count_rust_calls(out)
        ok = (engine == "rust" and calls > 0) or (engine == "perl" and calls == 0)
        status = "OK" if ok else "INVALID"
        lines = [l for l in open(mpath) if not l.startswith("status=")]
        lines += ["transcripts=%d\n" % ntx, "total_bases=%d\n" % nbases,
                  "rust_calls_in_logs=%d\n" % calls, "status=%s\n" % status,
                  "status_repaired=recheck\n"]
        with open(mpath, "w") as fh:
            fh.writelines(lines)
        print("%s: transcripts=%d rust_calls=%d -> %s" % (name, ntx, calls, status))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    recheck(sys.argv[1])
