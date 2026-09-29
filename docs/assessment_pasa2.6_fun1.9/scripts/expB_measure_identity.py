"""Read identity + map rate for experiment B genomes, using funannotate's train gate code.

Input: experiment_B/genomes.tsv. Reads: <training_dir>/normalize/left.norm.fq.gz
(the file funannotate train's gate samples). Genome: experiment_B/<name>/genome.fa.
Output: one TSV row per genome.
"""
import argparse, csv, os, sys, time
from funannotate import library as lib

ap = argparse.ArgumentParser()
ap.add_argument("--genomes", required=True)
ap.add_argument("--expdir", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--n_reads", type=int, default=200000)
ap.add_argument("--cpus", type=int, default=16)
ap.add_argument("--tmpdir", required=True)
a = ap.parse_args()

cols = ["name", "group", "stratum", "reads", "sampled", "mapped", "map_rate_pct",
        "median_identity_pct", "p10_identity_pct", "identity_reads", "seconds", "status"]
with open(a.genomes) as f, open(a.out, "w") as out:
    out.write("\t".join(cols) + "\n")
    for g in csv.DictReader(f, delimiter="\t"):
        row = dict.fromkeys(cols, "")
        row.update(name=g["name"], group=g["group"], stratum=g["stratum"])
        reads = os.path.join(g["training_dir"], "normalize", "left.norm.fq.gz")
        genome = os.path.join(a.expdir, g["name"], "genome.fa")
        row["reads"] = reads
        if not os.path.isfile(reads):
            row["status"] = "no_reads"
        else:
            t = time.time()
            try:
                s, m, ids = lib.sample_read_map_rate(reads, genome, a.n_reads, a.cpus, a.tmpdir)
                med, p10 = lib.summarize_identity(ids)
                row.update(sampled=s, mapped=m, map_rate_pct=round(100.0 * m / s, 2) if s else "",
                           median_identity_pct="" if med is None else med,
                           p10_identity_pct="" if p10 is None else p10,
                           identity_reads=len(ids), status="ok")
            except Exception as e:
                row["status"] = "error: " + str(e).replace("\t", " ")
            row["seconds"] = int(time.time() - t)
        out.write("\t".join(str(row[c]) for c in cols) + "\n")
        out.flush()
        print(row["name"], row["status"], row["map_rate_pct"], row["median_identity_pct"], flush=True)
