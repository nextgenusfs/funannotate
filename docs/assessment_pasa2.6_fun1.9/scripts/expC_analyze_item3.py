#!/usr/bin/python3.12
"""Item 3 (D121) analysis: does an own-genome Trinity/PASA set rescue PASA training?

Arms per genome (3 repeats each, masked scoring over all arms, scores_item3.tsv):
  pasa / busco         : shared-Trinity PASA set (production path); training from PASA / BUSCO
  pasa_own / busco_own : own-genome full train PASA set; training from PASA / BUSCO
Writes item3_summary.tsv and item3_analysis.txt next to this script.
"""
import csv
import gzip
import os
import statistics as st
import sys
from collections import defaultdict

C = os.path.dirname(os.path.abspath(__file__))
R = os.path.dirname(os.path.dirname(C))  # Fungi_BFD_runs
TRIAGE = os.path.join(C, "..", "wave1_triage", "wave1_triage.tsv")
T975_DF4 = 2.776  # t quantile, 4 df (two arms x 3 repeats)
ARMS = ["pasa", "busco", "pasa_own", "busco_own"]
LEVELS = ["locus", "exon", "intron_chain"]


def f1(sn, pr):
    return 0.0 if sn + pr == 0 else 2 * sn * pr / (sn + pr)


def diff_ci(a, b):
    d = st.mean(a) - st.mean(b)
    se = (st.variance(a) / len(a) + st.variance(b) / len(b)) ** 0.5
    return d, d - T975_DF4 * se, d + T975_DF4 * se


def decisions(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["stage"], row["decision"])] = row["value"]
    return out


def median_len(fasta):
    if not os.path.exists(fasta):
        return None
    op = gzip.open if fasta.endswith(".gz") else open
    lens, n = [], 0
    with op(fasta, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lens.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n:
        lens.append(n)
    return (len(lens), st.median(lens)) if lens else None


def main():
    scores = defaultdict(lambda: defaultdict(list))  # genome -> (arm, level) -> [f1]
    with open(os.path.join(C, "scores_item3.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            arm = row["variant"].rsplit("_r", 1)[0]
            for lv in LEVELS:
                scores[row["genome"]][(arm, lv)].append(
                    f1(float(row[f"{lv}_sn"]), float(row[f"{lv}_pr"]))
                )
            scores[row["genome"]]["scored"] = (row["ref_mrnas_scored"], row["ref_mrnas_total"])

    shared_med = {}
    with open(TRIAGE) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["shared_trinity_median_len"]:
                shared_med[row["out"]] = (row["shared_trinity_n"], row["shared_trinity_median_len"])

    rows = []
    skipped = [g for g in scores if not all(scores[g][(a, "locus")] for a in ARMS)]
    for g in scores:
        if g in skipped:
            continue
        rec = {"genome": g}
        rec["ref_scored"], rec["ref_total"] = scores[g]["scored"]
        for arm, tag in (("pasa", "shared"), ("pasa_own", "own")):
            d = decisions(os.path.join(C, g, f"{arm}_r1", "out", "logfiles", "training_decisions.tsv"))
            rec[f"{tag}_complete"] = d.get(("pasa_gate", "complete-ORF models in the PASA GFF3"), "")
            rec[f"{tag}_final"] = d.get(("select_final", "PASA training models"), "")
        sm = shared_med.get(g)
        rec["shared_trinity_n"], rec["shared_trinity_median"] = sm if sm else ("", "")
        om = median_len(os.path.join(C, g, "train_own", "out", "training", "funannotate_train.trinity-GG.fasta"))
        rec["own_trinity_n"], rec["own_trinity_median"] = (om[0], int(om[1])) if om else ("", "")
        for lv in LEVELS:
            for arm in ARMS:
                v = scores[g][(arm, lv)]
                rec[f"{arm}_{lv}"] = round(st.mean(v), 2)
                rec[f"{arm}_{lv}_sd"] = round(st.stdev(v), 2)
            for a, b in (("pasa", "busco"), ("pasa_own", "busco_own"), ("pasa_own", "pasa"),
                         ("busco_own", "busco")):
                d, lo, hi = diff_ci(scores[g][(a, lv)], scores[g][(b, lv)])
                rec[f"{a}-{b}_{lv}"] = f"{d:+.2f} [{lo:+.2f}, {hi:+.2f}]"
            best_shared = max(st.mean(scores[g][(a, lv)]) for a in ("pasa", "busco"))
            best_own = max(st.mean(scores[g][(a, lv)]) for a in ("pasa_own", "busco_own"))
            rec[f"best_own-best_shared_{lv}"] = round(best_own - best_shared, 2)
        rows.append(rec)
    rows.sort(key=lambda r: int(r["shared_complete"] or 0))

    with open(os.path.join(C, "item3_summary.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    out = []
    p = out.append
    p("Item 3 (D121): masked F1, mean of 3 repeats; 95% CI from repeat SDs (t, 4 df), run noise only.")
    for lv in LEVELS:
        p(f"\n== {lv} F1 ==")
        p("genome\tshared_med\town_med\tshared_cplt\town_cplt\tpasa\tbusco\tpasa_own\tbusco_own\t"
          "pasa-busco\tpasa_own-busco_own\tpasa_own-pasa\tbusco_own-busco\tbest_own-best_shared")
        for r in rows:
            p("\t".join(str(x) for x in [
                r["genome"], r["shared_trinity_median"], r["own_trinity_median"],
                r["shared_complete"], r["own_complete"],
                r[f"pasa_{lv}"], r[f"busco_{lv}"], r[f"pasa_own_{lv}"], r[f"busco_own_{lv}"],
                r[f"pasa-busco_{lv}"], r[f"pasa_own-busco_own_{lv}"],
                r[f"pasa_own-pasa_{lv}"], r[f"busco_own-busco_{lv}"],
                r[f"best_own-best_shared_{lv}"]]))
    p("\ngenomes without all four arms (not in item 3, shared arms only): " + ", ".join(skipped))
    p("\nmax repeat SD per arm (locus): " + ", ".join(
        f"{a} {max(r[f'{a}_locus_sd'] for r in rows):.2f}" for a in ARMS))
    txt = "\n".join(out)
    with open(os.path.join(C, "item3_analysis.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    sys.exit(main())
