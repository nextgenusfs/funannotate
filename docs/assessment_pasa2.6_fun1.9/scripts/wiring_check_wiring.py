#!/usr/bin/env python3
"""Gate wiring test checks (D118). Prints one PASS/FAIL/INFO line per check; writes checks.tsv.

Run after all wiring jobs finish: /usr/bin/python3.12 check_wiring.py
"""
import csv, glob, hashlib, os

BD = os.path.dirname(os.path.abspath(__file__))
EXPB = os.path.join(BD, "..", "experiment_B", "read_identity", "read_identity.tsv")
rows = []


def rec(status, genome, check, detail):
    rows.append((status, genome, check, detail)); print("%-4s %-34s %-44s %s" % (status, genome, check, detail))


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest() if os.path.isfile(p) else None


def decisions(run, command="predict"):
    p = os.path.join(BD, run, "out", "logfiles", "training_decisions.tsv")
    if not os.path.isfile(p):
        return None
    out = []
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["command"] == command:
            out.append((r["stage"], r["decision"], r["value"], r["threshold"], r["outcome"]))
    return out


def gate_report(train):
    p = os.path.join(BD, train, "out", "logfiles", "train_rnaseq_gate.tsv")
    return next(csv.DictReader(open(p), delimiter="\t")) if os.path.isfile(p) else None


def status(run):
    p = os.path.join(BD, run, "run_status.tsv")
    return dict(l.rstrip("\n").split("\t") for l in open(p)) if os.path.isfile(p) else {}


def gene_gff(run):
    g = [f for f in glob.glob(os.path.join(BD, run, "out", "predict_results", "*.gff3")) if ".tbl" not in f]
    return g[0] if g else None


def strip_gff(p):
    return hashlib.md5("".join(l for l in open(p) if not l.startswith("#")).encode()).hexdigest()


expb = {r["name"]: r for r in csv.DictReader(open(EXPB), delimiter="\t")}
EP = os.path.join(BD, "expect.tsv")
EXPECT = dict(l.rstrip("\n").split("\t") for l in open(EP)) if os.path.isfile(EP) else {}
NEW_STAGES = {"rnaseq_identity", "rnaseq_identity_gate"}

for g in sorted(d for d in os.listdir(BD) if os.path.isdir(os.path.join(BD, d, "train_new"))):
    st = status(g + "/train_new")
    rep = gate_report(g + "/train_new")
    if g.startswith("Colletotrichum_siamense"):
        ok = st.get("exit_status") == "3" and rep is not None and rep.get("passed") == "False"
        rec("PASS" if ok else "FAIL", g, "map-rate gate stops train (exit 3)", "exit %s; report %s" % (st.get("exit_status"), rep and {k: rep[k] for k in ("map_rate_pct", "median_identity_pct")}))
        continue
    if st.get("exit_status") != "0":
        rec("FAIL", g, "train_new finished", "exit %s" % st.get("exit_status")); continue
    has = rep is not None and rep.get("median_identity_pct") not in (None, "")
    rec("PASS" if has else "FAIL", g, "train report has read identity", rep and "median %s%%, p10 %s%%, map %s%%" % (rep.get("median_identity_pct"), rep.get("p10_identity_pct"), rep.get("map_rate_pct")))
    copy = os.path.join(BD, g, "train_new", "out", "training", "funannotate_train.rnaseq_gate.tsv")
    rec("PASS" if os.path.isfile(copy) else "FAIL", g, "training/ copy of report written", copy if os.path.isfile(copy) else "missing")
    if has and g in expb and expb[g]["median_identity_pct"]:
        d = abs(float(rep["median_identity_pct"]) - float(expb[g]["median_identity_pct"]))
        rec("PASS" if d <= 0.05 else "INFO", g, "identity matches experiment B (<= 0.05 pt)", "train %s vs expB %s" % (rep["median_identity_pct"], expb[g]["median_identity_pct"]))
    if os.path.isdir(os.path.join(BD, g, "train_old")):
        a = md5(os.path.join(BD, g, "train_old", "out", "training", "funannotate_train.pasa.gff3"))
        b = md5(os.path.join(BD, g, "train_new", "out", "training", "funannotate_train.pasa.gff3"))
        rec("PASS" if a and a == b else "INFO", g, "PASA GFF3 identical, train old vs new", "identical" if a == b else "differ (check determinism)")
        do = [x for x in decisions(g + "/train_old", "train") or [] if x[0] not in NEW_STAGES]
        dn = [x for x in decisions(g + "/train_new", "train") or [] if x[0] not in NEW_STAGES]
        rec("PASS" if do == dn else "FAIL", g, "train decisions identical (minus new rows)", "%d rows" % len(dn) if do == dn else "old %s / new %s" % (do, dn))
    if os.path.isdir(os.path.join(BD, g, "pred_old")) and os.path.isdir(os.path.join(BD, g, "pred_new")):
        do = [x for x in decisions(g + "/pred_old") or [] if x[0] not in NEW_STAGES]
        dn = decisions(g + "/pred_new") or []
        ig = [x for x in dn if x[0] == "rnaseq_identity_gate"]
        dn = [x for x in dn if x[0] not in NEW_STAGES]
        rec("PASS" if do and do == dn else "FAIL", g, "predict decisions identical (minus new rows)", "%d rows" % len(dn) if do == dn else "differ")
        rec("PASS" if ig and "disabled" in ig[0][3] + ig[0][2] + ig[0][4] or (ig and ig[0][4] == "PASA training allowed") else "FAIL", g, "identity gate off by default", str(ig))
        for f in ("predict_misc/final_training_models.gff3",):
            a, b = md5(os.path.join(BD, g, "pred_old", "out", f)), md5(os.path.join(BD, g, "pred_new", "out", f))
            rec("PASS" if a and a == b else "FAIL", g, "identical " + os.path.basename(f), "identical" if a == b else "differ")
        a, b = md5(os.path.join(BD, g, "pred_old", "out", "logfiles", "predict_training_gate.tsv")), md5(os.path.join(BD, g, "pred_new", "out", "logfiles", "predict_training_gate.tsv"))
        rec("PASS" if a == b else "FAIL", g, "identical predict_training_gate.tsv", "identical" if a == b else "differ")
        ga, gb = gene_gff(g + "/pred_old"), gene_gff(g + "/pred_new")
        if ga and gb:
            rec("PASS" if strip_gff(ga) == strip_gff(gb) else "INFO", g, "gene models identical (bonus)", "identical" if strip_gff(ga) == strip_gff(gb) else "differ; compare F1")
    for run in sorted(glob.glob(os.path.join(BD, g, "gate*"))) + sorted(glob.glob(os.path.join(BD, g, "decide*"))) + sorted(glob.glob(os.path.join(BD, g, "pred_default"))):
        r = os.path.relpath(run, BD); d = decisions(r) or []
        info = dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(run, "run_info.tsv")))
        exp = EXPECT.get(r, "?")
        if os.path.basename(run) == "pred_default":  # production default: BUSCO iff complete < 500
            n = [x for x in d if x[0] == "pasa_gate"]
            exp = ("busco" if int(n[0][2]) < 500 else "pasa") if n and n[0][2].isdigit() else "?"
        switched = any(x[0] == "training_mode_switch" for x in d)
        why = [x for x in d if x[0] in ("rnaseq_identity_gate", "pasa_gate", "training_mode_switch")]
        got = "busco" if switched else "pasa"
        rec("PASS" if got == exp else "FAIL", g, "%s: expect %s training" % (os.path.basename(run), exp), "got %s; %s; %s" % (got, info.get("extra", ""), why))
        if switched and os.path.isfile(os.path.join(run, "predict.capture.log")):
            log = open(os.path.join(run, "predict.capture.log")).read()
            kept = "pasa" in log.lower() and ("--pasa_gff" in log or "PASA:" in log)
            rec("INFO", g, "%s: PASA evidence still listed in log" % os.path.basename(run), "yes" if kept else "check EVM weights in log")

with open(os.path.join(BD, "checks.tsv"), "w") as o:
    o.write("status\tgenome\tcheck\tdetail\n")
    for r in rows:
        o.write("\t".join(str(x) for x in r) + "\n")
print("\n%d checks: %d PASS, %d FAIL, %d INFO" % (len(rows), sum(r[0] == "PASS" for r in rows), sum(r[0] == "FAIL" for r in rows), sum(r[0] == "INFO" for r in rows)))
