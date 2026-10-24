import glob, os, re, sys
root = sys.argv[1]
print("\t".join(["genome","arm","nuc_sn","nuc_sp","exon_sn","exon_sp","gene_sn","gene_sp"]))
def num(x):
    try: return float(x)
    except ValueError: return float("nan")
for path in sorted(glob.glob(os.path.join(root, "*", "*.A", "out", "predict_misc", "augustus.initial.training.txt"))):
    parts = path[len(root):].strip("/").split("/")
    genome, arm = parts[0], parts[1]
    vals = {}
    for line in open(path, errors="replace"):
        c = [x.strip() for x in line.rstrip("\n").split("|")]
        if line.startswith("nucleotide level") and len(c) > 2: vals["nuc"] = (num(c[1]), num(c[2]))
        elif line.startswith("exon level") and len(c) > 7: vals["exon"] = (num(c[6]), num(c[7]))
        elif line.startswith("gene level") and len(c) > 7: vals["gene"] = (num(c[6]), num(c[7]))
    if len(vals) == 3:
        print("\t".join([genome, arm] + [str(v) for k in ("nuc","exon","gene") for v in vals[k]]))
    else:
        sys.stderr.write("unparsed: %s\n" % path)
