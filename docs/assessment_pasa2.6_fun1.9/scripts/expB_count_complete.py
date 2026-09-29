"""Whole-genome and train-chromosome complete-ORF PASA model counts per experiment B genome.

Uses funannotate's gate function lib.count_complete_orf_models (the PASA gate's own count).
Also reads the number of final PASA training models from pasa.A training_decisions.tsv.
"""
import csv, os, sys
from funannotate import library as lib

E = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B"
out = open(os.path.join(E, "read_identity", "complete_models.tsv"), "w")
out.write("name\tcomplete_genome\ttotal_genome\tcomplete_train\ttotal_train\tfinal_training_models\n")
for g in csv.DictReader(open(os.path.join(E, "genomes.tsv")), delimiter="\t"):
    n = g["name"]; gd = os.path.join(E, n)
    if not os.path.isfile(os.path.join(gd, "pasa.genome.gff3")):
        continue
    w = lib.count_complete_orf_models(os.path.join(gd, "pasa.genome.gff3"), os.path.join(gd, "genome.fa"))
    t = lib.count_complete_orf_models(os.path.join(gd, "pasa.genome_train.gff3"), os.path.join(gd, "genome_train.fa"))
    final = ""
    dec = os.path.join(gd, "pasa.A", "out", "logfiles", "training_decisions.tsv")
    if os.path.isfile(dec):
        for line in open(dec):
            c = line.rstrip("\n").split("\t")
            if len(c) > 3 and c[1] == "select_final":
                final = c[3]
    out.write("\t".join(map(str, [n, w["complete"], w["total"], t["complete"], t["total"], final])) + "\n")
    out.flush(); print(n, w["complete"], t["complete"], final, flush=True)
