#!/bin/bash

#unit tests for funannotate

#run through tests for each genome
cmd='funannotate predict -i genome1.fasta -s "Genome one" --name GN1_ -o genome1 --protein_evidence proteins.fa --transcript_evidence transcripts.fa --augustus_species botrytis_cinerea --cpus 6'
echo $cmd; eval $cmd

cmd='funannotate predict -i genome2.fasta -s "Genome two" --name GN2_ -o genome2 --protein_evidence proteins.fa --transcript_evidence transcripts.fa --augustus_species botrytis_cinerea --cpus 6'
echo $cmd; eval $cmd

cmd='funannotate predict -i genome3.fasta -s "Genome three" --name GN3_ -o genome3 --protein_evidence proteins.fa --transcript_evidence transcripts.fa --augustus_species botrytis_cinerea --cpus 6'
echo $cmd; eval $cmd

#test BUSCO2 mediated training
cmd='funannotate predict -i genome4.fasta -s "Genome four" -o genome4 --name GN4_ --protein_evidence proteins.fa --transcript_evidence transcripts.fa --busco_seed_species botrytis_cinerea --cpus 6'
echo $cmd; eval $cmd

#test maker GFF
cmd='funannotate predict -i genome1.fasta -s "Genome one" --name GN12_ -o genome1_maker --protein_evidence proteins.fa --transcript_evidence transcripts.fa --augustus_species botrytis_cinerea --cpus 6 --maker_gff maker_genome1.all.gff'
echo $cmd; eval $cmd

#now annotate each genome
cmd='funannotate annotate -i genome1 -e palmer3@wisc.edu --cpus 6 --iprscan iprscan_results'
echo $cmd; eval $cmd

cmd='funannotate annotate -i genome2 -e palmer3@wisc.edu --cpus 6 --iprscan iprscan_results'
echo $cmd; eval $cmd

cmd='funannotate annotate -i genome3 -e palmer3@wisc.edu --cpus 6 --iprscan iprscan_results'
echo $cmd; eval $cmd

#now run compare
cmd='funannotate compare -i genome1 genome2 genome3 --cpus 6 --outgroup botrytis_cinerea.dikarya --run_dnds estimate'
echo $cmd; eval $cmd

#clean up augustus training
rm -r /opt/augustus-3.2.1/config/species/genome_four/


