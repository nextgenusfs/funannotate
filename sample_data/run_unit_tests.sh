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

#test BRAKER1 training
cmd='funannotate predict -i genome5.fasta -s "Genome five" -o genome5 --protein_evidence proteins.fa --transcript_evidence genome5.transcripts.fa --rna_bam genome5.bam --cpus 6'
echo $cmd; eval $cmd

#test maker GFF
cmd='funannotate predict -i genome1.fasta -s "Genome one" --name GN12_ -o genome1_maker --protein_evidence proteins.fa --transcript_evidence transcripts.fa --augustus_species botrytis_cinerea --cpus 6 --maker_gff maker_genome1.all.gff'
echo $cmd; eval $cmd

#now annotate each genome
cmd='funannotate annotate -i genome1 --cpus 6 --iprscan test_data.iprscan.xml --eggnog genome1.emapper.annotations'
echo $cmd; eval $cmd

cmd='funannotate annotate -i genome2 --cpus 6 --iprscan test_data.iprscan.xml --eggnog genome2.emapper.annotations'
echo $cmd; eval $cmd

cmd='funannotate annotate -i genome3 --cpus 6 --iprscan test_data.iprscan.xml --eggnog genome3.emapper.annotations'
echo $cmd; eval $cmd

#test annotation using direct input
cmd='funannotate annotate --gff genome1/predict_results/genome_one.gff3 --fasta genome1/predict_results/genome_one.scaffolds.fa --proteins genome1/predict_results/genome_one.proteins.fa --iprscan test_data.iprscan.xml --eggnog genome1.emapper.annotations -o direct -s "Aspergillus fumigatus"'
echo $cmd; eval $cmd

#now run compare
cmd='funannotate compare -i genome1 genome2 genome3 --cpus 6 --outgroup botrytis_cinerea.dikarya --run_dnds estimate'
echo $cmd; eval $cmd

#test RNAseq modules
cmd='funannotate train -i genome6.fasta -l genome6_R1.fq.gz -r genome6_R2.fq.gz --stranded RF --species "Rubeus macgubis" --cpus 6 -o genome6'
echo $cmd; eval $cmd
cmd='funannotate predict -i genome6.fasta --transcript_evidence genome6/training/funannotate_train.trinity-GG.fasta --rna_bam genome6/training/funannotate_train.coordSorted.bam --pasa_gff genome6/training/funannotate_train.pasa.gff3 -o genome6 -s "Rubeus macgubis" --cpus 6'
echo $cmd; eval $cmd
cmd='funannotate update -i genome6 --cpus 6'
echo $cmd; eval $cmd

#clean up augustus training
rm -r $AUGUSTUS_CONFIG_PATH/species/genome_four/
rm -r $AUGUSTUS_CONFIG_PATH/species/genome_five/
rm -r $AUGUSTUS_CONFIG_PATH/species/genome_six/

