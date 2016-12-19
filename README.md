# funannotate

funannotate is a pipeline for genome annotation (built specifically for fungi, but theoretically should work with other eukaryotes).  Genome annotation is a complicated process that uses software from many sources, thus the hardest part about using funannotate will be getting all of the dependencies installed.  After that, funannotate requires only a few simple commands to go from genome assembly all the way to a functional annotated genome (InterPro, PFAM, MEROPS, CAZymes, GO ontology, etc) that is ready for submission to NCBI.  Moreover, funannotate incorporates a light-weight comparative genomics package that can get you started looking at differences between funannotated genomes.

###Installation

funannotate will likely run on any POSIX system, although it has only been tested on Mac OSX and Ubuntu.

* [Mac OSX install instructions](docs/mac_install.md)
* [Ubuntu install instructions](docs/ubuntu_install.md)
* [FAQS](docs/faqs.md)

###Setup

See installation instructions, but funannotate will also download and format the databases that are required. The databases will occupy quite a bit of space, currently working (uncompressed) is ~ 24 GB.  Whenever possible, funannotate is configured to check external dependencies at runtime, however, `funannotate setup` will help you during the initial setup.

To run the setup script, type:
```
funannotate setup --all
```

Most problems that people have are with dependencies and installation of funannotate.  Here are some Frequently Asked Questions: [FAQ](cods/faqs.md)

###Funannotate help menu

To see the help menu, simply type `funannotate` in the terminal window.  Similarly, e.g `funannotate predict` without any arguments will give you the options available to pass to each script, this is consistent for all of the funannotate commands.
```
$  funannotate

Usage:       funannotate <command> <arguments>
version:     0.3.7

Description: Funannotate is a genome prediction, annotation, and comparison pipeline.
    
Command:     clean          Find/remove small repetitive contigs
             sort           Sort by size and rename contig headers (recommended)
             species        list pre-trained Augustus species
             
             predict        Run gene prediction pipeline
             annotate       Assign functional annotation to gene predictions
             compare        Compare funannotated genomes
             
             fix            Remove adapter/primer contamination from NCBI error report
             check          Check Python module versions installed
             
             setup          Setup/Install databases and check dependencies    
```

###Using funannotate: a simple walkthrough

move into the `sample_data` directory of funannotate.

```
#for example, funannotate on Mac OSX with HomeBrew
$ cd /usr/local/opt/funannotate/libexec/sample_data

#run funannotate predict on genome 1
$ funannotate predict -i genome1.fasta -o genome1 -s "Genome one" \
    --isolate fun1 --name GN01_ --augustus_species botrytis_cinerea \
    --protein_evidence proteins.fa --transcript_evidence transcripts.fa --cpus 6
```
This command should complete in ~ 5 minutes, will produce an output folder named `genome1` which contains the results.  To save time, here we are using pre-trained botrytis_cinerea to run AUGUSTUS - normally funannotate will train AUGUSTUS for you depending on which input you give it.  

```
#generate functional annotation for genome 1
$ funannotate annotate -i genome1 -e youremail@domain.edu --cpus 6
```
The second command, will add functional annotation to your protein models.  It should complete in ~ 15 minutes - 30 minutes (depending on how long remote query to InterProScan server takes.  Note you could also run InterProScan locally, funannotate requires the results to be in XML format one file per protein.  The results are in the `annotate_results` folder and have all the necessary files for NCBI WGS submission (.tbl, .sqn, .contigs.fsa, .agp).  A GBK flatfile is also provided.

You can now run similar commands for genome2.fasta and genome3.fasta
```
#predict genome2
$ funannotate predict -i genome2.fasta -o genome2 -s "Genome two" \
    --isolate fun2 --name GN02_ --augustus_species botrytis_cinerea \
    --protein_evidence proteins.fa --transcript_evidence transcripts.fa --cpus 6
    
#annotate genome2
$ funannotate annotate -i genome2 -e youremail@domain.edu --cpus 6

#predict genome3
$ funannotate predict -i genome3.fasta -o genome3 -s "Genome three" \
    --isolate fun3 --name GN03_ --augustus_species botrytis_cinerea \
    --protein_evidence proteins.fa --transcript_evidence transcripts.fa --cpus 6

#annotate genome3
$ funannotate annotate -i genome3 -e youremail@domain.edu --cpus 6
```

You can now run some "lightweight" comparative genomics on these funannotated genomes:

```
$ funannotate compare -i genome1 genome2 genome3 --outgroup Botrytis_cinerea
```

You can now visualize the results by opening up the `index.html` file produced in the `funannotate_compare` folder.  A phylogeny inferred from RAxML, genome stats, orthologs, InterPro summary, PFAM summary, MEROPS, CAZymes, and GO ontology enrichment results are all included in the browser-based output.  Additionally, the raw data is available in appropriate files in the output directory.

###Using funannotate: a more realistic walkthrough

Here is a step by step tutorial for annotating a genome using funannotate with a genome assembly and RNA-seq data.  This is a list of the data that I have available:
```
#my genome assembly
genome.scaffolds.fa

#my RNA Seq Reads
forward_R1.fastq
reverse_R2.fastq
```

1) Find/remove small repetitive contigs.  Assumption here is haploid fungal genome. (Optional) 
```
funannotate clean -i genome.scaffolds.fa -o genome.cleaned.fa
```
2) Now sort contigs by size and relabel header (Optional)
```
funannotate sort -i genome.cleaned.fa -o genome.final.fa
```
3) Align RNA seq reads to cleaned genome (I'll use hisat2 here, you could use a different aligner), but the BAM file needs to be sorted, e.g. `samtools sort`.
```
#build reference
hisat2-build genome.final.fa genome

#now align reads to reference, using 12 cpus
hisat2 --max_intronlen 3000 -p 12 -x genome -1 forward_R1.fastq -2 reverse_R2.fastq \
       | samtools view -bS - | samtools sort -o RNA_seq.bam - 
```
3b) You can also run something like Trinity/PASA/TransDecoder with the RNA-seq reads, for instructions see [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running%20Trinity) and [PASA](http://pasapipeline.github.io).
```
#run Trinity de novo
Trinity.pl --left forward_R1.fastq --right reverse_R2.fastq --max_memory 50G --CPU 6 \
    --jaccard_clip --trimmomatic --normalize_reads

#run Trinity genome guided
Trinity.pl --genome_guided_bam RNA_seq.bam --max_memory 50G --genome_guided_max_intron 3000 --CPU 6

#Concatenate the output of each
cat Trinity.fasta Trinity-GG.fasta > transcripts.fasta

#create transcript accessions
$PASA_HOME/misc_utilities/accession_extractor.pl < Trinity.fasta > tdn.accs

#now run PASA
$PASA_HOME/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g genome.final.fa \
    -t transcripts.fasta --TDN tdn.accs --ALIGNERS blat,gmap

#run PASA mediated TransDecoder
$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta genome.assemblies.fasta \
    --pasa_transcripts_gff3 genome.pasa_assemblies.gff3
```
4) Now you have an assembly, RNA_seq.bam, PASA_assemblies.gff3, and transcripts you can run funannotate like so:
```
funannotate predict -i genome.final.fa -o fun_out --species "Fungus specious" \
    --pasa_gff genome.fasta.transdecoder.genome.gff3 --rna_bam RNA_seq.bam \
    --transcript_evidence transcripts.fasta --cpus 12
```
This command will first run RepeatModeler on your genome, soft-mask repeats using RepeatMasker, align UniProtKB proteins to genome using tblastn/exonerate, align transcripts.fasta to genome using GMAP, launch BRAKER to train/run AUGUSTUS and GeneMark-ET, combine all predictions and evidence into gene models using Evidence Modeler, predict tRNAs, filter bad gene models, rename gene models, and finally convert to GenBank format.

5) You should now examine the output in the `fun_out/predict_results` folder, be sure to look at the NCBI discrepency report to identify any gene models that need to be adjusted manually.

6) If you are interested in secondary metabolism gene clusters, submit your genome.gbk file to the [antiSMASH](http://antismash.secondarymetabolites.org) web server. You can then download the results in genbank format and funannotate can parse them.

7) To add functional annotation to your genome, you would run the following command:
```
funannotate annotate -i fun_out -e youremail@domain.edu --antismash scaffold_1.final.gbk \
     --sbt my_ncbi_template.sbt --cpus 12
```
Your results will be located in the `fun_out/annotate_results` folder.  It contains the necessary files to submit to NCBI WGS submission system (assuming you have passed in an appropriate submission template), otherwise the template is a generic one used by funannotate.

####What if I've already run Maker2, can I use funannotate?

Yes, you can.  As of `v0.1.5` you can pass your Maker2 GFF file to the `--maker_gff` option of `funannotate predict` which will parse the alignment evidence and ab initio gene predictions from Maker into a format for EVidence Modeler.  So the `--maker_gff` will bypass gene predictions and evidence alignments done by funannotate and proceed to EVM - and then the rest of the script will run normally (filtering gene models and converting to GenBank).  For example:
```
#simple example
funannotate predict -i genome1.fasta --species "Genome one" -o test_output \
    --maker_gff maker_genome1.all.gff --cpus 6

#maker + pasa?
funannotate predict -i genome1.fasta --species "Genome one" -o test_output \
    --maker_gff maker_genome1.all.gff --pasa_gff my_pasa.gff --cpus 6

```
