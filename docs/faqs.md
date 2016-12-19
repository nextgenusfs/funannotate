##Frequently Asked Questions

###1) Homebrew/LinuxBrew installation failed?
Due to the large numer of dependencies that funannotate uses, I've tried to use Homebrew/Linuxbrew to make this easier to install.  On most systems, this works well.  On other systems, this does not work very well (notably large compute clusters, shared servers, etc).  Funannotate was built with a single user in mind, i.e. processes are spawned over several CPUs, but not across several nodes.  At any rate, you don't have to use Homebrew/Linuxbrew to install funannotate, you can do it manually and `git clone` the repository and then setup the external tools manually as well.  There are several scripts that can help you, notably the `funannotate check` will get you version of several python modules are installed.  Additionally `funannotate setup` will aid you in downloaded databases and checking external dependencies.

###2) I got an error, what should I do?
You can file an issue on GitHub.  Be sure to include the command you issued and the traceback error that you got in the script.  Attaching logfiles is also helpful to diagnose the problem.  Most of the time, problems arise because a particular external tool was not installed correctly.  There are of course bugs in the scripts, when/if you find them, please report them so I can get them fixed.

###3) Can I use a non-fungal genome?
Prokaryotes are not supported -> use [Prokka](https://github.com/tseemann/prokka).  If it is a Eukaryotic genome, then yes, at least I think so.  Options that you need to be conscientious about are setting the `--organism other` option, the `--busco_db` option, `--eggnog_db` option, as well you should increase the `--max_intronlen` option to at least 10000. The `--max_intronlen` is an important distinction between fungi and other higher eukaryotes, fungi typically have much higher gene density, so if you set the the max_intronlen to a higher value, you will see many gene models merge together.  The opposite then is likely true, running a higher eukaryote with the default setting will result in splitting gene models apart.

###4) How does funannotate train Augustus?
Training Augustus is not very easy.  There are several ways to do it in funannotate and the script will automatically pick a training path based on your input data.  For all of these training steps, the more evidence you can provide the better your training will be (`--protein_evidence` and `--transcript_evidence`). This is how the "logic" in the script is setup.  
1) If you pass a valid pre-trained species to `--augustus_species` or there is already one trained (`--species "Aspergillus nidulans"` will essentially be turned into `--augustus_species aspergillus_nidulans`) then the scripts will NOT train Augustus and will use the pre-trained parameters.  Note you can check which species have been pretrained with `funannotate species`. 
2) If you provide a coordinate sorted BAM file via `--rna_bam`, Augustus and GeneMark will be trained using BRAKER1.  
3) If you provide a PASA GFF file via `--pasa_gff` then Augustus will be trained using these PASA gene models. 
4) If you don't have PASA or a RNAseq BAM file, then Augustus will be trained using BUSCO2.  The `--busco_seed_species` option is for passing the most closely related pre-trained Augustus species parameter to BUSCO2 to improve its de novo prediction.  Funannotate uses a modified training regime where it takes BUSCO2 'Complete' models, de novo GeneMarkES models, and evidence in those regions and runs EvidenceModeler to predict gene models.  The models are then confirmed using BUSCO2 and a subset are used for training Augustus.

###5) Funannotate said I should manually fix problematic gene models, how???
In the 'predict_results' folder you will find the output from `funannotate predict` which is composed of a GenBank flatfile, feature table file, GFF3, proteins, transcripts, as well as 3 error reports from tbl2asn.  Gene models that show up as ERROR in the error.summary.txt file MUST be fixed prior to submission to NCBI.  All errors listed as FATAL in the discrepency.report.txt must also be fixed (with the exception of FATAL: DISC_BACTERIAL_PARTIAL_NONEXTENDABLE_PROBLEMS).  I try to parse the errors where I can automatically provide fixes or removing the gene models, however there are lots of tbl2asn errors I've either never seen before or don't know how to fix automatically.  Here is how you can fix those problematic gene models:
```
#Modify the GFF file in 'predict_results', then run it through GAG
gag.py -f $genome.scaffolds.fa -g $genome.fixed.gff -o gag

#move into output folder (monitor any errors you see) and change a file name
cd gag
mv genome.fasta genome.fsa

#run tbl2asn on this folder
tbl2asn -p . -t /path/to/funannotate/libexec/lib/test.sbt -M n -Z discrep -a r10u -l paired-ends -V b -c fx -j "[organism=Genus awesomenus] [isolate=GA1]"
```
You will now have to look again at the discrpenecy and error reports (discrep, genome.val, errorsummary.val).  This routine can be run until you are satisfied.  The last thing to do is to replace the files in the 'predict_results' folder.  I typically keep the old files and just rename them like so:
```
mv $genome.scaffolds.fa $genome.scaffolds.fa.bak
mv $genome.gff3 $genome.gff3.bak
mv $genome.gbk $genome.gbk.bak
etc....
```
Finally then you can move the updated files into the folder
cp gag/genome.gbf $genome.gbk
cp gag/genome.gff $genome.gff3
cp gag/genome.fsa $genome.scaffolds.fa



