# funannotate

###README is under construction....

funannotate is a pipeline for genome annotation (built specifically for fungi).  Genome annotation is a complicated process that uses software from many sources, thus the hardest part about using funannotate will be getting all of the dependencies installed.

####Python Dependencies:
* Python 2
* Biopython
* psutil
* natsort
* goatools

####Software Dependencies:
* Perl
* BioPerl
* Blast+
* Hmmer3
* RepeatModeler
* RepeatMasker
* GMAP
* Blat - if using PASA results to train Augustus
* pslCDnaFilter (kent tools) - if using PASA results to train Augustus
* BedTools
* Augustus
* GeneMark-ES/ET (gmes_petap.pl)
* BamTools
* Genome Annotation Generator (gag.py)
* tbl2asn
* BRAKER1 (optional if training Augustus with RNA-seq data BAM file)


####Environmental variables:
EVM_HOME, GENEMARK_PATH, BAMTOOLS_PATH, AUGUSTUS_CONFIG_PATH

