# funannotate

funannotate is a pipeline for genome annotation (built specifically for fungi).  Genome annotation is a complicated process that uses software from many sources, thus the hardest part about using funannotate will be getting all of the dependencies installed.

####Python Dependencies:
* Python 2
* Biopython
* psutil
* natsort

####Software Dependencies:
* Blast+
* Hmmer3
* RepeatModeler
* RepeatMasker
* GMAP
* Blat
* pslCDnaFilter (kent tools)
* BedTools
* Augustus
* GeneMark-ES/ET (gmes_petap.pl)
* BamTools
* Genome Annotation Generator (gag.py)
* tbl2asn
* BRAKER1 (optional if training with RNA-seq data)


####Environmental variables:
EVM_HOME, GENEMARK_PATH, BAMTOOLS_PATH, AUGUSTUS_CONFIG_PATH
