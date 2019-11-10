
.. _annotate:

Functional annotation
================================
 
After your genome has gone through the gene prediction module and you have gene models that pass NCBI specs the next step is to add functional annotate to the protein-coding genes. Funannotate accomplishes this using several curated databases and is run using the :code:`funannotate annotate` command. 

Funannotate will parse the protein-coding models from the annotation and identify Pfam domains, CAZYmes, secreted proteins, proteases (MEROPS), and BUSCO groups.  If you provide the script with InterProScan5 data :code:`--iprscan`, funannotate will also generate additional annotation: InterPro terms, GO ontology, and fungal transcription factors. If Eggnog-mapper is installed locally or you pass eggnog results via :code:`--eggnog`, then Eggnog annotations and COGs will be added to the functional annotation.  The scripts will also parse UniProtKb/SwissProt searches with Eggnog-mapper searches (optional) to generate gene names and product descriptions. 

InterProScan5 and Eggnog-Mapper are two functional annotation pipelines that can be parsed by funannotate, however due to the large database sizes they are not run directly.  If :code:`emapper.py` (Eggnog-mapper) is installed, then it will be run automatically during the functional annotation process. Because InterProScan5 is Linux only, it must be run outside funannotate and the results passed to the script. If you are on Mac, I've included a method to run InterProScan5 using Docker and the :code:`funannotate predict` output will let the user know how to run this script.  Alternatively, you can run the InterProScan5 search remotely using the :code:`funannotate remote` command.

Phobius and SignalP will be run automatically if they are installed (i.e. in the PATH), however, Phobius will not run on Mac.  If you are on Mac you can run Phobius with the :code:`funannotate remote` script. 

If you are annotating a fungal genome, you can run Secondary Metabolite Gene Cluster prediction using antiSMASH.  This can be done on the webserver, submit your GBK file from predict (predict_results/yourGenome.gbk) or alternatively you can submit from the command line using :code:`funannotate remote`.  Of course, if you are on Linux you can install the antiSMASH program locally and run that way as well.  The annotated GBK file is fed back to this script with the :code:`--antismash` option.

Similarily to :code:`funannotate predict`, the output from :code:`funannotate annotate` will be populated in the output/annotate_results folder. The output files are:

+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| **File Name**                      | **Description**                                                                                                                  |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.gbk                       | Annotated Genome in GenBank Flat File format                                                                                     |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.contigs.fsa               | Multi-fasta file of contigs, split at gaps (use for NCBI submission)                                                             |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.agp                       | AGP file; showing linkage/location of contigs (use for NCBI submission)                                                          |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.tbl                       | NCBI tbl annotation file (use for NCBI submission)                                                                               |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.sqn                       | NCBI Sequin genome file (use for NCBI submission)                                                                                |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.scaffolds.fa              | Multi-fasta file of scaffolds                                                                                                    |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.proteins.fa               | Multi-fasta file of protein coding genes                                                                                         |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.transcripts.fa            | Multi-fasta file of transcripts (mRNA)                                                                                           |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.discrepency.report.txt    | tbl2asn summary report of annotated genome                                                                                       |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Basename.annotations.txt           | TSV file of all annotations added to genome. (i.e. import into excel)                                                            |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Gene2Products.must-fix.txt         | TSV file of Gene Name/Product deflines that failed to pass tbl2asn checks and must be fixed                                      |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Gene2Products.need-curating.txt    | TSV file of Gene Name/Product defines that need to be curated                                                                    |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Gene2Products.new-names-passed.txt | TSV file of Gene Name/Product deflines that passed tbl2asn but are not in Gene2Products database. Please submit a PR with these. |
+------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+

.. code-block:: none

	$ funannotate annotate

	Usage:       funannotate annotate <arguments>
	version:     1.7.0

	Description: Script functionally annotates the results from funannotate predict.  It pulls
				 annotation from PFAM, InterPro, EggNog, UniProtKB, MEROPS, CAZyme, and GO ontology.
	
	Required:    
	  -i, --input        Folder from funannotate predict
		or
	  --genbank          Genome in GenBank format
	  -o, --out          Output folder for results
		or   
	  --gff              Genome GFF3 annotation file
	  --fasta            Genome in multi-fasta format
	  -s, --species      Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
	  -o, --out          Output folder for results

	Optional:    
	  --sbt              NCBI submission template file. (Recommended)
	  -a, --annotations  Custom annotations (3 column tsv file)
	  --eggnog           Eggnog-mapper annotations file (if NOT installed)
	  --antismash        antiSMASH secondary metabolism results (GBK file from output)
	  --iprscan          InterProScan5 XML file
	  --phobius          Phobius pre-computed results (if phobius NOT installed)
	  --isolate          Isolate name
	  --strain           Strain name
	  --rename           Rename GFF gene models with locus_tag from NCBI.
	  --fix              Gene/Product names fixed (TSV: GeneID	Name	Product)
	  --remove           Gene/Product names to remove (TSV: Gene	Product)
	  --busco_db         BUSCO models. Default: dikarya
	  -t, --tbl2asn      Additional parameters for tbl2asn. Default: "-l paired-ends"
	  -d, --database     Path to funannotate database. Default: $FUNANNOTATE_DB
	  --force            Force over-write of output folder
	  --cpus             Number of CPUs to use. Default: 2



