
.. _compare:

Comparative genomics
================================
A typical workflow in a genomics project would be to compare your newly sequenced/assembled/annotated genome to other organisms. The impetus behind :code:`funannotate compare` was that there was previously no way to easily compare multiple genomes. Funannotate stores all annotation in GenBank flat file format, while some people don't like this format as it is difficult to parse with standard unix tools, the main advantage is that the annotation can be stored in a standardized format and retrieved in the same way for each genome. GFF3 is the common output of many annotation tools, however, this doesn't work well for functional annotation as all of the "information" is stored in a single column.  At any rate, :code:`funannotate compare` can take either folders containing "funannotated" genomes or GBK files --> the output is stats, graphs, CSV files, phylogeny, etc all summarized in HTML format.

.. code-block:: none
    
	Usage:       funannotate compare <arguments>
	version:     1.7.0

	Description: Script does light-weight comparative genomics between funannotated genomes.  Output
				 is graphs, phylogeny, CSV files, etc --> visualized in web-browser.  
	
	Required:    
	  -i, --input         List of funannotate genome folders or GBK files

	Optional:    
	  -o, --out           Output folder name. Default: funannotate_compare
	  -d, --database      Path to funannotate database. Default: $FUNANNOTATE_DB
	  --cpus              Number of CPUs to use. Default: 2
	  --run_dnds          Calculate dN/dS ratio on all orthologs. [estimate,full]
	  --go_fdr            P-value for FDR GO-enrichment. Default: 0.05
	  --heatmap_stdev     Cut-off for heatmap. Default: 1.0
	  --num_orthos        Number of Single-copy orthologs to use for ML. Default: 500
	  --bootstrap         Number of boostrap replicates to run with RAxML. Default: 100
	  --outgroup          Name of species to use for ML outgroup. Default: no outgroup
	  --proteinortho      ProteinOrtho5 POFF results.
	  --ml_method         Maxmimum Liklihood method: Default: raxml [raxml,iqtree]

