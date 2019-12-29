.. Funannotate documentation master file, created by
   sphinx-quickstart on Sat Nov 18 22:41:39 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Funannotate documentation
=======================================

.. toctree::
   :hidden:
  
   install
   prepare
   predict
   evidence
   update
   annotate
   compare
   databases
   tutorials
   commands
   utilities


Funannotate is a genome prediction, annotation, and comparison software package. It was originally written to annotate fungal genomes (small eukaryotes ~ 30 Mb genomes), but has evolved over time to accomodate larger genomes. The impetus for this software package was to be able to accurately and easily annotate a genome for submission to NCBI GenBank. Existing tools (such as Maker) require significant manually editing to comply with GenBank submission rules, thus funannotate is aimed at simplifying the genome submission process.

Funannotate is also a lightweight comparative genomics platform. Genomes that have had functional annotation added via the :code:`funannotate annotate` command can be run through the :code:`funannotate compare` script that outputs html based whole genome comparisons. The software can run orthologous clustering, construct whole-genome phylogenies, run Gene Ontology enrichment analysis, as well as calculate dN/dS ratios for orthologous clusters under positive selection.



* :ref:`install`
* :ref:`prepare`
* :ref:`predict`
* :ref:`update`
* :ref:`annotate`
* :ref:`compare`
* :ref:`tutorials`
* :ref:`utilities`
