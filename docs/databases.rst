
.. _databases:

Annotation Databases
================================
 
Funannotate uses several publicly available databases, they can be installed with the :code:`funannotate setup` command.  The currently installed databases and version numbers can be displayed with the :code:`funannotate database` command.

Initial setup is simple and requires only a path to a database location, this can (should) be set using the $FUNANNOTATE_DB environmental variable.  If $FUNANNOTATE_DB is set, then the script will use that location by default, otherwise you will need to specify a location to the script i.e.:

.. code-block:: none

    funannotate setup -d $HOME/funannotate_db
    
    
You could then update the databases if $FUNANNOTATE_DB is set like this:

.. code-block:: none

    funannotate setup -i all --update
    
    #or force update of just one database
    funannotate setup -i uniprot --force
    

This will download and format the databases, they can be displayed like so:

.. code-block:: none

    $ funannotate database

	Funannotate Databases currently installed:

	  Database          Type        Version      Date         Num_Records   Md5checksum                     
          merops            diamond     12.5         2023-01-19          5098   6cd3c3dd85650394ce4e3dacb591f2a5
          uniprot           diamond     2024_01      2024-01-24        570830   c7507ea16b3c4807971c663994cad329
          dbCAN             hmmer3      11.0         2022-08-09           699   fb112af319a5001fbf547eac29e7c3b5
          pfam              hmmer3      36.0         2023-07            20795   0725495ccf049a4f198fcc0a92f7f38c
          repeats           diamond     1.0          2022-03-13         11950   4e8cafc3eea47ec7ba505bb1e3465d21
          go                text        2024-01-17   2024-01-17         47729   7e6b9974184dda306e6e07631f1783af
          mibig             diamond     1.4          2022-03-13         31023   118f2c11edde36c81bdea030a0228492
          interpro          xml         98.0         2024-01-25         40768   502ea05009761b893dedb56d5ea89c48
          gene2product      text        1.92         2023-10-02         34459   32a4a80987720e0872377de3207dc0f5

	To update a database type:
		funannotate setup -i DBNAME -d $HOME/funannotate_db --force

	To see BUSCO tree type:
		funannotate database --show-buscos



Similarly, database sources can be updated with the :code:`funannotate setup` command, for example to update the gene2product database to its most recent version you would run:

.. code-block:: none

    $ funannotate setup -d $HOME/funannotate_db -i gene2product --update
    
    


BUSCO lineages
--------------

BUSCO models are installed with :code:`funannotate setup -i busco -b LINEAGE` into :code:`$FUNANNOTATE_DB/LINEAGE` and selected with
:code:`--busco_db LINEAGE` in :code:`predict` and :code:`annotate`. The lineage name is the folder name. The dataset inside the download
has its own name, which is kept in :code:`LINEAGE/dataset.cfg`: it is the OrthoDB 9 dataset (:code:`_odb9`) for the lineages below, except
for :code:`sordariomycetes` and :code:`saccharomycetes`, whose odb9 datasets are named by taxonomic rank (:code:`sordariomyceta_odb9`,
:code:`saccharomyceta_odb9`), and :code:`protists` and :code:`alveolata_stramenophiles`, which are Ensembl-derived datasets
(:code:`_ensembl`). :code:`setup` renames the dataset folder to the lineage name. The BUSCO ids written to the annotation (:code:`BUSCO:...`)
come from these datasets, so they are not comparable with ids from BUSCO v5 or later datasets.
:code:`predict` writes the dataset, its creation date and its number of BUSCOs to :code:`logfiles/training_decisions.tsv` (stage
:code:`busco_training`, decision :code:`BUSCO lineage`), and :code:`annotate` logs it next to the BUSCO step.

.. list-table::
   :header-rows: 1

   * - Lineage (``--busco_db``)
     - Dataset
     - BUSCOs
   * - ``fungi``
     - ``fungi_odb9``
     - 290
   * - ``microsporidia``
     - ``microsporidia_odb9``
     - 518
   * - ``dikarya``
     - ``dikarya_odb9``
     - 1312
   * - ``ascomycota``
     - ``ascomycota_odb9``
     - 1315
   * - ``pezizomycotina``
     - ``pezizomycotina_odb9``
     - 3156
   * - ``eurotiomycetes``
     - ``eurotiomycetes_odb9``
     - 4046
   * - ``sordariomycetes``
     - ``sordariomyceta_odb9``
     - 3725
   * - ``saccharomycetes``
     - ``saccharomyceta_odb9``
     - 1759
   * - ``saccharomycetales``
     - ``saccharomycetales_odb9``
     - 1711
   * - ``basidiomycota``
     - ``basidiomycota_odb9``
     - 1335
   * - ``eukaryota``
     - ``eukaryota_odb9``
     - 303
   * - ``protists``
     - ``protists_ensembl``
     - 215
   * - ``alveolata_stramenophiles``
     - ``alveolata_stramenophiles_ensembl``
     - 234
   * - ``metazoa``
     - ``metazoa_odb9``
     - 978
   * - ``nematoda``
     - ``nematoda_odb9``
     - 982
   * - ``arthropoda``
     - ``arthropoda_odb9``
     - 1066
   * - ``insecta``
     - ``insecta_odb9``
     - 1658
   * - ``endopterygota``
     - ``endopterygota_odb9``
     - 2442
   * - ``hymenoptera``
     - ``hymenoptera_odb9``
     - 4415
   * - ``diptera``
     - ``diptera_odb9``
     - 2799
   * - ``vertebrata``
     - ``vertebrata_odb9``
     - 2586
   * - ``actinopterygii``
     - ``actinopterygii_odb9``
     - 4584
   * - ``tetrapoda``
     - ``tetrapoda_odb9``
     - 3950
   * - ``aves``
     - ``aves_odb9``
     - 4915
   * - ``mammalia``
     - ``mammalia_odb9``
     - 4104
   * - ``euarchontoglires``
     - ``euarchontoglires_odb9``
     - 6192
   * - ``laurasiatheria``
     - ``laurasiatheria_odb9``
     - 6253
   * - ``embryophyta``
     - ``embryophyta_odb9``
     - 1440


OrthoDB 10 BUSCO lineages
-------------------------

Newer lineages (OrthoDB 10, for example :code:`fungi_odb10` with 758 BUSCOs, :code:`ascomycota_odb10` with 1,706) can be installed next to
the odb9 ones, from the BUSCO download server rather than the funannotate mirror:

.. code-block:: none

    funannotate setup -d $FUNANNOTATE_DB -i busco -b fungi_odb10
    funannotate predict -i genome.fa -s "Genus species" --busco_db fungi_odb10 ...
    funannotate annotate -i predict_dir --busco_db fungi_odb10

- Any eukaryotic :code:`<lineage>_odb10` name in the BUSCO dataset list works (:code:`busco --list-datasets`); :code:`database --show-buscos` lists the
  fungal ones. :code:`-b all` still means the 28 odb9 lineages above. The download is checked against the md5 in the BUSCO dataset list.
- OrthoDB 12 datasets (:code:`_odb12`, :code:`_odb12.2`) are refused: the BUSCO that funannotate runs cannot read them.
- There is no :code:`dikarya`, :code:`pezizomycotina` or :code:`saccharomycetales` lineage in odb10; :code:`fungi_odb10` and :code:`ascomycota_odb10` are the closest.
- The BUSCO ids of odb9 (:code:`EOG09260004`) and odb10 (:code:`100957at4751`) are different orthogroup definitions with no mapping between them.
  Annotations made with different lineages must not be compared by BUSCO id, and the same applies to the BUSCO column of :code:`funannotate-compare`.
- Each installed lineage is recorded in :code:`funannotate-db-info.txt` (dataset, creation date, number of BUSCOs, md5), and the dataset used by a run is
  written to :code:`logfiles/training_decisions.tsv` by :code:`predict`.
- :code:`fungi_odb10` has 758 BUSCOs against 1,312 for :code:`dikarya`, so BUSCO-based training starts from about half as many gene models. On a small
  or partial assembly this can fall below :code:`--min_training_models` (default 200); use a larger lineage such as :code:`ascomycota_odb10` (1,706)
  or :code:`basidiomycota_odb10` (1,764), or lower the minimum. :code:`funannotate test -t busco --busco_db fungi_odb10` passes
  :code:`--min_training_models 100` for this reason.

