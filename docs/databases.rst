
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
          busco_outgroups   outgroups   1.0          2024-03-04             8   6795b1d4545850a4226829c7ae8ef058
          gene2product      text        1.92         2023-10-02         34459   32a4a80987720e0872377de3207dc0f5

	To update a database type:
		funannotate setup -i DBNAME -d $HOME/funannotate_db --force

	To see install BUSCO outgroups type:
		funannotate database --show-outgroups

	To see BUSCO tree type:
		funannotate database --show-buscos



Similarly, database sources can be updated with the :code:`funannotate setup` command, for example to update the gene2product database to its most recent version you would run:

.. code-block:: none

    $ funannotate setup -d $HOME/funannotate_db -i gene2product --update
    
    
