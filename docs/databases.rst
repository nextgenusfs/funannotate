
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
	  pfam              hmmer3      32.0         2018-08            17929   de7496fad69c1040fd74db1cb5eef0fc
	  gene2product      text        1.45         2019-07-31         30103   657bb30cf3247fcb74ca4f51a4ab7c18
	  interpro          xml         76.0         2019-09-18         37113   328f66a791f9866783764f24a74a5aa3
	  dbCAN             hmmer3      8.0          2019-08-08           607   51c724c1f9ac45687f08d0faa689ed58
	  busco_outgroups   outgroups   1.0          2019-10-20             7   6795b1d4545850a4226829c7ae8ef058
	  merops            diamond     12.0         2017-10-04          5009   a6dd76907896708f3ca5335f58560356
	  mibig             diamond     1.4          2019-10-20         31023   118f2c11edde36c81bdea030a0228492
	  uniprot           diamond     2019_09      2019-10-16        561176   9fc7871b8c4e3b755fe2086d77ed0645
	  go                text        2019-10-07   2019-10-07         47375   3bc9ba43a98bf8fcd01db6e7e7813dd2
	  repeats           diamond     1.0          2019-10-20         11950   4e8cafc3eea47ec7ba505bb1e3465d21

	To update a database type:
		funannotate setup -i DBNAME -d $HOME/funannotate_db --force

	To see install BUSCO outgroups type:
		funannotate database --show-outgroups

	To see BUSCO tree type:
		funannotate database --show-buscos



Similarly, database sources can be updated with the :code:`funannotate setup` command, for example to update the gene2product database to its most recent version you would run:

.. code-block:: none

    $ funannotate setup -d $HOME/funannotate_db -i gene2product --update
    
    
