
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

    --------------------------------------------------------------
    Funannotate Databases currently installed:
    --------------------------------------------------------------
    Database       Type      Version      Date         Num_Records   Md5checksum                       
    pfam           hmmer3    31.0         2017-02      16712         b6fda5fdc90d24fbc1484d3b641d4e32  
    gene2product   text      1.2          2017-12-11   24728         a75679565d42bb93d0a343d18c631ff6  
    interpro       xml       66.0         2017-11-23   32568         c230e27471cad4a352cbe45beada0c52  
    dbCAN          hmmer3    6.0          2017-09-12   585           3cb06f6f93c72a56c9fa12a6294b41d5  
    merops         diamond   12.0         2017-10-04   4968          d923f0177c6d27c3d2886c705347adc0  
    mibig          diamond   1.3          2017-12-02   24085         84b3cd16e0b3b074e4b7ee18c6aa31fd  
    uniprot        diamond   2017_11      2017-11-22   556196        90c8910ef07e3ac601cdb43462f82c45  
    go             text      2017-12-01   2017-12-01   47071         f561c5193bcc9048fb337a0d7bb44b24  
    repeats        diamond   1.0          2017-12-02   11950         4e8cafc3eea47ec7ba505bb1e3465d21  
    --------------------------------------------------------------

Similarly, database sources can be updated with the :code:`funannotate setup` command, for example to update the gene2product database to its most recent version you would run:

.. code-block:: none

    $ funannotate setup -d $HOME/funannotate_db -i gene2product --update
    
    
