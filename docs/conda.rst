.. _conda:

Conda mediated Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Funannotate 1.9 is installed in two parts:

1. A conda (or mamba) environment that provides the external tools
   (Augustus, BLAST, HMMER, HISAT2, samtools, Perl and R libraries, and so on).
2. The funannotate Python package, which you install into that environment
   from PyPI with ``pip``.

Three components are not taken from conda. EvidenceModeler (EVM), PASA and
Trinity are built from source with the Rust-accelerated forks that 1.9 uses.
The repository provides scripts that do this. GeneMark-ES/ET and
RepeatMasker libraries (Repbase/Dfam) are licensed software that you install
yourself.

The conda environment is generated from the ``[dependencies]`` table of
``pixi.toml``. This is the same dependency list that the 1.9 container images
use. Linux x86-64 is the supported platform. macOS is experimental; see
:ref:`conda-macos`.

.. note::

    If you only want to run funannotate and do not need to modify the
    environment, the container image is simpler. See :ref:`containers`.
    If you prefer a locked, reproducible environment, you can use
    `pixi <https://pixi.sh>`_ in the repository directly: ``pixi install``,
    then ``pixi run install-externals``. The steps below do the same
    thing with plain conda.

1. Create the conda environment
--------------------------------

Get the environment file and the install scripts by cloning the repository at
the release tag. The install scripts are not part of the PyPI package.

.. code-block:: none

    git clone --branch v1.9.0 https://github.com/nextgenusfs/funannotate.git
    cd funannotate

    # use mamba (or conda) with the conda-forge and bioconda channels
    mamba env create -f environment.yml
    mamba activate funannotate

``environment.yml`` pins ``python>=3.6,<3.9`` and sets the channel order to
``conda-forge`` then ``bioconda``. Do not add the ``defaults`` channel.

2. Install funannotate from PyPI
--------------------------------

.. code-block:: none

    python -m pip install funannotate==1.9.0

The Python dependencies that conda already installed (biopython, numpy,
matplotlib, scikit-learn, and so on) are reused by ``pip``.
Check the installed version:

.. code-block:: none

    funannotate --version

3. Set the environment variables
--------------------------------

The Rust builds of EVM, PASA and Trinity install under ``$CONDA_PREFIX/opt``.
Funannotate finds them through these variables. Put them in an activation
script so that they are set each time you activate the environment:

.. code-block:: none

    mkdir -p $CONDA_PREFIX/etc/conda/activate.d $CONDA_PREFIX/etc/conda/deactivate.d

    cat > $CONDA_PREFIX/etc/conda/activate.d/funannotate.sh << 'EOF'
    export TRINITY_HOME=$CONDA_PREFIX/opt/trinityrnaseq
    export TRINITYHOME=$CONDA_PREFIX/opt/trinityrnaseq
    export PASAHOME=$CONDA_PREFIX/opt/pasa/src
    export EVM_HOME=$CONDA_PREFIX/opt/evm
    export PERL5LIB=$PASAHOME/SAMPLE_HOOKS:$PASAHOME/PerlLib${PERL5LIB:+:$PERL5LIB}
    export FUNANNOTATE_EVM_ENGINE=rust
    # fasta36 is installed under that name; funannotate calls it "fasta"
    [ -e $CONDA_PREFIX/bin/fasta ] || ln -s $CONDA_PREFIX/bin/fasta36 $CONDA_PREFIX/bin/fasta
    EOF

    cat > $CONDA_PREFIX/etc/conda/deactivate.d/funannotate.sh << 'EOF'
    unset TRINITY_HOME TRINITYHOME PASAHOME EVM_HOME FUNANNOTATE_EVM_ENGINE
    EOF

    # reload the environment so the variables take effect
    mamba deactivate && mamba activate funannotate

``$PERL5LIB`` is needed so that PASA can load its annotation hooks during
``funannotate update``. The scripts do not set ``FUNANNOTATE_DB``. Set it in your
shell startup file (for example ``~/.bashrc``) so that the location of the
databases does not change when you activate or deactivate the environment.

4. Build EVM, PASA, Trinity and bowtie2
---------------------------------------

Run these from the cloned repository, with the environment active. The build
needs the compilers and Rust toolchain that ``environment.yml`` installs.
Each script exits at once if its component is already built.

.. code-block:: none

    bash install_scripts/pixi_install_evm.sh
    bash install_scripts/pixi_install_pasa.sh
    bash install_scripts/pixi_install_trinity.sh
    bash install_scripts/pixi_install_bowtie2.sh

The scripts fetch the Rust forks from GitHub, so the build node needs
internet access. They need ``$CONDA_PREFIX`` to be set, which activating the
environment does. ``pixi_install_bowtie2.sh`` rebuilds bowtie2 from source. It fixes a problem
where the conda bowtie2 package does not start its AVX2-optimized binaries
correctly.

5. Install the licensed software and the databases
--------------------------------------------------

GeneMark-ES/ET is optional but recommended for ``funannotate predict``.
Download and license it from http://topaz.gatech.edu/GeneMark/license_download.cgi.
Then:

* Change the shebang line of each GeneMark Perl script to
  ``#!/usr/bin/env perl``.
* Add the directory that contains ``gmes_petap.pl`` to ``$PATH``, or set
  ``GENEMARK_PATH`` to that directory.

If you cannot install GeneMark, run it separately and pass the result to
``funannotate predict --genemark_gtf``.

RepeatMasker and RepeatModeler are installed by conda. Add Dfam or Repbase
libraries as described at http://www.repeatmasker.org if you need them.

EggNOG-mapper is installed by conda. Download its data with
``download_eggnog_data.py`` and set ``EGGNOG_DATA_DIR`` to that location
before you run ``funannotate annotate``. If the data are missing,
``funannotate annotate`` stops when it calls ``emapper.py``.

Then download the funannotate databases:

.. code-block:: none

    export FUNANNOTATE_DB=$HOME/funannotate_db
    funannotate setup -d $FUNANNOTATE_DB

6. Check the installation
-------------------------

.. code-block:: none

    # list every dependency and its version; fix any that are missing
    funannotate check --show-versions

    # run the test suite (needs internet access to download test data)
    funannotate test -t all --cpus 8

On a shared cluster, run the tests on a compute node, not on a login node.

.. _conda-macos:

macOS (experimental)
--------------------

macOS support is experimental. Use Linux or the container for production
annotation. The steps below were tested on one Apple Silicon Mac (macOS 26.7,
arm64, Xcode installed) with a native arm64 conda. They covered the environment
build, ``funannotate check``, the unit tests (``pytest tests``: 237 passed), the
EVM, PASA and Trinity builds, and the Trinity self-test. A full annotation run
(``funannotate train``, ``predict``, ``update``, ``annotate``) was **not** run on
macOS. Intel Macs were not tested.

**Use a native arm64 conda on Apple Silicon.** Check with
``conda info | grep platform``. It must say ``osx-arm64``. If you have an
x86-64 conda (it runs under Rosetta), install an arm64 Miniforge or Miniconda.
An x86-64 environment on an Apple Silicon Mac builds EVM and PASA, but the Trinity
build fails in ``bamsifter``: that Makefile calls the system ``g++`` (arm64) and
links against an x86-64 ``htslib``. Do not use it for Trinity.

**What is not available on macOS.** These have no working macOS package, so the
commands that need them do not run:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Missing
     - Effect
   * - ``tbl2asn``
     - ``funannotate annotate`` cannot make GenBank submission files.
   * - CodingQuarry
     - The CodingQuarry predictor in ``funannotate predict`` is not available.
   * - ``fasta`` (``fasta3``)
     - Not installed.
   * - RepeatModeler
     - Not installed. RepeatMasker is installed.
   * - R and Bioconductor packages that Trinity can use for differential expression
     - Not requested by ``environment-macos.yml``. Some R packages still come in as
       dependencies of other packages.
   * - GeneMark-ES/ET, SignalP
     - Licensed software, as on Linux.

1. Create the environment
~~~~~~~~~~~~~~~~~~~~~~~~~

Install the Xcode command line tools or Xcode (``xcode-select --install``). The
Trinity and PASA builds need the compilers. Then clone the repository and use the
macOS environment file, not ``environment.yml``:

.. code-block:: none

    git clone --branch v1.9.0 https://github.com/nextgenusfs/funannotate.git
    cd funannotate
    mamba env create -f environment-macos.yml
    mamba activate funannotate

``environment-macos.yml`` differs from ``environment.yml``. It leaves out the
packages above. It uses **Python 3.9** (Linux uses 3.8) because the bioconda
``kmer-jellyfish`` 2.x package for ``osx-arm64`` exists only for Python 3.9 and
newer, and Trinity needs jellyfish 2. If you change the Python pin, conda may
silently install jellyfish 1.1.12 and Trinity stops with
``need jellyfish v2``.

2. Install funannotate and set the variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

    python -m pip install funannotate==1.9.0

    mkdir -p $CONDA_PREFIX/etc/conda/activate.d $CONDA_PREFIX/etc/conda/deactivate.d

    cat > $CONDA_PREFIX/etc/conda/activate.d/funannotate.sh << 'EOF'
    export TRINITY_HOME=$CONDA_PREFIX/opt/trinityrnaseq
    export TRINITYHOME=$CONDA_PREFIX/opt/trinityrnaseq
    export PASAHOME=$CONDA_PREFIX/opt/pasa/src
    export EVM_HOME=$CONDA_PREFIX/opt/evm
    export FUNANNOTATE_EVM_ENGINE=rust
    # The conda perl-db_file package for osx-arm64 crashes on "use DB_File" (it
    # does not link libdb), and its copy hides perl's working core copy. Put perl's
    # core library first so that Trinity's insilico_read_normalization.pl runs.
    export PERL5LIB=$CONDA_PREFIX/lib/perl5/5.32/core_perl:$PASAHOME/SAMPLE_HOOKS:$PASAHOME/PerlLib${PERL5LIB:+:$PERL5LIB}
    EOF

    cat > $CONDA_PREFIX/etc/conda/deactivate.d/funannotate.sh << 'EOF'
    unset TRINITY_HOME TRINITYHOME PASAHOME EVM_HOME FUNANNOTATE_EVM_ENGINE
    EOF

    mamba deactivate && mamba activate funannotate

The openjdk package sets ``JAVA_HOME`` itself, so you do not need to. It has no
``java`` in ``$CONDA_PREFIX/bin``. ``which java`` shows ``/usr/bin/java``, the
macOS launcher, which uses ``JAVA_HOME``. Check with ``java -version``: it should
print OpenJDK 25.

The ``perl5/5.32`` directory name belongs to the Perl version in the
environment. If ``ls $CONDA_PREFIX/lib/perl5`` shows a different version,
change the path. The check for the problem is
``perl -e 'use DB_File'``: it must print nothing. Without the ``PERL5LIB`` line
it ends with a segmentation fault on ``osx-arm64``. The segfault does not affect
the Linux or ``osx-64`` packages.

3. Build EVM, PASA and Trinity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run from the cloned repository, with the environment active. Skip
``pixi_install_bowtie2.sh``. The conda bowtie2 for arm64 works, and that script
builds x86 (AVX2) variants.

.. code-block:: none

    bash install_scripts/pixi_install_evm.sh
    bash install_scripts/pixi_install_pasa.sh
    bash install_scripts/pixi_install_trinity.sh

The scripts detect macOS and use the conda clang compilers with OpenMP. The
tested build times on Apple Silicon were about 1 minute for EVM and PASA together
and about 3 minutes for Trinity. EVM prints a warning that three converter
programs (``convert_EVM_outputs_to_GFF3``, ``augustus_to_evm_gff3``,
``gff3_gene_prediction_file_validator``) were not found. The Perl versions are in
``$EVM_HOME/EvmUtils``. The warning did not stop the build.

4. Fix blat
~~~~~~~~~~~

The ``ucsc-blat`` package for ``osx-arm64`` has an invalid code signature. macOS
stops the program (exit code 137), and ``funannotate check`` reports ``blat`` as
not installed. Sign it again, locally:

.. code-block:: none

    codesign -f -s - $CONDA_PREFIX/bin/blat
    blat        # prints the usage text

5. Check the installation
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

    funannotate check --show-versions
    python -m pip install pytest && python -m pytest tests -q

In the tested environment, ``funannotate check`` reported only these items as
missing: ``tbl2asn``, CodingQuarry, ``fasta``, ``gmes_petap.pl``, ``signalp``,
``GENEMARK_PATH`` and ``DB_File``. The ``DB_File`` message is expected until
you confirm that the ``PERL5LIB`` line above fixes it: the ``check`` run in the test
was before ``PERL5LIB`` was tested, and the Trinity self-test
(``$TRINITY_HOME/sample_data/test_Trinity_Assembly/runMe.sh``, with ``TRINITY_HOME``
set) passed only with it. It finished in 57 seconds and made 80 transcripts.

Updating
--------

To move to a newer release, update the tag in the clone, recreate or update
the environment from ``environment.yml``, and install the matching package:

.. code-block:: none

    git fetch --tags && git checkout v1.9.1
    mamba env update -n funannotate -f environment.yml --prune
    python -m pip install --upgrade funannotate==1.9.1
    bash install_scripts/pixi_install_evm.sh   # repeat for pasa and trinity
