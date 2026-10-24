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

``environment.yml`` pins ``python>=3.10,<3.12`` (Biopython 1.87 and later, which carry the Bio.Entrez XXE fix, need Python 3.10) and sets the channel order to
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
    # keep the user's PERL5LIB so that deactivating restores it (and so that repeated
    # activation does not grow the path)
    export FUNANNOTATE_PERL5LIB_BACKUP="${PERL5LIB:-}"
    export PERL5LIB=$PASAHOME/SAMPLE_HOOKS:$PASAHOME/PerlLib${PERL5LIB:+:$PERL5LIB}
    export FUNANNOTATE_EVM_ENGINE=rust
    # fasta36 is installed under that name; funannotate calls it "fasta"
    [ -e $CONDA_PREFIX/bin/fasta ] || ln -s $CONDA_PREFIX/bin/fasta36 $CONDA_PREFIX/bin/fasta
    EOF

    cat > $CONDA_PREFIX/etc/conda/deactivate.d/funannotate.sh << 'EOF'
    if [ -n "${FUNANNOTATE_PERL5LIB_BACKUP:-}" ]; then
        export PERL5LIB="$FUNANNOTATE_PERL5LIB_BACKUP"
    else
        unset PERL5LIB
    fi
    unset FUNANNOTATE_PERL5LIB_BACKUP TRINITY_HOME TRINITYHOME PASAHOME EVM_HOME FUNANNOTATE_EVM_ENGINE
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
annotation. The steps below were run on one Apple Silicon Mac (macOS 26.7,
arm64) with a native arm64 conda (Miniforge) and Python 3.9. They covered:

* the environment build and ``funannotate check``,
* the unit tests (``pytest tests``),
* the EVM, PASA and Trinity builds and the Trinity self-test,
* ``funannotate test -t annotate``, which builds a GenBank file with NCBI's
  arm64 ``table2asn``,
* ``funannotate test -t predict``, which made 1,690 genes (the test expects
  1,500 to 1,800) with the Rust EVM engine, with GeneMark and CodingQuarry absent,
* the RNA-seq pipeline of ``funannotate test -t rna-seq``: ``train`` (Trinity and PASA)
  completed, ``predict`` made 1,696 genes and ``update`` made 1,707 (the test expects
  1,630 to 1,830 for each). ``train`` was resumed from its Trinity output after the
  ``pblat`` fix in step 5, so this was not one uninterrupted run.

Intel Macs were not tested. ``compare`` was not run on macOS.

**Installer script (Apple Silicon).** Each release has an ``install-funannotate-macos.sh``
asset (and a ``.sha256``) that builds the supported environment on your Mac. It installs
``pixi`` if needed, clones the release tag, runs ``pixi install --locked``, builds EVM, PASA and
Trinity, applies the macOS fixes below and runs the smoke checks. It takes about 30 minutes and
needs Homebrew (for the GCC that builds Augustus) and ``git``. CI runs the same script on an
Apple Silicon runner before attaching it to the release. Nothing is copied from another machine,
so there are no relocation or signing problems::

    gh release download v1.9.1 -R nextgenusfs/funannotate -p 'install-funannotate-macos.sh*'
    shasum -a 256 -c install-funannotate-macos.sh.sha256
    bash install-funannotate-macos.sh           # clones into ./funannotate-1.9.1
    cd funannotate-1.9.1 && pixi shell          # sets PASAHOME, EVM_HOME, TRINITY_HOME, ...
    funannotate check --show-versions

Use ``pixi shell`` (or ``pixi run``), not just ``PATH``: the environment variables the pipeline
needs are set on activation. You still need ``FUNANNOTATE_DB`` (``funannotate setup``). To do the
steps by hand, use ``pixi install --locked`` and ``pixi run install-externals``, or the conda
steps below.

**Use a native arm64 conda on Apple Silicon.** Check with
``conda info | grep platform``. It must say ``osx-arm64``. An x86-64 conda (it
runs under Rosetta) gives an Intel environment, and the Trinity build fails in
``bamsifter`` there: that Makefile calls the system ``g++`` (arm64) and links
against an x86-64 ``htslib``.

**What is not available on macOS arm64**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Missing
     - Effect
   * - CodingQuarry
     - The CodingQuarry predictor in ``funannotate predict`` is not available.
   * - GeneMark-ES/ET
     - Licensed software, as on Linux. ``funannotate predict`` skips GeneMark
       with a warning. Run it elsewhere and pass ``--genemark_gtf`` if you need it.
   * - SignalP
     - Licensed software, as on Linux.
   * - RepeatModeler
     - It needs ``perl-devel-size``, which has no arm64 build. RepeatMasker is
       installed.
   * - legacy ``tbl2asn``
     - bioconda has no arm64 build. Use NCBI's ``table2asn`` (step 7).
   * - R and Bioconductor packages for Trinity differential expression
     - Not requested by ``environment-macos.yml``.

1. Prerequisites and environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install the Xcode command line tools (``xcode-select --install``) and Homebrew,
then Homebrew's GCC (``brew install gcc``). Step 4 needs GCC. Then clone the
repository and use the macOS environment file, not ``environment.yml``:

.. code-block:: none

    git clone --branch v1.9.0 https://github.com/nextgenusfs/funannotate.git
    cd funannotate
    mamba env create -f environment-macos.yml
    mamba activate funannotate

``environment-macos.yml`` leaves out the packages above. It uses the same Python
range as Linux (3.10 or 3.11). Keep it at 3.9 or newer because the bioconda ``kmer-jellyfish`` 2.x package for
``osx-arm64`` exists only for Python 3.9 and newer, and Trinity needs
jellyfish 2. If you change the Python pin, conda may install jellyfish 1.1.12
and Trinity stops with ``need jellyfish v2``. It includes ``fasta3``, which
installs ``fasta36``. funannotate and PASA accept that name.

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
    export FUNANNOTATE_PERL5LIB_BACKUP="${PERL5LIB:-}"
    export PERL5LIB=$CONDA_PREFIX/lib/perl5/5.32/core_perl:$PASAHOME/SAMPLE_HOOKS:$PASAHOME/PerlLib${PERL5LIB:+:$PERL5LIB}
    EOF

    cat > $CONDA_PREFIX/etc/conda/deactivate.d/funannotate.sh << 'EOF'
    if [ -n "${FUNANNOTATE_PERL5LIB_BACKUP:-}" ]; then
        export PERL5LIB="$FUNANNOTATE_PERL5LIB_BACKUP"
    else
        unset PERL5LIB
    fi
    unset FUNANNOTATE_PERL5LIB_BACKUP TRINITY_HOME TRINITYHOME PASAHOME EVM_HOME FUNANNOTATE_EVM_ENGINE
    EOF

    mamba deactivate && mamba activate funannotate

Set ``FUNANNOTATE_DB`` in your shell startup file (step 8).

The openjdk package sets ``JAVA_HOME`` itself. It has no ``java`` in
``$CONDA_PREFIX/bin``, and ``which java`` shows ``/usr/bin/java``, the macOS
launcher, which uses ``JAVA_HOME``. ``java -version`` should print OpenJDK 25.

The ``perl5/5.32`` directory name belongs to the Perl version in the
environment. If ``ls $CONDA_PREFIX/lib/perl5`` shows a different version, change
the path. The segfault itself is ``perl -e 'use DB_File'`` without the
``PERL5LIB`` line. The Linux and ``osx-64`` packages are not affected.

3. Build EVM, PASA and Trinity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run from the cloned repository, with the environment active. Skip
``pixi_install_bowtie2.sh``: the conda bowtie2 for arm64 works, and the script
needs x86 SIMD headers.

.. code-block:: none

    bash install_scripts/pixi_install_evm.sh
    bash install_scripts/pixi_install_pasa.sh
    bash install_scripts/pixi_install_trinity.sh

The scripts build pinned release tags (EVM ``v2.5.0_rust``, PASA ``v2.6.1-rc.4``,
Trinity ``v2.16.2_rust``). They use a checkout in ``../EVidenceModeler`` or
``../PASApipeline`` if one exists, and print its commit and warn when it is not
the pinned tag. A stale checkout there gives an old build. On macOS the Trinity
script uses the conda clang compilers with OpenMP. Build times on Apple Silicon:
under one minute each for EVM and PASA, and about 3 minutes for Trinity.

On macOS PASA's C++ assembler (``pasa``) needs an ``rpath`` to ``$CONDA_PREFIX/lib``;
PASA ``v2.6.1-rc.4`` adds it in its Makefile, and the PASA script still checks and adds it for
builds from older tags. Without it the binary builds without error and then stops at the
first assembly with ``Library not loaded: @rpath/libomp.dylib``. The script does this
automatically, and also repairs a PASA build made earlier when you run it again.

4. Build Augustus with GCC
~~~~~~~~~~~~~~~~~~~~~~~~~~

``augustus --proteinprofile`` crashes (exit 139) in every Augustus built with
Apple clang, including the bioconda packages. ``funannotate predict`` runs that
test before BUSCO training, so it stops there. The cause is in Augustus's
``src/pp_scoring.cc``: a ``std::map`` iterator is decremented past ``begin()``,
which libstdc++ (GCC) tolerates and libc++ (clang) does not. Build Augustus with
GCC and replace the conda binary:

.. code-block:: none

    bash install_scripts/install_augustus_gcc_macos.sh

The script builds Augustus v3.5.0 with Homebrew GCC (about one minute), links the
GCC runtime statically, keeps the conda binary as ``$CONDA_PREFIX/bin/augustus.clang``,
and ends by running funannotate's own ``--proteinprofile`` test. Linux is not
affected.

5. Run pblat single-threaded
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bioconda's ``pblat`` 2.5.1 crashes with a bus error whenever it uses two or more
threads on Apple Silicon. The older arm64 build and builds from source with GCC and
clang crash too; the Intel build under Rosetta works. PASA runs ``pblat -threads=$CPU``
in its blat alignment step, so ``funannotate train`` and ``update`` stop there. One
thread gives byte-identical output to the Intel build and takes seconds for this step.
Install a wrapper that drops ``-threads`` and keeps the real binary as ``pblat.real``:

.. code-block:: none

    bash install_scripts/install_pblat_shim_macos.sh

Run it again after ``conda install`` or an update replaces ``pblat``. It does nothing on
Linux or on Intel Macs.

6. Fix blat
~~~~~~~~~~~

The ``ucsc-blat`` package for ``osx-arm64`` has an invalid code signature. macOS
stops the program (exit code 137), and ``funannotate check`` reports ``blat`` as
not installed. Sign it again, locally:

.. code-block:: none

    codesign -f -s - $CONDA_PREFIX/bin/blat
    blat        # prints the usage text

7. Install NCBI table2asn
~~~~~~~~~~~~~~~~~~~~~~~~~

bioconda has no ``table2asn`` for arm64, but NCBI's own macOS build is native
arm64. funannotate prefers ``table2asn`` and falls back to ``tbl2asn``
(``$FUN_TBL2ASN`` overrides the choice).

.. code-block:: none

    curl -fsSL -o table2asn.gz https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/mac.table2asn.gz
    gunzip table2asn.gz && chmod +x table2asn
    mv table2asn $CONDA_PREFIX/bin/
    table2asn -version

``curl`` adds no quarantine flag. A browser download does; clear it with
``xattr -d com.apple.quarantine table2asn``. The version prints as
``table2asn: 1.29.324``, with a note that the copy is more than a year old; this
is informational.

8. Databases
~~~~~~~~~~~~

Set ``FUNANNOTATE_DB`` in your shell startup file, then install the databases.
Do not use ``-i all``. ``-b`` only applies together with ``-i busco``.

.. code-block:: none

    export FUNANNOTATE_DB=$HOME/db/funannotate
    funannotate setup -d $FUNANNOTATE_DB -i repeats busco merops dbCAN go gene2product uniprot pfam interpro -b dikarya

This took about 20 minutes and 6.1 GB on disk, mostly the Pfam files.

9. Check the installation
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none

    funannotate check --show-versions
    python -m pip install pytest && python -m pytest tests -q
    funannotate test -t annotate --cpus 4
    funannotate test -t predict --cpus 4
    funannotate test -t rna-seq --cpus 4

The predict test took about 25 minutes on 4 CPUs. The RNA-seq test runs ``train``,
``predict`` and ``update`` and takes about an hour; Trinity's genome-guided assembly is
about half of it. It needs steps 3, 4 and 8. Without
the Augustus build in step 4 it stops at ``augustus --proteinprofile test failed``.

``funannotate check`` is expected to report only CodingQuarry, ``gmes_petap.pl``,
``GENEMARK_PATH``, ``signalp`` and ``DB_File`` as missing. The ``DB_File`` message
remains because the broken conda package is still installed; the ``PERL5LIB`` line
makes Trinity work anyway (the Trinity self-test,
``$TRINITY_HOME/sample_data/test_Trinity_Assembly/runMe.sh``, finished in 57
seconds with 80 transcripts).

``table2asn`` in the annotate test may log an error with exit code 2. The test data
use an organism name that is not in NCBI's taxonomy ("Genome one"). With a real
species name ``table2asn`` returns 0.

Updating
--------

To move to a newer release, update the tag in the clone, recreate or update
the environment from ``environment.yml``, and install the matching package:

.. code-block:: none

    git fetch --tags && git checkout v1.9.1
    mamba env update -n funannotate -f environment.yml --prune
    python -m pip install --upgrade funannotate==1.9.1
    bash install_scripts/pixi_install_evm.sh   # repeat for pasa and trinity
