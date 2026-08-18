
.. _containers:

Funannotate Containers (rust-optimized fork)
================================

This page documents the container build in *this* repository -- a fork of
funannotate that packages the ``rust_EVM_trinity_PASA`` branch (Rust-optimized
PASA, EVidenceModeler, and Trinity forks) into a self-contained image via
conda-pack. It is separate from the upstream ``nextgenusfs/funannotate``
Docker Hub image documented in :ref:`docker`, which does not include the Rust
components and predates this fork's build process.

Building the image
-------------------

.. code-block:: none

    # release build
    docker build -t funannotate-live -f Dockerfile .

    # iterative dev build (faster rebuilds while developing against this fork)
    docker build -t funannotate-live-dev -f Dockerfile.dev .

For HPC use, convert to Singularity/Apptainer:

.. code-block:: none

    singularity build funannotate-live.sif docker-daemon://funannotate-live:latest

Environment variables
-----------------------

The image bakes in sane defaults for everything funannotate/PASA/EVM need to
find each other, but several variables are *placeholders* you are expected to
override at runtime by bind-mounting your own data onto the same path (rather
than baking any single site's storage layout into the image). ``FUNANNOTATE_DB``
and ``EGGNOG_DATA_DIR`` are both subpaths under a common ``/opt/databases``
parent, so a single ``-B /your/dbs:/opt/databases``-style bind (with your own
``funannotate_db/`` and ``eggnog_db/`` subdirectories laid out to match) covers
both at once, or bind/override each independently:

.. list-table::
   :header-rows: 1

   * - Variable
     - Image default
     - Notes
   * - ``FUNANNOTATE_DB``
     - ``/opt/databases/funannotate_db``
     - Core reference DBs (UniProt/SwissProt, taxonomy, BUSCO lineages,
       etc). Bind your own directory onto this path, or override the
       variable to point elsewhere.
   * - ``EGGNOG_DATA_DIR``
     - ``/opt/databases/eggnog_db``
     - EggNOG-mapper database (~50GB: ``eggnog.db``, ``eggnog_proteins.dmnd``,
       taxonomy files). Deliberately not baked into the image -- same
       rationale as the upstream image (see :ref:`docker`). **See the
       "EggNOG-mapper" section below -- an unmounted default here doesn't
       fail gracefully.**
   * - ``PASAHOME``
     - ``/venv/opt/pasa/src``
     - PASA install root. Not normally overridden.
   * - ``PERL5LIB``
     - ``/venv/opt/pasa/src/SAMPLE_HOOKS:/venv/opt/pasa/src/PerlLib``
     - Required for PASA's hook loader to resolve
       ``GFF3::GFF3_annot_retriever`` during ``funannotate update``'s PASA
       annotation-comparison step. See "Known gotchas" below.
   * - ``EVM_HOME``
     - ``/venv/opt/evm``
     - EVidenceModeler install root. Not normally overridden.

GeneMark-ES/ET is deliberately not included
-----------------------------------------------

Unlike ``FUNANNOTATE_DB``/``EGGNOG_DATA_DIR``, there is no ``GENEMARK_PATH``
placeholder in the image and none is planned. GeneMark is separately
licensed software (see :ref:`docker` for the same caveat on the upstream
image), and integrating it cleanly is more than a bind mount: funannotate
invokes ``gmes_petap.pl`` from a fresh temp run directory rather than from
GeneMark's own install directory, and ``gmes_petap.pl`` resolves its bundled
``LineFit.pm`` via a *relative* ``use lib './lib'`` -- so a plain
``GENEMARK_PATH`` bind mount plus an env var is not sufficient on its own;
correctly wiring it in also requires ``PERL5LIB`` to include GeneMark's own
``lib/`` directory, and depending on your install, other bundled non-CPAN
dependencies besides ``LineFit.pm``.

The recommended pattern -- and the one this project's own downstream
Nextflow pipeline uses -- is to **run GeneMark-ES/ET outside the container
entirely** (via a locally installed, licensed copy) and feed its output GTF
into ``funannotate predict`` via ``--genemark_gtf``, which needs no GeneMark
runtime dependency inside the container at all. If you specifically want
GeneMark integrated *inside* a container/environment rather than run
separately, build your own local conda or pixi environment with GeneMark
installed and configured there (mirroring how this image itself is built),
rather than expecting this image to bundle it.

Running with Singularity on shared HPC storage
-------------------------------------------------

A worked example, binding real host paths 1:1 onto themselves (so relative
paths inside output files stay meaningful) rather than relying on the
image's internal placeholder paths:

.. code-block:: none

    SIF=/path/to/funannotate-live.sif
    FUNANNOTATE_DB=/shared/lib/funannotate_db
    EGGNOG_DB=/shared/lib/eggnog_db   # your own copy, or your site's shared install

    singularity exec \
        --bind ${FUNANNOTATE_DB}:${FUNANNOTATE_DB} \
        --bind ${EGGNOG_DB}:${EGGNOG_DB} \
        --env FUNANNOTATE_DB=${FUNANNOTATE_DB} \
        --env EGGNOG_DATA_DIR=${EGGNOG_DB} \
        ${SIF} funannotate annotate -i /path/to/output --cpus 16 ...

Note that ``singularity exec`` (without ``--cleanenv``) passes the host
shell's environment through to the container by default -- so exporting
``FUNANNOTATE_DB``/``EGGNOG_DATA_DIR`` in the calling shell *before* the
``singularity exec`` call reaches the container even without an explicit
``--env`` flag. The explicit form above is shown for clarity and portability
to plain Docker, where this passthrough does not happen automatically.

Known gotchas
--------------

**EggNOG-mapper crashes instead of skipping cleanly when its DB is missing.**
``funannotate annotate`` auto-runs ``emapper.py`` whenever the binary is on
``$PATH`` (it is, in this image) -- there is no flag to disable this when the
binary is present but its DB isn't. If ``EGGNOG_DATA_DIR`` points at a path
with no real database, ``emapper.py --version`` prints a "database not
found" warning to *stderr* ahead of its actual version line. funannotate's
``get_emapper_version()`` (``annotate.py``) does a naive stderr-prefix regex
match for ``emapper-\\S+``; when that fails it returns ``False``, and the
caller does ``LooseVersion(False) >= LooseVersion("2.1.0")``, which raises
``AttributeError: 'LooseVersion' object has no attribute 'version'`` instead
of a clean "eggnog-mapper not available, skipping" message. **You must bind a
real EggNOG-mapper database onto** ``EGGNOG_DATA_DIR`` **before running
``funannotate annotate``**, or the run will fail outright partway through
functional annotation (after Pfam/Diamond/UniProt steps have already run and
written output). Reproduced against a real run, 2026-08-18.

**PASA's ``funannotate update`` annotation-comparison step needs PERL5LIB
set.** PASA's hook system (``Pasa_conf.pm::_get_hook``, used by
``Load_Current_Gene_Annotations.dbi`` during the ``-A`` annotation-compare
step that only ``funannotate update`` triggers -- ``funannotate train``'s
PASA alignment-assembly step does not hit this code path, which is why this
went uncaught until update was first exercised) resolves
``HOOK_EXISTING_GENE_ANNOTATION_LOADER`` targets like
``GFF3::GFF3_annot_retriever`` by searching Perl's ``@INC`` for
``GFF3/GFF3_annot_retriever.pm``. That file ships in the image under
``$PASAHOME/SAMPLE_HOOKS/GFF3/GFF3_annot_retriever.pm``, but neither
``SAMPLE_HOOKS`` nor ``PerlLib`` is on Perl's default ``@INC`` without
``PERL5LIB`` set -- so without it, every real ``funannotate update`` run
crashes with ``Error, couldn't resolve path for GFF3::GFF3_annot_retriever``.
Fixed in the image by baking the ``PERL5LIB`` value above into ``ENV``
(2026-08-18); if you're running an older image tag built before this fix,
export ``PERL5LIB`` yourself before invoking ``funannotate update``.

**A crashed ``funannotate update`` can permanently break a later
``funannotate annotate`` retry on the same output directory.**
``funannotate annotate``'s input-directory detection is unconditional:
``if os.path.isdir(update_results): use it`` (``annotate.py``), with no
check that the directory actually contains anything and no fallback to
``predict_results`` when it doesn't. ``funannotate update`` creates
``update_results/`` near the start of a run, before it writes anything into
it -- so if ``update`` crashes (e.g. the PASA hook issue above, or any other
mid-run failure), it leaves behind an *empty* ``update_results/`` directory
that permanently poisons every subsequent ``annotate`` attempt on that
output dir with ``Properly formatted 'funannotate predict' files do no
exist in this directory``, even though ``predict_results/`` is completely
intact. If you hit this, ``rmdir`` (or ``rm -rf``, if it has other partial
content) the empty/partial ``update_results/`` directory before retrying.

**A host-side ``which`` shell function can leak into the container and look
like a broken image.** Some ``bash -c '...'`` invocations under this image
report a ``which: command not found``-style failure for tools that are
clearly installed. This is not an image defect -- it's Rocky Linux's
``BASH_FUNC_which%%`` bash function leaking in from the *host* environment
via Singularity's default env passthrough (the same passthrough this page
recommends using for ``FUNANNOTATE_DB``/``EGGNOG_DATA_DIR`` above). It does
not affect funannotate/PASA's own Perl ``system()``/backtick calls, which
spawn ``dash`` and never import ``BASH_FUNC_*%%``. Defensively
``unset -f which`` and ``unset which_declare`` before invoking funannotate
commands from a wrapper shell if you see this.

See also
---------

* :ref:`docker` -- the legacy upstream ``nextgenusfs/funannotate`` Docker Hub
  image workflow (no Rust components, no ``PERL5LIB``/``EGGNOG_DATA_DIR``
  fixes above -- predates this fork).
