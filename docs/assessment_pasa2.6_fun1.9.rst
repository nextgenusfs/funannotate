.. _assessment_pasa2.6_fun1.9:

Training and evidence assessment (v1.9.0)
=========================================

The assessment of gene-prediction training and transcript evidence in **funannotate v1.9.0** with **PASApipeline v2.6.1**,
measured against RefSeq annotations in September 2026, has moved out of this repository, so that funannotate holds code and its
unit tests only. The results, figures (PNG and PDF), data tables, analysis scripts and the methods draft are in the
`funannotate-manuscript <https://github.com/hyphaltip/funannotate-manuscript>`_ repository, in ``assessment_pasa2.6_fun1.9/``.

The behaviour that these experiments set is described in :doc:`training`; the tests that pin that behaviour are
``tests/test_training_source_decision.py`` and ``tests/test_augustus_parameters.py`` in this repository.
