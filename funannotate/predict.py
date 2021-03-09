#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import sys
import os
import subprocess
import shutil
import argparse
import json
import datetime
import funannotate.library as lib
from natsort import natsorted


def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(
        prog='funannotate-predict.py', usage="%(prog)s [options] -i genome.fasta",
        description='''Script that does it all.''',
        epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
        formatter_class=MyFormatter)
    parser.add_argument('-i', '--input', help='Genome in FASTA format')
    parser.add_argument('-o', '--out', required=True,
                        help='Basename of output files')
    parser.add_argument('-s', '--species', required=True,
                        help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
    parser.add_argument('-w', '--weights', nargs='+',
                        help='Gene predictors and weights')
    parser.add_argument('-p', '--parameters',
                        help='Training parameters JSON file')
    parser.add_argument('--isolate', help='Isolate name (e.g. Af293)')
    parser.add_argument('--strain', help='Strain name (e.g. CEA10)')
    parser.add_argument('--header_length', default=16,
                        type=int, help='Max length for fasta headers')
    parser.add_argument('--name', default="FUN_",
                        help='Shortname for genes, perhaps assigned by NCBI, eg. VC83')
    parser.add_argument('--numbering', default=1, type=int,
                        help='Specify start of gene numbering')
    parser.add_argument('--augustus_species', '--trained_species', '--species_parameters',
                        dest='augustus_species', help='Specify species for Augustus')
    parser.add_argument('--genemark_mod',
                        help='Use pre-existing Genemark training file (e.g. gmhmm.mod)')
    parser.add_argument('--protein_evidence', nargs='+',
                        help='Specify protein evidence (multiple files can be separaed by a space)')
    parser.add_argument('--protein_alignments', dest='exonerate_proteins',
                        help='Pre-computed Exonerate protein alignments (see README for how to run exonerate)')
    parser.add_argument('--transcript_evidence', nargs='+',
                        help='Transcript evidence (map to genome with minimap2)')
    parser.add_argument('--transcript_alignments',
                        help='Transcript evidence in GFF3 format')
    parser.add_argument('-gm', '--genemark_mode', default='ES',
                        choices=['ES', 'ET'], help='Mode to run genemark in')
    parser.add_argument('--pasa_gff',
                        help='Pre-computed PASA/TransDecoder high quality models')
    parser.add_argument('--other_gff', nargs='+',
                        help='GFF gene prediction pass-through to EVM')
    parser.add_argument('--augustus_gff',
                        help='Pre-computed Augustus gene models (GFF3)')
    parser.add_argument('--genemark_gtf',
                        help='Pre-computed GeneMark gene models (GTF)')
    parser.add_argument('--soft_mask', type=int, default=2000,
                        help='Threshold used in GeneMark for use of softmasked regions')
    parser.add_argument('--maker_gff', help='MAKER2 GFF output')
    parser.add_argument('--repeats2evm', action='store_true',
                        help='Pass repeat GFF3 to EVM')
    parser.add_argument('--repeat_filter', default=['overlap', 'blast'], nargs='+', choices=[
                        'overlap', 'blast', 'none'], help='Repeat filters to apply')
    parser.add_argument('--rna_bam',
                        help='BAM (sorted) of RNAseq aligned to reference')
    parser.add_argument('--stringtie', help='StringTie GTF')
    parser.add_argument('--min_intronlen', default=10,
                        type=int, help='Minimum intron length for gene models')
    parser.add_argument('--max_intronlen', default=3000,
                        type=int, help='Maximum intron length for gene models')
    parser.add_argument('--min_protlen', default=50, type=int,
                        help='Minimum amino acid length for valid gene model')
    parser.add_argument('--keep_no_stops', action='store_true',
                        help='Keep gene models without valid stop codons')
    parser.add_argument('--ploidy', default=1, type=int,
                        help='Ploidy of assembly')
    parser.add_argument('--cpus', default=2, type=int,
                        help='Number of CPUs to use')
    parser.add_argument('--busco_seed_species', default='anidulans',
                        help='Augustus species to use as initial training point for BUSCO')
    parser.add_argument('--optimize_augustus', action='store_true',
                        help='Run "long" training of Augustus')
    parser.add_argument('--force', action='store_true',
                        help='Annotated if genome not masked and skip bad contigs')
    parser.add_argument('--busco_db', default='dikarya',
                        help='BUSCO model database')
    parser.add_argument('-t', '--tbl2asn', default='-l paired-ends',
                        help='Parameters for tbl2asn, linkage and gap info')
    parser.add_argument('--organism', default='fungus',
                        choices=['fungus', 'other'],
                        help='Fungal specific settings')
    parser.add_argument('--SeqCenter', default='CFMR',
                        help='Sequencing center for GenBank tbl file')
    parser.add_argument('--SeqAccession', default='12345',
                        help='Sequencing accession number')
    parser.add_argument('-d', '--database',
                        help='Path to funannotate database, $FUNANNOTATE_DB')
    parser.add_argument('--keep_evm', action='store_true',
                        help='dont rerun EVM')
    parser.add_argument('--no-evm-partitions', action='store_false',
                        dest='evm_partitions',
                        help='do not split contigs in EVM')
    parser.add_argument('--evm-partition-interval', default=1500,
                        dest='evm_interval',
                        help='min space between genes to use for partition')
    parser.add_argument('--aligners', default=['minimap2'], nargs='+', choices=[
                        'minimap2', 'gmap', 'blat'], help='transcript alignment programs')
    parser.add_argument('--EVM_HOME',
                        help='Path to Evidence Modeler home directory, $EVM_HOME')
    parser.add_argument('--AUGUSTUS_CONFIG_PATH',
                        help='Path to Augustus config directory, $AUGUSTUS_CONFIG_PATH')
    parser.add_argument('--GENEMARK_PATH',
                        help='Path to GeneMark exe (gmes_petap.pl) directory, $GENEMARK_PATH')
    parser.add_argument('--min_training_models', default=200, type=int,
                        help='Minimum number of BUSCO or BUSCO_EVM gene models to train Augustus')
    parser.add_argument('--p2g_pident', default=80, help='Exonerate pct identity')
    parser.add_argument('--p2g_diamond_db', help='Premade diamond genome database')
    parser.add_argument('--p2g_prefilter', default='diamond', choices=['diamond', 'tblastn'],
                        help='Option for p2g on which prefilter')
    parser.add_argument('--no-progress', dest='progress', action='store_false',
                        help='no progress on multiprocessing')
    parser.add_argument('--trnascan',
                        help='Pre-computed tRNAScan results')
    args = parser.parse_args(args)

    parentdir = os.path.join(os.path.dirname(__file__))
    TODAY = datetime.date.today().strftime("%Y-%m-%d")

    # create folder structure
    if not os.path.isdir(args.out):
        os.makedirs(args.out)
        os.makedirs(os.path.join(args.out, 'predict_misc'))
        os.makedirs(os.path.join(args.out, 'predict_results'))
        os.makedirs(os.path.join(args.out, 'logfiles'))
    else:
        if os.path.isdir(os.path.join(args.out, 'predict_results')):
            shutil.rmtree(os.path.join(args.out, 'predict_results'))
            os.makedirs(os.path.join(args.out, 'predict_results'))
        # make sure subdirectories exist
        dirs = [os.path.join(args.out, 'predict_misc'), os.path.join(
            args.out, 'logfiles'), os.path.join(args.out, 'predict_results')]
        for d in dirs:
            if not os.path.isdir(d):
                os.makedirs(d)

    # create log file
    log_name = os.path.join(args.out, 'logfiles', 'funannotate-predict.log')
    if os.path.isfile(log_name):
        os.remove(log_name)

    # initialize script, log system info and cmd issue at runtime
    lib.setupLogging(log_name)
    FNULL = open(os.devnull, 'w')
    cmd_args = " ".join(sys.argv)+'\n'
    lib.log.debug(cmd_args)
    sys.stderr.write(
        "-------------------------------------------------------\n")
    lib.SystemInfo()

    # get version of funannotate
    version = lib.get_version()
    lib.log.info("Running funannotate v{:}".format(version))

    # check for conflicting folder names to avoid problems
    conflict = ['busco', 'busco_proteins', 'genemark', 'EVM_tmp']
    if args.out in conflict:
        lib.log.error(
            "%s output folder conflicts with a hard coded tmp folder, please change -o parameter" % args.out)
        sys.exit(1)

    # setup funannotate DB path
    if args.database:
        FUNDB = args.database.strip()
    else:
        try:
            FUNDB = os.environ["FUNANNOTATE_DB"].strip()
        except KeyError:
            lib.log.error(
                'Funannotate database not properly configured, run funannotate setup.')
            sys.exit(1)

    # check if database setup
    blastdb = os.path.join(FUNDB, 'repeats.dmnd')
    if not os.path.isfile(blastdb):
        lib.log.error("Can't find Repeat Database at {:}, you may need to re-run funannotate setup".format(
            os.path.join(FUNDB, 'repeats.dmnd')))
        sys.exit(1)
    # check buscos, download if necessary
    if not os.path.isdir(os.path.join(FUNDB, args.busco_db)):
        lib.log.error("ERROR: %s busco database is not found, install with funannotate setup -b %s" %
                      (args.busco_db, args.busco_db))
        sys.exit(1)

    # do some checks and balances
    if args.EVM_HOME:
        EVM = args.EVM_HOME.strip()
    else:
        try:
            EVM = os.environ["EVM_HOME"].strip()
        except KeyError:
            lib.log.error(
                "$EVM_HOME environmental variable not found, Evidence Modeler is not properly configured.  You can use the --EVM_HOME argument to specifiy a path at runtime")
            sys.exit(1)

    if args.AUGUSTUS_CONFIG_PATH:
        AUGUSTUS = args.AUGUSTUS_CONFIG_PATH.strip()
    else:
        try:
            AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"].strip()
        except KeyError:
            lib.log.error("$AUGUSTUS_CONFIG_PATH environmental variable not found, Augustus is not properly configured. You can use the --AUGUSTUS_CONFIG_PATH argument to specify a path at runtime.")
            sys.exit(1)

    # if you want to use BRAKER1, you also need some additional config paths
    if args.GENEMARK_PATH:
        GENEMARK_PATH = args.GENEMARK_PATH.strip()
    else:
        try:
            GENEMARK_PATH = os.environ["GENEMARK_PATH"].strip()
        except KeyError:
            gmes_path = which_path('gmes_petap.pl')
            if not gmes_path:
                lib.log.error(
                    "GeneMark not found and $GENEMARK_PATH environmental variable missing. Will skip GeneMark ab-initio prediction.")
                GENEMARK_PATH = False
            else:
                GENEMARK_PATH = os.path.dirname(gmes_path)

    # check for some Augustus scripts
    scripts_missing = []
    if os.path.basename(os.path.normcase(os.path.abspath(AUGUSTUS))) == 'config':
        AUGUSTUS_BASE = os.path.dirname(os.path.abspath(AUGUSTUS))
    if lib.which('bam2hints'):
        BAM2HINTS = 'bam2hints'
    else:
        BAM2HINTS = os.path.join(AUGUSTUS_BASE, 'bin', 'bam2hints')
        if not lib.which(BAM2HINTS):
            scripts_missing.append(BAM2HINTS)
    if lib.which_path('join_mult_hints.pl'):
        JOINHINTS = 'join_mult_hints.pl'
    else:
        JOINHINTS = os.path.join(AUGUSTUS_BASE, 'scripts', 'join_mult_hints.pl')
        if not lib.which_path(JOINHINTS):
            scripts_missing.append(JOINHINTS)
    if lib.which('gff2gbSmallDNA.pl'):
        GFF2GB = 'gff2gbSmallDNA.pl'
    else:
        GFF2GB = os.path.join(AUGUSTUS_BASE, 'scripts', 'gff2gbSmallDNA.pl')
        if not lib.which(GFF2GB):
            scripts_missing.append(GFF2GB)
    if len(scripts_missing) > 0:
        lib.log.error('ERROR: unable to locate the following Augustus scripts:\n{}'.format('\n'.join(scripts_missing)))
        lib.log.error('Either add these scripts to your PATH or ensure that $AUGUSTUS_CONFIG_PATH/../scripts/ exists.')
        sys.exit(1)
    # get genemark scripts and determine if installed
    GeneMark2GFF = os.path.join(parentdir, 'aux_scripts', 'genemark_gtf2gff3.pl')
    if GENEMARK_PATH:
        try:
            GENEMARKCMD = os.path.join(GENEMARK_PATH, 'gmes_petap.pl')
            lib.log.debug('GeneMark path: {:}'.format(GENEMARK_PATH))
            genemarkcheck = lib.which(GENEMARKCMD)
            lib.log.debug('Full path to gmes_petap.pl: {:}'.format(GENEMARKCMD))
            lib.log.debug('GeneMark appears to be functional? {:}'.format(genemarkcheck))
        except NameError:
            GENEMARKCMD = ''
            genemarkcheck = False
    else:
        GENEMARKCMD = ''
        genemarkcheck = False

    # setup dictionary to store weights
    # default=['genemark:1', 'pasa:6', 'codingquarry:2', 'snap:1', 'glimmerhmm:1']
    StartWeights = {'augustus': 1, 'hiq': 2, 'genemark': 1, 'pasa': 6,
                    'codingquarry': 0, 'snap': 1, 'glimmerhmm': 1,
                    'proteins': 1, 'transcripts': 1}
    if not genemarkcheck:
        StartWeights['genemark'] = 0

    EVMBase = {'augustus': 'ABINITIO_PREDICTION',
               'genemark': 'ABINITIO_PREDICTION',
               'snap': 'ABINITIO_PREDICTION',
               'glimmerhmm': 'ABINITIO_PREDICTION',
               'codingquarry': 'OTHER_PREDICTION',
               'pasa': 'OTHER_PREDICTION',
               'hiq': 'OTHER_PREDICTION',
               'proteins': 'PROTEIN',
               'transcripts': 'TRANSCRIPT'}


    programs = ['exonerate', 'diamond', 'tbl2asn', 'bedtools',
                'augustus', 'etraining', 'tRNAscan-SE', BAM2HINTS]
    programs = programs + args.aligners
    if 'blat' in args.aligners:
        programs.append('pslCDnaFilter')
    lib.CheckDependencies(programs)

    # check if diamond version matches database version
    if not lib.CheckDiamondDB(blastdb):
        lib.log.error(
            'Diamond repeat database was created with different version of diamond, please re-run funannotate setup')
        sys.exit(1)

    # check that variables are correct, i.e. EVM should point to correct folder
    if not os.path.isfile(os.path.join(EVM, 'EvmUtils', 'partition_EVM_inputs.pl')):
        lib.log.error(
            'EvidenceModeler $EVM_HOME variable is not correct\nEVM scripts not found in $EVM_HOME: {:}'.format(EVM))
        sys.exit(1)

    # look for pre-existing data in training folder
    # look for pre-existing training data to use
    pre_existing = []
    if os.path.isdir(os.path.join(args.out, 'training')):
        traindir = os.path.join(args.out, 'training')
        if os.path.isfile(os.path.join(traindir, 'funannotate_train.coordSorted.bam')):
            if not args.rna_bam:
                args.rna_bam = os.path.join(
                    traindir, 'funannotate_train.coordSorted.bam')
                pre_existing.append(
                    '  --rna_bam '+os.path.join(traindir, 'funannotate_train.coordSorted.bam'))
        if os.path.isfile(os.path.join(traindir, 'funannotate_train.pasa.gff3')):
            if not args.pasa_gff:
                args.pasa_gff = os.path.join(
                    traindir, 'funannotate_train.pasa.gff3')
                pre_existing.append(
                    '  --pasa_gff '+os.path.join(traindir, 'funannotate_train.pasa.gff3'))
        if os.path.isfile(os.path.join(traindir, 'funannotate_train.stringtie.gtf')):
            if not args.stringtie:
                args.stringtie = os.path.join(
                    traindir, 'funannotate_train.stringtie.gtf')
                pre_existing.append(
                    '  --stringtie '+os.path.join(traindir, 'funannotate_train.stringtie.gtf'))
        if os.path.isfile(os.path.join(traindir, 'funannotate_train.transcripts.gff3')):
            if not args.transcript_alignments:
                args.transcript_alignments = os.path.join(
                    traindir, 'funannotate_train.transcripts.gff3')
                pre_existing.append(
                    '  --transcript_alignments '+os.path.join(traindir, 'funannotate_train.transcripts.gff3'))
        else:
            if os.path.isfile(os.path.join(traindir, 'funannotate_train.trinity-GG.fasta')):
                if not args.transcript_evidence:
                    args.transcript_evidence = [os.path.join(
                        traindir, 'funannotate_train.trinity-GG.fasta')]
                    pre_existing.append(
                        '  --transcript_evidence '+os.path.join(traindir, 'funannotate_train.trinity-GG.fasta'))
                else:  # maybe passed a different one? then append to the list
                    if not os.path.join(traindir, 'funannotate_train.trinity-GG.fasta') in args.transcript_evidence:
                        args.transcript_evidence.append(os.path.join(
                            traindir, 'funannotate_train.trinity-GG.fasta'))
                        pre_existing.append(
                            '  --transcript_evidence '+' '.join(args.transcript_evidence))
            if os.path.isfile(os.path.join(traindir, 'funannotate_long-reads.fasta')):
                if not args.transcript_evidence:
                    args.transcript_evidence = [os.path.join(
                        traindir, 'funannotate_long-reads.fasta')]
                    pre_existing.append(
                        '  --transcript_evidence '+os.path.join(traindir, 'funannotate_long-reads.fasta'))
                else:  # maybe passed a different one? then append to the list
                    if not os.path.join(traindir, 'funannotate_long-reads.fasta') in args.transcript_evidence:
                        args.transcript_evidence.append(os.path.join(
                            traindir, 'funannotate_long-reads.fasta'))
                        pre_existing.append(
                            '  --transcript_evidence '+' '.join(args.transcript_evidence))

    if len(pre_existing) > 0:
        lib.log.info("Found training files, will re-use these files:\n%s" %
                     '\n'.join(pre_existing))

    # see if organism/species/isolate was passed at command line, build PASA naming scheme
    organism = None
    if args.species:
        organism = args.species
    else:
        organism = os.path.basename(args.input).split('.fa')[0]
    if args.strain:
        organism_name = organism+'_'+args.strain
    elif args.isolate:
        organism_name = organism+'_'+args.isolate
    else:
        organism_name = organism
    organism_name = organism_name.replace(' ', '_')

    # check augustus species now, so that you don't get through script and then find out it is already in DB
    LOCALPARAMETERS = os.path.join(
        args.out, 'predict_misc', 'ab_initio_parameters')
    LOCALAUGUSTUS = os.path.join(LOCALPARAMETERS, 'augustus')
    lib.copyDirectory(os.path.join(FUNDB, 'trained_species', args.busco_seed_species), os.path.join(
        LOCALAUGUSTUS, 'species', args.busco_seed_species), overwrite=True)
    lib.copyDirectory(os.path.join(FUNDB, 'trained_species', 'anidulans'), os.path.join(
        LOCALAUGUSTUS, 'species', 'anidulans'), overwrite=True)
    lib.copyDirectory(os.path.join(AUGUSTUS, 'species', 'generic'), os.path.join(
        LOCALAUGUSTUS, 'species', 'generic'), overwrite=True)
    lib.copyDirectory(os.path.join(AUGUSTUS, 'extrinsic'), os.path.join(
        LOCALAUGUSTUS, 'extrinsic'), overwrite=True)
    lib.copyDirectory(os.path.join(AUGUSTUS, 'model'), os.path.join(
        LOCALAUGUSTUS, 'model'), overwrite=True)
    lib.copyDirectory(os.path.join(AUGUSTUS, 'profile'), os.path.join(
        LOCALAUGUSTUS, 'profile'), overwrite=True)
    if not args.augustus_species:
        aug_species = organism_name.lower()
    else:
        aug_species = args.augustus_species

    # copy the necessary config files to local dir to help/prevent permissions issues
    # --AUGUSTUS_CONFIG_PATH=LOCALAUGUSTUS
    if args.parameters:
        augspeciescheck = True
        lib.log.info(
            'Ab initio training parameters file passed: {:}'.format(args.parameters))
        with open(args.parameters) as infile:
            trainingData = json.load(infile)
    else:
        augspeciescheck = lib.CheckFunannotateSpecies(aug_species, FUNDB)
        if augspeciescheck:
            # load in the trainingData
            with open(os.path.join(FUNDB, 'trained_species', aug_species, 'info.json')) as infile:
                trainingData = json.load(infile)
        else:
            trainingData = {
                'augustus': [{}],
                'genemark': [{}],
                'codingquarry': [{}],
                'snap': [{}],
                'glimmerhmm': [{}],
            }

    # if fungus and RNA-bam then make codingquarry run
    if args.rna_bam and args.organism == 'fungus':
        StartWeights['codingquarry'] = 2

    # parse input programs/weights then cross ref with what is installed
    # respect user input here, ie codingquary:0 should turn it off
    if args.weights:
        for x in args.weights:
            if ':' in x:
                predictor, weight = x.split(':')
                weight = int(weight)
            else:
                predictor = x
                weight = 1
            if not predictor.lower() in StartWeights:
                lib.log.info(
                    '{:} is unknown source for evidence modeler'.format(predictor))
            else:
                StartWeights[predictor.lower()] = weight

    lib.log.debug(StartWeights)

    # move files around appropriately
    RunModes = {}
    if 'path' in trainingData['augustus'][0]:
        RunModes['augustus'] = 'pretrained'
        lib.copyDirectory(trainingData['augustus'][0]['path'], os.path.join(
            LOCALAUGUSTUS, 'species', aug_species), overwrite=True)
        for f in os.listdir(os.path.join(LOCALAUGUSTUS, 'species', aug_species)):
            if not f.startswith(aug_species):
                ff = os.path.join(os.path.join(
                    LOCALAUGUSTUS, 'species', aug_species, f))
                if ff.endswith('_exon_probs.pbl'):
                    shutil.copyfile(ff, os.path.join(os.path.join(
                        LOCALAUGUSTUS, 'species', aug_species, aug_species+'_exon_probs.pbl')))
                elif ff.endswith('_igenic_probs.pbl'):
                    shutil.copyfile(ff, os.path.join(os.path.join(
                        LOCALAUGUSTUS, 'species', aug_species, aug_species+'_igenic_probs.pbl')))
                elif ff.endswith('_intron_probs.pbl'):
                    shutil.copyfile(ff, os.path.join(os.path.join(
                        LOCALAUGUSTUS, 'species', aug_species, aug_species+'_intron_probs.pbl')))
                elif ff.endswith('_metapars.cfg'):
                    shutil.copyfile(ff, os.path.join(os.path.join(
                        LOCALAUGUSTUS, 'species', aug_species, aug_species+'_metapars.cfg')))
                elif ff.endswith('_parameters.cfg'):
                    shutil.copyfile(ff, os.path.join(os.path.join(
                        LOCALAUGUSTUS, 'species', aug_species, aug_species+'_parameters.cfg')))
                elif ff.endswith('_weightmatrix.txt'):
                    shutil.copyfile(ff, os.path.join(os.path.join(
                        LOCALAUGUSTUS, 'species', aug_species, aug_species+'_weightmatrix.txt')))
                elif ff.endswith('_metapars.cgp.cfg'):
                    shutil.copyfile(ff, os.path.join(os.path.join(
                        LOCALAUGUSTUS, 'species', aug_species, aug_species+'_metapars.cgp.cfg')))
                elif ff.endswith('_metapars.utr.cfg'):
                    shutil.copyfile(ff, os.path.join(os.path.join(
                        LOCALAUGUSTUS, 'species', aug_species, aug_species+'_metapars.utr.cfg')))
    else:
        if args.pasa_gff:
            RunModes['augustus'] = 'pasa'
            sourceText = 'PASA: {:}'.format(os.path.abspath(args.pasa_gff))
        else:
            RunModes['augustus'] = 'busco'
            # need to check that args.busoco_seed_species
            if not os.path.isdir(os.path.join(LOCALAUGUSTUS, 'species', args.busco_seed_species)):
                lib.log.error('ERROR: --busco_seed_species {} is not valid as it is not in database. Try `funannotate species` command and check spelling.')
                sys.exit(1)
            sourceText = 'BUCSCO {:}'.format(args.busco_db)
        trainingData['augustus'] = [{'version': 'funannotate v{:}'.format(version), 'source': sourceText,
                                     'date': TODAY, 'path': os.path.abspath(os.path.join(LOCALAUGUSTUS, 'species', aug_species))}]
    if genemarkcheck:
        if 'path' in trainingData['genemark'][0]:
            RunModes['genemark'] = 'pretrained'
            shutil.copyfile(trainingData['genemark'][0]['path'], os.path.join(
                LOCALPARAMETERS, aug_species+'.genemark.mod'))
            args.genemark_mod = os.path.join(
                LOCALPARAMETERS, aug_species+'.genemark.mod')
        else:
            RunModes['genemark'] = 'selftraining'
            trainingData['genemark'] = [{'version': 'funannotate v{:}'.format(version), 'source': RunModes['genemark']+' ' + args.genemark_mode,
                                         'date': TODAY, 'path': os.path.abspath(os.path.join(LOCALPARAMETERS, aug_species+'.genemark.mod'))}]
    if 'path' in trainingData['snap'][0]:
        RunModes['snap'] = 'pretrained'
        shutil.copyfile(trainingData['snap'][0]['path'], os.path.join(
            LOCALPARAMETERS, aug_species+'.snap.hmm'))
    else:
        if args.pasa_gff:
            RunModes['snap'] = 'pasa'
            sourceText = 'PASA: {:}'.format(os.path.abspath(args.pasa_gff))
        else:
            RunModes['snap'] = 'busco'
            sourceText = 'BUCSCO {:}'.format(args.busco_db)
        trainingData['snap'] = [{'version': 'funannotate v{:}'.format(version), 'source': sourceText,
                                 'date': TODAY, 'path': os.path.abspath(os.path.join(LOCALPARAMETERS, aug_species+'.snap.hmm'))}]
    if 'path' in trainingData['glimmerhmm'][0]:
        RunModes['glimmerhmm'] = 'pretrained'
        lib.copyDirectory(trainingData['glimmerhmm'][0]['path'], os.path.join(
            LOCALPARAMETERS, 'glimmerhmm'), overwrite=True)
    else:
        if args.pasa_gff:
            RunModes['glimmerhmm'] = 'pasa'
            sourceText = 'PASA: {:}'.format(os.path.abspath(args.pasa_gff))
        else:
            RunModes['glimmerhmm'] = 'busco'
            sourceText = 'BUCSCO {:}'.format(args.busco_db)
        trainingData['glimmerhmm'] = [{'version': 'funannotate v{:}'.format(version), 'source': sourceText,
                                       'date': TODAY, 'path': os.path.abspath(os.path.join(LOCALPARAMETERS, 'glimmerhmm'))}]
    if StartWeights['codingquarry'] > 0:
        if 'path' in trainingData['codingquarry'][0] and 'QUARRY_PATH' in os.environ:
            RunModes['codingquarry'] = 'pretrained'
            lib.copyDirectory(os.environ['QUARRY_PATH'], os.path.join(
                args.out, 'predict_misc', 'CodingQuarry', 'QuarryFiles'), overwrite=True)
            lib.copyDirectory(trainingData['codingquarry'][0]['path'], os.path.join(
                args.out, 'predict_misc', 'CodingQuarry', 'QuarryFiles', 'species', aug_species), overwrite=True)
        else:
            if args.rna_bam and 'QUARRY_PATH' in os.environ:
                lib.copyDirectory(os.environ['QUARRY_PATH'], os.path.join(
                    args.out, 'predict_misc', 'CodingQuarry', 'QuarryFiles'), overwrite=True)
                RunModes['codingquarry'] = 'rna-bam'
                sourceText = 'BAM: {:}'.format(os.path.abspath(args.rna_bam))
                trainingData['codingquarry'] = [{'version': 'funannotate v{:}'.format(version), 'source': sourceText,
                                                 'date': TODAY, 'path': os.path.abspath(os.path.join(LOCALPARAMETERS, 'codingquarry'))}]
    elif args.rna_bam and 'QUARRY_PATH' in os.environ:
        lib.log.info('Skipping CodingQuarry as --organism=other. Pass a weight larger than 0 to run CQ, ie --weights codingquarry:1')
    elif not args.rna_bam and 'QUARRY_PATH' in os.environ:
        lib.log.info('Skipping CodingQuarry as no --rna_bam passed')
    else:
        lib.log.info('Skipping CodingQuarry as $QUARRY_PATH not found as ENV')

    # let user know what will be run and from what source
    lib.log.debug(RunModes)
    RunBusco = False
    lib.log.info(
        'Parsed training data, run ab-initio gene predictors as follows:')
    AbInitio = [['Program', 'Training-Method']]
    for k, v in natsorted(list(RunModes.items())):
        AbInitio.append([k, v])
        if 'busco' == v:
            RunBusco = True
    abinitio_table = lib.print_table(AbInitio, return_str=True)
    sys.stderr.write(abinitio_table)
    if 'QUARRY_PATH' in os.environ and not 'codingquarry' in RunModes and StartWeights['codingquarry'] > 0:
        lib.log.info(
            'CodingQuarry will be skipped --> --rna_bam required for training')
        StartWeights['codingquarry'] = 0

    # check augustus functionality
    augustus_version, augustus_functional = lib.checkAugustusFunc()
    system_os = lib.systemOS()
    if not augustus_functional and RunBusco:
        lib.log.error(
            'ERROR: augustus --proteinprofile test failed, likely a compilation error. This is required to run BUSCO, exiting.')
        if 'MacOSX' in system_os:
            lib.log.error(
                'OS={:}, try to install augustus v3.2.1 from https://github.com/nextgenusfs/augustus or from HomeBrew'.format(system_os))
        sys.exit(1)
    # check bam2hints, which often also not compiled correctly, but only need if passing rna_bam
    if args.rna_bam:
        bam2hints_functional = False
        for line in lib.execute([BAM2HINTS, '-h']):
            if line.startswith('bam2hints --'):
                bam2hints_functional = True
        if not bam2hints_functional:
            lib.log.error(
                "ERROR: {:} is not installed properly check compilation".format(BAM2HINTS))
            sys.exit(1)

    # check input files to make sure they are not empty, first check if multiple files passed to transcript/protein evidence
    input_checks = [args.input, args.genemark_mod,
                    args.exonerate_proteins, args.pasa_gff, args.rna_bam]
    if not args.protein_evidence:
        args.protein_evidence = [os.path.join(FUNDB, 'uniprot_sprot.fasta')]
    input_checks = input_checks + args.protein_evidence
    if args.transcript_evidence:  # if transcripts passed, otherwise ignore
        input_checks = input_checks + args.transcript_evidence
    if args.other_gff:
        input_checks = input_checks + args.other_gff
    # now check the inputs
    for i in input_checks:
        if i:
            if ':' in i:
                i = i.split(':')[0]
            lib.checkinputs(i)

    # convert PASA GFF and/or GFF pass-through
    # convert PASA to have 'pasa' in second column to make sure weights work with EVM
    PASA_GFF = os.path.join(args.out, 'predict_misc', 'pasa_predictions.gff3')
    if args.pasa_gff:
        if ':' in args.pasa_gff:
            args.pasa_gff, PASA_weight = args.pasa_gff.split(':')
            StartWeights['pasa'] = int(PASA_weight)
        lib.renameGFF(os.path.abspath(args.pasa_gff), 'pasa', PASA_GFF)
        # validate it will work with EVM
        if not lib.evmGFFvalidate(PASA_GFF, EVM, lib.log):
            lib.log.error(
                "ERROR: %s is not a properly formatted PASA GFF file, please consult EvidenceModeler docs" % args.pasa_gff)
            sys.exit(1)

    # parse and convert other GFF files for pass through to EVM
    OTHER_GFFs = []
    other_weights = []
    other_files = []
    other_contigs = set()
    if args.other_gff:
        if any(':' in s for s in args.other_gff):
            for x in args.other_gff:
                if ':' in x:
                    other_weights.append(x.split(':')[-1])
                    other_files.append(x.split(':')[0])
                else:
                    other_weights.append('1')
                    other_files.append(x)
        else:
            other_weights = ['1', ]*len(args.other_gff)
            other_files = args.other_gff

    if len(other_files) > 0:
        for i, file in enumerate(other_files):
            featurename = 'other_pred'+str(i+1)
            lib.log.info('Parsing GFF pass-through: {:} --> setting source to {:}'.format(file, featurename))
            outputGFF = os.path.join(args.out, 'predict_misc', 'other'+str(i+1)+'_predictions.gff3')
            contig_names = lib.renameGFF(os.path.abspath(file), featurename, outputGFF)
            other_contigs.update(contig_names)
            # validate output with EVM
            if not lib.evmGFFvalidate(outputGFF, EVM, lib.log):
                lib.log.error("ERROR: {} is not proper GFF file, please consult EvidenceModeler docs".format(file))
                sys.exit(1)
            OTHER_GFFs.append(outputGFF)
            if not featurename in StartWeights:
                StartWeights[featurename] = other_weights[i]
    lib.log.debug(StartWeights)

    # setup the genome fasta file, need either args.input or need to have args.masked_genome + args.repeatmasker_gff3
    # declare output location
    MaskGenome = os.path.join(args.out, 'predict_misc', 'genome.softmasked.fa')
    RepeatMasker = os.path.join(args.out, 'predict_misc', 'repeatmasker.bed')
    AssemblyGaps = os.path.join(args.out, 'predict_misc', 'assembly-gaps.bed')
    Scaffoldsort = os.path.join(
        args.out, 'predict_misc', 'scaffold.sort.order.txt')
    Renamingsort = os.path.join(
        args.out, 'predict_misc', 'scaffold.sort.rename.txt')
    # check inputs
    if args.input:
        # check fasta header length and fasta alphabet/characters
        #header_test = lib.checkFastaHeaders(args.input, args.header_length)
        bad_headers, bad_contigs, suspect_contigs = lib.analyzeAssembly(
            args.input, header_max=args.header_length)
        if len(bad_headers) > 0:
            lib.log.error("Genome assembly error: headers contain more characters than the max ({}), reformat headers to continue.".format(
                args.header_length))
            lib.log.error("First {:} headers that failed names:\n{}".format(len(bad_headers[:5]),
                                                                           '\n'.join(bad_headers[:5])))
            sys.exit(1)
        elif len(bad_contigs) > 0:
            lib.log.error('Found {:,} contigs contain non-IUPAC characters:'.format(len(bad_contigs)))
            for k, v in natsorted(bad_contigs.items()):
                print(k)
                for x in v:
                    print('  {}\t{}'.format(x[0], x[1]))
                lib.log.debug('{} {}'.format(k, v))
            sys.exit(1)
        elif len(suspect_contigs) > 0 and not args.force:
            lib.log.error('Found {:,} bad contigs, where alphabet is less than 4 [this should not happen]'.format(
                len(suspect_contigs)))
            for k, v in natsorted(suspect_contigs.items()):
                lib.log.debug('{} {}'.format(k, v))
                print(k)
                total = 0
                for nuc, num in natsorted(v.items()):
                    print('  {:}: {:,}'.format(nuc, num))
                    total += int(num)
                print('len: {:,}'.format(total))
                print('-----------------------')
            lib.log.info('If you really want to keep and annotate these contigs (not recommended), pass --force')
            sys.exit(1)
        else:
            # not sure why this is in the code...
            '''
            with open(Scaffoldsort, 'w') as contigsout:
                sortedHeaders = natsorted(header_test[1])
                contigsout.write('%s' % '\n'.join(sortedHeaders))
            with open(Renamingsort, 'w') as renameout:
                counter = 0
                with open(Scaffoldsort, 'r') as contigsin:
                    for line in contigsin:
                        counter += 1
                        line = line.replace('\n', '')
                        renameout.write('%s\t%i\n' % (line, counter))
            '''

        # if BAM file passed, check if headers are same as input
        if args.rna_bam:
            if not lib.BamHeaderTest(args.input, args.rna_bam):
                lib.log.error(
                    "Fasta headers in BAM file do not match genome, exiting.")
                sys.exit(1)
        # check that the genome is soft-masked
        lib.log.info(
            'Loading genome assembly and parsing soft-masked repetitive sequences')
        ContigSizes, GenomeLength, maskedSize, percentMask = lib.checkMasklowMem(
            args.input, RepeatMasker, AssemblyGaps, args.cpus)
        if maskedSize == 0 and not args.force:
            lib.log.error(
                'Error: Genome is not repeat-masked, to ignore use --force. Or soft-mask using `funannotate mask` command or suitable external program.')
            sys.exit(1)
        else:
            lib.log.info('Genome loaded: {:,} scaffolds; {:,} bp; {:.2%} repeats masked'.format(
                len(ContigSizes), GenomeLength, percentMask))

        # if other_gff passed, check contigs
        if len(other_contigs) > 0:
            extra_contigs = []
            for contig in other_contigs:
                if contig not in ContigSizes:
                    extra_contigs.append(contig)
            if len(extra_contigs) > 0:
                lib.log.error('ERROR: found {:,} contigs in --other_gff that are not found in the genome assembly.'.format(len(extra_contigs)))
                lib.log.error('  --other_gff should contain gene models from. Sample of offending contigs: {}\n  {}'.format(args.input, ', '.join(extra_contigs[:10])))
                sys.exit(1)

        # just copy the input fasta to the misc folder and move on.
        shutil.copyfile(args.input, MaskGenome)
    else:
        lib.log.error('Error: Please provide a genome file, -i or --input')
        sys.exit(1)

    # setup augustus parallel command
    AUGUSTUS_PARALELL = os.path.join(
        parentdir, 'aux_scripts', 'augustus_parallel.py')

    # EVM command line scripts
    Converter = os.path.join(EVM, 'EvmUtils', 'misc',
                             'augustus_GFF3_to_EVM_GFF3.pl')
    Converter2 = os.path.join(EVM, 'EvmUtils', 'misc',
                              'augustus_GTF_to_EVM_GFF3.pl')
    EVM2proteins = os.path.join(EVM, 'EvmUtils', 'gff3_file_to_proteins.pl')

    # make sure absolute path
    RepeatMasker = os.path.abspath(RepeatMasker)
    MaskGenome = os.path.abspath(MaskGenome)

    # final output for augustus hints, declare ahead of time for checking portion of script
    hintsE = os.path.join(args.out, 'predict_misc', 'hints.E.gff')
    hintsP = os.path.join(args.out, 'predict_misc', 'hints.P.gff')
    hintsBAM = os.path.join(args.out, 'predict_misc', 'hints.BAM.gff')
    hints_all = os.path.join(args.out, 'predict_misc', 'hints.ALL.gff')
    hintsM = os.path.join(args.out, 'predict_misc', 'hints.M.gff')

    # check longest 10 contigs
    longest10 = natsorted(list(ContigSizes.values()), reverse=True)[:10]

    # check for previous files and setup output files
    Predictions = os.path.join(
        args.out, 'predict_misc', 'gene_predictions.gff3')
    Exonerate = os.path.join(args.out, 'predict_misc',
                             'protein_alignments.gff3')
    Transcripts = os.path.join(
        args.out, 'predict_misc', 'transcript_alignments.gff3')
    Weights = os.path.join(args.out, 'predict_misc', 'weights.evm.txt')
    EVM_out = os.path.join(args.out, 'predict_misc', 'evm.round1.gff3')
    EVMWeights = {}  # this will store the ones that are actually used

    # if maker_gff passed, use that info and move on, if pasa present than run EVM.
    if args.maker_gff:
        lib.log.info("Parsing Maker2 GFF for use in EVidence Modeler")
        lib.maker2evm(args.maker_gff, os.path.join(args.out, 'predict_misc'))
        # append PASA data if exists
        if args.pasa_gff:
            with open(Predictions, 'a') as output:
                with open(PASA_GFF) as input:
                    output.write(input.read())
        if OTHER_GFFs:
            for y in OTHER_GFFs:
                with open(Predictions, 'a') as output:
                    with open(y) as input:
                        output.write(input.read())
        # setup weights file for EVM
        with open(Weights, 'w') as output:
            genesources = []
            with open(Predictions, 'r') as preds:
                for line in preds:
                    if line.startswith('\n'):
                        continue
                    source = line.split('\t')[1]
                    if not source in genesources:
                        genesources.append(source)
            if not genesources:
                lib.log.error(
                    "Maker2 GFF not parsed correctly, no gene models found, exiting.")
                sys.exit(1)
            for i in genesources:
                if i == 'maker':
                    output.write("ABINITIO_PREDICTION\t{:}\t1\n".format(i))
                    if not 'maker' in EVMWeights:
                        EVMWeights['MAKER'] = '1'
                elif i == 'pasa':
                    if 'pasa' in StartWeights:
                        output.write("OTHER_PREDICTION\t{:}\t{:}\n".format(
                            i, StartWeights.get('pasa')))
                    else:
                        output.write(
                            "OTHER_PREDICTION\t{:}\t{:}\n".format(i, 6))
                        EVMWeights['pasa'] = '6'
                elif i.startswith('other_pred'):
                    output.write("OTHER_PREDICTION\t{:}\t{:}\n".format(
                        i, StartWeights.get(i)))
                    if not i in EVMWeights:
                        EVMWeights[i] = str(StartWeights.get(i))
                else:
                    output.write("OTHER_PREDICTION\t{:}\t1\n".format(i))
                    if not i in EVMWeights:
                        EVMWeights[i] = '1'
            tr_sources = []
            with open(Transcripts, 'r') as trns:
                for line in trns:
                    source = line.split('\t')[1]
                    if source not in tr_sources:
                        tr_sources.append(source)
            for i in tr_sources:
                output.write("TRANSCRIPT\t{:}\t{:}\n".format(
                    i, StartWeights.get('transcripts')))
                EVMWeights['transcripts'] = str(
                    StartWeights.get('transcripts'))
            output.write("PROTEIN\tprotein2genome\t{:}\n".format(
                StartWeights.get('proteins')))
            EVMWeights['proteins'] = str(StartWeights.get('proteins'))
        Exonerate = os.path.abspath(Exonerate)
        Transcripts = os.path.abspath(Transcripts)
    else:
        # no maker_gff, so let funannotate handle gene prediction
        # check for transcript evidence/format as needed
        trans_out = os.path.join(args.out, 'predict_misc', 'transcript_alignments.gff3')
        trans_temp = os.path.join(args.out, 'predict_misc', 'transcripts.combined.fa')
        minimapGFF3 = os.path.join(args.out, 'predict_misc', 'transcript_minimap2.gff3')
        gmapGFF3 = os.path.join(args.out, 'predict_misc', 'transcript_gmap.gff3')
        blat_out = os.path.join(args.out, 'predict_misc', 'blat.psl')
        blat_filt = os.path.join(args.out, 'predict_misc', 'blat.filt.psl')
        blat_sort1 = os.path.join(args.out, 'predict_misc', 'blat.sort.tmp.psl')
        blat_sort2 = os.path.join(args.out, 'predict_misc', 'blat.sort.psl')
        maxINT = '-maxIntron='+str(args.max_intronlen)
        b2h_input = '--in='+blat_sort2
        b2h_output = '--out='+hintsE
        FinalTrainingModels = os.path.join(args.out, 'predict_misc', 'final_training_models.gff3')
        if args.transcript_alignments:
            lib.harmonize_transcripts(
                MaskGenome, args.transcript_alignments, trans_out,
                hintsM, evidence=args.transcript_evidence,
                tmpdir=os.path.join(args.out, 'predict_misc'), cpus=args.cpus,
                maxintron=args.max_intronlen)
        if not lib.checkannotations(trans_out):
            # combine transcript evidence into a single file
            if args.transcript_evidence:
                if os.path.isfile(trans_temp):
                    lib.SafeRemove(trans_temp)
                with open(trans_temp, 'w') as output:
                    for f in args.transcript_evidence:
                        with open(f) as input:
                            output.write(input.read())
                if 'minimap2' in args.aligners:
                    minimapBAM = os.path.join(args.out, 'predict_misc',
                                              'transcripts.minimap2.bam')
                    if not lib.checkannotations(minimapGFF3) or not lib.checkannotations(hintsM):
                        lib.log.info("Aligning transcript evidence to genome with minimap2")
                        lib.minimap2Align(trans_temp, MaskGenome, args.cpus,
                                          args.max_intronlen, minimapBAM)
                        minimapCount = lib.bam2ExonsHints(minimapBAM,
                                                          minimapGFF3, hintsM)
                        lib.log.info(
                            "Found {:,} alignments, wrote GFF3 and Augustus hints to file".format(minimapCount))
                    else:
                        lib.log.info('Existing minimap2 alignments found: {:} and {:}'.format(
                            minimapGFF3, hintsM))
                if 'gmap' in args.aligners:
                    # run Gmap of transcripts to genome
                    if not lib.checkannotations(gmapGFF3):
                        lib.log.info(
                            "Aligning transcript evidence to genome with GMAP")
                        lib.runGMAP(trans_temp, MaskGenome, args.cpus,
                                    args.max_intronlen,
                                    os.path.join(args.out, 'predict_misc'),
                                    gmapGFF3)
                        gmapCount = lib.countGMAPtranscripts(gmapGFF3)
                        lib.log.info("Found {:,} alignments, wrote GFF3 to file".format(gmapCount))
                    else:
                        lib.log.info(
                            'Existing gmap alignments found: {:}'.format(gmapGFF3))
                if 'blat' in args.aligners:
                    if not lib.checkannotations(hintsE):
                        # now run BLAT for Augustus hints
                        lib.log.info(
                            "Aligning transcript evidence to genome with BLAT")
                        cmd = ['blat', '-noHead', '-minIdentity=80',
                               maxINT, MaskGenome, trans_temp, blat_out]
                        lib.runSubprocess(cmd, '.', lib.log)
                        cmd = ['pslCDnaFilter', '-minId=0.9', '-localNearBest=0.005',
                               '-ignoreNs', '-bestOverlap', blat_out, blat_filt]
                        lib.runSubprocess(cmd, '.', lib.log)
                        cmd = ['sort', '-n', '-k', '16,16', blat_filt]
                        lib.runSubprocess2(cmd, '.', lib.log, blat_sort1)
                        cmd = ['sort', '-s', '-k', '14,14', blat_sort1]
                        lib.runSubprocess2(cmd, '.', lib.log, blat_sort2)
                        # run blat2hints
                        if lib.which('blat2hints.pl'):
                            blat2hints = 'blat2hints.pl'
                        else:
                            blat2hints = os.path.join(
                                AUGUSTUS_BASE, 'scripts', 'blat2hints.pl')
                        cmd = [blat2hints, b2h_input, b2h_output,
                               '--minintronlen=20', '--trunkSS']
                        lib.runSubprocess(cmd, '.', lib.log)
                        total = lib.line_count(blat_sort2)
                        lib.log.info('{0:,}'.format(total) +
                                     ' filtered BLAT alignments')
                    else:
                        lib.log.info(
                            'Existing blat hintsfile found {:}'.format(hintsE))

                # combine transcripts for EVM (need to process GMAP ones here)
                if lib.checkannotations(minimapGFF3) and lib.checkannotations(gmapGFF3):
                    # write function to rename/gmap and combine with minimap data
                    lib.combineTranscripts(minimapGFF3, gmapGFF3, trans_out)
                elif lib.checkannotations(minimapGFF3):
                    shutil.copyfile(minimapGFF3, trans_out)
                elif lib.checkannotations(gmapGFF3):
                    lib.combineTranscripts(False, gmapGFF3, trans_out)
                Transcripts = os.path.abspath(trans_out)
            else:
                Transcripts = False
        else:
            if not args.transcript_alignments:
                lib.log.info('Existing transcript alignments found: {:}'.format(trans_out))
            Transcripts = os.path.abspath(trans_out)
        # check if BAM file passed, if so run bam2hints
        if args.rna_bam:
            if not lib.checkannotations(hintsBAM):
                lib.log.info(
                    "Extracting hints from RNA-seq BAM file using bam2hints")
                bamhintstmp = os.path.join(
                    args.out, 'predict_misc', 'bam_hints.tmp')
                cmd = [BAM2HINTS, '--intronsonly', '--in',
                       args.rna_bam, '--out', bamhintstmp]
                lib.runSubprocess(cmd, '.', lib.log)
                # sort the hints
                bamhintssorted = os.path.join(
                    args.out, 'predict_misc', 'bam_hints.sorted.tmp')
                lib.sortHints(bamhintstmp, bamhintssorted)
                # join hints
                bamjoinedhints = os.path.join(
                    args.out, 'predict_misc', 'bam_hints.joined.tmp')
                cmd = [JOINHINTS]
                lib.runSubprocess5(cmd, '.', lib.log,
                                   bamhintssorted, bamjoinedhints)
                # filter intron hints
                cmd = [os.path.join(
                    parentdir, 'aux_scripts', 'filterIntronsFindStrand.pl'), MaskGenome, bamjoinedhints, '--score']
                lib.runSubprocess2(cmd, '.', lib.log, hintsBAM)
            else:
                lib.log.info(
                    "Existing RNA-seq BAM hints found: {:}".format(hintsBAM))

        # check for protein evidence/format as needed
        Exonerate = os.path.join(
            args.out, 'predict_misc', 'protein_alignments.gff3')
        prot_temp = os.path.join(
            args.out, 'predict_misc', 'proteins.combined.fa')
        P2G = os.path.join(parentdir, 'aux_scripts', 'funannotate-p2g.py')
        # this is alignments variable name is confusing for historical reasons...
        if not args.exonerate_proteins:
            if args.protein_evidence:
                if lib.checkannotations(prot_temp):
                    lib.SafeRemove(prot_temp)
                # clean up headers, etc
                lib.cleanProteins(args.protein_evidence, prot_temp)
                # run funannotate-p2g to map to genome
                p2g_cmd = [sys.executable, P2G, '-p', prot_temp, '-g', MaskGenome,
                           '-o', Exonerate, '--maxintron', str(args.max_intronlen),
                           '--cpus', str(args.cpus),
                           '--exonerate_pident', str(args.p2g_pident),
                           '--ploidy', str(args.ploidy),
                           '-f', args.p2g_prefilter,
                           '--tblastn_out', os.path.join(args.out, 'predict_misc', 'p2g.diamond.out'),
                           '--logfile', os.path.join(args.out, 'logfiles', 'funannotate-p2g.log')]
                if args.p2g_diamond_db:
                    p2g_cmd += ['-d', args.p2g_diamond_db]
                if not args.progress:
                    p2g_cmd.append('--no-progress')
                # check if protein evidence is same as old evidence
                if not lib.checkannotations(Exonerate):
                    # lib.log.info("Mapping proteins to genome using Diamond blastx/Exonerate")
                    subprocess.call(p2g_cmd)
                else:
                    lib.log.info(
                        "Existing protein alignments found: {:}".format(Exonerate))
                Exonerate = os.path.abspath(Exonerate)
            else:
                Exonerate = False
        else:
            lib.log.info("Loading protein alignments {:}".format(
                args.exonerate_proteins))
            shutil.copyfile(args.exonerate_proteins, Exonerate)
            Exonerate = os.path.abspath(Exonerate)

        # generate Augustus hints file from protein_alignments
        if Exonerate:
            lib.exonerate2hints(Exonerate, hintsP)

        # combine hints for Augustus
        allhintstmp = os.path.join(args.out, 'predict_misc', 'hints.all.tmp')
        if lib.checkannotations(hintsP) or lib.checkannotations(hintsE) or lib.checkannotations(hintsBAM) or lib.checkannotations(hintsM):
            if lib.checkannotations(allhintstmp):
                os.remove(allhintstmp)
            with open(allhintstmp, 'a') as out:
                if lib.checkannotations(hintsP):
                    with open(hintsP) as input:
                        out.write(input.read())
                if lib.checkannotations(hintsE):
                    with open(hintsE) as input2:
                        out.write(input2.read())
                if lib.checkannotations(hintsBAM):
                    with open(hintsBAM) as input3:
                        out.write(input3.read())
                if lib.checkannotations(hintsM):
                    with open(hintsM) as input4:
                        out.write(input4.read())
        # now sort hints file, and join multiple hints_all
        allhintstmp_sort = os.path.join(
            args.out, 'predict_misc', 'hints.all.sort.tmp')
        lib.sortHints(allhintstmp, allhintstmp_sort)
        cmd = [JOINHINTS]
        lib.runSubprocess5(cmd, '.', lib.log, allhintstmp_sort, hints_all)

        Augustus, GeneMark = (None,)*2

        # Walk thru data available and determine best approach.
        if args.genemark_gtf:
            # convert the predictors to EVM format and merge
            # convert GeneMark
            GeneMarkGFF3 = os.path.join(
                args.out, 'predict_misc', 'genemark.gff')
            cmd = [GeneMark2GFF, args.genemark_gtf]
            lib.runSubprocess2(cmd, '.', lib.log, GeneMarkGFF3)
            GeneMarkTemp = os.path.join(
                args.out, 'predict_misc', 'genemark.temp.gff')
            cmd = ['perl', Converter, GeneMarkGFF3]
            lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp)
            GeneMark = os.path.join(
                args.out, 'predict_misc', 'genemark.evm.gff3')
            with open(GeneMark, 'w') as output:
                with open(GeneMarkTemp, 'r') as input:
                    lines = input.read().replace("Augustus", "GeneMark")
                    output.write(lines)

        ##############
        #  GeneMark  #
        ##############
        if not GeneMark and 'genemark' in RunModes and StartWeights['genemark'] > 0:
            GeneMarkGFF3 = os.path.join(
                args.out, 'predict_misc', 'genemark.gff')
            # count contigs
            num_contigs = lib.countfasta(MaskGenome)
            if longest10[0] < 50000:
                lib.log.error("GeneMark-ES may fail because this assembly appears to be highly fragmented:\n\
-------------------------------------------------------\n\
The longest %s scaffolds are: %s.\n\
If you can run GeneMark outside funannotate you can add with --genemark_gtf option.\n\
-------------------------------------------------------" % (len(longest10), ', '.join([str(x) for x in longest10])))
            # now run GeneMark, check for number of contigs and ini
            if num_contigs < 2 and not args.genemark_mod:
                lib.log.error(
                    "GeneMark-ES cannot run with only a single contig, you must provide --ini_mod file to run GeneMark")
            elif num_contigs < 2 and args.genemark_mod:
                with open(MaskGenome, 'r') as genome:
                    for line in genome:
                        if line.startswith('>'):
                            header = line.replace('>', '')
                            header = header.replace('\n', '')
                GeneMark = os.path.join(
                    args.out, 'predict_misc', 'genemark.evm.gff3')
                GeneMarkTemp = os.path.join(
                    args.out, 'predict_misc', 'genemark.temp.gff')
                if not os.path.isfile(GeneMarkGFF3):
                    lib.log.info("Running GeneMark on single-contig assembly")
                    cmd = ['gmhmme3', '-m', args.genemark_mod, '-o',
                           GeneMarkGFF3, '-f', 'gff3', MaskGenome]
                    lib.runSubprocess(cmd, '.', lib.log)
                # now open output and reformat
                # lib.log.info("Converting GeneMark GTF file to GFF3")
                with open(GeneMarkTemp, 'w') as geneout:
                    with open(GeneMarkGFF3, 'r') as genein:
                        for line in genein:
                            if not line.startswith('#'):
                                if not '\tIntron\t' in line:
                                    newline = line.replace('seq', header)
                                    newline = newline.replace('.hmm3', '')
                                    geneout.write(newline)
                GeneMarkTemp2 = os.path.join(
                    args.out, 'predict_misc', 'genemark.temp2.gff')
                cmd = ['perl', Converter, GeneMarkTemp]
                lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp2)
                with open(GeneMark, 'w') as output:
                    with open(GeneMarkTemp2, 'r') as input:
                        lines = input.read().replace("Augustus", "GeneMark")
                        output.write(lines)
            else:
                if not lib.checkannotations(GeneMarkGFF3):
                    if args.genemark_mode == 'ES':
                        lib.RunGeneMarkES(GENEMARKCMD, MaskGenome, args.genemark_mod, args.max_intronlen,
                                          args.soft_mask, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, args.organism)
                    else:
                        lib.RunGeneMarkET(GENEMARKCMD, MaskGenome, args.genemark_mod, hints_all, args.max_intronlen,
                                          args.soft_mask, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, args.organism)
                else:
                    lib.log.info(
                        "Existing GeneMark annotation found: {:}".format(GeneMarkGFF3))
                if lib.checkannotations(GeneMarkGFF3):
                    GeneMarkTemp = os.path.join(
                        args.out, 'predict_misc', 'genemark.temp.gff')
                    cmd = ['perl', Converter, GeneMarkGFF3]
                    lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp)
                    GeneMark = os.path.join(
                        args.out, 'predict_misc', 'genemark.evm.gff3')
                    with open(GeneMark, 'w') as output:
                        with open(GeneMarkTemp, 'r') as input:
                            lines = input.read().replace("Augustus", "GeneMark")
                            output.write(lines)
            # GeneMark has occasionally failed internally resulting in incomplete output, check that contig names are okay
            Contigsmissing = []
            if GeneMark:
                os.rename(GeneMark, GeneMark+'.bak')
                with open(GeneMark, 'w') as output:
                    with open(GeneMark+'.bak', 'r') as input:
                        for line in input:
                            if line.startswith('#') or line.startswith('\n'):
                                output.write(line)
                            else:
                                contig = line.split('\t')[0]
                                if not contig in ContigSizes:
                                    Contigsmissing.append(contig)
                                else:
                                    output.write(line)
                Contigsmissing = set(Contigsmissing)
                if len(Contigsmissing) > 0:
                    lib.log.error(
                        "Warning: GeneMark might have failed on at least one contig, double checking results")
                    fileList = []
                    genemark_folder = os.path.join(
                        args.out, 'predict_misc', 'genemark', 'output', 'gmhmm')
                    for file in os.listdir(genemark_folder):
                        if file.endswith('.out'):
                            fileList.append(os.path.join(
                                genemark_folder, file))
                    genemarkGTFtmp = os.path.join(
                        args.out, 'predict_misc', 'genemark', 'genemark.gtf.tmp')
                    genemarkGTF = os.path.join(
                        args.out, 'predict_misc', 'genemark', 'genemark.gtf')
                    lib.SafeRemove(genemarkGTFtmp)
                    lib.SafeRemove(genemarkGTF)
                    for x in fileList:
                        cmd = [os.path.join(GENEMARK_PATH, 'hmm_to_gtf.pl'), '--in',
                               x, '--app', '--out', genemarkGTFtmp, '--min', '300']
                        subprocess.call(cmd)
                    cmd = [os.path.join(GENEMARK_PATH, 'reformat_gff.pl'), '--out', genemarkGTF, '--trace', os.path.join(
                        args.out, 'predict_misc', 'genemark', 'info', 'dna.trace'), '--in', genemarkGTFtmp, '--back']
                    subprocess.call(cmd)
                    # lib.log.info("Converting GeneMark GTF file to GFF3")
                    with open(GeneMarkGFF3, 'w') as out:
                        subprocess.call(
                            [GeneMark2GFF, genemarkGTF], stdout=out)
                    GeneMarkTemp = os.path.join(
                        args.out, 'predict_misc', 'genemark.temp.gff')
                    cmd = ['perl', Converter, GeneMarkGFF3]
                    lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp)
                    GeneMark = os.path.join(
                        args.out, 'predict_misc', 'genemark.evm.gff3')
                    with open(GeneMark, 'w') as output:
                        with open(GeneMarkTemp, 'r') as input:
                            lines = input.read().replace("Augustus", "GeneMark")
                            output.write(lines)
                lib.log.info('{:,} predictions from GeneMark'.format(
                    lib.countGFFgenes(GeneMark)))
                if os.path.isfile(os.path.join(args.out, 'predict_misc', 'gmhmm.mod')):
                    shutil.copyfile(os.path.join(args.out, 'predict_misc', 'gmhmm.mod'), os.path.join(
                        LOCALPARAMETERS, aug_species+'.genemark.mod'))
                else:
                    lib.log.debug(
                        'Genemark mod file not found, removing genemark from training parameters')
                    trainingData['genemark'] = [{}]
        # GeneMark can fail if you try to pass a single contig, check file length
        if StartWeights['genemark'] > 0:
            if genemarkcheck:
                if GeneMark:
                    GM_check = lib.line_count(GeneMark)
                    if GM_check < 3:
                        lib.log.error(
                            "GeneMark predictions failed. If you can run GeneMark outside of funannotate, then pass the results to --genemark_gtf.")
                else:
                    lib.log.error(
                        "GeneMark predictions failed. If you can run GeneMark outside of funannotate, then pass the results to --genemark_gtf.")

        ##############
        #   BUSCO    #
        ##############
        if RunBusco:
            busco_final = os.path.join(
                args.out, 'predict_misc', 'busco.final.gff3')
            busco_log = os.path.join(args.out, 'logfiles', 'busco.log')
            if not lib.checkannotations(busco_final):
                if (sys.version_info > (3, 0)):
                    BUSCO = os.path.join(parentdir,
                                         'aux_scripts', 'funannotate-BUSCO2.py')
                else:
                    BUSCO = os.path.join(parentdir,
                                         'aux_scripts', 'funannotate-BUSCO2-py2.py')
                BUSCO_FUNGI = os.path.join(FUNDB, args.busco_db)
                busco_location = os.path.join(
                    args.out, 'predict_misc', 'busco')
                busco_fulltable = os.path.join(
                    busco_location, 'run_'+aug_species, 'full_table_'+aug_species+'.tsv')
                if lib.CheckAugustusSpecies(args.busco_seed_species):
                    busco_seed = args.busco_seed_species
                else:
                    busco_seed = 'anidulans'
                runbuscocheck = True
                # check if complete run
                if lib.checkannotations(busco_fulltable):
                    lib.log.info(
                        "BUSCO has already been run, using existing data")
                    runbuscocheck = False
                else:
                    if os.path.isdir(busco_location):
                        shutil.rmtree(busco_location)
                if runbuscocheck:
                    lib.log.info(
                        "Running BUSCO to find conserved gene models for training ab-initio predictors")
                    tblastn_version = lib.vers_tblastn()
                    if tblastn_version > '2.2.31':
                        lib.log.info(
                            "Multi-threading in tblastn v{:} is unstable, running in single threaded mode for BUSCO".format(tblastn_version))
                    if not os.path.isdir(busco_location):
                        os.makedirs(busco_location)
                    else:
                        shutil.rmtree(busco_location)  # delete if it is there
                        os.makedirs(busco_location)  # create fresh folder

                    busco_cmd = [sys.executable, BUSCO, '-i', MaskGenome,
                                 '-m', 'genome', '--lineage', BUSCO_FUNGI,
                                 '-o', aug_species, '-c', str(args.cpus),
                                 '--species', busco_seed, '-f',
                                 '--local_augustus', os.path.abspath(LOCALAUGUSTUS)]
                    lib.log.debug(' '.join(busco_cmd))

                    with open(busco_log, 'w') as logfile:
                        subprocess.call(busco_cmd, cwd=busco_location,
                                        stdout=logfile, stderr=logfile)

                    # check if BUSCO found models for training, if not error out and exit.
                    if not lib.checkannotations(busco_fulltable):
                        lib.log.error("BUSCO training of Augusus failed, check busco logs, exiting")
                        sys.exit(1)

                # open output and pull locations to make bed file
                busco_bed = os.path.join(
                    args.out, 'predict_misc', 'buscos.bed')
                busco_complete = lib.parseBUSCO2genome(busco_fulltable,
                                                       args.ploidy,
                                                       ContigSizes,
                                                       busco_bed)

                # proper training files exist, now run EVM on busco models to get high quality predictions.
                lib.log.info('{:,} valid BUSCO predictions found, validating protein sequences'.format(
                    len(busco_complete)))

                # now get BUSCO GFF models
                busco_augustus_tmp = os.path.join(
                    args.out, 'predict_misc', 'busco_augustus.tmp')
                with open(busco_augustus_tmp, 'w') as output:
                    for i in busco_complete:
                        file = os.path.join(busco_location, 'run_'+aug_species, 'augustus_output', 'gffs', i+'.gff')
                        subprocess.call(['perl', Converter2, file],
                                        stderr=FNULL, stdout=output)
                # finally rename models so they are not redundant
                busco_augustus = os.path.join(args.out, 'predict_misc', 'busco_augustus.gff3')
                busco_proteins = os.path.join(args.out, 'predict_misc', 'busco_augustus.proteins.fasta')
                lib.fix_busco_naming(busco_fulltable, MaskGenome,
                                     busco_augustus_tmp, busco_augustus,
                                     ploidy=args.ploidy,
                                     proteins=busco_proteins)
                # double check protein models with busco proteome
                buscoProtDir = os.path.join(args.out, 'predict_misc', 'busco_proteins')
                if not os.path.isdir(buscoProtDir):
                    os.makedirs(buscoProtDir)
                with open(busco_log, 'a') as logfile:
                    subprocess.call([sys.executable, BUSCO,
                                     '-i', os.path.abspath(busco_proteins),
                                     '-m', 'proteins',
                                     '--lineage', BUSCO_FUNGI,
                                     '-o', aug_species,
                                     '--cpu', str(args.cpus),
                                     '--species', busco_seed, '-f',
                                     '--local_augustus', os.path.abspath(LOCALAUGUSTUS)],
                                    cwd=buscoProtDir,
                                    stdout=logfile,
                                    stderr=logfile)
                # rename models
                buscoProtOutput = os.path.join(buscoProtDir, 'run_'+aug_species, 'full_table_'+aug_species+'.tsv')
                buscoProtComplete = lib.getCompleteBuscos(buscoProtOutput,
                                                          ploidy=args.ploidy)
                lib.filterGFF3(buscoProtComplete, MaskGenome, busco_augustus,
                               busco_final)
                total = lib.countGFFgenes(busco_final)
                lib.log.info('{:,} BUSCO predictions validated'.format(total))
            else:
                lib.log.info("Existing BUSCO results found: {:} containing {:,} predictions".format(
                    busco_final, lib.countGFFgenes(busco_final)))
            FinalTrainingModels = busco_final

        ######################
        #     Augustus       #
        ######################
        aug_out = os.path.join(args.out, 'predict_misc', 'augustus.gff3')
        if args.augustus_gff:
            aug_out = args.augustus_gff
        else:
            if not lib.checkannotations(aug_out):
                if args.pasa_gff and RunModes['augustus'] == 'pasa':
                    # check for training data, if no training data, then train using PASA
                    lib.log.info("Filtering PASA data for suitable training set")
                    trainingModels = os.path.join(args.out, 'predict_misc', 'pasa.training.tmp.gtf')
                    # convert PASA GFF to GTF format
                    lib.gff3_to_gtf(PASA_GFF, MaskGenome, trainingModels)
                    # now get best models by cross-ref with intron BAM hints
                    if lib.which('filterGenemark.pl'):
                        FILTERGENE = 'filterGenemark.pl'
                    else:
                        FILTERGENE = os.path.join(parentdir, 'aux_scripts', 'filterGenemark.pl')
                    cmd = [FILTERGENE, os.path.abspath(trainingModels),
                           os.path.abspath(hints_all)]
                    lib.runSubprocess4(cmd,
                                       os.path.join(args.out, 'predict_misc'),
                                       lib.log)
                    totalTrain = lib.selectTrainingModels(PASA_GFF,
                                                          MaskGenome,
                                                          os.path.join(args.out, 'predict_misc', 'pasa.training.tmp.f.good.gtf'),
                                                          FinalTrainingModels)
                elif RunModes['augustus'] == 'busco':
                    totalTrain = lib.countGFFgenes(FinalTrainingModels)

                if RunModes['augustus'] != 'pretrained':
                    # now train augustus
                    if totalTrain < int(args.min_training_models):
                        lib.log.error("Not enough gene models {:} to train Augustus ({:} required), exiting".format(
                            totalTrain, int(args.min_training_models)))
                        sys.exit(1)
                    if totalTrain > 1000:
                        numTrainingSet = round(totalTrain * 0.10)
                    elif totalTrain < 500:
                        numTrainingSet = round(totalTrain * 0.20)
                    else:
                        numTrainingSet = 100
                    lib.log.info("Training Augustus using {:} gene models".format(RunModes['augustus'].upper()))
                    trainingset = os.path.join(args.out, 'predict_misc', 'augustus.training.{:}.gb'.format(RunModes['augustus']))
                    cmd = [GFF2GB, FinalTrainingModels, MaskGenome, '600',
                           trainingset]
                    lib.runSubprocess(cmd, '.', lib.log)
                    lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset,
                                      MaskGenome, args.out, args.cpus,
                                      numTrainingSet, args.optimize_augustus,
                                      os.path.abspath(LOCALAUGUSTUS))

                # now run Augustus multithreaded...
                lib.log.info(
                    "Running Augustus gene prediction using {:} parameters".format(aug_species))
                if os.path.isfile(hints_all):
                    cmd = [AUGUSTUS_PARALELL, '--species', aug_species,
                           '--hints', hints_all,
                           '--AUGUSTUS_CONFIG_PATH', AUGUSTUS,
                           '--local_augustus', LOCALAUGUSTUS,
                           '-i', MaskGenome, '-o', aug_out,
                           '--cpus', str(args.cpus),
                           '-e', os.path.join(parentdir, 'config', 'extrinsic.E.XNT.RM.cfg'),
                           '--logfile', os.path.join(args.out, 'logfiles', 'augustus-parallel.log')]
                else:
                    cmd = [AUGUSTUS_PARALELL, '--species', aug_species,
                           '-i', MaskGenome, '--AUGUSTUS_CONFIG_PATH', AUGUSTUS,
                           '--local_augustus', LOCALAUGUSTUS,
                           '-o', aug_out, '--cpus', str(args.cpus),
                           '-e', os.path.join(parentdir, 'config', 'extrinsic.E.XNT.RM.cfg'),
                           '--logfile', os.path.join(args.out, 'logfiles', 'augustus-parallel.log')]
                if not args.progress:
                    cmd.append('--no-progress')
                subprocess.call(cmd)
            else:
                lib.log.info("Existing Augustus annotations found: {:}".format(aug_out))
        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        cmd = ['perl', Converter, aug_out]
        lib.runSubprocess2(cmd, '.', lib.log, Augustus)

        # make sure Augustus finished successfully
        if not lib.checkannotations(Augustus):
            lib.log.error(
                "Augustus prediction failed, check `logfiles/augustus-parallel.log`")
            sys.exit(1)

        # if hints used for Augustus, get high quality models > 90% coverage to pass to EVM
        if os.path.isfile(hints_all) or args.rna_bam:
            lib.log.info("Pulling out high quality Augustus predictions")
            hiQ_models = []
            with open(aug_out, 'r') as augustus:
                for pred in lib.readBlocks(augustus, '# start gene'):
                    values = []
                    geneID = ''
                    support = ''
                    if pred[0].startswith('# This output'):
                        continue
                    if pred[0].startswith('##gff-version 3'):
                        continue
                    for line in pred:
                        line = line.replace('\n', '')
                        if line.startswith('# start gene'):
                            geneID = line.split(' ')[-1]
                            values.append(geneID)
                        if not args.rna_bam:
                            if line.startswith('# % of transcript supported by hints'):
                                support = line.split(' ')[-1]
                                values.append(support)
                        else:  # if BRAKER is run then only intron CDS evidence is passed, so get models that fullfill that check
                            if line.startswith('# CDS introns:'):
                                intronMatch = line.split(' ')[-1]
                                try:
                                    support = int(intronMatch.split(
                                        '/')[0]) / float(intronMatch.split('/')[1]) * 100
                                except ZeroDivisionError:
                                    support = 0
                                values.append(support)
                    # greater than ~90% of exons supported, this is really stringent which is what we want here, as we are going to weight these models 5 to 1 over genemark
                    if float(values[1]) > 89:
                        hiQ_models.append(values[0])

            # now open evm augustus and rename models that are HiQ
            HiQ = set(hiQ_models)
            lib.log.info(
                "Found {:,} high quality predictions from Augustus (>90% exon evidence)".format(len(HiQ)))
            os.rename(Augustus, Augustus+'.bak')
            with open(Augustus, 'w') as HiQ_out:
                with open(Augustus+'.bak', 'r') as evm_aug:
                    for line in evm_aug:
                        if line.startswith('\n'):
                            HiQ_out.write(line)
                        else:
                            contig, source, feature, start, end, score, strand, phase, attributes = line.split(
                                '\t')
                            info = attributes.split(';')
                            ID = None
                            for x in info:
                                if x.startswith('ID='):
                                    ID = x.replace('-', '.')
                            if ID:
                                IDparts = ID.split('.')
                                for y in IDparts:
                                    if y.startswith('g'):
                                        ID = y
                            if feature == 'gene':
                                if ID in HiQ:
                                    HiQ_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
                                        contig, 'HiQ', feature, start, end, score, strand, phase, attributes))
                                else:
                                    HiQ_out.write(line)
                            elif feature == 'mRNA' or feature == 'exon' or feature == 'CDS':
                                ID = ID.split('.')[0]
                                if ID in HiQ:
                                    HiQ_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
                                        contig, 'HiQ', feature, start, end, score, strand, phase, attributes))
                                else:
                                    HiQ_out.write(line)

        # CodingQuarry installed and rna_bam and/or stringtie then run CodingQuarry and add to EVM
        Quarry = os.path.join(args.out, 'predict_misc', 'coding_quarry.gff3')
        if 'codingquarry' in RunModes:
            if not lib.checkannotations(Quarry):
                if RunModes['codingquarry'] == 'rna-bam':
                    if lib.cq_run_check(Quarry, args.rna_bam, args.stringtie, StartWeights['codingquarry']):
                        if not args.stringtie:
                            args.stringtie = os.path.join(
                                args.out, 'predict_misc', 'stringtie.gtf')
                            lib.log.info(
                                'Running stringie on RNA-seq alignments ')
                            lib.runStringtie(
                                args.rna_bam, args.cpus, args.stringtie)
                        if lib.checkannotations(args.stringtie):
                            lib.log.info(
                                'Running CodingQuarry prediction using stringtie alignments')
                            checkQuarry = lib.runCodingQuarry(
                                MaskGenome, args.stringtie, args.cpus, Quarry)
                            if not checkQuarry:
                                Quarry = False
                    else:
                        Quarry = False
                else:
                    # run CQ with pre-trained parameters
                    if StartWeights['codingquarry'] > 0:
                        checkQuarry = lib.runCodingQuarryTrained(MaskGenome, aug_species, os.path.join(
                            args.out, 'predict_misc', 'CodingQuarry'), args.cpus, Quarry)
                        if not checkQuarry:
                            Quarry = False
                    else:
                        Quarry = False
            else:
                lib.log.info(
                    'Using existing CodingQuarry results: {:}'.format(Quarry))
        else:
            Quarry = False

        # report gene numbers
        if Quarry:
            cqCount = lib.countGFFgenes(Quarry)
            lib.log.info('{:,} predictions from CodingQuarry'.format(cqCount))
            if os.path.isdir(os.path.join(args.out, 'predict_misc', 'CodingQuarry', 'ParameterFiles')):
                lib.copyDirectory(os.path.join(args.out, 'predict_misc', 'CodingQuarry',
                                               'ParameterFiles'), os.path.join(LOCALPARAMETERS, 'codingquarry'))
        else:
            if StartWeights['codingquarry'] > 0:
                lib.log.debug(
                    'CodingQuarry failed, removing from training parameters')
                trainingData['codingquarry'] = [{}]

        # run snap prediction
        SNAP = False
        SnapPredictions = os.path.join(
            args.out, 'predict_misc', 'snap-predictions.gff3')
        if 'snap' in RunModes:
            if not lib.checkannotations(SnapPredictions):
                if RunModes['snap'] == 'pretrained' and StartWeights['snap'] > 0 and os.path.isfile(os.path.join(LOCALPARAMETERS, aug_species+'.snap.hmm')):
                    lib.log.info(
                        'Running SNAP gene prediction, using pre-trained HMM profile')
                    lib.runSnapTrained(MaskGenome, os.path.join(
                        LOCALPARAMETERS, aug_species+'.snap.hmm'), os.path.join(args.out, 'predict_misc'), SnapPredictions)
                else:
                    if lib.snap_run_check(SnapPredictions, FinalTrainingModels, StartWeights['snap']):
                        lib.log.info('Running SNAP gene prediction, using training data: {:}'.format(
                            FinalTrainingModels))
                        lib.runSnap(MaskGenome, FinalTrainingModels, args.min_intronlen, args.max_intronlen, os.path.join(
                            args.out, 'predict_misc'), SnapPredictions)
            else:
                lib.log.info(
                    'Existing snap predictions found {:}'.format(SnapPredictions))
            if lib.checkannotations(SnapPredictions):
                snapCount = lib.countGFFgenes(SnapPredictions)
                lib.log.info('{:,} predictions from SNAP'.format(snapCount))
                if snapCount > 0:
                    SNAP = True
                else:
                    lib.log.info(
                        'SNAP prediction failed, moving on without result')
            if SNAP and os.path.isfile(os.path.join(args.out, 'predict_misc', 'snap-trained.hmm')):
                shutil.copyfile(os.path.join(args.out, 'predict_misc', 'snap-trained.hmm'),
                                os.path.join(LOCALPARAMETERS, aug_species+'.snap.hmm'))
            else:
                lib.log.debug('snap failed removing from training parameters')
                trainingData['snap'] = [{}]

        # run Glimmer predictions
        GLIMMER = False
        GlimmerPredictions = os.path.join(
            args.out, 'predict_misc', 'glimmerhmm-predictions.gff3')
        if 'glimmerhmm' in RunModes:
            if not lib.checkannotations(GlimmerPredictions):
                if RunModes['glimmerhmm'] == 'pretrained' and StartWeights['glimmerhmm'] > 0 and os.path.isdir(os.path.join(LOCALPARAMETERS, 'glimmerhmm')):
                    lib.log.info(
                        'Running GlimmerHMM gene prediction, using pretrained HMM profile')
                    lib.runGlimmerHMMTrained(MaskGenome, os.path.join(
                        LOCALPARAMETERS, 'glimmerhmm'), os.path.join(args.out, 'predict_misc'), GlimmerPredictions)
                else:
                    if lib.glimmer_run_check(GlimmerPredictions, FinalTrainingModels, StartWeights['glimmerhmm']):
                        lib.log.info('Running GlimmerHMM gene prediction, using training data: {:}'.format(
                            FinalTrainingModels))
                        lib.runGlimmerHMM(MaskGenome, FinalTrainingModels, os.path.join(
                            args.out, 'predict_misc'), GlimmerPredictions)
            else:
                lib.log.info('Existing GlimmerHMM predictions found: {:}'.format(
                    GlimmerPredictions))
            if lib.checkannotations(GlimmerPredictions):
                glimmerCount = lib.countGFFgenes(GlimmerPredictions)
                lib.log.info(
                    '{:,} predictions from GlimmerHMM'.format(glimmerCount))
                if glimmerCount > 0:
                    GLIMMER = True
                else:
                    lib.log.info(
                        'GlimmerHMM prediction failed, moving on without result')
            if GLIMMER and os.path.isdir(os.path.join(args.out, 'predict_misc', 'glimmerhmm')):
                lib.copyDirectory(os.path.join(args.out, 'predict_misc', 'glimmerhmm'), os.path.join(
                    LOCALPARAMETERS, 'glimmerhmm'))
            else:
                lib.log.debug(
                    'GlimmerHMM failed, removing from training parameters')
                trainingData = [{}]

        # EVM related input tasks, find all predictions and concatenate together
        pred_in = [Augustus]
        if GeneMark:
            pred_in.append(GeneMark)
        if args.pasa_gff:
            pred_in.append(PASA_GFF)
        if OTHER_GFFs:
            pred_in += OTHER_GFFs
        if Quarry:
            pred_in.append(Quarry)
        if SNAP:
            pred_in.append(SnapPredictions)
        if GLIMMER:
            pred_in.append(GlimmerPredictions)

        # write gene predictions file
        PredictionSources = []
        with open(Predictions+'.tmp', 'w') as output:
            for f in sorted(pred_in):
                if f:
                    with open(f, 'r') as input:
                        for line in input:
                            if not line.startswith('#'):
                                if line.count('\t') == 8:
                                    cols = line.split('\t')
                                    if not cols[1] in PredictionSources:
                                        PredictionSources.append(cols[1])
                                output.write(line)
        lib.log.debug('Prediction sources: {:}'.format(PredictionSources))
        # make sure spaces in between gene models for EVM
        with open(Predictions, 'w') as outfile:
            with open(Predictions+'.tmp', 'r') as infile:
                previous_line = None
                for line in infile:
                    if line.startswith('#'):
                        outfile.write(line)
                    if '\tgene\t' in line:
                        if not previous_line:
                            outfile.write(line)
                        elif previous_line == '\n':
                            outfile.write(line)
                            previous_line = line
                        else:
                            outfile.write('\n')
                            outfile.write(line)
                            previous_line = line
                    else:
                        outfile.write(line)
                        previous_line = line
        os.remove(Predictions+'.tmp')

        # set Weights file dependent on which data is present.
        Weights = os.path.join(args.out, 'predict_misc', 'weights.evm.txt')
        with open(Weights, 'w') as output:
            for y in PredictionSources:
                if y.lower() in EVMBase:
                    BASE = EVMBase.get(y.lower())
                else:
                    BASE = 'OTHER_PREDICTION'
                if y.lower() in StartWeights:
                    numWeight = StartWeights.get(y.lower())
                else:
                    print(('ERROR:{:} not in {:}'.format(
                        y.lower(), StartWeights)))
                    numWeight = 1
                output.write('{:}\t{:}\t{:}\n'.format(BASE, y, numWeight))
                EVMWeights[y] = numWeight
            if Exonerate:
                output.write("PROTEIN\texonerate\t{:}\n".format(
                    StartWeights.get('proteins')))
                EVMWeights['proteins'] = StartWeights.get('proteins')
            if Transcripts:
                output.write("TRANSCRIPT\tgenome\t{:}\n".format(
                    StartWeights.get('transcripts')))
                EVMWeights['transcripts'] = StartWeights.get('transcripts')

    # total up Predictions, get source counts
    EVMCounts = lib.countEVMpredictions(Predictions)
    lib.log.debug('Summary of gene models: {:}'.format(EVMCounts))
    lib.log.debug('EVM Weights: {:}'.format(EVMWeights))
    lib.log.info('Summary of gene models passed to EVM (weights):')
    lib.log.debug('Launching EVM via funannotate-runEVM.py')
    TableHeader = ['Source', 'Weight', 'Count']
    InputListCounts = []
    for k, v in list(EVMCounts.items()):
        if k in EVMWeights:
            eviweight = EVMWeights.get(k)
            if k == 'HiQ':
                k = 'Augustus HiQ'
            InputListCounts.append([k, eviweight, v])
    InputListCounts = natsorted(InputListCounts, key=lambda x: x[0])
    InputListCounts.append(['Total', '-', EVMCounts['total']])
    InputListCounts.insert(0, TableHeader)
    evm_table = lib.print_table(InputListCounts, return_str=True)
    sys.stderr.write(evm_table)

    if args.keep_evm and os.path.isfile(EVM_out):
        lib.log.info("Using existing EVM predictions: {:}".format(EVM_out))
    else:
        # setup EVM run
        EVM_script = os.path.join(
            parentdir, 'aux_scripts', 'funannotate-runEVM.py')

        # get absolute paths for everything
        Weights = os.path.abspath(Weights)
        EVM_out = os.path.abspath(EVM_out)
        Predictions = os.path.abspath(Predictions)
        EVMFolder = os.path.abspath(
            os.path.join(args.out, 'predict_misc', 'EVM'))

        # setup base evm command
        evm_cmd = [sys.executable, EVM_script, '-w', Weights,
                    '-c', str(args.cpus), '-g', Predictions,
                    '-d', EVMFolder, '-f', MaskGenome,
                    '-l', os.path.join(args.out,
                                       'logfiles',
                                       'funannotate-EVM.log'),
                    '-m', str(args.min_intronlen),
                    '-i', str(args.evm_interval),
                    '-o', EVM_out, '--EVM_HOME', EVM]

        if args.repeats2evm:
            RepeatGFF = os.path.join(
                args.out, 'predict_misc', 'repeatmasker.gff3')
            lib.bed2gff3(RepeatMasker, RepeatGFF)
            RepeatGFF = os.path.abspath(RepeatGFF)
            evm_cmd += ['-r', RepeatGFF]
        # parse entire EVM command to script
        if Exonerate:
            Exonerate = os.path.abspath(Exonerate)
            evm_cmd += ['-p', Exonerate]
        if Transcripts:
            Transcripts = os.path.abspath(Transcripts)
            evm_cmd += ['-t', Transcripts]
        if not args.evm_partitions:
            evm_cmd.append('--no-partitions')
        if not args.progress:
            evm_cmd.append('--no-progress')
        # run EVM
        lib.log.debug(' '.join(evm_cmd))
        subprocess.call(evm_cmd)

        # check output
        try:
            total = lib.countGFFgenes(EVM_out)
        except IOError:
            lib.log.error("EVM did not run correctly, output file missing")
            sys.exit(1)
        # check number of gene models, if 0 then failed, delete output file for re-running
        if total < 1:
            lib.log.error("Evidence modeler has failed, exiting")
            os.remove(EVM_out)
            sys.exit(1)
        else:
            lib.log.info('{0:,}'.format(total) + ' total gene models from EVM')

    # get protein fasta files
    evmCount = lib.countGFFgenes(EVM_out)
    lib.log.info(
        "Generating protein fasta files from {:,} EVM models".format(evmCount))
    EVM_proteins = os.path.join(
        args.out, 'predict_misc', 'evm.round1.proteins.fa')
    # translate GFF3 to proteins
    EVMGenes = {}
    EVMGenes = lib.gff2dict(EVM_out, MaskGenome, EVMGenes)
    with open(EVM_proteins, 'w') as evmprots:
        for k, v in natsorted(list(EVMGenes.items())):
            for i, x in enumerate(v['ids']):
                Prot = v['protein'][i]
                evmprots.write('>{:} {:}\n{:}\n'.format(x, k, Prot))

    # now filter bad models
    lib.log.info("now filtering out bad gene models (< %i aa in length, transposable elements, etc)." % args.min_protlen)
    Blast_rep_remove = os.path.join(
        args.out, 'predict_misc', 'repeat.gene.models.txt')
    # need to run this every time if gene models have changed from a re-run
    if os.path.isfile(Blast_rep_remove):
        os.remove(Blast_rep_remove)
    lib.RepeatBlast(EVM_proteins, args.cpus, 1e-10, FUNDB,
                    os.path.join(args.out, 'predict_misc'), Blast_rep_remove)
    EVMCleanGFF = os.path.join(args.out, 'predict_misc', 'evm.cleaned.gff3')
    if os.path.isfile(EVMCleanGFF):
        os.remove(EVMCleanGFF)
    lib.RemoveBadModels(EVM_proteins, EVM_out, args.min_protlen, RepeatMasker, Blast_rep_remove,
                        os.path.join(args.out, 'predict_misc'), args.repeat_filter, EVMCleanGFF)
    total = lib.countGFFgenes(EVMCleanGFF)
    lib.log.info('{:,} gene models remaining'.format(total))

    # run tRNAscan
    tRNAscan = os.path.join(args.out, 'predict_misc', 'trnascan.gff3')
    if args.trnascan:
        lib.log.info("Existing tRNAscan results passed: {}".format(args.trnascan))
        trna_result = lib.runtRNAscan(MaskGenome, os.path.join(args.out, 'predict_misc'),
                                      tRNAscan, cpus=args.cpus, precalc=args.trnascan)
    else:
        lib.log.info("Predicting tRNAs")
        if not os.path.isfile(tRNAscan):
            trna_result = lib.runtRNAscan(MaskGenome, os.path.join(args.out, 'predict_misc'),
                                          tRNAscan, cpus=args.cpus)
        else:
            trna_result = True

    # combine tRNAscan with EVM gff, dropping tRNA models if they overlap with EVM models
    cleanTRNA = os.path.join(args.out, 'predict_misc', 'trnascan.no-overlaps.gff3')
    if trna_result:
        lib.validate_tRNA(tRNAscan, EVMCleanGFF, AssemblyGaps, cleanTRNA)
        lib.log.info("{:,} tRNAscan models are valid (non-overlapping)".format(lib.countGFFgenes(cleanTRNA)))
    else:
        with open(cleanTRNA, 'w') as outfile:
            outfile.write('##gff-version 3\n')

    # load EVM models and tRNAscan models, output tbl annotation file
    lib.log.info("Generating GenBank tbl annotation file")
    prefix = args.name.replace('_', '')
    gag3dir = os.path.join(args.out, 'predict_misc', 'tbl2asn')
    if os.path.isdir(gag3dir):
        lib.SafeRemove(gag3dir)
    os.makedirs(gag3dir)
    tbl_file = os.path.join(gag3dir, 'genome.tbl')
    lib.GFF2tbl(EVMCleanGFF, cleanTRNA, MaskGenome, ContigSizes, prefix,
                args.numbering, args.SeqCenter, args.SeqAccession, tbl_file)

    # setup final output files
    final_fasta = os.path.join(
        args.out, 'predict_results', organism_name + '.scaffolds.fa')
    final_gff = os.path.join(
        args.out, 'predict_results', organism_name + '.gff3')
    final_gbk = os.path.join(
        args.out, 'predict_results', organism_name + '.gbk')
    final_tbl = os.path.join(
        args.out, 'predict_results', organism_name + '.tbl')
    final_proteins = os.path.join(
        args.out, 'predict_results', organism_name + '.proteins.fa')
    final_transcripts = os.path.join(
        args.out, 'predict_results', organism_name + '.mrna-transcripts.fa')
    final_cds_transcripts = os.path.join(
        args.out, 'predict_results', organism_name + '.cds-transcripts.fa')
    final_validation = os.path.join(
        args.out, 'predict_results', organism_name+'.validation.txt')
    final_error = os.path.join(
        args.out, 'predict_results', organism_name+'.error.summary.txt')
    final_fixes = os.path.join(
        args.out, 'predict_results', organism_name+'.models-need-fixing.txt')
    final_stats = os.path.join(args.out, 'predict_results', organism_name+'.stats.json')

    # run tbl2asn in new directory directory
    # setup SBT file
    SBT = os.path.join(parentdir, 'config', 'test.sbt')
    discrep = os.path.join(args.out, 'predict_results',
                           organism_name + '.discrepency.report.txt')
    lib.log.info("Converting to final Genbank format")
    # have to run as subprocess because of multiprocessing issues
    cmd = [sys.executable, os.path.join(parentdir, 'aux_scripts', 'tbl2asn_parallel.py'),
           '-i', tbl_file, '-f', MaskGenome, '-o', gag3dir, '--sbt', SBT, '-d', discrep,
           '-s', args.species, '-t', args.tbl2asn, '-v', '1', '-c', str(args.cpus)]
    if args.isolate:
        cmd += ['--isolate', args.isolate]
    if args.strain:
        cmd += ['--strain', args.strain]
    lib.log.debug(' '.join(cmd))
    subprocess.call(cmd)
    # check if completed successfully
    if not lib.checkannotations(os.path.join(gag3dir, 'genome.gbf')):
        lib.log.info('ERROR: GBK file conversion failed, tbl2asn parallel script has died')
        sys.exit(1)

    # retrieve files/reorganize
    shutil.copyfile(os.path.join(gag3dir, 'genome.gbf'), final_gbk)
    shutil.copyfile(os.path.join(gag3dir, 'genome.tbl'), final_tbl)
    shutil.copyfile(os.path.join(gag3dir, 'genome.val'), final_validation)
    shutil.copyfile(os.path.join(gag3dir, 'errorsummary.val'), final_error)
    lib.tbl2allout(final_tbl, MaskGenome, final_gff, final_proteins,
                   final_transcripts, final_cds_transcripts, final_fasta)
    lib.annotation_summary(MaskGenome, final_stats, tbl=final_tbl,
                           transcripts=Transcripts, proteins=Exonerate,
                           database=FUNDB, command=' '.join(sys.argv),
                           organism=organism_name)
    total = lib.countGFFgenes(final_gff)
    lib.log.info("Collecting final annotation files for {:,} total gene models".format(total))

    lib.log.info("Funannotate predict is finished, output files are in the %s/predict_results folder" % (args.out))

    # check if there are error that need to be fixed
    errors = lib.ncbiCheckErrors(
        final_error, final_validation, prefix, final_fixes)
    if errors > 0:
        sys.stderr.write(
            '-------------------------------------------------------\n')
        lib.log.info("Manually edit the tbl file %s, then run:\n\nfunannotate fix -i %s -t %s\n" %
                     (final_tbl, final_gbk, final_tbl))
        lib.log.info(
            "After the problematic gene models are fixed, you can proceed with functional annotation.")

    # give a suggested command
    if args.rna_bam and args.pasa_gff and os.path.isdir(os.path.join(args.out, 'training')):
        lib.log.info("Your next step to capture UTRs and update annotation using PASA:\n\n\
  funannotate update -i {:} --cpus {:}\n".format(args.out, args.cpus))
    elif args.rna_bam:  # means you have RNA-seq, but did not use funannotate train
        lib.log.info("Your next step to capture UTRs and update annotation using PASA:\n\n\
  funannotate update -i {:} --cpus {:} \\\n\
            --left illumina_forward_RNAseq_R1.fastq.gz \\\n\
            --right illumina_forward_RNAseq_R2.fastq.gz \\\n\
            --jaccard_clip\n".format(args.out, args.cpus))
    else:
        lib.log.info("Your next step might be functional annotation, suggested commands:\n\
-------------------------------------------------------\n\
Run InterProScan (Docker required): \nfunannotate iprscan -i {:} -m docker -c {:}\n\n\
Run antiSMASH: \nfunannotate remote -i {:} -m antismash -e youremail@server.edu\n\n\
Annotate Genome: \nfunannotate annotate -i {:} --cpus {:} --sbt yourSBTfile.txt\n\
-------------------------------------------------------\n\
                ".format(args.out,
                         args.cpus,
                         args.out,
                         args.out,
                         args.cpus))
    # write parameters file
    with open(os.path.join(args.out, 'predict_results', aug_species+'.parameters.json'), 'w') as outfile:
        json.dump(trainingData, outfile)
    lib.log.info('Training parameters file saved: {:}'.format(
        os.path.join(args.out, 'predict_results', aug_species+'.parameters.json')))
    lib.log.info('Add species parameters to database:\n\n\
  funannotate species -s {:} -a {:}\n'.format(aug_species, os.path.join(args.out, 'predict_results', aug_species+'.parameters.json')))
    # clean up intermediate folders
    if os.path.exists('discrepency.report.txt') and os.path.isfile('discrepency.report.txt'):
        os.rename('discrepency.report.txt', os.path.join(
            gag3dir, 'discrepency.report.txt'))
    if os.path.isfile('funannotate-EVM.log'):
        os.rename('funannotate-EVM.log', os.path.join(args.out,
                                                      'logfiles', 'funannotate-EVM.log'))
    if os.path.isfile('funannotate-p2g.log'):
        os.rename('funannotate-p2g.log', os.path.join(args.out,
                                                      'logfiles', 'funannotate-p2g.log'))


if __name__ == "__main__":
    main(sys.argv[1:])
