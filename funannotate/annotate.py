#!/usr/bin/env python
# -*- coding: utf-8 -*-


import funannotate.library as lib
import sys
import os
import subprocess
import shutil
import argparse
import re
from natsort import natsorted
import warnings
from Bio import SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO


def MEROPSBlast(input, cpus, evalue, tmpdir, output, diamond=True):
    # run blastp against merops
    blast_tmp = os.path.join(tmpdir, 'merops.xml')
    if diamond:
        blastdb = os.path.join(FUNDB, 'merops.dmnd')
        cmd = ['diamond', 'blastp', '--sensitive', '--query', input, '--threads', str(cpus),
               '--out', blast_tmp, '--db', blastdb, '--evalue', str(
                   evalue), '--max-target-seqs', '1',
               '--outfmt', '5']
    else:
        blastdb = os.path.join(FUNDB, 'MEROPS')
        cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus),
               '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input]
    if not os.path.isfile(blast_tmp):
        lib.runSubprocess4(cmd, '.', lib.log)
    # parse results
    with open(output, 'w') as out:
        with open(blast_tmp, 'r') as results:
            for qresult in SearchIO.parse(results, "blast-xml"):
                hits = qresult.hits
                ID = qresult.id
                num_hits = len(hits)
                if num_hits > 0:
                    if hits[0].hsps[0].evalue > evalue:
                        continue
                    sseqid = hits[0].id
                    out.write("%s\tnote\tMEROPS:%s\n" % (ID, sseqid))


def SwissProtBlast(input, cpus, evalue, tmpdir, GeneDict, diamond=True):
    # run blastp against uniprot
    blast_tmp = os.path.join(tmpdir, 'uniprot.xml')
    if diamond:
        blastdb = os.path.join(FUNDB, 'uniprot.dmnd')
        cmd = ['diamond', 'blastp', '--sensitive', '--query', input,
               '--threads', str(cpus), '--out', blast_tmp, '--db', blastdb,
               '--evalue', str(evalue), '--max-target-seqs',
               '1', '--outfmt', '5']
    else:
        blastdb = os.path.join(FUNDB, 'uniprot')
        cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp,
               '-num_threads', str(cpus), '-max_target_seqs', '1',
               '-evalue', str(evalue), '-query', input]
    if not lib.checkannotations(blast_tmp):
        lib.runSubprocess4(cmd, '.', lib.log)
    # parse results
    counter = 0
    total = 0
    with open(blast_tmp, 'r') as results:
        for qresult in SearchIO.parse(results, "blast-xml"):
            hits = qresult.hits
            qlen = qresult.seq_len
            ID = qresult.id
            num_hits = len(hits)
            if num_hits > 0:
                length = hits[0].hsps[0].aln_span
                pident = hits[0].hsps[0].ident_num / float(length)
                if pident < 0.6:
                    continue
                diff = length / float(qlen)
                if diff < 0.6:
                    continue
                hdescript = hits[0].description.split(' OS=')[0]
                name = hits[0].description.split('GN=')[-1]
                name = name.split(' ')[0].upper()
                name = name.replace('-', '')
                passname = None
                if not '_' in name and not ' ' in name and not '.' in name and number_present(name) and len(name) > 2 and not morethanXnumbers(name, 3):
                    passname = name
                # need to do some filtering here of certain words
                bad_words = ['(Fragment)', 'homolog', 'homolog,', 'AltName:']
                # turn string into array, splitting on spaces
                descript = hdescript.split(' ')
                final_desc = [x for x in descript if x not in bad_words]
                final_desc = ' '.join(final_desc)
                total += 1
                # add to GeneDict
                if passname:
                    counter += 1
                    if not ID in GeneDict:
                        GeneDict[ID] = [
                            {'name': passname, 'product': final_desc,
                             'source': 'UniProtKB'}]
                    else:
                        GeneDict[ID].append(
                            {'name': passname, 'product': final_desc,
                             'source': 'UniProtKB'})
    lib.log.info(
        '{:,} valid gene/product annotations from {:,} total'.format(counter, total))


def number_present(s):
    return any(i.isdigit() for i in s)


def morethanXnumbers(s, num):
    count = 0
    for i in s:
        if number_present(i):
            count += 1
    if count >= num:
        return True
    else:
        return False


def capfirst(x):
    return x[0].upper() + x[1:]


def item2index(inputList, item):
    # return the index of an item in the input list
    item_index = None
    for x in inputList:
        if item.lower() in x.lower():
            item_index = inputList.index(x)
    return item_index


def getEggNogHeaders(input):
    '''
    function to get the headers from eggnog mapper annotations
    web-based eggnog mapper has no header....
    #web based 'guess'
    0   query_name
    1   seed_eggNOG_ortholog
    2   seed_ortholog_evalue
    3   seed_ortholog_score
    4   predicted_gene_name
    5   GO_terms
    6   KEGG_KOs
    7   BiGG_reactions
    8   Annotation_tax_scope
    9   OGs
    10  bestOG|evalue|score
    11  COG cat
    12  eggNOG annot
    '''
    IDi, DBi, OGi, Genei, COGi, Desci = (None,)*6
    with open(input, 'r') as infile:
        for line in infile:
            if line.startswith('#query_name'):  # this is HEADER
                line = line.rstrip()
                headerCols = line.split('\t')
                IDi = item2index(headerCols, 'query_name')
                Genei = item2index(headerCols, 'predicted_gene_name')
                DBi = item2index(headerCols, 'Annotation_tax_scope')
                OGi = item2index(headerCols, 'OGs')
                COGi = item2index(headerCols, 'COG cat')
                Desci = item2index(headerCols, 'eggNOG annot')
                break
    if not IDi:  # then no header file, so have to guess
        IDi, DBi, OGi, Genei, COGi, Desci = (0, 8, 9, 4, 11, 12)
    return IDi, DBi, OGi, Genei, COGi, Desci, None


def getEggNogHeadersv2(input):
    '''
    function to get the headers from eggnog mapper annotations
    web-based eggnog mapper has no header....
    '''
    IDi, DBi, OGi, Genei, COGi, Desci, ECi = (None,)*7
    with open(input, 'r') as infile:
        for line in infile:
            if line.startswith('#query'):  # this is HEADER
                line = line.rstrip()
                headerCols = line.split('\t')
                IDi = 0
                Genei = item2index(headerCols, 'Preferred_name')
                DBi = item2index(headerCols, 'eggNOG OGs')
                OGi = item2index(headerCols, 'best_og_name')
                COGi = item2index(headerCols, 'best_og_cat')
                Desci = item2index(headerCols, 'best_og_desc')
                ECi = item2index(headerCols, 'EC')
                break
    if not IDi:  # then no header file, so have to guess
        IDi, DBi, OGi, Genei, COGi, Desci, ECi = (0, 4, 8, 11, 9, 10, 13)
    return IDi, DBi, OGi, Genei, COGi, Desci, ECi

def parseEggNoggMapper(input, output, GeneDict):
    # try to parse header
    version, prefix = getEggnogVersion(input)
    if version and version > ('2.0.0') and version < ('2.0.5'):
        lib.log.error('Unable to parse emapper results from v{}, please use either v1.0.3 or >=v2.0.5'.format(version))
        return {}
    if not prefix:  # we have to guess here, sorry
        prefix = 'ENOG50'
    if not version:  # also then we guess
        version = '2.1.0'
    lib.log.debug('EggNog annotation detected as emapper v{} and DB prefix {}'.format(version, prefix))
    Definitions = {}
    # indexes from header file
    if version < ('2.0.0'):
        IDi, DBi, OGi, Genei, COGi, Desci, ECi = getEggNogHeaders(input)
    else:
        IDi, DBi, OGi, Genei, COGi, Desci, ECi = getEggNogHeadersv2(input)
    # take annotations file from eggnog-mapper and create annotations
    with open(output, 'w') as out:
        with open(input, 'r') as infile:
            for line in infile:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                cols = ['' if x=='-' else x for x in cols]
                ID = cols[IDi]
                Description = cols[Desci].split('. ')[0]
                Gene = ''
                if cols[Genei] != '':
                    if not '_' in cols[Genei] and not '.' in cols[Genei] and number_present(cols[Genei]) and len(cols[Genei]) > 2 and not morethanXnumbers(cols[Genei], 3):
                        Gene = cols[Genei]
                if version < ('2.0.0'):
                    EC = None
                    DB = cols[DBi].split('[')[0]
                    OGs = cols[OGi].split(',')
                    NOG = ''
                    for x in OGs:
                        if DB in x:
                            NOG = prefix + x.split('@')[0]
                    COGs = cols[COGi].replace(' ', '')
                else:  # means we have v2 or great
                    NOG, DB = cols[OGi].split('@')
                    OGs = cols[DBi].split(',')
                    if NOG == 'seed_ortholog': # not sure if this is bug, but get second to last OG from all
                        NOG, DB = OGs[-2].split('@')
                    DB = DB.split('|')[-1]
                    NOG = prefix+NOG
                    EC = cols[ECi]
                    if ',' in EC: # this is least common ancestor approach
                        EC = os.path.commonprefix(EC.split(',')).rstrip('.')
                    COGs = cols[COGi].replace(' ', '')
                    if len(COGs) > 1:
                        COGs = ''.join([c + ',' for c in COGs]).rstrip(',')

                if EC and EC != '':
                    out.write("%s\tEC_number\t%s\n" % (ID, EC))
                if NOG == '':
                    continue
                if not NOG in Definitions:
                    Definitions[NOG] = Description
                out.write("%s\tnote\tEggNog:%s\n" % (ID, NOG))
                if COGs != '':
                    out.write("%s\tnote\tCOG:%s\n" % (ID, COGs))
                if Gene != '':
                    product = Gene.lower()+'p'
                    product = capfirst(product)
                    GeneID = ID
                    if not GeneID in GeneDict:
                        GeneDict[GeneID] = [
                            {'name': Gene, 'product': Description, 'source': 'EggNog-Mapper'}]
                    else:
                        GeneDict[GeneID].append(
                            {'name': Gene, 'product': Description, 'source': 'EggNog-Mapper'})
    return Definitions


def getEggnogVersion(annotfile):
    # try to parse the version of eggnog mapper used
    # caveat here is web eggnog has no header!
    vers = None
    prefix = None
    with open(annotfile, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if not line.startswith('#'):
                return vers, prefix
            else:
                if line.startswith('# emapper version:'):
                    vers = line.split('emapper-')[-1].split()[0]
                    prefix = 'ENOG41'
                if line.startswith('## emapper-'):
                    vers = line.split('## emapper-')[-1]
                    prefix = 'ENOG50'
    return vers, prefix


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(
        prog='funannotate-functional.py',
        usage="%(prog)s [options] -i folder --eggnog emapper.annotations --iprscan proteins.xml --cpus 12",
        description='''Script that adds functional annotation to a genome.''',
        epilog="""Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com""",
        formatter_class=MyFormatter)
    parser.add_argument('-i', '--input',
                        help='Folder from funannotate predict.')
    parser.add_argument('--genbank', help='Annotated genome in GenBank format')
    parser.add_argument('--fasta', help='Genome in FASTA format')
    parser.add_argument('--gff', help='GFF3 annotation file')
    parser.add_argument('-o', '--out', help='Basename of output files')
    parser.add_argument('--sbt', default='SBT',
                        help='Basename of output files')
    parser.add_argument('-s', '--species',
                        help='Species name (e.g. "Aspergillus fumigatus")')
    parser.add_argument('-t', '--tbl2asn', default='-l paired-ends',
                        help='Custom parameters for tbl2asn, example: linkage and gap info')
    parser.add_argument('-a', '--annotations',
                        help='Custom annotations, tsv 3 column file')
    parser.add_argument('--isolate', help='Isolate name (e.g. Af293)')
    parser.add_argument('--strain', help='Strain name (e.g. CEA10)')
    parser.add_argument('--cpus', default=2, type=int,
                        help='Number of CPUs to use')
    parser.add_argument('--iprscan',
                        help='IPR5 XML file or folder of pre-computed InterProScan results')
    parser.add_argument('--antismash',
                        help='antiSMASH results in genbank format')
    parser.add_argument('--signalp',
                        help='signalp results caculted elsewhere')
    parser.add_argument('--force', action='store_true',
                        help='Over-write output folder')
    parser.add_argument('--phobius', help='Phobius results')
    parser.add_argument('--eggnog', help='EggNog Mapper annotations')
    parser.add_argument('--busco_db', default='dikarya',
                        help='BUSCO model database')
    parser.add_argument('--p2g', help='NCBI p2g file from previous annotation')
    parser.add_argument('-d', '--database',
                        help='Path to funannotate database, $FUNANNOTATE_DB')
    parser.add_argument('--fix',
                        help='TSV ID GeneName Product file to over-ride automated process')
    parser.add_argument('--remove',
                        help='TSV ID GeneName Product file to remove from annotation')
    parser.add_argument('--rename', help='Rename locus tag')
    parser.add_argument('--no-progress', dest='progress', action='store_false',
                        help='no progress on multiprocessing')
    parser.add_argument('--header_length', default=16,
                        type=int, help='Max length for fasta headers')
    args = parser.parse_args(args)

    global parentdir, IPR2ANNOTATE, FUNDB
    parentdir = os.path.join(os.path.dirname(__file__))
    IPR2ANNOTATE = os.path.join(
        parentdir, 'aux_scripts', 'iprscan2annotations.py')

    # start here rest of script
    # create log file
    log_name = 'funannotate-annotate.'+str(os.getpid())+'.log'
    if os.path.isfile(log_name):
        os.remove(log_name)

    # initialize script, log system info and cmd issue at runtime
    lib.setupLogging(log_name)
    cmd_args = " ".join(sys.argv)+'\n'
    lib.log.debug(cmd_args)
    print("-------------------------------------------------------")
    lib.SystemInfo()

    # get version of funannotate
    version = lib.get_version()
    lib.log.info("Running %s" % version)

    # check dependencies
    if args.antismash:
        programs = ['hmmscan', 'hmmsearch', 'diamond', 'bedtools']
    else:
        programs = ['hmmscan', 'hmmsearch', 'diamond']
    lib.CheckDependencies(programs)

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

    # check database sources, so no problems later
    sources = [os.path.join(FUNDB, 'Pfam-A.hmm.h3p'), os.path.join(FUNDB, 'dbCAN.hmm.h3p'),
               os.path.join(FUNDB, 'merops.dmnd'), os.path.join(FUNDB, 'uniprot.dmnd')]
    if not all([os.path.isfile(f) for f in sources]):
        lib.log.error(
            'Database files not found in %s, run funannotate database and/or funannotate setup' % FUNDB)
        sys.exit(1)

    # check if diamond version matches database version
    if not lib.CheckDiamondDB(os.path.join(FUNDB, 'merops.dmnd')):
        lib.log.error(
            'Diamond merops database was created with different version of diamond, please re-run funannotate setup')
        sys.exit(1)
    if not lib.CheckDiamondDB(os.path.join(FUNDB, 'uniprot.dmnd')):
        lib.log.error(
            'Diamond uniprot database was created with different version of diamond, please re-run funannotate setup')
        sys.exit(1)

    # write versions of Databases used to logfile
    versDB = {}
    if not lib.checkannotations(os.path.join(FUNDB, 'funannotate-db-info.txt')):
        lib.log.error('Database not properly configured, %s missing. Run funannotate database and/or funannotate setup.' %
                      os.path.join(FUNDB, 'funannotate-db-info.txt'))
        sys.exit(1)
    with open(os.path.join(FUNDB, 'funannotate-db-info.txt'), 'r') as dbfile:
        for line in dbfile:
            line = line.strip()
            name, type, file, version, date, num_records, mdchecksum = line.split(
                '\t')
            versDB[name] = version

    # take care of some preliminary checks
    if args.sbt == 'SBT':
        SBT = os.path.join(parentdir, 'config', 'test.sbt')
        lib.log.info(
            "No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'")
    else:
        SBT = args.sbt

    # check other input files
    if not os.path.isfile(SBT):
        lib.log.error("SBT file not found, exiting")
        sys.exit(1)
    if args.antismash:
        if not os.path.isfile(args.antismash):
            lib.log.error("Antismash GBK file not found, exiting")
            sys.exit(1)

    # check buscos, download if necessary
    if not os.path.isdir(os.path.join(FUNDB, args.busco_db)):
        lib.log.error("ERROR: %s busco database is not found, install with funannotate setup -b %s" %
                      (args.busco_db, args.busco_db))
        sys.exit(1)

    # need to do some checks here of the input
    genbank, Scaffolds, Protein, Transcripts, GFF, TBL = (None,)*6
    existingStats = False
    GeneCounts = 0
    GeneDB = {}
    if not args.input:
        # did not parse folder of funannotate results, so need either gb + gff or fasta + proteins, + gff and also need to have args.out for output folder
        if not args.out:
            lib.log.error(
                "If you are not providing funannotate predict input folder, then you need to provide an output folder (--out)")
            sys.exit(1)
        else:
            outputdir = args.out
            if os.path.isdir(outputdir):
                lib.log.error(
                    "Found existing output directory %s. Warning, will re-use any intermediate files found." % (outputdir))
            # create outputdir and subdirs if not already present
            lib.createdir(outputdir)
            lib.createdir(os.path.join(outputdir, 'annotate_misc'))
            lib.createdir(os.path.join(outputdir, 'annotate_results'))
            lib.createdir(os.path.join(outputdir, 'logfiles'))

        if not args.genbank:
            if not args.fasta or not args.gff:
                lib.log.error(
                    "You did not specifiy the apropriate input files, either: \n1) GenBank \n2) Genome FASTA + GFF3")
                sys.exit(1)
            else:
                Scaffolds = args.fasta
                GFF = args.gff
                Proteins = os.path.join(outputdir, 'annotate_misc', 'genome.proteins.fa')
                Transcripts = os.path.join(outputdir, 'annotate_misc', 'genome.transcripts.fasta')
                annotTBL = os.path.join(outputdir, 'annotate_misc', 'genome.tbl')
                prefix = None
                if args.rename:
                    prefix = args.rename.replace('_', '')
                lib.log.info(
                    "Parsing annotation and preparing annotation files.")
                GeneCounts, GeneDB = lib.convertgff2tbl(
                    GFF, prefix, Scaffolds, Proteins, Transcripts, annotTBL, external=True)
        else:
            genbank = args.genbank
            Scaffolds = os.path.join(outputdir, 'annotate_misc', 'genome.scaffolds.fasta')
            Proteins = os.path.join(outputdir, 'annotate_misc', 'genome.proteins.fasta')
            Transcripts = os.path.join(outputdir, 'annotate_misc', 'genome.transcripts.fasta')
            GFF = os.path.join(outputdir, 'annotate_misc', 'genome.gff3')
            annotTBL = os.path.join(outputdir, 'annotate_misc', 'genome.tbl')
            lib.log.info("Checking GenBank file for annotation")
            if not lib.checkGenBank(genbank):
                lib.log.error("Found no annotation in GenBank file, exiting")
                sys.exit(1)
            GeneCounts = lib.gb2parts(
                genbank, annotTBL, GFF, Proteins, Transcripts, Scaffolds)
    else:
        # should be a folder, with funannotate files, thus store results there, no need to create output folder
        if not os.path.isdir(args.input):
            lib.log.error("%s directory does not exist" % args.input)
            sys.exit(1)
        # funannotate results 1) in update folder or 2) in predict folder
        if os.path.isdir(os.path.join(args.input, 'update_results')):
            inputdir = os.path.join(args.input, 'update_results')
            outputdir = args.input
        elif os.path.isdir(os.path.join(args.input, 'predict_results')):
            inputdir = os.path.join(args.input, 'predict_results')
            outputdir = args.input
        else:
            # here user specified the predict_results folder, or it is a custom folder
            inputdir = os.path.abspath(args.input)
            if '_results' in inputdir:  # then it is the _results dir, so move up one directory
                outputdir = os.path.dirname(inputdir)
            else:
                lib.log.error(
                    'Unable to detect funannotate folder as input, please provide -o,--out directory')
                sys.exit(1)

        annotTBL = os.path.join(outputdir, 'annotate_misc', 'genome.tbl')

        # get files that you need
        for file in os.listdir(inputdir):
            if file.endswith('.gbk'):
                genbank = os.path.join(inputdir, file)
            if file.endswith('.gff3'):
                GFF = os.path.join(inputdir, file)
            if file.endswith('.tbl'):
                TBL = os.path.join(inputdir, file)
            if file.endswith('.stats.json'):
                existingStats = os.path.join(inputdir, file)

        # now create the files from genbank input file for consistency in gene naming, etc
        if not genbank or not GFF:
            lib.log.error(
                "Properly formatted 'funannotate predict' files do no exist in this directory")
            sys.exit(1)
        else:
            # if user gave predict_results folder, then set output to up one directory
            if 'predict_results' in inputdir or 'update_results' in inputdir:
                outputdir = lib.get_parent_dir(inputdir)
            else:
                if not args.out:
                    outputdir = inputdir  # output the results in the input directory
                else:
                    outputdir = args.out
                    if not os.path.isdir(outputdir):
                        os.makedirs(outputdir)
            # create output directories
            if not os.path.isdir(os.path.join(outputdir, 'annotate_misc')):
                os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            if not os.path.isdir(os.path.join(outputdir, 'annotate_results')):
                os.makedirs(os.path.join(outputdir, 'annotate_results'))
            else:
                lib.log.error(
                    "Found existing output directory %s. Warning, will re-use any intermediate files found." % (outputdir))
            lib.log.info("Parsing input files")
            Scaffolds = os.path.join(outputdir, 'annotate_misc', 'genome.scaffolds.fasta')
            Proteins = os.path.join(outputdir, 'annotate_misc', 'genome.proteins.fasta')
            Transcripts = os.path.join(outputdir, 'annotate_misc', 'genome.transcripts.fasta')
            if TBL:
                lib.log.info('Existing tbl found: {:}'.format(TBL))
                shutil.copyfile(TBL, annotTBL)
                if not lib.checkannotations(GFF):
                    GFF = os.path.join(
                        outputdir, 'annotate_misc', 'genome.gff3')
                    GeneCounts = lib.gb2gffnuc(
                        genbank, GFF, Proteins, Transcripts, Scaffolds)
                else:
                    GeneCounts = lib.gb2nucleotides(
                        genbank, Proteins, Transcripts, Scaffolds)
            else:
                GFF = os.path.join(outputdir, 'annotate_misc', 'genome.gff3')
                GeneCounts = lib.gb2parts(
                    genbank, annotTBL, GFF, Proteins, Transcripts, Scaffolds)

    # double check that you have a TBL file, otherwise will have nothing to append to.
    if not lib.checkannotations(annotTBL):
        lib.log.error("NCBI tbl file not found, exiting")
        sys.exit(1)
    lib.log.debug('TBL file: {}'.format(annotTBL))
    if not lib.checkannotations(GFF):
        lib.log.error("GFF file not found, exiting")
        sys.exit(1)
    lib.log.debug('GFF3 file: {}'.format(GFF))
    if not lib.checkannotations(Proteins):
        lib.log.error("Protein FASTA file not found, exiting")
        sys.exit(1)
    lib.log.debug('Proteins file: {}'.format(Proteins))

    # parse prefix from tbl file for existing
    locusTagPrefix = None
    with open(annotTBL, 'r') as infile:
        for line in infile:
            if line.startswith('\t\t\tlocus_tag\t'):
                prelimTag = line.split('\t')[-1].rstrip()
                if '_' in prelimTag:
                    locusTagPrefix = prelimTag.split('_')[0]
                break
    if args.rename and not locusTagPrefix:
        lib.log.error('Error parsing existing locus_tag, expecting underscore "_" in locus_tag')
        sys.exit(1)
    # make sure logfiles directory is present, will need later
    if not os.path.isdir(os.path.join(outputdir, 'logfiles')):
        os.makedirs(os.path.join(outputdir, 'logfiles'))
    if not os.path.isdir(os.path.join(outputdir, 'annotate_results')):
        os.makedirs(os.path.join(outputdir, 'annotate_results'))

    # get absolute path for all input so there are no problems later, not using Transcripts yet could be error? so take out here
    Scaffolds, Proteins, GFF = [os.path.abspath(i) for i in [Scaffolds, Proteins, GFF]]

    # check the genome fasta for any potential errors
    bad_headers, bad_contigs, suspect_contigs = lib.analyzeAssembly(
        Scaffolds, header_max=args.header_length)
    if len(bad_headers) > 0 and not args.force:
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

    # get organism and isolate from GBK file
    organism, strain, isolate, accession, WGS_accession, gb_gi, version = (None,)*7
    if genbank:
        organism, strain, isolate, accession, WGS_accession, gb_gi, version = lib.getGBKinfo(
            genbank)
        # since can't find a way to propage the WGS_accession, writing to a file and then parse here
        if os.path.isfile(os.path.join(outputdir, 'update_results', 'WGS_accession.txt')):
            with open(os.path.join(outputdir, 'update_results', 'WGS_accession.txt'), 'r') as infile:
                for line in infile:
                    line = line.replace('\n', '')
                    if line == 'None':
                        WGS_accession = None
                    else:
                        WGS_accession = line

    # if command line species/strain/isolate passed, over-write detected
    # check if organism/species/isolate passed at command line, if so, overwrite what you detected.
    if args.species:
        organism = args.species
    if args.strain:
        strain = args.strain
    if args.isolate:
        isolate = args.isolate
    if not organism:
        lib.log.error(
            "No GenBank species and no species name given will cause problems downstream, please pass a name to -s,--species")
        sys.exit(1)
    if strain:
        organism_name = organism+'_'+strain
    elif isolate:
        organism_name = organism+'_'+isolate
    else:
        organism_name = organism
    organism_name = organism_name.replace(' ', '_')

    lib.log.info("Adding Functional Annotation to %s, NCBI accession: %s" % (
        organism, WGS_accession))
    lib.log.info(
        "Annotation consists of: {:,} gene models".format(int(GeneCounts)))

    ############################################################################
    # start workflow here
    ProtCount = lib.countfasta(Proteins)
    lib.log.info('{0:,}'.format(ProtCount) + ' protein records loaded')
    if ProtCount < 1:
        lib.log.error("There are no gene models in this genbank file")
        sys.exit(1)

    # create tmpdir folder and split proteins into X CPUs to run with HMMER3 searches
    protDir = os.path.join(outputdir, 'annotate_misc', 'split_prots')
    if not os.path.isdir(protDir):
        os.makedirs(protDir)
    lib.fasta2chunks(Proteins, args.cpus, os.path.join(
        outputdir, 'annotate_misc'), 'split_prots')

    # run PFAM-A search
    pfam_results = os.path.join(
        outputdir, 'annotate_misc', 'annotations.pfam.txt')
    if not lib.checkannotations(pfam_results):
        lib.log.info("Running HMMer search of PFAM version %s" %
                     versDB.get('pfam'))
        cmd = [sys.executable, os.path.join(parentdir, 'aux_scripts', 'hmmer_parallel.py'),
               '-c', str(args.cpus), '-d', FUNDB, '-i', protDir,
               '-o', pfam_results, '-m', 'pfam']
        subprocess.call(cmd)
    else:
        lib.log.info('Existing Pfam-A results found: {:}'.format(pfam_results))
    num_annotations = lib.line_count(pfam_results)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')

    # initiate Gene Name/Product dictionary
    GeneProducts = {}

    # run SwissProt Blast search
    lib.log.info("Running Diamond blastp search of UniProt DB version %s" %
                 versDB.get('uniprot'))
    blast_out = os.path.join(outputdir, 'annotate_misc',
                             'annotations.swissprot.txt')
    SwissProtBlast(Proteins, args.cpus, 1e-5,
                   os.path.join(outputdir, 'annotate_misc'), GeneProducts)

    # Check for EggNog annotations, parse if present
    eggnog_out = os.path.join(
        outputdir, 'annotate_misc', 'annotations.eggnog.txt')
    eggnog_result = os.path.join(
        outputdir, 'annotate_misc', 'eggnog.emapper.annotations')
    if args.eggnog:
        if os.path.isfile(eggnog_result):
            os.remove(eggnog_result)
        shutil.copyfile(args.eggnog, eggnog_result)
    if not lib.checkannotations(eggnog_result):
        if lib.which('emapper.py'):  # eggnog installed, so run it
            lib.log.info("Running Eggnog-mapper")
            cmd = ['emapper.py', '-m', 'diamond', '-i', Proteins,
                   '-o', 'eggnog', '--cpu', str(args.cpus)]
            lib.runSubprocess(cmd, os.path.join(
                outputdir, 'annotate_misc'), lib.log)
        else:
            lib.log.info(
                "Install eggnog-mapper or use webserver to improve functional annotation: https://github.com/jhcepas/eggnog-mapper")
    else:
        lib.log.info(
            'Existing Eggnog-mapper results found: {:}'.format(eggnog_result))
    if lib.checkannotations(eggnog_result):
        lib.log.info("Parsing EggNog Annotations")
        EggNog = parseEggNoggMapper(eggnog_result, eggnog_out, GeneProducts)
        if lib.checkannotations(eggnog_out):
            num_annotations = lib.line_count(eggnog_out)
            lib.log.info('{0:,}'.format(num_annotations) +
                        ' COG and EggNog annotations added')
    else:
        lib.log.error("No Eggnog-mapper results found.")
        EggNog = {}

    RawProductNames = os.path.join(
        outputdir, 'annotate_misc', 'uniprot_eggnog_raw_names.txt')
    # GeneDict[ID] = [{'name': passname, 'product': final_desc}]
    with open(RawProductNames, 'w') as uniprottmp:
        for k, v in natsorted(list(GeneProducts.items())):
            for x in v:  # v is list of dictionaries
                uniprottmp.write('{:}\t{:}\t{:}\t{:}\n'.format(
                    k, x['name'], x['product'], x['source']))

    # combine the results from UniProt and Eggnog to parse Gene names and product descriptions
    # load curated list
    lib.log.info("Combining UniProt/EggNog gene and product names using Gene2Product version %s" %
                 versDB.get('gene2product'))
    CuratedNames = {}
    with open(os.path.join(FUNDB, 'ncbi_cleaned_gene_products.txt'), 'r') as input:
        for line in input:
            line = line.strip()
            if line.startswith('#'):
                continue
            ID, product = line.split('\t')
            if not ID in CuratedNames:
                CuratedNames[ID] = product

    GeneSeen = {}
    NeedCurating = {}
    NotInCurated = {}
    thenots = []
    for k, v in natsorted(list(GeneProducts.items())):
        GeneName = None
        GeneProduct = None
        for x in v:
            if x['name'] in CuratedNames:
                GeneProduct = CuratedNames.get(x['name'])
                GeneName = x['name']
            elif x['name'].lower() in CuratedNames:
                GeneProduct = CuratedNames.get(x['name'].lower())
                GeneName = x['name']
        if not GeneName:  # taking first one will default to swissprot if products for both
            GeneName = v[0]['name']
            GeneProduct = v[0]['product']
            OriginalProd = GeneProduct
            thenots.append(GeneName)
            # if not GeneName in NotInCurated:
            #    NotInCurated[GeneName] = GeneProduct
        # now attempt to clean the product name
        rep = {'potential': 'putative', 'possible': 'putative', 'probable': 'putative', 'predicted': 'putative',
               'uncharacterized': 'putative', 'uncharacterised': 'putative', 'homolog': '', 'EC': '', 'COG': '',
               'inactivated': '', 'related': '', 'family': '', 'gene': 'protein', 'homologue': '', 'open reading frame': '',
               'frame': '', 'yeast': '', 'Drosophila': '', 'Yeast': '', 'drosophila': ''}
        # replace words in dictionary, from https://stackoverflow.com/questions/6116978/python-replace-multiple-strings
        rep = dict((re.escape(k), v) for k, v in rep.items())
        pattern = re.compile("|".join(list(rep.keys())))
        GeneProduct = pattern.sub(
            lambda m: rep[re.escape(m.group(0))], GeneProduct)
        # if gene name in product, convert to lowercase
        if GeneName in GeneProduct:
            GeneProduct = GeneProduct.replace(GeneName, GeneName.lower())
        # check for some obvious errors, then change product description to gene name + p
        if not GeneName in CuratedNames:
            # some eggnog descriptions are paragraphs....
            if 'By similarity' in GeneProduct or 'Required for' in GeneProduct or 'nvolved in' in GeneProduct or 'protein '+GeneName == GeneProduct or 'nherit from' in GeneProduct or len(GeneProduct) > 100:
                OriginalProd = GeneProduct
                GeneProduct = GeneName.lower()+'p'
                GeneProduct = capfirst(GeneProduct)
                if not GeneName in NeedCurating:
                    NeedCurating[GeneName] = [(OriginalProd, GeneProduct)]
                else:
                    NeedCurating[GeneName].append((OriginalProd, GeneProduct))
        # make sure not multiple spaces
        GeneProduct = ' '.join(GeneProduct.split())
        GeneProduct = GeneProduct.replace('()', '')
        if '(' in GeneProduct and not ')' in GeneProduct:
            GeneProduct = GeneProduct.split('(')[0].rstrip()
        GeneProduct = GeneProduct.replace(' ,', ',')
        # populate dictionary of NotInCurated
        if GeneName in thenots:
            if not GeneName in NotInCurated:
                NotInCurated[GeneName] = [(OriginalProd, GeneProduct)]
            else:
                NotInCurated[GeneName].append((OriginalProd, GeneProduct))
        if not GeneName in GeneSeen:
            GeneSeen[GeneName] = [(k, GeneProduct)]
        else:
            GeneSeen[GeneName].append((k, GeneProduct))

    # finally output the annotations
    # which genes are duplicates, need to append numbers to those gene names and then finally output annotations
    Gene2ProdFinal = {}
    with open(os.path.join(outputdir, 'annotate_misc', 'annotations.genes-products.txt'), 'w') as gene_annotations:
        for key, value in natsorted(list(GeneSeen.items())):
            if len(value) > 1:
                try:
                    testMultiple = len(set([x[0].split('-T')[0] for x in value]))
                except:
                    testMultiple = len(value)
                for i in range(0, len(value)):
                    if testMultiple > 1:
                        gene_annotations.write(
                            "%s\tname\t%s_%i\n" % (value[i][0], key, i+1))
                    else:
                        gene_annotations.write(
                            "%s\tname\t%s\n" % (value[i][0], key))
                    gene_annotations.write(
                        "%s\tproduct\t%s\n" % (value[i][0], value[i][1]))
                    Gene2ProdFinal[value[i][0]] = (
                        key+'_'+str(i+1), value[i][1])
            else:
                gene_annotations.write("%s\tname\t%s\n" % (value[0][0], key))
                gene_annotations.write("%s\tproduct\t%s\n" %
                                       (value[0][0], value[0][1]))
                Gene2ProdFinal[value[0][0]] = (key, value[0][1])
    num_annotations = int(lib.line_count(os.path.join(
        outputdir, 'annotate_misc', 'annotations.genes-products.txt')) / 2)
    lib.log.info('{:,} gene name and product description annotations added'.format(
        num_annotations))

    # run MEROPS Blast search
    blast_out = os.path.join(outputdir, 'annotate_misc',
                             'annotations.merops.txt')
    if not lib.checkannotations(blast_out):
        lib.log.info(
            "Running Diamond blastp search of MEROPS version %s" % versDB.get('merops'))
        MEROPSBlast(Proteins, args.cpus, 1e-5,
                    os.path.join(outputdir, 'annotate_misc'), blast_out)
    else:
        lib.log.info('Existing MEROPS results found: {:}'.format(blast_out))
    num_annotations = lib.line_count(blast_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')

    # run dbCAN search
    dbCAN_out = os.path.join(outputdir, 'annotate_misc',
                             'annotations.dbCAN.txt')
    if not lib.checkannotations(dbCAN_out):
        lib.log.info(
            "Annotating CAZYmes using HMMer search of dbCAN version %s" % versDB.get('dbCAN'))
        cmd = [sys.executable,
               os.path.join(parentdir, 'aux_scripts', 'hmmer_parallel.py'),
               '-c', str(args.cpus), '-d', FUNDB, '-i', protDir,
               '-o', dbCAN_out, '-m', 'cazy']
        subprocess.call(cmd)
    else:
        lib.log.info('Existing CAZYme results found: {:}'.format(dbCAN_out))
    num_annotations = lib.line_count(dbCAN_out)
    lib.log.info('{:,} annotations added'.format(num_annotations))

    # run BUSCO OGS search
    busco_out = os.path.join(
        outputdir, 'annotate_misc', 'annotations.busco.txt')
    buscoDB = os.path.join(FUNDB, args.busco_db)
    if not lib.checkannotations(busco_out):
        lib.log.info("Annotating proteins with BUSCO %s models" %
                     args.busco_db)
        lib.runBUSCO(Proteins, buscoDB, args.cpus, os.path.join(
            outputdir, 'annotate_misc'), busco_out)
    else:
        lib.log.info('Existing BUSCO2 results found: {:}'.format(busco_out))
    num_annotations = lib.line_count(busco_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')

    # run Phobius if local is installed, otherwise you will have to use funannotate remote
    phobius_out = os.path.join(
        outputdir, 'annotate_misc', 'phobius.results.txt')
    phobiusLog = os.path.join(outputdir, 'logfiles', 'phobius.log')
    if args.phobius:
        if os.path.isfile(phobius_out):
            os.remove(phobius_out)
        shutil.copyfile(args.phobius, phobius_out)
    if not lib.checkannotations(phobius_out):
        if lib.which('phobius.pl'):
            if not lib.checkannotations(phobius_out):
                lib.log.info(
                    "Predicting secreted and transmembrane proteins using Phobius")
                subprocess.call([os.path.join(parentdir, 'aux_scripts', 'phobius-multiproc.py'),
                                 '-i', Proteins, '-o', phobius_out, '-l', phobiusLog])
        else:
            lib.log.info(
                "Skipping phobius predictions, try funannotate remote -m phobius")
    else:
        lib.log.info('Existing Phobius results found: {:}'.format(phobius_out))

    # run signalP if installed, have to manually install, so test if exists first, then run it if it does, parse results
    signalp_out = os.path.join(
        outputdir, 'annotate_misc', 'signalp.results.txt')
    secreted_out = os.path.join(
        outputdir, 'annotate_misc', 'annotations.secretome.txt')
    membrane_out = os.path.join(
        outputdir, 'annotate_misc', 'annotations.transmembrane.txt')
    if args.signalp:
        shutil.copyfile(args.signalp, signalp_out)
    if lib.which('signalp') or lib.checkannotations(signalp_out):
        if not lib.checkannotations(signalp_out):
            lib.log.info("Predicting secreted proteins with SignalP")
            lib.signalP(Proteins, os.path.join(
                outputdir, 'annotate_misc'), signalp_out)
        else:
            lib.log.info(
                'Existing SignalP results found: {:}'.format(signalp_out))
        if lib.checkannotations(phobius_out):
            lib.parsePhobiusSignalP(
                phobius_out, signalp_out, membrane_out, secreted_out)
        else:
            lib.parseSignalP(signalp_out, secreted_out)
    else:
        if not lib.checkannotations(phobius_out):
            lib.log.info(
                "Skipping secretome: neither SignalP nor Phobius searches were run")
        else:
            lib.log.info(
                "SignalP not installed, secretome prediction less accurate using only Phobius")
            lib.parsePhobiusSignalP(
                phobius_out, False, membrane_out, secreted_out)
    if lib.checkannotations(secreted_out):
        num_secreted = lib.line_count(secreted_out)
    else:
        num_secreted = 0
    if lib.checkannotations(membrane_out):
        num_mem = lib.line_count(membrane_out)
    else:
        num_mem = 0
    lib.log.info('{0:,}'.format(num_secreted) + ' secretome and ' +
                 '{0:,}'.format(num_mem) + ' transmembane annotations added')

    # interproscan
    IPRCombined = os.path.join(outputdir, 'annotate_misc', 'iprscan.xml')
    IPR_terms = os.path.join(outputdir, 'annotate_misc',
                             'annotations.iprscan.txt')
    if args.iprscan and args.iprscan != IPRCombined:
        if os.path.isfile(IPRCombined):
            os.remove(IPRCombined)
        shutil.copyfile(args.iprscan, IPRCombined)
    if not lib.checkannotations(IPRCombined):
        lib.log.error(
            "InterProScan error, %s is empty, or no XML file passed via --iprscan. Functional annotation will be lacking." % IPRCombined)
    else:
        if os.path.isfile(IPR_terms):
            if os.path.getmtime(IPR_terms) < os.path.getmtime(IPRCombined):
                os.remove(IPR_terms)
        if not lib.checkannotations(IPR_terms):
            lib.log.info("Parsing InterProScan5 XML file")
            cmd = [sys.executable, IPR2ANNOTATE, IPRCombined, IPR_terms]
            lib.runSubprocess(cmd, '.', lib.log)

    # check if antiSMASH data is given, if so parse and reformat for annotations and cluster textual output
    antismash_input = os.path.join(
        outputdir, 'annotate_misc', 'antiSMASH.results.gbk')
    if args.antismash:
        if os.path.isfile(antismash_input):
            os.remove(antismash_input)
        shutil.copyfile(args.antismash, antismash_input)
    if lib.checkannotations(antismash_input):  # result found
        AntiSmashFolder = os.path.join(outputdir, 'annotate_misc', 'antismash')
        AntiSmashBed = os.path.join(AntiSmashFolder, 'clusters.bed')
        GFF2clusters = os.path.join(AntiSmashFolder, 'secmet.clusters.txt')
        AntiSmash_annotations = os.path.join(
            outputdir, 'annotate_misc', 'annotations.antismash.txt')
        Cluster_annotations = os.path.join(
            outputdir, 'annotate_misc', 'annotations.antismash.clusters.txt')
        if os.path.isdir(AntiSmashFolder):
            shutil.rmtree(AntiSmashFolder)
        os.makedirs(AntiSmashFolder)
        # results in several dictionaries
        bbDomains, bbSubType, BackBone = lib.ParseAntiSmash(antismash_input,
                                                            AntiSmashFolder,
                                                            AntiSmashBed,
                                                            AntiSmash_annotations)
        # results in dictClusters dictionary
        dictClusters = lib.GetClusterGenes(AntiSmashBed, GFF, Scaffolds,
                                           Cluster_annotations)

    # if custom annotations passed, parse here
    '''
    if args.annotations:
        lib.log.info("Parsing custom annotations from %s" % args.annotations)
        shutil.copyfile(args.annotations, os.path.join(
            outputdir, 'annotate_misc', 'annotations.custom.txt'))
        num_annotations = lib.line_count(os.path.join(
            outputdir, 'annotate_misc', 'annotations.custom.txt'))
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    '''
    # now bring all annotations together and annotated genome using gag, remove any duplicate annotations
    ANNOTS = os.path.join(outputdir, 'annotate_misc', 'all.annotations.txt')
    GeneNames = lib.getGeneBasename(Proteins)
    total_annotations = 0
    filtered_annotations = 0
    lines_seen = set()
    with open(ANNOTS, 'w') as output:
        for file in os.listdir(os.path.join(outputdir, 'annotate_misc')):
            if file.startswith('annotations'):
                file = os.path.join(outputdir, 'annotate_misc', file)
                with open(file) as input:
                    for line in input:
                        total_annotations += 1
                        if not line.startswith(tuple(GeneNames)):
                            continue
                        if line.count('\t') != 2:  # make sure it is 3 columns
                            continue
                        if line not in lines_seen:
                            output.write(line)
                            lines_seen.add(line)
                            filtered_annotations += 1
    ANNOTS = os.path.abspath(ANNOTS)
    diff_annotations = total_annotations - filtered_annotations
    lib.log.info("Found " + '{0:,}'.format(diff_annotations) + " duplicated annotations, adding " +
                 '{0:,}'.format(filtered_annotations) + ' valid annotations')

    # setup tbl2asn folder
    if os.path.isdir(os.path.join(outputdir, 'annotate_misc', 'tbl2asn')):
        lib.SafeRemove(os.path.join(outputdir, 'annotate_misc', 'tbl2asn'))
    os.makedirs(os.path.join(outputdir, 'annotate_misc', 'tbl2asn'))
    TBLOUT = os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.tbl')
    shutil.copyfile(Scaffolds, os.path.join(
        outputdir, 'annotate_misc', 'tbl2asn', 'genome.fsa'))

    # add annotation to tbl annotation file, generate dictionary of dictionaries with values as a list
    # need to keep multiple transcripts annotations separate, so this approach may have to modified
    # custom annotations take precedence so parse differently
    Annotations = lib.annotations2dict(ANNOTS, geneDB=GeneDB,
                                       custom=args.annotations)

    # to update annotations, user can pass --fix or --remove, update Annotations here
    if args.fix:
        with open(args.fix, 'r') as fixfile:
            for line in fixfile:
                line = line.strip()
                if line.startswith('#'):
                    continue
                # ID Name Description Error (could be more columns i guess)
                cols = line.split('\t')
                if len(cols) < 3:  # skip if number of columns isn't correct
                    continue
                if cols[0] in Annotations:
                    Annotations[cols[0]]['name'] = [cols[1]]
                    Annotations[cols[0]]['product'] = [cols[2]]
                if cols[1] in NotInCurated:
                    NotInCurated[cols[1]] = [cols[2]]
                if cols[1] in NeedCurating:
                    old = NeedCurating.get(cols[1])
                    NeedCurating[cols[1]] = (old[0], cols[2])
                if cols[0] in Gene2ProdFinal:
                    Gene2ProdFinal[cols[0]] = (cols[1], cols[2])

    if args.remove:
        with open(args.remove, 'r') as removefile:
            for line in removefile:
                line = line.strip()
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                if cols[0] in Annotations:
                    if 'name' in Annotations[cols[0]]:
                        del Annotations[cols[0]]['name']
                    if 'product' in Annotations[cols[0]]:
                        del Annotations[cols[0]]['product']
                if cols[0] in Gene2ProdFinal:
                    del Gene2ProdFinal[cols[0]]

    # grab some info from the annotation dictionary
    IPRterms = []
    NoteHeaders = []
    for k, v in natsorted(Annotations.items()):
        if 'note' in v:
            for x in v['note']:
                if ':' in x:
                    h = x.split(':', 1)[0]
                    if h.startswith('SMCOG'):
                        continue
                    if h not in NoteHeaders:
                        NoteHeaders.append(h)
        elif 'db_xref' in v:
            for y in v['db_xref']:
                if y.startswith('InterPro'):
                    g = y.split(':', 1)[1]
                    if not g in IPRterms:
                        IPRterms.append(g)
    NoteHeaders = natsorted(NoteHeaders)


    # now parse tbl file and add annotations
    if args.rename and '_' in args.rename:
        args.rename = args.rename.split('_')[0]
    lib.updateTBL(annotTBL, Annotations, TBLOUT, prefix=locusTagPrefix,
                  newtag=args.rename)

    # if this is reannotation, then need to fix tbl file to track gene changes
    if WGS_accession:
        shutil.copyfile(os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.tbl'),
                        os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.tbl.bak'))
        p2g = {}
        # see if p2g file is present
        p2gfile = None
        if os.path.isfile(os.path.join(outputdir, 'update_results', 'ncbi.p2g')):
            p2gfile = os.path.join(outputdir, 'update_results', 'ncbi.p2g')
        else:
            if args.p2g:
                p2gfile = args.p2g
        if p2gfile:
            with open(p2gfile, 'r') as input:
                for line in input:
                    cols = line.split('\t')
                    if not cols[0] in p2g:
                        p2g[cols[0]] = cols[1]
            with open(os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.tbl'), 'w') as outfile:
                with open(os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.tbl.bak'), 'r') as infile:
                    for line in infile:
                        line = line.replace('\n', '')
                        if line.startswith('\t\t\tprotein_id') or line.startswith('\t\t\ttranscript_id'):
                            ID = line.rsplit('|', 1)[-1].replace('_mrna', '')
                            type = 'prot'
                            if 'transcript_id' in line:
                                type = 'transcript'
                            if not ID in p2g:
                                if type == 'prot':
                                    outfile.write(
                                        '\t\t\tprotein_id\tgnl|%s|%s\n' % (WGS_accession, ID))
                                elif type == 'transcript':
                                    outfile.write(
                                        '\t\t\ttranscript_id\tgnl|%s|%s_mrna\n' % (WGS_accession, ID))
                            else:
                                p2gID = p2g.get(ID)
                                if type == 'prot':
                                    outfile.write('\t\t\tprotein_id\tgnl|%s|%s|gb|%s\n' % (
                                        WGS_accession, ID, p2gID))
                                elif type == 'transcript':
                                    outfile.write(
                                        '\t\t\ttranscript_id\tgnl|%s|%s_mrna\n' % (WGS_accession, ID))
                        else:
                            outfile.write('%s\n' % line)
        else:
            lib.log.error(
                "Detected NCBI reannotation, but couldn't locate p2g file, please pass via --p2g")
            shutil.copyfile(os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.tbl.bak'),
                            os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.tbl'))

    # launch tbl2asn to create genbank submission files
    discrep = 'discrepency.report.txt'
    lib.log.info("Converting to final Genbank format, good luck!")
    if not version:
        annot_version = 1
    else:
        annot_version = version
    # have to run as subprocess because of multiprocessing issues
    cmd = [sys.executable,
           os.path.join(parentdir, 'aux_scripts', 'tbl2asn_parallel.py'),
           '-i', TBLOUT,
           '-f', os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.fsa'),
           '-o', os.path.join(outputdir, 'annotate_misc', 'tbl2asn'),
           '--sbt', SBT, '-d', discrep,
           '-s', organism, '-t', args.tbl2asn,
           '-v', str(annot_version), '-c', str(args.cpus)]
    if args.isolate:
        cmd += ['--isolate', args.isolate]
    if args.strain:
        cmd += ['--strain', args.strain]
    lib.log.debug(' '.join(cmd))
    subprocess.call(cmd)
    # check if completed succesfully
    if not lib.checkannotations(os.path.join(outputdir, 'annotate_misc', 'tbl2asn', 'genome.gbf')):
        lib.log.info('ERROR: GBK file conversion failed, tbl2asn parallel script has died')
        sys.exit(1)

    # parse discrepancy report to see which names/product descriptions failed/passed
    # return dict containing tuples of (GeneName, GeneProduct, [reason])
    BadProducts = []
    if os.path.isfile(discrep) and os.path.exists(discrep):
        BadProducts = lib.getFailedProductNames(discrep, Gene2ProdFinal)

    Gene2ProductPassed = os.path.join(
        outputdir, 'annotate_results', 'Gene2Products.new-names-passed.txt')
    PassedCounts = 0
    with open(Gene2ProductPassed, 'w') as prodpassed:
        prodpassed.write('#Name\tPassed Description\n')
        for key, value in natsorted(list(NotInCurated.items())):
            if not key in BadProducts and not key in NeedCurating:
                PassedCounts += 1
                prodpassed.write('%s\t%s\n' % (key, value[0][1]))
    Gene2ProductHelp = os.path.join(
        outputdir, 'annotate_results', 'Gene2Products.need-curating.txt')
    MustFixHelp = os.path.join(
        outputdir, 'annotate_results', 'Gene2Products.must-fix.txt')
    CurateCount = 0
    MustFixCount = 0
    with open(Gene2ProductHelp, 'w') as needhelp:
        needhelp.write(
            '#Name\tOriginal Description\tCleaned Description\tError-message\n')
        for key, value in natsorted(list(NeedCurating.items())):
            CurateCount += 1
            needhelp.write('%s\t%s\t%s\tProduct defline failed funannotate checks\n' % (
                key, value[0][0], value[0][1]))
    with open(MustFixHelp, 'w') as musthelp:
        musthelp.write('#GeneID\tName\tProduct Description\ttbl2asn Error\n')
        if BadProducts:
            for key, value in natsorted(list(BadProducts.items())):
                MustFixCount += 1
                musthelp.write('%s\t%s\t%s\t%s\n' %
                               (value[1], key, value[0], ', '.join(value[2])))

    # collected output files and rename accordingly
    ResultsFolder = os.path.join(outputdir, 'annotate_results')
    if os.path.exists(discrep) and os.path.isfile(discrep):
        shutil.copyfile(discrep, os.path.join(ResultsFolder,
                    organism_name+'.discrepency.report.txt'))
        os.remove(discrep)
    else:
        lib.log.error('no discrepency file %s found'%(discrep))

    final_tbl = os.path.join(ResultsFolder, organism_name+'.tbl')
    final_gbk = os.path.join(ResultsFolder, organism_name+'.gbk')
    final_gff = os.path.join(ResultsFolder, organism_name+'.gff3')
    final_proteins = os.path.join(ResultsFolder, organism_name+'.proteins.fa')
    final_transcripts = os.path.join(
        ResultsFolder, organism_name+'.mrna-transcripts.fa')
    final_cds_transcripts = os.path.join(
        ResultsFolder, organism_name+'.cds-transcripts.fa')
    final_fasta = os.path.join(ResultsFolder, organism_name+'.scaffolds.fa')
    final_annotation = os.path.join(
        ResultsFolder, organism_name+'.annotations.txt')
    final_stats = os.path.join(ResultsFolder, organism_name+'.stats.json')
    shutil.copyfile(os.path.join(outputdir, 'annotate_misc',
                    'tbl2asn', 'genome.gbf'), final_gbk)
    shutil.copyfile(os.path.join(outputdir, 'annotate_misc',
                    'tbl2asn', 'genome.tbl'), final_tbl)
    # because of possible splitting tbl2asn output, loop through and get sqn and tbl parts
    for file in os.listdir(os.path.join(outputdir, 'annotate_misc', 'tbl2asn')):
        if file.endswith('.sqn') or file.endswith('.tbl'):
            if 'genome.' in file:
                updatedName = file.replace('genome', organism_name)
            else:
                updatedName = file.replace('genome', organism_name+'.part_')
            shutil.copyfile(os.path.join(outputdir, 'annotate_misc',
                                         'tbl2asn', file), os.path.join(ResultsFolder, updatedName))
    lib.tbl2allout(final_tbl, Scaffolds, final_gff, final_proteins,
                   final_transcripts, final_cds_transcripts, final_fasta)

    lib.annotation_summary(Scaffolds, final_stats, tbl=final_tbl,
                           previous=existingStats, database=FUNDB,
                           command=' '.join(sys.argv),
                           organism=organism_name)

    # write AGP output so all files in correct directory
    lib.log.info("Creating AGP file and corresponding contigs file")
    agp2fasta = os.path.join(parentdir, 'aux_scripts', 'fasta2agp.pl')
    AGP = os.path.join(ResultsFolder, organism_name+'.agp')
    cmd = ['perl', agp2fasta, organism_name+'.scaffolds.fa']
    lib.runSubprocess2(cmd, ResultsFolder, lib.log, AGP)

    # write secondary metabolite clusters output using the final genome in gbk format
    if lib.checkannotations(antismash_input):
        lib.log.info(
            "Cross referencing SM cluster hits with MIBiG database version %s" % versDB.get('mibig'))
        # do a blast best hit search against MIBiG database for cluster annotation, but looping through gene cluster hits
        AllProts = []
        SMgenes = []
        for k, v in list(dictClusters.items()):
            for i in v:
                if '-T' in i:
                    ID = i.split('-T')[0]
                else:
                    ID = i
                if not i in AllProts:
                    AllProts.append(i)
                if not ID in SMgenes:
                    SMgenes.append(ID)
        AllProts = set(AllProts)
        mibig_fasta = os.path.join(AntiSmashFolder, 'smcluster.proteins.fasta')
        mibig_blast = os.path.join(
            AntiSmashFolder, 'smcluster.MIBiG.blast.txt')
        mibig_db = os.path.join(FUNDB, 'mibig.dmnd')
        with open(mibig_fasta, 'w') as output:
            with open(Proteins, 'r') as input:
                SeqRecords = SeqIO.parse(Proteins, 'fasta')
                for record in SeqRecords:
                    genename = record.id
                    if genename in AllProts:
                        SeqIO.write(record, output, 'fasta')
        cmd = ['diamond', 'blastp', '--sensitive', '--query', mibig_fasta,
               '--threads', str(args.cpus), '--out', mibig_blast,
               '--db', mibig_db, '--max-hsps', '1',
               '--evalue', '0.001', '--max-target-seqs', '1',
               '--outfmt', '6']
        lib.runSubprocess4(cmd, '.', lib.log)
        # now parse blast results to get {qseqid: hit}
        MIBiGBlast = {}
        with open(mibig_blast, 'r') as input:
            for line in input:
                cols = line.split('\t')
                if '-T' in cols[0]:
                    ID = cols[0].split('-T')[0]
                else:
                    ID = cols[0]
                hit = cols[1].split('|')
                desc = hit[5]
                cluster = hit[0]
                db_ref = hit[6]
                evalue = cols[10]
                pident = cols[2]
                result = (desc, cluster, db_ref, pident, evalue)
                MIBiGBlast[ID] = result

        lib.log.info("Creating tab-delimited SM cluster output")

        # load in antismash cluster bed file to slice record
        slicing = []
        with open(AntiSmashBed, 'r') as antibed:
            for line in antibed:
                cols = line.split('\t')
                # chr, cluster, start, stop in a tuple
                cluster = (cols[0], cols[3], cols[1], cols[2])
                slicing.append(cluster)
        Offset = {}
        # Get each cluster + 15 Kb in each direction to make sure you can see the context of the cluster
        with open(os.path.join(ResultsFolder, organism_name+'.gbk'), 'r') as gbk:
            SeqRecords = SeqIO.parse(gbk, 'genbank')
            for record in SeqRecords:
                for f in record.features:
                    if f.type == "source":
                        record_end = f.location.end
                for slice in slicing:
                    if record.id == slice[0]:
                        sub_start = int(slice[2]) - 15000
                        sub_stop = int(slice[3]) + 15000
                        if sub_start < 1:
                            sub_start = 1
                        if sub_stop > record_end:
                            sub_stop = record_end
                        sub_record = record[sub_start:sub_stop]
                        # this seems to be either py3 requirement or required in newer biopython
                        sub_record.annotations = record.annotations
                        cluster_name = slice[1]
                        sub_record_name = os.path.join(
                            AntiSmashFolder, cluster_name+'.gbk')
                        Offset[cluster_name] = sub_start
                        with open(sub_record_name, 'w') as clusterout:
                            try:
                                SeqIO.write(sub_record, clusterout, 'genbank')
                            except ValueError:
                                print(slice)
                                print(subrecord.id)
                                print(sub_record.annotations)
                                sys.exit(1)

        # okay, now loop through each cluster
        for file in os.listdir(AntiSmashFolder):
            if file.endswith('.gbk'):
                base = file.replace('.gbk', '')
                outputName = os.path.join(
                    AntiSmashFolder, base+'.secmet.cluster.txt')
                file = os.path.join(AntiSmashFolder, file)
                with open(outputName, 'w') as output:
                    output.write("#%s\n" % base)
                    output.write(
                        "#GeneID\tChromosome:start-stop\tStrand\tClusterPred\tBackbone Enzyme\tBackbone Domains\tProduct\tsmCOGs\tEggNog\tInterPro\tPFAM\tGO terms\tNotes\tMIBiG Blast\tProtein Seq\tDNA Seq\n")
                    with open(file, 'r') as input:
                        SeqRecords = SeqIO.parse(input, 'genbank')
                        for record in SeqRecords:
                            for f in record.features:
                                if f.type == "CDS":
                                    name = f.qualifiers["locus_tag"][0]
                                    prot_seq = f.qualifiers['translation'][0]
                                    start = f.location.nofuzzy_start
                                    # account for python numbering shift?
                                    actualStart = int(
                                        start) + int(Offset.get(base)) + 1
                                    end = f.location.nofuzzy_end
                                    actualEnd = int(end) + \
                                        int(Offset.get(base))
                                    strand = f.location.strand
                                    if strand == 1:
                                        strand = '+'
                                        DNA_seq = record.seq[start:end]
                                    elif strand == -1:
                                        strand = '-'
                                        DNA_seq = record.seq[start:end].reverse_complement(
                                        )
                                    chr = record.id
                                    product = f.qualifiers["product"][0]
                                    # now get the info out of the note and db_xref fields, need to clear each field for each record
                                    note = []
                                    goTerms = []
                                    pFAM = []
                                    IPR = []
                                    eggnogDesc = 'NA'
                                    if name in SMgenes:
                                        location = 'cluster'
                                    else:
                                        location = 'flanking'
                                    cog = '.'
                                    for k, v in list(f.qualifiers.items()):
                                        if k == 'note':
                                            # multiple notes are split with a semi colon
                                            items = v[0].split('; ')
                                            for i in items:
                                                if i.startswith('EggNog:'):
                                                    eggnogID = i.replace(
                                                        'EggNog:', '')
                                                    eggnogDesc = EggNog.get(
                                                        eggnogID)
                                                elif i.startswith('GO_'):
                                                    goterm = i.split(
                                                        ': ', 1)[-1]
                                                    goTerms.append(goterm)
                                                elif i.startswith('SMCOG'):
                                                    cog = i
                                                else:
                                                    note.append(i)
                                        if k == 'db_xref':
                                            for i in v:
                                                if i.startswith('InterPro:'):
                                                    r = i.replace(
                                                        'InterPro:', '')
                                                    IPR.append(r)
                                                if i.startswith('PFAM:'):
                                                    p = i.replace('PFAM:', '')
                                                    pFAM.append(p)
                                    if name in bbDomains:
                                        domains = ";".join(bbDomains.get(name))
                                    else:
                                        domains = '.'
                                    if name in bbSubType:
                                        enzyme = bbSubType.get(name)
                                    else:
                                        if name in BackBone:
                                            enzyme = BackBone.get(name)
                                        else:
                                            enzyme = '.'
                                    if name in MIBiGBlast:
                                        mibigTup = MIBiGBlast.get(name)
                                        mibig = mibigTup[0]+' from '+mibigTup[1] + \
                                            ' ('+mibigTup[2]+':pident=' + \
                                            mibigTup[3]+', evalue=' + \
                                            mibigTup[4]+')'
                                        mibig = str(mibig)
                                    else:
                                        mibig = '.'
                                    if IPR:
                                        IP = ";".join(IPR)
                                    else:
                                        IP = '.'
                                    if pFAM:
                                        PF = ";".join(pFAM)
                                    else:
                                        PF = '.'
                                    if goTerms:
                                        GO = ";".join(goTerms)
                                    else:
                                        GO = '.'
                                    if note:
                                        No = ";".join(note)
                                    else:
                                        No = '.'
                                    output.write("%s\t%s:%i-%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, chr, actualStart,
                                                                                                                             actualEnd, strand, location, enzyme, domains, product, cog, eggnogDesc, IP, PF, GO, No, mibig, prot_seq, DNA_seq))

        # now put together into a single file
        finallist = []
        ClustersOut = os.path.join(
            ResultsFolder, organism_name+'.clusters.txt')
        for file in os.listdir(AntiSmashFolder):
            if file.endswith('secmet.cluster.txt'):
                file = os.path.join(AntiSmashFolder, file)
                finallist.append(file)
        with open(ClustersOut, 'w') as output:
            for file in natsorted(finallist):
                with open(file, 'r') as input:
                    output.write(input.read())
                    output.write('\n\n')

    # write tsv annotation table
    lib.log.info("Writing genome annotation table.")
    #lib.annotationtable(final_gbk, FUNDB, final_annotation)
    INTERPRO = lib.iprTSV2dict(os.path.join(FUNDB, 'interpro.tsv'), IPRterms)
    lib.annotationtable(final_gbk, FUNDB, NoteHeaders, INTERPRO,
                        final_annotation)

    # final wrap up message
    if MustFixCount == 0 and PassedCounts == 0 and CurateCount == 0:
        lib.log.info("Funannotate annotate has completed successfully!")
    else:
        lib.log.info("Funannotate annotate has completed successfully!\n\n\
        We need YOUR help to improve gene names/product descriptions:\n\
           {:,} gene/products names MUST be fixed, see {:}\n\
           {:,} gene/product names need to be curated, see {:}\n\
           {:,} gene/product names passed but are not in Database, see {:}\n\n\
        Please consider contributing a PR at https://github.com/nextgenusfs/gene2product\n".format(MustFixCount, MustFixHelp, CurateCount, Gene2ProductHelp, PassedCounts, Gene2ProductPassed))

    if MustFixCount > 0:  # show user how to update
        lib.log.info("To fix gene names/product deflines, manually fix or can remove in {:}\n\n\
  funannotate annotate -i {:} --fix fixed_file.txt --remove delete.txt\n".format(MustFixHelp, args.input))
    print("-------------------------------------------------------")
    # move logfile to logfiles directory
    if os.path.isfile(log_name):
        if not os.path.isdir(os.path.join(outputdir, 'logfiles')):
            os.makedirs(os.path.join(outputdir, 'logfiles'))
        shutil.copyfile(log_name, os.path.join(
            outputdir, 'logfiles', 'funannotate-annotate.log'))
        os.remove(log_name)


if __name__ == "__main__":
    main(sys.argv[1:])
