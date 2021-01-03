#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import uuid
import subprocess
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
import socket
import argparse
import shutil
import errno


def checkFile(input):
    def _getSize(filename):
        st = os.stat(filename)
        return st.st_size
    if os.path.isfile(input):
        filesize = _getSize(input)
        if int(filesize) < 1:
            return False
        else:
            return True
    elif os.path.islink(input):
        return True
    else:
        return False


def countfasta(input):
    count = 0
    with open(input, 'r') as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def countGFFgenes(input):
    count = 0
    with open(input, 'r') as f:
        for line in f:
            if "\tgene\t" in line:
                count += 1
    return count


def runCMD(cmd, dir):
    print(('CMD: {:}'.format(' '.join(cmd))))
    print("#########################################################")
    subprocess.call(cmd, cwd=dir)


def download(url, name):
    file_name = name
    try:
        u = urlopen(url)
        f = open(file_name, 'wb')
        meta = u.info()
        file_size = 0
        for x in meta.items():
            if x[0].lower() == 'content-length':
                file_size = int(x[1])
        print("Downloading: {} Bytes: {}".format(url, file_size))
        file_size_dl = 0
        block_sz = 8192
        while True:
            buffer = u.read(block_sz)
            if not buffer:
                break
            file_size_dl += len(buffer)
            f.write(buffer)
            p = float(file_size_dl) / file_size
            status = r"{}  [{:.2%}]".format(file_size_dl, p)
            status = status + chr(8)*(len(status)+1)
            sys.stdout.write(status)
        sys.stdout.flush()
        f.close()
    except socket.error as e:
        if e.errno != errno.ECONNRESET:
            raise
        pass


def runMaskTest(args):
    print("#########################################################")
    print('Running `funannotate mask` unit testing: RepeatModeler --> RepeatMasker')
    tmpdir = 'test-mask_'+pid
    os.makedirs(tmpdir)
    inputFasta = 'test.fa'
    if not os.path.isfile(inputFasta):
        if not os.path.isfile('test-mask.tar.gz'):
            download(download_links.get('mask'), 'test-mask.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-mask.tar.gz'])
    shutil.copyfile(inputFasta, os.path.join(tmpdir, inputFasta))
    runCMD(['funannotate', 'mask', '-i', inputFasta, '-o',
            'test.masked.fa', '--cpus', str(args.cpus)], tmpdir)
    # check that everything worked
    print("#########################################################")
    try:
        assert checkFile(os.path.join(tmpdir, 'test.masked.fa'))
        print('SUCCESS: `funannotate mask` test complete.')
        shutil.rmtree(tmpdir)
    except AssertionError:
        print('ERROR: `funannotate mask` test failed.')
    print("#########################################################\n")


def runCleanTest(args):
    print("#########################################################")
    print('Running `funannotate clean` unit testing: minimap2 mediated assembly duplications')
    tmpdir = 'test-clean_'+pid
    os.makedirs(tmpdir)
    inputFasta = 'test.clean.fa'
    if not os.path.isfile(inputFasta):
        if not os.path.isfile('test-clean.tar.gz'):
            download(download_links.get('clean'), 'test-clean.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-clean.tar.gz'])
    shutil.copyfile(inputFasta, os.path.join(tmpdir, inputFasta))
    assert countfasta(os.path.join(tmpdir, inputFasta)) == 6
    # run exhaustive
    runCMD(['funannotate', 'clean', '-i', inputFasta, '-o',
            'test.exhaustive.fa', '--exhaustive'], tmpdir)
    print("#########################################################")
    try:
        assert countfasta(os.path.join(tmpdir, 'test.exhaustive.fa')) == 3
        print('SUCCESS: `funannotate clean` test complete.')
        if not args.debug:
            shutil.rmtree(tmpdir)
    except AssertionError:
        print('ERROR: `funannotate clean` test failed.')
    print("#########################################################\n")


def runPredictTest(args):
    print("#########################################################")
    print('Running `funannotate predict` unit testing')
    tmpdir = 'test-predict_'+pid
    os.makedirs(tmpdir)
    inputFasta = 'test.softmasked.fa'
    protEvidence = 'protein.evidence.fasta'
    if not checkFile(inputFasta) or not checkFile(protEvidence):
        if not os.path.isfile('test-predict.tar.gz'):
            download(download_links.get('predict'), 'test-predict.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-predict.tar.gz'])
    shutil.copyfile(inputFasta, os.path.join(tmpdir, inputFasta))
    shutil.copyfile(protEvidence, os.path.join(tmpdir, protEvidence))
    # run predict
    runCMD(['funannotate', 'predict', '-i', inputFasta,
            '--protein_evidence', protEvidence,
            '-o', 'annotate', '--augustus_species', 'saccharomyces',
            '--cpus', str(args.cpus), '--species', "Awesome testicus"], tmpdir)
    print("#########################################################")
    # check results
    try:
        assert 1500 <= countGFFgenes(os.path.join(
            tmpdir, 'annotate', 'predict_results', 'Awesome_testicus.gff3')) <= 1800
        print('SUCCESS: `funannotate predict` test complete.')
        if not args.debug:
            shutil.rmtree(tmpdir)
    except AssertionError:
        print('ERROR: `funannotate predict` test failed - check logfiles')
    print("#########################################################\n")


def runBuscoTest(args):
    print("#########################################################")
    print('Running `funannotate predict` BUSCO-mediated training unit testing')
    # need to delete any pre-existing Augustus training data
    try:
        FUNDB = os.environ["FUNANNOTATE_DB"]
    except KeyError:
        print(
            'Funannotate database not properly configured, run funannotate setup.')
        return
    if os.path.isdir(os.path.join(FUNDB, 'trained_species', 'awesome_busco')):
        shutil.rmtree(os.path.join(FUNDB, 'trained_species', 'awesome_busco'))
    tmpdir = 'test-busco_'+pid
    os.makedirs(tmpdir)
    inputFasta = 'test.softmasked.fa'
    protEvidence = 'protein.evidence.fasta'
    if not checkFile(inputFasta) or not checkFile(protEvidence):
        if not os.path.isfile('test-busco.tar.gz'):
            download(download_links.get('predict'), 'test-busco.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-busco.tar.gz'])
    shutil.copyfile(inputFasta, os.path.join(tmpdir, inputFasta))
    shutil.copyfile(protEvidence, os.path.join(tmpdir, protEvidence))
    # run predict
    runCMD(['funannotate', 'predict', '-i', inputFasta,
            '--protein_evidence', protEvidence,
            '-o', 'annotate', '--cpus', str(args.cpus),
            '--species', "Awesome busco"], tmpdir)
    print("#########################################################")
    # check results
    try:
        assert 1500 <= countGFFgenes(os.path.join(
            tmpdir, 'annotate', 'predict_results', 'Awesome_busco.gff3')) <= 1800
        print('SUCCESS: `funannotate predict` BUSCO-mediated training test complete.')
        print("#########################################################")
        # print('Adding training parameters to database')
        # now lets try to add it to the database as new species
        # cmd2 = ['funannotate', 'species', '-s', 'awesome_busco', '-p', 'annotate/predict_results/awesome_busco.parameters.json']
        # runCMD(cmd2, tmpdir)
        # now lets try to run this again, using the parameters
        # print("#########################################################")
        print('Now running predict using all pre-trained ab-initio predictors')
        runCMD(['funannotate', 'predict', '-i', inputFasta,
                '--protein_evidence', protEvidence,
                '-o', 'annotate2', '--cpus', str(args.cpus),
                '--species', "Awesome busco", '-p', 'annotate/predict_results/awesome_busco.parameters.json'],
               tmpdir)
        print("#########################################################")
        try:
            assert 1500 <= countGFFgenes(os.path.join(
                tmpdir, 'annotate2', 'predict_results', 'Awesome_busco.gff3')) <= 1800
            print(
                'SUCCESS: `funannotate predict` using existing parameters test complete.')
            if not args.debug:
                shutil.rmtree(tmpdir)
        except AssertionError:
            print(
                'ERROR: `funannotate predict` using existing parameters test failed - check logfiles')
        # delete the training data from database
        # shutil.rmtree(os.path.join(FUNDB, 'trained_species', 'awesome_busco'))
    except AssertionError:
        print('ERROR: `funannotate predict` BUSCO-mediated training test failed - check logfiles')
    print("#########################################################\n")


def runAnnotateTest(args):
    print("#########################################################")
    print('Running `funannotate annotate` unit testing')
    tmpdir = 'test-annotate_'+pid
    os.makedirs(tmpdir)
    input = 'Genome_one.gbk'
    iprscan = 'genome_one.iprscan.xml'
    emapper = 'genome_one.emapper.annotations'
    if not checkFile(input) or not checkFile(iprscan) or not checkFile(emapper):
        if not os.path.isfile('test-annotate.tar.gz'):
            download(download_links.get('annotate'), 'test-annotate.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-annotate.tar.gz'])
    shutil.copyfile(input, os.path.join(tmpdir, input))
    shutil.copyfile(iprscan, os.path.join(tmpdir, iprscan))
    shutil.copyfile(emapper, os.path.join(tmpdir, emapper))
    # run predict
    runCMD(['funannotate', 'annotate', '--genbank', input,
            '-o', 'annotate', '--cpus', str(args.cpus),
            '--iprscan', iprscan,
            '--eggnog', emapper], tmpdir)
    print("#########################################################")
    # check results
    try:
        assert checkFile(os.path.join(tmpdir, 'annotate',
                                      'annotate_results', 'Genome_one.gbk'))
        assert checkFile(os.path.join(tmpdir, 'annotate',
                                      'annotate_results', 'Genome_one.sqn'))
        assert checkFile(os.path.join(tmpdir, 'annotate',
                                      'annotate_results', 'Genome_one.agp'))
        assert checkFile(os.path.join(tmpdir, 'annotate',
                                      'annotate_results', 'Genome_one.tbl'))
        assert checkFile(os.path.join(tmpdir, 'annotate',
                                      'annotate_results', 'Genome_one.annotations.txt'))
        print('SUCCESS: `funannotate annotate` test complete.')
        if not args.debug:
            shutil.rmtree(tmpdir)
    except AssertionError:
        print('ERROR: `funannotate annotate` test failed - check logfiles')
    print("#########################################################\n")


def runCompareTest(args):
    print("#########################################################")
    print('Running `funannotate compare` unit testing')
    tmpdir = 'test-compare_'+pid
    os.makedirs(tmpdir)
    input1 = 'Genome_one.gbk'
    input2 = 'Genome_two.gbk'
    input3 = 'Genome_three.gbk'
    if not checkFile(input1) or not checkFile(input2) or not checkFile(input3):
        if not os.path.isfile('test-compare.tar.gz'):
            download(download_links.get('compare'), 'test-compare.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-compare.tar.gz'])
    shutil.copyfile(input1, os.path.join(tmpdir, input1))
    shutil.copyfile(input2, os.path.join(tmpdir, input2))
    shutil.copyfile(input3, os.path.join(tmpdir, input3))
    # run predict
    runCMD(['funannotate', 'compare',
            '-i', input1, input2, input3,
            '-o', 'compare', '--cpus', str(args.cpus),
            '--outgroup', 'botrytis_cinerea.dikarya'], tmpdir)
    print("#########################################################")
    # check results
    try:
        assert checkFile(os.path.join(tmpdir, 'compare', 'index.html'))
        assert checkFile(os.path.join(tmpdir, 'compare', 'phylogeny.html'))
        assert checkFile(os.path.join(tmpdir, 'compare.tar.gz'))
        print('SUCCESS: `funannotate compare` test complete.')
        if not args.debug:
            shutil.rmtree(tmpdir)
    except AssertionError:
        print('ERROR: `funannotate compare` test failed - check logfiles')
    print("#########################################################\n")


def runRNAseqTest(args):
    print("#########################################################")
    print('Running funannotate RNA-seq training/prediction unit testing')
    # need to delete any pre-existing Augustus training data
    try:
        AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
    except KeyError:
        print(
            "$AUGUSTUS_CONFIG_PATH environmental variable not found, set to continue.")
        return
    if os.path.isdir(os.path.join(AUGUSTUS, 'species', 'awesome_rna')):
        shutil.rmtree(os.path.join(AUGUSTUS, 'species', 'awesome_rna'))
    tmpdir = 'test-rna_seq_'+pid
    os.makedirs(tmpdir)
    inputFasta = 'test.softmasked.fa'
    protEvidence = 'protein.evidence.fasta'
    illumina = 'rna-seq.illumina.fastq.gz'
    nanopore = 'rna-seq.nanopore.fastq.gz'
    if not checkFile(inputFasta) or not checkFile(protEvidence) or not checkFile(illumina) or not checkFile(nanopore):
        if not os.path.isfile('test-rna_seq.tar.gz'):
            download(download_links.get('rna-seq'), 'test-rna_seq.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-rna_seq.tar.gz'])
    for f in [inputFasta, protEvidence, illumina, nanopore]:
        shutil.copyfile(f, os.path.join(tmpdir, f))
    # run train
    runCMD(['funannotate', 'train', '-i', inputFasta,
            '--single', illumina, '--nanopore_mrna', nanopore,
            '-o', 'rna-seq', '--cpus', str(args.cpus), '--jaccard_clip',
            '--species', "Awesome rna"], tmpdir)
    # run predict
    print("#########################################################")
    print('Now running `funannotate predict` using RNA-seq training data')
    runCMD(['funannotate', 'predict', '-i', inputFasta,
            '--protein_evidence', protEvidence,
            '-o', 'rna-seq', '--cpus', str(
                args.cpus), '--min_training_models', '150',
            '--species', "Awesome rna"], tmpdir)
    try:
        assert 1630 <= countGFFgenes(os.path.join(
            tmpdir, 'rna-seq', 'predict_results', 'Awesome_rna.gff3')) <= 1830
        # run update
        print("#########################################################")
        print('Now running `funannotate update` to run PASA-mediated UTR addition and multiple transcripts')
        runCMD(['funannotate', 'update', '-i', 'rna-seq',
                '--cpus', str(args.cpus)], tmpdir)
        print("#########################################################")
        # check results
        try:
            assert 1630 <= countGFFgenes(os.path.join(
                tmpdir, 'rna-seq', 'update_results', 'Awesome_rna.gff3')) <= 1830
            print('SUCCESS: funannotate RNA-seq training/prediction test complete.')
            if not args.debug:
                shutil.rmtree(tmpdir)
        except AssertionError:
            print(
                'ERROR: funannotate RNA-seq training/prediction test failed - check logfiles')
        print("#########################################################\n")
    except AssertionError:
        print('ERROR: funannotate RNA-seq prediction test failed - check logfiles')
    print("#########################################################\n")


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='funannotate-test.py',
                                     description='''Script to download and then test funannotate installation''',
                                     epilog="""Written by Jon Palmer (2016-2018) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-t', '--tests', required=True, nargs='+',
                        choices=['all', 'clean', 'mask', 'predict',
                                 'busco', 'annotate', 'rna-seq', 'compare'],
                        help='select which tests to run')
    parser.add_argument('--debug', action='store_true',
                        help='keep folder')
    parser.add_argument('--cpus', default=2, type=int,
                        help='Number of CPUs to use')
    args = parser.parse_args(args)

    global download_links
    download_links = {'mask': 'https://osf.io/hbryz/download?version=1',
                      'clean': 'https://osf.io/8pjbe/download?version=1',
                      'predict': 'https://osf.io/te2pf/download?version=1',
                      'busco': 'https://osf.io/kyrd9/download?version=1',
                      'rna-seq': 'https://osf.io/t7j83/download?version=1',
                      'annotate': 'https://osf.io/97pyn/download?version=1',
                      'compare': 'https://osf.io/7s9xh/download?version=1'}
    global pid
    pid = str(uuid.uuid4())
    if 'clean' in args.tests or 'all' in args.tests:
        runCleanTest(args)
    if 'mask' in args.tests or 'all' in args.tests:
        runMaskTest(args)
    if 'predict' in args.tests or 'all' in args.tests:
        runPredictTest(args)
    if 'busco' in args.tests or 'all' in args.tests:
        runBuscoTest(args)
    if 'rna-seq' in args.tests or 'all' in args.tests:
        runRNAseqTest(args)
    if 'annotate' in args.tests or 'all' in args.tests:
        runAnnotateTest(args)
    if 'compare' in args.tests or 'all' in args.tests:
        runCompareTest(args)


if __name__ == "__main__":
    main(sys.argv[1:])
