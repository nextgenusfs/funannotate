#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import hashlib
import os
import uuid
import subprocess
import requests
import argparse
import shutil


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
    if not os.path.isfile(input):
        # The command exited before it wrote its output; say so, and let the
        # caller's count check fail instead of raising FileNotFoundError.
        print('ERROR: expected output file not found: {}'.format(input))
        return 0
    with open(input, 'r') as f:
        for line in f:
            if "\tgene\t" in line:
                count += 1
    return count


def runCMD(cmd, dir):
    """Run a command in dir and return its exit status.

    A failed command is reported here, so that a later check does not have to explain a
    missing output file.
    """
    print(('CMD: {:}'.format(' '.join(cmd))))
    print("#########################################################")
    rc = subprocess.call(cmd, cwd=dir)
    if rc != 0:
        print('ERROR: command exited with status {:}: {:}'.format(rc, ' '.join(cmd[:3])))
    return rc


# Test data is mirrored in the `test-data-v1` GitHub release (byte-identical to the original
# OSF files, https://osf.io/bj7v4/). The OSF link is the fallback. The sha256 values are in
# the release's SHA256SUMS file.
TEST_DATA_URL = 'https://github.com/nextgenusfs/funannotate/releases/download/test-data-v1/{}'
TEST_DATA_LINKS = {'mask': 'test-mask.tar.gz',
                   'clean': 'test-clean.tar.gz',
                   'predict': 'test-predict.tar.gz',
                   'busco': 'test-busco.tar.gz',
                   'rna-seq': 'test-rna-seq.tar.gz',
                   'annotate': 'test-annotate.tar.gz',
                   'compare': 'test-compare.tar.gz'}
TEST_DATA_OSF = {'mask': 'https://osf.io/hbryz/download?version=1',
                 'clean': 'https://osf.io/8pjbe/download?version=1',
                 'predict': 'https://osf.io/te2pf/download?version=1',
                 'busco': 'https://osf.io/kyrd9/download?version=1',
                 'rna-seq': 'https://osf.io/t7j83/download?version=1',
                 'annotate': 'https://osf.io/97pyn/download?version=1',
                 'compare': 'https://osf.io/7s9xh/download?version=1'}
TEST_DATA_SHA256 = {'test-mask.tar.gz': 'fc0c465ebe08822b2764cb3077c6ee19d878903719cb07979bb6a78f42693824',
                    'test-clean.tar.gz': 'ec92f67010b2720d0181c2fc0410666b2982bf87aaa6338eb25ebc2ff47b9893',
                    'test-predict.tar.gz': 'caaba07e407520928fcc6d2bbcdecb60d3309970997912b671acded770742261',
                    'test-busco.tar.gz': '22fd0032c4a309dca2f9332e4de0917314f54ae8183ffcb71ef0e9ec06e26d91',
                    'test-rna-seq.tar.gz': '65e8f3c3cacc6e684294fbfccff7f3c0eeb3bfbd499fd02785dc435085c45baa',
                    'test-annotate.tar.gz': '2f2765b6def5466d7b6355b24cfecc8ecec1bb658d1982887cafcfd032cce45c',
                    'test-compare.tar.gz': '3b23491f06d03c1570d266f223421b93477275afca3eb8fa8c8068f96e24c902'}


def sha256sum(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def download_test_data(key, name):
    """Download a test archive from the GitHub release, else from OSF, and verify its sha256.

    `key` selects the entry in TEST_DATA_LINKS and TEST_DATA_OSF; `name` is the local file.
    The archive's checksum is looked up by its release file name.
    """
    release_name = TEST_DATA_LINKS[key]
    expected = TEST_DATA_SHA256.get(release_name)
    for url in (TEST_DATA_URL.format(release_name), TEST_DATA_OSF[key]):
        try:
            download(url, name)
        except requests.exceptions.RequestException:
            continue
        if expected and sha256sum(name) != expected:
            print('ERROR: {} from {} does not match the expected sha256'.format(name, url))
            os.remove(name)
            continue
        return
    raise SystemExit('ERROR: could not download test data {}'.format(release_name))


def download(url, name):
    try:
        r = requests.get(url, allow_redirects=True, stream=True)
        r.raise_for_status()
        file_size = int(r.headers.get('content-length', 0))
        print("Downloading: {} Bytes: {}".format(url, file_size))
        file_size_dl = 0
        with open(name, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    file_size_dl += len(chunk)
                    if file_size > 0:
                        p = float(file_size_dl) / file_size
                        status = r"{}  [{:.2%}]".format(file_size_dl, p)
                        status = status + chr(8)*(len(status)+1)
                        sys.stdout.write(status)
        sys.stdout.flush()
    except requests.exceptions.RequestException as e:
        print("Download failed: {}".format(e))
        raise


def runMaskTest(args):
    print("#########################################################")
    print('Running `funannotate mask` unit testing: RepeatModeler --> RepeatMasker')
    tmpdir = 'test-mask_'+pid
    os.makedirs(tmpdir)
    inputFasta = 'test.fa'
    if not os.path.isfile(inputFasta):
        if not os.path.isfile('test-mask.tar.gz'):
            download_test_data('mask', 'test-mask.tar.gz')
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
            download_test_data('clean', 'test-clean.tar.gz')
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
            download_test_data('predict', 'test-predict.tar.gz')
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
            download_test_data('busco', 'test-busco.tar.gz')
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
            download_test_data('annotate', 'test-annotate.tar.gz')
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
            download_test_data('compare', 'test-compare.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-compare.tar.gz'])
    shutil.copyfile(input1, os.path.join(tmpdir, input1))
    shutil.copyfile(input2, os.path.join(tmpdir, input2))
    shutil.copyfile(input3, os.path.join(tmpdir, input3))
    # run predict
    runCMD(['funannotate', 'compare',
            '-i', input1, input2, input3,
            '-o', 'compare', '--cpus', str(args.cpus), '--ml_model', 'LG+G4',
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
            download_test_data('rna-seq', 'test-rna_seq.tar.gz')
        subprocess.call(['tar', '-zxf', 'test-rna_seq.tar.gz'])
    for f in [inputFasta, protEvidence, illumina, nanopore]:
        shutil.copyfile(f, os.path.join(tmpdir, f))
    # run train
    rc = runCMD(['funannotate', 'train', '-i', inputFasta,
                 '--single', illumina, '--nanopore_mrna', nanopore,
                 '-o', 'rna-seq', '--cpus', str(args.cpus), '--jaccard_clip',
                 '--species', "Awesome rna"], tmpdir)
    if rc != 0:
        # predict and update need the training data; running them would only hide this error
        print('ERROR: funannotate RNA-seq test failed - `funannotate train` exited with status {:}, '
              'so predict and update were skipped. See rna-seq/logfiles in {:}'.format(rc, tmpdir))
        print("#########################################################\n")
        return
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
    parser.add_argument('-c','--cpus', default=2, type=int,
                        help='Number of CPUs to use')
    args = parser.parse_args(args)

    global download_links
    download_links = TEST_DATA_LINKS
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
