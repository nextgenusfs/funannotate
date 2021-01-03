#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import uuid
import argparse
import multiprocessing
import subprocess
import time
import shutil
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
import socket
import errno
import funannotate.library as lib
from xml.etree import ElementTree as et


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=50)


parser = argparse.ArgumentParser(
    prog='funannotate-iprscan.py',
    description='''Script to run InterProscan locally or with Docker''',
    epilog="""Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i', '--input', required=True,
                    help='FASTA file or funannotate folder')
parser.add_argument('-n', '--num', type=int, default=1000,
                    help='Number of files in each chunk')
parser.add_argument('-c', '--cpus', default=12, type=int,
                    help='Total number of CPUs')
parser.add_argument('-m', '--method', required=True,
                    choices=['local', 'docker'], help='Method to use')
parser.add_argument('--cpus_per_chunk', default=4, type=int,
                    help='Number of CPUs per IPR instance')
parser.add_argument('--iprscan_path', default='interproscan.sh',
                    help='Local Path to interproscan.sh')
parser.add_argument('-o', '--out', help='Final output XML file')
parser.add_argument('--debug', action='store_true', help='Keep intermediate files')
parser.add_argument('--no-progress', dest='progress', action='store_false',
                    help='no progress on multiprocessing')
args = parser.parse_args()


def combine_xml(files, output):
    first = None
    for filename in files:
        data = et.parse(filename).getroot()
        if first is None:
            first = data
        else:
            first.extend(data.getchildren())
    #generate new tree
    tree = et.ElementTree()
    tree._setroot(first)
    et.register_namespace("", "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5")
    tree.write(output, encoding='utf-8', xml_declaration=True)


def checkDocker():
    try:
        proc = subprocess.Popen(
            ['docker', 'images'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print('Docker is not installed, exiting.')
            sys.exit(1)
    stdout, stderr = proc.communicate()
    if 'Cannot connect' in stderr:
        print('Docker is not running, please launch Docker app and re-run script.')
        sys.exit(1)
    # now check for interproscan images
    check = False
    for line in stdout.split('\n'):
        if line.startswith('blaxterlab/interproscan'):
            check = True
    if not check:
        print('Downloading InterProScan Docker images:')
        subprocess.call(['docker', 'pull', 'blaxterlab/interproscan'])
    print('Docker InterProScan container is ready.')


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


def countfasta(input):
    count = 0
    with open(input, 'r') as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with open(path) as inf:
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos


def return_lines(path, linepos, nstart, nstop):
    """return nsamp lines from path where line offsets are in linepos"""
    offsets = linepos[nstart:nstop]
    lines = []
    with open(path) as inf:
        for offset in offsets:
            inf.seek(offset)
            lines.append(inf.readline())
    return lines


def split_fasta(input, outputdir, chunks):
    # function to return line positions of fasta files for chunking
    fastapos = []
    position = 0
    numseqs = 0
    basename = os.path.basename(input).split('.fa', -1)[0]
    with open(input, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                numseqs += 1
                fastapos.append(position)
            position += 1
    splits = []
    n = int(numseqs / chunks)
    num = 0
    for i in range(chunks):
        if i == 0:
            start = 0
            num = n
            lastpos = fastapos[n+1]
        else:
            start = lastpos
            num = num + n
            try:
                lastpos = fastapos[num+1]  # find the n+1 seq
            except IndexError:
                lastpos = fastapos[-1]
        splits.append((start, lastpos))
    # check if output folder exists, if not create it
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    # get line positions from file
    linepos = scan_linepos(input)
    # loop through the positions and write output
    for i, x in enumerate(splits):
        num = i+1
        with open(os.path.join(outputdir, basename+'_'+str(num)+'.fasta'), 'w') as output:
            lines = return_lines(input, linepos, x[0], x[1])
            output.write('%s' % ''.join(lines))


def downloadIPRproperties(name, cpus):
    download('https://raw.githubusercontent.com/ebi-pf-team/interproscan/5.22-61.0/core/jms-implementation/support-mini-x86-32/interproscan.properties', 'ipr.tmp')
    with open(name, 'w') as outfile:
        with open('ipr.tmp', 'r') as infile:
            for line in infile:
                if line.startswith('number.of.embedded.workers='):
                    outfile.write('number.of.embedded.workers=1\n')
                elif line.startswith('maxnumber.of.embedded.workers='):
                    outfile.write('maxnumber.of.embedded.workers=%i\n' % cpus)
                else:
                    outfile.write(line)
    os.remove('ipr.tmp')


def runDocker(input):
    current = os.getcwd()
    UID = str(os.getuid())
    GROUP = str(os.getgid())
    prop = os.path.join(current, tmpdir, 'interproscan.properties')
    cmd = ['docker', 'run', '-u', UID+':'+GROUP, '--rm', '-v', os.path.join(current, tmpdir)+':/dir',
           '-v', os.path.join(current, tmpdir)+':/in', '-v', prop +
           ':/interproscan-5.22-61.0/interproscan.properties',
           'blaxterlab/interproscan:latest', 'interproscan.sh', '-i', '/in/'+input, '-d', '/dir',
           '-dp', '-f', 'XML', '-goterms', '-dra']
    logfile = os.path.join(tmpdir, input.split('.fasta')[0])
    logfile = logfile + '.log'
    with open(logfile, 'w') as log:
        log.write('%s\n' % ' '.join(cmd))
        subprocess.call(cmd, cwd=tmpdir, stdout=log, stderr=log)


def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        runDocker(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def runLocal(input):
    cmd = [iprpath, '-i', input, '-d', '.', '-f', 'XML', '-goterms', '-pa']
    logfile = os.path.join(tmpdir, input.split('.fasta')[0])
    logfile = logfile + '.log'
    with open(logfile, 'w') as log:
        subprocess.call(cmd, cwd=tmpdir, stdout=log, stderr=log)


def safe_run2(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        runLocal(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def runMultiProgress(function, inputList, cpus, progress=True):
    # setup pool
    p = multiprocessing.Pool(cpus)
    # setup results and split over cpus
    tasks = len(inputList)
    results = []
    for i in inputList:
        results.append(p.apply_async(function, [i]))
    # refresh pbar every 5 seconds
    if progress:
        while True:
            incomplete_count = sum(1 for x in results if not x.ready())
            if incomplete_count == 0:
                break
            sys.stdout.write("Progress: %.2f%% \r" %
                            (float(tasks - incomplete_count) / tasks * 100))
            sys.stdout.flush()
            time.sleep(1)
    p.close()
    p.join()


# make temp directory
tmpdir = 'iprscan_' + str(uuid.uuid4())
os.makedirs(tmpdir)

# check input, if folder
input = None
finalOut = None
if os.path.isfile(args.input):
    input = args.input
    if not args.out:
        print('Please specify an output file, -o,--out.')
        sys.exit(1)
    finalOut = args.out
elif os.path.isdir(args.input):  # now run through funannotate folders
    # funannotate results 1) in update folder or 2) in predict folder
    if os.path.isdir(os.path.join(args.input, 'update_results')):
        inputdir = os.path.join(args.input, 'update_results')
        outputdir = args.input
    elif os.path.isdir(os.path.join(args.input, 'predict_results')):
        inputdir = os.path.join(args.input, 'predict_results')
        outputdir = args.input
    else:
        # here user specified the predict_results folder, or it is a custom folder
        inputdir = os.path.join(args.input)
        if not args.out:
            print('Please specify an output file, -o,--out.')
            sys.exit(1)
    # get files that you need
    for file in os.listdir(inputdir):
        if file.endswith('.proteins.fa'):
            input = os.path.join(inputdir, file)
    if not args.out:
        # setup output if funannotate outputdir
        if not os.path.isdir(os.path.join(outputdir, 'annotate_misc')):
            os.makedirs(os.path.join(outputdir, 'annotate_misc'))
        finalOut = os.path.join(outputdir, 'annotate_misc', 'iprscan.xml')
    else:
        finalOut = args.out
else:
    print(('%s input does not exist' % args.input))
    sys.exit(1)
if not input:
    print('Error: could not parse input. Should be base funannotate folder or protein fasta file.')
    sys.exit(1)
if not finalOut:
    print('Error: could not parse output, specify')

# figure out number of chunks
count = countfasta(input)
print(('Running InterProScan5 on %i proteins' % count))
if args.num > count:
    chunks = 1
else:
    chunks = int(round(count / float(args.num)))
    if chunks == 1:  # means we rounded down to 1, but we actually want to split this into 2
        chunks = 2
threads = int(round(args.cpus / float(args.cpus_per_chunk)))
if threads < 1:
    threads = 1
# split protein fasta files into chunks
if chunks > 1:
    split_fasta(input, tmpdir, chunks)
else:
    tmpfile = os.path.basename(input).split('.fa')[0]
    shutil.copyfile(input, os.path.join(tmpdir, tmpfile+'.fasta'))

# get list of inputs
file_list = []
for file in os.listdir(tmpdir):
    if file.endswith('.fasta'):
        file_list.append(file)

if args.method == 'docker':
    checkDocker()
    ipr_properties = os.path.join(tmpdir, 'interproscan.properties')
    downloadIPRproperties(ipr_properties, args.cpus_per_chunk)
    if chunks > 1:
        runMultiProgress(safe_run, file_list, threads,
                         progress=args.progress)
    else:
        runDocker(file_list[0])

elif args.method == 'local':
    if not args.iprscan_path:
        print('You must specify location of interproscan.sh to --iprscan_path, exiting.')
        sys.exit(1)
    else:
        if args.iprscan_path.endswith('interproscan.sh'):
            iprpath = lib.which_path(args.iprscan_path)
        else:
            iprpath = os.path.join(args.iprscan_path, 'interproscan.sh')
    if not os.path.isfile(iprpath):
        print(('%s is not a valid path to interproscan.sh' % iprpath))
        sys.exit(1)
    print('Important: you need to manually configure your interproscan.properties file for embedded workers.')
    print(('Will try to launch %i interproscan processes, adjust -c,--cpus for your system' % args.cpus))
    if chunks > 1:
        runMultiProgress(safe_run2, file_list, args.cpus,
                         progress=args.progress)
    else:
        runLocal(file_list[0])

final_list = []
logfiles = []
for file in os.listdir(tmpdir):
    if file.endswith('.xml'):
        final_list.append(os.path.join(tmpdir, file))
    elif file.endswith('.log'):
        logfiles.append(os.path.join(tmpdir, file))

# apparently IPRscan XML has changed the header format in newest version [accidental?]
with open(finalOut, 'w') as output:
    for i, x in enumerate(final_list):
        with open(x, 'r') as infile:
            lines = infile.readlines()
            if i == 0:
                if '<protein-matches xml' in lines[0]:
                    linestart = 1
                elif '<protein-matches xml' in lines[1]:
                    linestart = 2
                for line in lines[:-1]:
                    output.write(line)
            else:
                for line in lines[linestart:-1]:
                    output.write(line)
    output.write('</protein-matches>\n')

# sometimes docker fails because can't mount from this directory, i.e. if not in docker preferences, check logfile
doublecheck = True
with open(logfiles[0], 'r') as logcheck:
    for line in logcheck:
        if line.startswith('docker:'):
            if 'Error' in line:
                print(line)
                doublecheck = False
if doublecheck:
    if not args.debug:
        if os.path.isfile(finalOut):
            shutil.rmtree(tmpdir)
    print('InterProScan5 search has completed successfully!')
    print(('Results are here: %s' % finalOut))
else:
    print(('Docker IPRscan run has failed, see log file: %s' % logfiles[0]))
