#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import sys
import os
import subprocess
import argparse
import signal
import uuid
import shutil
from multiprocessing import Pool
from tempfile import NamedTemporaryFile as NTF
from Bio.SeqIO.FastaIO import SimpleFastaParser
from funannotate.library import CheckDependencies, softwrap, countfasta


def calcN50(input):
    lengths = []
    with open(input, 'r') as infile:
        for id, sequence in SimpleFastaParser(infile):
            lengths.append(len(sequence))
    # now get N50
    lengths.sort()
    nlist = []
    for x in lengths:
        nlist += [x]*x
    if len(nlist) % 2 == 0:
        medianpos = int(len(nlist) / 2)
        N50 = int((nlist[medianpos] + nlist[medianpos-1]) / 2)
    else:
        medianpos = int(len(nlist) / 2)
        N50 = int(nlist[medianpos])
    return N50


def Sortbysize(input, n50, minlen=500):
    IUPAC = {'A', 'C', 'G', 'T', 'R', 'Y',
             'S', 'W', 'K', 'M', 'B', 'D',
             'H', 'V', 'N'}
    contigs = []
    keep = []
    Seqs = []
    bad = []
    with open(input, 'r') as infile:
        for header, sequence in SimpleFastaParser(infile):
            characters = {}
            for nuc in sequence:
                nuc = nuc.upper()
                if nuc not in characters:
                    characters[nuc] = 1
                else:
                    characters[nuc] += 1
            for c in characters.keys():
                if c not in IUPAC:
                    bad.append(header)
                    print('ERROR: {} contains non-IUPAC [{}] character, removing'.format(header, c))
            if len(characters) < 4:
                bad.append(header)
                print('ERROR: {} contains only {} nucleotides, removing'.format(header, len(characters)))
            Seqs.append((header, len(sequence)))
    # sort by length
    sortedSeqs = sorted(Seqs, key=lambda x: x[1], reverse=True)
    # loop through and return contigs and keepers
    for name, length in sortedSeqs:
        if name in bad:
            continue
        if length >= minlen:
            if n50:
                if length >= n50:
                    keep.append(name)
                else:
                    contigs.append(name)
            else:
                contigs.append(name)
    return contigs, keep


def generateFastas(input, qFasta, rFasta, index, Contigs, query):
    # loop through fasta once, generating query and reference
    contiglist = Contigs[index+1:] + keepers
    with open(input, 'r') as infile:
        for Id, Sequence in SimpleFastaParser(infile):
            if Id == query:
                qFasta.write('>%s\n%s\n' % (Id, softwrap(Sequence)))
            elif Id in contiglist:
                rFasta.write('>%s\n%s\n' % (Id, softwrap(Sequence)))
    # flush changes to disk before calling minimap2
    # note: keeping file open for writing while trying to read it will not work on Windows
    qFasta.flush()
    rFasta.flush()
    return


def runMinimap2(query, reference, output, index, min_pident=95, min_cov=95):
    '''
    I have not found parameters that mirror mummer yet, do not use minimap method
    '''
    FNULL = open(os.devnull, 'w')
    minitmp = ""
    with NTF(mode='w', suffix=".tmp", prefix="minmap_", dir=tmpdir, delete=False) as out:
        minitmp = out.name
        subprocess.call(['minimap2', '-x', 'asm5', '-N5',
                         reference, query], stdout=out, stderr=FNULL)
    # now load in results and filter
    garbage = False  # assume this is a good contig
    with open(minitmp, 'r') as data:
        for line in data:
            line = line.replace('\n', '')
            qID, qLen, qStart, qEnd, strand, tID, tLen, tStart, tEnd, matches, alnLen, mapQ = line.split('\t')[
                :12]
            pident = float(matches) / int(alnLen) * 100
            coverage = float(alnLen) / int(qLen) * 100
            # print qID, str(qLen), tID, matches, alnLen, str(pident), str(coverage)
            if pident > min_pident and coverage > min_cov:
                print(("{} appears duplicated: {:.0f}% identity over {:.0f}% of the contig. contig length: {}".format(
                    output, pident, coverage, qLen)))
                garbage = True
                break
    os.remove(minitmp)
    return (output, garbage)


def align_contigs(mp_args):
    scaffolds, i = mp_args
    # Create/open named temporary files for query.fa and reference.fa; they will be deleted when closed.
    with NTF(mode='w', suffix=".fa", prefix="query_", dir=tmpdir, delete=True) as qFasta, \
         NTF(mode='w', suffix=".fa", prefix="reference_", dir=tmpdir, delete=True) as rFasta:
        generateFastas(GENOME, qFasta, rFasta, i, scaffolds, scaffolds[i])
        out = runMinimap2(qFasta.name, rFasta.name,
                          scaffolds[i], i, min_pident=PIDENT, min_cov=COV)
    return out


def multithread_aligning(scaffolds):
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    p = Pool(CPUS)
    signal.signal(signal.SIGINT, original_sigint_handler)
    mp_args = [(scaffolds, i) for i in range(0, len(scaffolds))]
    try:
        out = p.map_async(align_contigs, mp_args)
        result = out.get(999999999)
    except KeyboardInterrupt:
        p.terminate()
    else:
        p.close()
    p.join()
    return result


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='contig_cleaner.py', usage="%(prog)s [options] -i genome.fa -o cleaned.fa",
                                     description='''Script that removes short scaffolds that are duplicated elsewhere.''',
                                     epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help='Multi-fasta genome file')
    parser.add_argument('-o', '--out', required=True,
                        help='Cleaned output (FASTA)')
    parser.add_argument('-p', '--pident', type=int,
                        default=95, help='percent identity of contig')
    parser.add_argument('-c', '--cov', type=int,
                        default=95, help='coverage of contig')
    parser.add_argument('-m', '--minlen', type=int,
                        default=500, help='Minimum length of contig')
    parser.add_argument('--cpus', default=2, type=int,
                        help='Number of CPUs to use')
    parser.add_argument('--tmpdir',
                        help='TMP directory to use')
    parser.add_argument('--exhaustive', action='store_true',
                        help='Compute every contig, else stop at N50')
    parser.add_argument('--debug', action='store_true',
                        help='Debug the output')
    args = parser.parse_args(args)

    # setup some global variables used in functions above
    global GENOME, CPUS, PIDENT, COV, keepers, repeats, tmpdir
    GENOME = args.input
    CPUS = args.cpus
    PIDENT = args.pident
    COV = args.cov
    keepers, repeats = ([],)*2

    # write in a tmpdir
    if args.tmpdir:
        tmpdir = args.tmpdir
    else:
        tmpdir = 'clean_{}'.format(str(uuid.uuid4()))
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

    # run some checks of dependencies first
    programs = ['minimap2']
    CheckDependencies(programs)

    # calculate N50 of assembly
    n50 = calcN50(args.input)

    # now get list of scaffolds, shortest->largest
    if args.exhaustive:
        scaffolds, keepers = Sortbysize(args.input, False, minlen=args.minlen)
    else:
        scaffolds, keepers = Sortbysize(args.input, n50, minlen=args.minlen)

    print("-----------------------------------------------")
    PassSize = len(scaffolds)+len(keepers)
    print(("{:,} input contigs, {:,} larger than {:,} bp, N50 is {:,} bp".format(
        countfasta(args.input), PassSize, args.minlen, n50)))
    if args.exhaustive:
        print(("Checking duplication of {:,} contigs".format(len(scaffolds))))
    else:
        print(("Checking duplication of {:,} contigs shorter than N50".format(
            len(scaffolds))))
    print("-----------------------------------------------")

    # now generate pool and parallel process the list
    mp_output = multithread_aligning(scaffolds)

    for output, garbage in mp_output:
        if not garbage:
            keepers.append(output)
        else:
            repeats.append(output)

    print("-----------------------------------------------")
    print(("{:,} input contigs; {:,} larger than {:} bp; {:,} duplicated; {:,} written to file".format(
        countfasta(args.input), PassSize, args.minlen, len(repeats), len(keepers))))
    if args.debug:
        print(("\nDuplicated contigs are:\n{:}\n".format(', '.join(repeats))))
        print(("Contigs to keep are:\n{:}\n".format(', '.join(keepers))))

    # finally write a new reference based on list of keepers
    with open(args.out, 'w') as output:
        with open(args.input, 'r') as input:
            for title, sequence in SimpleFastaParser(input):
                if title in keepers:
                    output.write('>{}\n{}\n'.format(title, softwrap(sequence)))

    shutil.rmtree(tmpdir)


if __name__ == "__main__":
    main(sys.argv[1:])
