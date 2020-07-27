#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import subprocess
import funannotate.library as lib


def runGOenrichment(input):
    basename = os.path.basename(input).replace('.txt', '')
    goa_out = os.path.join(args.out, basename+'.go.enrichment.txt')
    go_log = os.path.join(args.out, basename+'.go.enrichment.log')
    if not lib.checkannotations(goa_out):
        cmd = ['find_enrichment.py', '--obo', os.path.join(FUNDB, 'go.obo'),
               '--pval', '0.001', '--alpha', '0.001', '--method', 'fdr',
               '--outfile', goa_out, input, os.path.join(args.input, 'population.txt'),
               os.path.join(args.input, 'associations.txt')]
        with open(go_log, 'w') as outfile:
            outfile.write('{}\n'.format(' '.join(cmd)))
        with open(go_log, 'a') as outfile:
            subprocess.call(cmd, stdout=outfile, stderr=outfile)


def GO_safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        runGOenrichment(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))

# setup menu with argparse


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)


parser = argparse.ArgumentParser(prog='enrichment_parallel.py',
                                 description='''Run goatools enrichment in parallel.''',
                                 epilog="""Written by Jon Palmer (2019) nextgenusfs@gmail.com""",
                                 formatter_class=MyFormatter)
parser.add_argument('-i', '--input', required=True,
                    help='folder of protein fasta files')
parser.add_argument('-d', '--db', required=True,
                    help='location of HMM database')
parser.add_argument('-c', '--cpus', default=1, type=int,
                    help='location of HMM database')
parser.add_argument('-o', '--out', required=True, help='output file')
args = parser.parse_args()

global FUNDB, FNULL
FUNDB = args.db
FNULL = open(os.devnull, 'w')

# now loop through each genome comparing to population
file_list = []
for f in os.listdir(args.input):
    if f.startswith('associations'):
        continue
    if f.startswith('population'):
        continue
    file = os.path.join(args.input, f)
    if lib.checkannotations(file):
        file_list.append(file)
    else:
        print('  WARNING: skipping {} as no GO terms'.format(f))

# run over multiple CPUs
if len(file_list) > args.cpus:
    procs = args.cpus
else:
    procs = len(file_list)

lib.runMultiNoProgress(GO_safe_run, file_list, procs)
