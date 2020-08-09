#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import warnings
import subprocess
from natsort import natsorted
import funannotate.library as lib
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO


def PfamHmmer(input):
    HMM = os.path.join(FUNDB, 'Pfam-A.hmm')
    base = os.path.basename(input).split('.fa')[0]
    pfam_out = os.path.join(os.path.dirname(input), base+'.pfam.txt')
    cmd = ['hmmsearch', '--domtblout', pfam_out,
           '--cpu', '1', '--cut_ga', HMM, input]
    subprocess.call(cmd, stdout=FNULL, stderr=FNULL)


def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        PfamHmmer(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def combineHmmerOutputs(inputList, output):
    # function to combine multiple HMMER runs with proper header/footer so biopython can read
    allHeadFoot = []
    with open(inputList[0], 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                allHeadFoot.append(line)
    with open(output, 'w') as out:
        for x in allHeadFoot[:3]:
            out.write(x)
        for file in inputList:
            with open(file, 'r') as resultin:
                for line in resultin:
                    if line.startswith('#') or line.startswith('\n'):
                        continue
                    out.write(line)
        for y in allHeadFoot[3:]:
            out.write(y)


def multiPFAMsearch(inputList, cpus, tmpdir, output):
    # run hmmerscan multithreaded by running at same time
    # input is a list of files, run multiprocessing on them
    pfam_results = os.path.join(os.path.dirname(tmpdir), 'pfam.txt')
    pfam_filtered = os.path.join(os.path.dirname(tmpdir), 'pfam.filtered.txt')
    lib.runMultiNoProgress(safe_run, inputList, cpus)

    # now grab results and combine, kind of tricky as there are header and footers for each
    resultList = [os.path.join(tmpdir, f) for f in os.listdir(
        tmpdir) if os.path.isfile(os.path.join(tmpdir, f)) and f.endswith('.pfam.txt')]
    combineHmmerOutputs(resultList, pfam_results)

    # now parse results
    with open(output, 'w') as out:
        with open(pfam_filtered, 'w') as filtered:
            with open(pfam_results, 'r') as results:
                for qresult in SearchIO.parse(results, "hmmsearch3-domtab"):
                    hits = qresult.hits
                    num_hits = len(hits)
                    if num_hits > 0:
                        for i in range(0, num_hits):
                            hit_evalue = hits[i].evalue
                            query = hits[i].id
                            pfam = qresult.accession.split('.')[0]
                            hmmLen = qresult.seq_len
                            hmm_aln = int(hits[i].hsps[0].hit_end) - \
                                int(hits[i].hsps[0].hit_start)
                            coverage = hmm_aln / float(hmmLen)
                            if coverage < 0.50:  # coverage needs to be at least 50%
                                continue
                            filtered.write("%s\t%s\t%s\t%f\n" %
                                           (query, pfam, hit_evalue, coverage))
                            out.write("%s\tdb_xref\tPFAM:%s\n" % (query, pfam))


def dbCANHmmer(input):
    HMM = os.path.join(FUNDB, 'dbCAN.hmm')
    base = os.path.basename(input).split('.fa')[0]
    outfiles = os.path.join(os.path.dirname(input), base+'.dbcan.txt')
    cmd = ['hmmscan', '--domtblout', outfiles,
           '--cpu', '1', '-E', '1e-17', HMM, input]
    subprocess.call(cmd, stdout=FNULL, stderr=FNULL)


def safe_run2(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        dbCANHmmer(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def dbCANsearch(inputList, cpus, evalue, tmpdir, output):
    # run hmmerscan
    dbCAN_out = os.path.join(tmpdir, 'dbCAN.txt')
    dbCAN_filtered = os.path.join(tmpdir, 'dbCAN.filtered.txt')
    lib.runMultiNoProgress(safe_run2, inputList, cpus)
    # now grab results
    resultList = [os.path.join(tmpdir, f) for f in os.listdir(
        tmpdir) if os.path.isfile(os.path.join(tmpdir, f)) and f.endswith('.dbcan.txt')]
    combineHmmerOutputs(resultList, dbCAN_out)

    # now parse results
    Results = {}
    with open(dbCAN_filtered, 'w') as filtered:
        filtered.write(
            "#HMM_family\tHMM_len\tQuery_ID\tQuery_len\tE-value\tHMM_start\tHMM_end\tQuery_start\tQuery_end\tCoverage\n")
        with open(dbCAN_out, 'r') as results:
            for qresult in SearchIO.parse(results, "hmmscan3-domtab"):
                query_length = qresult.seq_len
                hits = qresult.hits
                num_hits = len(hits)
                if num_hits > 0:
                    for i in range(0, num_hits):
                        hit_evalue = hits[i].evalue
                        if hit_evalue > evalue:
                            continue
                        hit = hits[i].id
                        hmmLen = hits[i].seq_len
                        hmm_aln = int(hits[i].hsps[0].hit_end) - \
                            int(hits[i].hsps[0].hit_start)
                        coverage = hmm_aln / float(hmmLen)
                        if coverage < 0.45:
                            continue
                        query = hits[i].query_id
                        filtered.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\n" % (
                            hit, hmmLen, query, query_length, hit_evalue,
                            hits[i].hsps[0].hit_start,
                            hits[i].hsps[0].hit_end,
                            hits[i].hsps[0].query_start,
                            hits[i].hsps[0].query_end,
                            coverage))
                        if query not in Results:
                            Results[query] = [hit]
                        else:
                            Results[query].append(hit)
    # run through results and simplify subdomain hits
    with open(output, 'w') as out:
        for k, v in natsorted(Results.items()):
            simplified = []
            for x in v:
                if '_' in x:
                    cazy, subdomain = x.rsplit('_', 1)
                    if cazy not in simplified:
                        simplified.append(cazy)
                else:
                    if not x in simplified:
                        simplified.append(x)
            for hit in simplified:
                out.write("{}\tnote\tCAZy:{}\n".format(k, hit))


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)


parser = argparse.ArgumentParser(
    prog='hmmer_parallel.py',
    description='''Run hmmer3 multipthreaded.''',
    epilog="""Written by Jon Palmer (2019) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)
parser.add_argument('-i', '--input', required=True,
                    help='folder of protein fasta files')
parser.add_argument('-m', '--method', default='pfam',
                    choices=['pfam', 'cazy'], help='database to search')
parser.add_argument('-d', '--db', required=True,
                    help='location of HMM database')
parser.add_argument('-c', '--cpus', default=1, type=int,
                    help='location of HMM database')
parser.add_argument('-o', '--out', required=True, help='output file')
args = parser.parse_args()

global FUNDB, FNULL
FUNDB = args.db
FNULL = open(os.devnull, 'w')
splitProts = [os.path.join(args.input, f) for f in os.listdir(
    args.input) if os.path.isfile(os.path.join(args.input, f))]
if args.method == 'pfam':
    multiPFAMsearch(splitProts, args.cpus, args.input, args.out)
elif args.method == 'cazy':
    dbCANsearch(splitProts, args.cpus, 1e-17, args.input, args.out)
