#!/usr/bin/env python

import sys, subprocess, os, itertools, argparse
from Bio import SeqIO

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='contig_cleaner.py', usage="%(prog)s [options] -i genome.fa -o cleaned.fa",
    description='''Script that removes short scaffolds that are duplicated elsewhere.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Multi-fasta genome file')
parser.add_argument('-o','--out', required=True, help='Cleaned output (FASTA)')
parser.add_argument('-p','--pident', type=int, default=90, help='percent identity of contig')
parser.add_argument('-c','--cov', type=int, default=90, help='coverage of contig')
args=parser.parse_args()

def which(name):
    try:
        with open(os.devnull) as devnull:
            subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def CheckDependencies(input):
    missing = []
    for p in input:
        if which(p) == False:
            missing.append(p)
    if missing != []:
        error = ", ".join(missing)
        log.error("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        sys.exit(1)

def Sortbysize(input):
    #sort records and return a list of scaffolds in descending size order
    contigs = []
    with open(input, 'rU') as input:
        records = list(SeqIO.parse(input, 'fasta'))
        records.sort(cmp=lambda x,y: cmp(len(y),len(x)), reverse=True)
        for rec in records:
            contigs.append(rec.id)
        return contigs

def getFasta(sequences, header):
    with open('query.fa', 'w') as fasta:
        with open(sequences, 'rU') as input:
            SeqRecords = SeqIO.parse(input, 'fasta')
            for rec in SeqRecords:
                if rec.id == header:
                    SeqIO.write(rec, fasta, 'fasta')

def getReference(sequences, index, Contigs):
    #pull out rest of list plus any keepers
    contiglist = Contigs[index+1:] + keepers
    with open('reference.fa', 'w') as output:
        with open(sequences, 'rU') as input:
            SeqRecords = SeqIO.parse(input, 'fasta')
            for rec in SeqRecords:
                if rec.id in contiglist:
                    SeqIO.write(rec, output, 'fasta')

def runNucmer(query, reference, output):
    FNULL = open(os.devnull, 'w')
    subprocess.call(['nucmer', '-p', output, query, reference], stdout = FNULL, stderr = FNULL)
    input = output + '.delta'
    coord_out = output + '.coords'
    with open(coord_out, 'w') as coords:
        subprocess.call(['show-coords', '-r', '-c', '-l', '-T', '-o', '-I', '75', input], stdout = coords, stderr = FNULL)
    #now load in results and filter
    garbage = False #assume this is a good contig
    with open(coord_out, 'rU') as c:
        for line in itertools.islice(c, 4, None):
            cols = line.split('\t')
            match = (float(cols[6]), float(cols[9]))
            if match[0] > args.pident and match[1] > args.cov:
                print "%s appears duplicated: %i%% identity over %i%% of the contig. contig length: %s " % (output, match[0], match[1], cols[7])
                #print match
                garbage = True
                break
        if not garbage:
            keepers.append(output)
    os.remove(input)
    os.remove(coord_out)

#run some checks of dependencies first
programs = ['nucmer', 'show-coords']
CheckDependencies(programs)

#now get list of scaffolds, shortest->largest
scaffolds = Sortbysize(args.input)

#make a list of contigs that are to keep from the analysis
global keepers
keepers = []

#now loop through the list
for i in range(0, len(scaffolds)):
    getFasta(args.input, scaffolds[i])
    getReference(args.input, i, scaffolds)
    runNucmer('query.fa', 'reference.fa', scaffolds[i])
    os.remove('query.fa')
    os.remove('reference.fa')

print"------------------------------------"
print"%i input contigs, %i duplicated, %i written to file" % (len(scaffolds), (len(scaffolds) - len(keepers)), len(keepers))

#finally write a new reference based on list of keepers
with open(args.out, 'w') as output:
    with open(args.input, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fasta')
        for rec in SeqRecords:
            if rec.id in keepers:
                SeqIO.write(rec, output, 'fasta')
