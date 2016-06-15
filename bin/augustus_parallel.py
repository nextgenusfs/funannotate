#!/usr/bin/env python

import sys, multiprocessing, subprocess, os, shutil, argparse, time
from Bio import SeqIO

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='augustus_parallel.py', usage="%(prog)s [options] -i genome.fasta -s botrytis_cinera -o new_genome",
    description='''Script that does it all...''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('-s','--species', required=True, help='Augustus species name')
parser.add_argument('--hints', help='Hints file (PE)')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to run')
args=parser.parse_args()

#check for augustus installation
try:
    AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
except KeyError:
    if not args.AUGUSTUS_CONFIG_PATH:
        print("$AUGUSTUS_CONFIG_PATH environmental variable not found, Augustus is not properly configured")
        os._exit(1)
if AUGUSTUS.endswith('config'):
    AUGUSTUS_BASE = AUGUSTUS.replace('config', '')
elif AUGUSTUS.endswith('config'+os.sep):
    AUGUSTUS_BASE = AUGUSTUS.replace('config'+os.sep, '')

#setup hints and extrinic input, hard coded for protein and transcript alignments from funannotate
extrinsic = '--extrinsicCfgFile='+os.path.join(AUGUSTUS_BASE, 'config', 'extrinsic', 'extrinsic.E.XNT.cfg')

def countGFFgenes(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if "\tgene\t" in line:
                count += 1
    return count

def runAugustus(Input):
    FNULL = open(os.devnull, 'w')
    species='--species='+args.species
    hints_input = '--hintsfile='+os.path.join(tmpdir, Input+'.hints.gff')
    aug_out = os.path.join(tmpdir, Input+'.augustus.gff3')   
    with open(aug_out, 'w') as output:
        if args.hints:
            subprocess.call(['augustus', species, hints_input, extrinsic, '--gff3=on', os.path.join(tmpdir, Input+'.fa')], stdout = output, stderr= FNULL)
        else:
            subprocess.call(['augustus', species, '--gff3=on', os.path.join(tmpdir, Input+'.fa')], stdout = output, stderr = FNULL)


#first step is to split input fasta file into individual files in tmp folder
print("Splitting contigs and hints files")
tmpdir = 'augustus_tmp_'+str(os.getpid())
os.makedirs(tmpdir)
scaffolds = []
with open(args.input, 'rU') as InputFasta:
    for record in SeqIO.parse(InputFasta, 'fasta'):
        name = str(record.id)
        scaffolds.append(name)
        outputfile = os.path.join(tmpdir, name+'.fa')
        with open(outputfile, 'w') as output:
            SeqIO.write(record, output, 'fasta')

#if hints file passed, split it up by scaffold
if args.hints:
    for i in scaffolds:
        with open(os.path.join(tmpdir, i+'.hints.gff'), 'w') as output:
            with open(args.hints, 'rU') as hintsfile:
                for line in hintsfile:
                    cols = line.split('\t')
                    if cols[0] == i:
                        output.write(line)

#now loop through each scaffold running augustus
if args.cpus > len(scaffolds):
    num = len(scaffolds)
else:
    num = args.cpus
print("Running augustus on %i scaffolds, using %i CPUs" % (len(scaffolds), num))
p = multiprocessing.Pool(num)
rs = p.map_async(runAugustus, scaffolds)
p.close()
while (True):
    if (rs.ready()): break
    remaining = rs._number_left
    print "Waiting for", remaining, "augustus jobs to complete..."
    time.sleep(30)
print("Augustus prediction is finished, now concatenating results")
with open(os.path.join(tmpdir, 'augustus_all.gff3'), 'w') as output:
    for file in scaffolds:
        file = os.path.join(tmpdir, file+'.augustus.gff3')
        with open(file) as input:
            output.write(input.read())

join_script = os.path.join(AUGUSTUS_BASE, 'scripts', 'join_aug_pred.pl')
with open(args.out, 'w') as finalout:
    with open(os.path.join(tmpdir, 'augustus_all.gff3'), 'rU') as input:
        subprocess.call([join_script],stdin = input, stdout = finalout)
shutil.rmtree(tmpdir)
print("Found %i total gene models" % countGFFgenes(args.out))
