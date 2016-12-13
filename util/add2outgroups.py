#!/usr/bin/env python

import sys, os, subprocess, shutil, argparse, inspect
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib


#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='funannotate-predict.py', usage="%(prog)s [options] -i genome.fasta",
    description = '''Script that adds a proteome to the outgroups.''',
    epilog = """Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--input', required=True, help='Proteome in FASTA format')
parser.add_argument('--species', required=True, help='Species name "binomial in quotes"')
parser.add_argument('--busco_db', default='dikarya', required=True, help='BUSCO database to use')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to use')
args=parser.parse_args()

#get base name
species = args.species.replace(' ', '_').lower()+'.'+args.busco_db
OUTGROUPS = os.path.join(parentdir, 'DB', 'outgroups')

#create log file
log_name = species+'-add2outgroups.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.SystemInfo()

#get version of funannotate
version = lib.get_version()
lib.log.info("Running %s" % version)

#check buscos, download if necessary
if not os.path.isdir(os.path.join(parentdir, 'DB', args.busco_db)):
    lib.download_buscos(args.busco_db)

ProtCount = lib.countfasta(args.input)
lib.log.info('{0:,}'.format(ProtCount) + ' protein records loaded')  

#convert to proteins and screen with busco
lib.log.info("Looking for BUSCO models with %s DB" % args.busco_db)
BUSCODB = os.path.join(parentdir, 'DB', args.busco_db)
BUSCO = os.path.join(parentdir, 'util', 'funannotate-BUSCO2.py')
cmd = [sys.executable, BUSCO, '-i', os.path.abspath(args.input), '-m', 'proteins', '--lineage', BUSCODB, '-o', species, '--cpu', str(args.cpus), '-f']
lib.runSubprocess(cmd, '.', lib.log)

#check that it ran correctly
busco_results = os.path.join('run_'+species, 'full_table_'+species+'.tsv')
if not lib.checkannotations(busco_results):
    lib.log.error("BUSCO failed, check logfile")
    sys.exit(1)
nameChange = {}
with open(busco_results, 'rU') as input:
    for line in input:
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        if cols[1] == 'Complete':
            if not cols[2] in nameChange:
                nameChange[cols[2]] = cols[0]
            else:
                lib.log.error("Duplicate ID found: %s %s. Removing from results" % (cols[2], cols[0]))
                del nameChange[cols[2]]

#output counts
lib.log.info('{0:,}'.format(len(nameChange)) + ' BUSCO models found')  

#index the proteome for parsing
SeqRecords = SeqIO.index(args.input, 'fasta')

#setup output proteome
busco_out = os.path.join(OUTGROUPS, species+'_buscos.fa')
with open(busco_out, 'w') as output:
    for k,v in nameChange.items():
        rec = SeqRecords[k]
        output.write('>%s\n%s\n' % (v, rec.seq))
lib.log.info("Results written to: %s" % busco_out)

#clean up your mess
shutil.rmtree('run_'+species)
shutil.rmtree('tmp')


