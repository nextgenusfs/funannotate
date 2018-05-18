#!/usr/bin/env python

import sys, os, inspect, argparse, shutil, subprocess
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='gbk2parts.py', 
    description = '''Script to convert GBK file to its components.''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--tbl', required=True, help='Genome annotation in tbl format')
parser.add_argument('-f', '--fasta', required=True, help='Genome in FASTA format')
parser.add_argument('-s', '--species', required=True, help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('--isolate', help='Isolate name (e.g. Af293)')
parser.add_argument('--strain', help='Strain name (e.g. CEA10)')
parser.add_argument('-t','--tbl2asn', help='Custom parameters for tbl2asn, example: linkage and gap info')
parser.add_argument('--sbt', help='tbl2asn template file')
parser.add_argument('-o', '--output', help='Output basename')
args=parser.parse_args()

def runSubprocess(cmd, dir):
    proc = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stdout:
        print(stdout)


def runtbl2asn(folder, template, discrepency, organism, isolate, strain, parameters, version):
    '''
    function to run NCBI tbl2asn
    '''
    #get funannotate version
    fun_version = lib.get_version()
    #input should be a folder
    if not os.path.isdir(folder):
        log.error("tbl2asn error: %s is not a directory, exiting" % folder)
        sys.exit(1)
    #based on organism, isolate, strain, construct meta info for -j flag
    if not organism:
        log.error("tbl2asn error: organism not specified")
        sys.exit(1)       
    meta = "[organism=" + organism + "]"
    if isolate:
        isolate_meta = "[isolate=" + isolate + "]"
        meta = meta + " " + isolate_meta
    if strain:
        strain_meta = "[strain=" + strain + "]"
        meta = meta + " " + strain_meta
    cmd = ['tbl2asn', '-y', '"Annotated using '+fun_version+'"', '-N', str(version), '-p', folder, '-t', template, '-M', 'n', '-Z', discrepency, '-j', '"'+meta+'"', '-V', 'b', '-c', 'fx', '-T', '-a', 'r10u']
    #check for custom parameters
    if parameters:
        params = parameters.split(' ')
        cmd = cmd + params
    runSubprocess(cmd, '.')
    return ' '.join(cmd)


#see if organism/species/isolate was passed at command line
organism = None
if args.species:
    organism = args.species
else:
    organism = os.path.basename(args.tbl).split('.t')[0]
if args.strain:
    organism_name = organism+'_'+args.strain
elif args.isolate:
    organism_name = organism+'_'+args.isolate
else:
    organism_name = organism
organism_name = organism_name.replace(' ', '_')
if args.output:
	outputname = args.output
else:
	outputname = organism_name

#create tmp folder to run tbl2asn from
#make tmp folder
tmp = outputname + '_tmp'
if not os.path.exists(tmp):
    os.makedirs(tmp)

#now move files into proper location
if not lib.checkannotations(args.fasta):
	print('FASTA genome file not found: {:}'.format(args.fasta))
	sys.exit(1)
if not lib.checkannotations(args.tbl):
	print('TBL annotations file not found: {:}'.format(args.tbl))
	sys.exit(1)
shutil.copyfile(args.fasta, os.path.join(tmp, 'genome.fsa'))
shutil.copyfile(args.tbl, os.path.join(tmp, 'genome.tbl'))

#now we can run tbl2asn
if args.sbt:
	SBT = args.sbt
else:
	SBT = os.path.join(parentdir, 'lib', 'test.sbt')
discrep = outputname+'.discrepency.txt'
version = 1
tbl2asn_cmd = runtbl2asn(tmp, SBT, discrep, organism, args.isolate, args.strain, args.tbl2asn, version)

#get output files
gbkout = outputname+'.gbk'
shutil.copyfile(os.path.join(tmp, 'genome.gbf'), gbkout)
lib.SafeRemove(tmp)