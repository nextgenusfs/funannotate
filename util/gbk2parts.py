#!/usr/bin/env python

import sys, os, inspect, argparse
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
parser.add_argument('-g', '--gbk', required=True, help='Genome in GenBank format')
parser.add_argument('-o', '--output', required=True, help='Output basename')
args=parser.parse_args()

#setup output files
tblout = args.output+'.tbl'
protout = args.output+'.proteins.fasta'
transout = args.output+'.transcripts.fasta'
dnaout = args.output+'.scaffolds.fasta'
lib.gb2parts(args.gbk, tblout, protout, transout, dnaout)