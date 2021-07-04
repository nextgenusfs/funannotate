#!/usr/bin/env python
# -*- coding: utf-8 -*-

# based on fasta2agp.pl from david.studholme@tsl.ac.uk
# rewritten in python by Jason Stajich @hyphaltip

import os
import sys
import re
import csv
import argparse
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_scaffolds_makeagp(scaffolds,agpout,ctgsout):
    x = 0
    i = 0
    spadesnamepat = re.compile(r'^NODE_(\d+)_length_\d+_cov_\d+')
    numnamepat = re.compile(r'^(\d+)$')
    validSeq   = re.compile(r'^[ACGTRYSWKMBDHVN]+$',flags=re.IGNORECASE)
    with open(agpout, 'w') as agpoutfh:
        csvout = csv.writer(agpoutfh,delimiter="\t",lineterminator="\n")
        with open(ctgsout,"w") as ctgoutfh:
            with open(scaffolds, 'r') as scaff_in:
                for seq in SeqIO.parse(scaff_in, "fasta"):
                    supercontig_id = seq.id
                    supercontig_seq = seq.seq
                    supercontig_desc = seq.description
                    supercontig_length = len(seq);
                    x = 0
                    m = spadesnamepat.match(supercontig_id) or spadesnamepat.match(supercontig_id)
                    if m:
                        supercontig_id = "scf_%s"%(m.match(1))
                    start_pos = 1 # keep track of whereabouts in this supercontig we are
                    substring_sequences = {}
                    for substring_sequence in re.split(r'(N{10,})',str(supercontig_seq),maxsplit=0,flags=re.IGNORECASE):
                        if len(substring_sequence) == 0:
                            continue
                        object1         = supercontig_id
                        object_beg2     = start_pos
                        object_end3     = start_pos + len(substring_sequence) - 1
                        part_number4    =  x
                        x += 1
                        component_type5 = None
                        component_id6a  = None
                        gap_length6b    = None
                        component_beg7a = None
                        gap_type7b      = None
                        component_end8a = None
                        linkage8b       = None
                        orientation9a   = None
                        filler9b        = None
                        if re.match(r'^N+$',substring_sequence):
                            ### This is poly-N gap between contigs
                            component_type5 = 'N'
                            gap_length6b    = len(substring_sequence)
                            gap_type7b      = 'scaffold'
                            linkage8b       = 'yes'
                            filler9b        = 'paired-ends'
                        elif validSeq.match(substring_sequence):
                            ### This is a contig
                            i+=1 # a counter, used for generating unique contig names
                            component_type5 = 'W'
                            component_id6a = "contig_%d"%(i)
                            component_beg7a = 1
                            component_end8a = len(substring_sequence)
                            orientation9a = '+'
                            ### Print FastA formatted contig
                            record = SeqRecord( Seq(substring_sequence),
                                                id=component_id6a,
                                                description="")
                            SeqIO.write(record, ctgoutfh, "fasta")
                        else:
                            print("Illegal characters in sequence")
                            print(substring_sequence)
                            return

                        start_pos += len (substring_sequence)
                        part_number4 += 1
                        if component_type5 == 'N':
                            ### print AGP line for gap
                            csvout.writerow([object1,object_beg2,object_end3, part_number4,component_type5,gap_length6b,gap_type7b,linkage8b,filler9b])
                        else:
                            ### print AGP line for contig
                            csvout.writerow([object1,object_beg2,object_end3, part_number4,component_type5,component_id6a,component_beg7a,component_end8a,orientation9a])


def main(args):
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)


    parser = argparse.ArgumentParser(
        prog='fasta2agp.py',
        description='''Convert FastA format scaffolds file into contigs file and print the AGP based on parsing gaps (N runs).''',
        #usage='''fasta2agp.py scaffolds.fa > scaffolds.agp''',
        epilog="""Written by Jason Stajich @hyphaltip (2021) jasonstajich.phd@gmail.com""",
        formatter_class=MyFormatter)
    parser.add_argument('--ext', default='contigs.fsa',
                        help='Default extensions for output contigs file')
    parser.add_argument('scaffoldfile', nargs='?',help='Scaffolds FastA file')
    parser.add_argument('agpfile', nargs='?',type=argparse.FileType('w'), default=sys.stdout,
                        help='AGP output file (defaults to STDOUT)')
    args = parser.parse_args(args)
    ctgfile = args.scaffoldfile + "." + args.ext
    m = re.match(r'^(\S+)\.(fa|fasta|fsa)$',args.scaffoldfile)
    if m:
        ctgfile = m.group(1)
        m = re.match(r'^(\S+)\.scaffolds?$',ctgfile)
        if m:
            ctgfile = "{}.{}".format(m.group(1),args.ext)
    # run cmd
    parse_scaffolds_makeagp(args.scaffoldfile,args.agpfile,ctgfile)

if __name__ == "__main__":
    main(sys.argv[1:])