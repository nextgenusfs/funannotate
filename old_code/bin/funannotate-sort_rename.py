#!/usr/bin/env python

import sys
import argparse
import os
import inspect
from Bio import SeqIO
from funannotate.utilities import countfasta

def SortRenameHeaders(input, basename, output):
    #sort records and write temp file
    with open(output, 'w') as output:
        with open(input, 'rU') as input:
            records = list(SeqIO.parse(input, 'fasta'))
            records.sort(cmp=lambda x,y: cmp(len(y),len(x)))
            counter = 1
            for rec in records:
                rec.name = ''
                rec.description = ''
                rec.id = basename + '_' + str(counter)
                if len(rec.id) > 16:
                    print "Error. Fasta header too long %s.  Choose a different --base name. NCBI/GenBank max is 16 characters" % rec.id
                    os._exit(1)
                if args.minlen:
                    if len(rec.seq) < int(args.minlen):
                        continue
                counter +=1
                SeqIO.write(rec, output, 'fasta')

def main(args):
	#setup menu with argparse
	class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
		def __init__(self,prog):
			super(MyFormatter,self).__init__(prog,max_help_position=48)
	parser=argparse.ArgumentParser(prog='sort_rename.py', usage="%(prog)s [options] -i genome.fa -o sorted.fa",
		description='''Script that sorts input by length and then renames contig headers.''',
		epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
		formatter_class = MyFormatter)
	parser.add_argument('-i','--input', required=True, help='Multi-fasta genome file')
	parser.add_argument('-o','--out', required=True, help='Cleaned output (FASTA)')
	parser.add_argument('-b','--base', default='scaffold', help='Basename of contig header')
	parser.add_argument('--minlen', help='Contigs shorter than threshold are discarded')
	args=parser.parse_args(args)

	Count = countfasta(args.input)
	print('{0:,}'.format(Count) + ' contigs records loaded')
	print("Sorting and renaming contig headers")
	if args.minlen:
		print("Removing contigs less than %s bp" % args.minlen)
	SortRenameHeaders(args.input, args.base, args.out)
	FCount = countfasta(args.out)
	print('{0:,}'.format(FCount) + ' contigs saved to file')
	
if __name__ == "__main__":
	main(sys.argv[1:])