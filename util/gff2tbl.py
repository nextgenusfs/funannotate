#!/usr/bin/env python

import sys, os, inspect, argparse
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib
from Bio import SeqIO
from collections import OrderedDict

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='gff2prot.py', 
    description = '''Script to convert GFF3 and FASTA to tbl, proteins, transcripts.''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-g', '--gff3', required=True, help='Genome annotation GFF3 format')
parser.add_argument('-f', '--fasta', required=True, help='Genome in FASTA format')
args=parser.parse_args()

def scaffold2Dict(input):
    #get scaffold names/lengths
    scaffLen = {}
    with open(input, 'rU') as seqin:
        for record in SeqIO.parse(seqin, 'fasta'):
            if not record.id in scaffLen:
                scaffLen[record.id] = len(record.seq)
    return scaffLen

def dict2tbl(geneDict, scaffDict):
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    #get ordered dict to get genes in proper location
    sGenes = sorted(Genes.iteritems(), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    scaff2genes = {}
    for k,v in sortedGenes.items():
        if not v['contig'] in scaff2genes:
            scaff2genes[v['contig']] = [k]
        else:
            scaff2genes[v['contig']].append(k)
    #now have scaffolds dict and gene dict, loop through scaff dict printing tbl
    for k,v in natsorted(scaff2genes.items()):
        sys.stdout.write('>Feature %s\n' % k)
        sys.stdout.write('1\t%s\tREFERENCE\n' % scaffDict.get(k))
        sys.stdout.write('\t\t\t%s\t%s\n' % ('CFMR', '12345'))
        for genes in v: #now loop through each gene on the scaffold
            geneInfo = geneDict.get(genes)
            #check for partial models
            if True in geneInfo['partialStart']:
                partialStart = '<'
            else:
                partialStart = ''
            if True in geneInfo['partialStop']:
                partialStop = '>'
            else:
                partialStop = ''
            #now write gene model
            if geneInfo['strand'] == '+':
                sys.stdout.write('%s%i\t%s%i\tgene\n' % (partialStart, geneInfo['location'][0], partialStop, geneInfo['location'][1]))
            else:
                sys.stdout.write('%s%i\t%s%i\tgene\n' % (partialStart, geneInfo['location'][1], partialStop, geneInfo['location'][0]))
            if geneInfo['name']:
            	sys.stdout.write('\t\t\tgene\t%s\n' % geneInfo['name'])
            sys.stdout.write('\t\t\tlocus_tag\t%s\n' % genes)                                 
            #support multiple transcripts
            for i in range(0,len(geneInfo['ids'])):
                if geneInfo['type'] == 'mRNA':
                    if geneInfo['partialStart'][i] == False:
                        partialStart = ''
                    else:
                        partialStart = '<'
                    if geneInfo['partialStop'][i] == False:
                        partialStop = ''
                    else:
                        partialStop = '>'
                    if geneInfo['strand'] == '+':
                        for num, exon in enumerate(geneInfo['mRNA'][i]):
                            if num == 0 and num == len(geneInfo['mRNA'][i]) - 1: #single exon, so slightly differnt method
                                sys.stdout.write('%s%s\t%s%s\tmRNA\n' % (partialStart, exon[0], partialStop, exon[1]))
                            elif num == 0:
                                sys.stdout.write('%s%s\t%s\tmRNA\n' % (partialStart, exon[0], exon[1]))
                            elif num == len(geneInfo['mRNA'][i]) - 1: #this is last one
                                sys.stdout.write('%s\t%s%s\n' % (exon[0], partialStop, exon[1]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (exon[0], exon[1]))
                        sys.stdout.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                        sys.stdout.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,i+1))
                        sys.stdout.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes, i+1)) 
                        for num, cds in enumerate(geneInfo['CDS'][i]):
                            if num == 0 and num == len(geneInfo['CDS'][i]) - 1: #single exon, so slightly differnt method
                                sys.stdout.write('%s%s\t%s%s\tCDS\n' % (partialStart, cds[0], partialStop, cds[1]))
                            elif num == 0:
                                sys.stdout.write('%s%s\t%s\tCDS\n' % (partialStart, cds[0], cds[1]))
                            elif num == len(geneInfo['CDS'][i]) - 1: #this is last one
                                sys.stdout.write('%s\t%s%s\n' % (cds[0], partialStop, cds[1]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (cds[0], cds[1]))
                        sys.stdout.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'][i])
                        if geneInfo['db_xref'][i]:
                        	for x in geneInfo['db_xref'][i]:
                        		sys.stdout.write('\t\t\tdb_xref\t%s\n' % x)
                        if geneInfo['note'][i]:
                        	for x in geneInfo['note'][i]:
                        		sys.stdout.write('\t\t\tnote\t%s\n' % x)
                        sys.stdout.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                        sys.stdout.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,i+1))
                        sys.stdout.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,i+1))                                       
                    else:
                        for num, exon in enumerate(geneInfo['mRNA'][i]):
                            if num == 0 and num == len(geneInfo['mRNA'][i]) - 1: #single exon, so slightly differnt method
                                sys.stdout.write('%s%s\t%s%s\tmRNA\n' % (partialStart, exon[1], partialStop, exon[0]))
                            elif num == 0:
                                sys.stdout.write('%s%s\t%s\tmRNA\n' % (partialStart, exon[1], exon[0]))
                            elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                sys.stdout.write('%s\t%s%s\n' % (exon[1], partialStop, exon[0]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (exon[1], exon[0]))                 
                        sys.stdout.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                        sys.stdout.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,i+1))
                        sys.stdout.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,i+1))
                        for num, cds in enumerate(geneInfo['CDS'][i]):
                            if num == 0 and num == len(geneInfo['CDS'][i]) - 1: #single exon, so slightly differnt method
                                sys.stdout.write('%s%s\t%s%s\tCDS\n' % (partialStart, cds[1], partialStop, cds[0]))
                            elif num == 0:
                                sys.stdout.write('%s%s\t%s\tCDS\n' % (partialStart, cds[1], cds[0]))
                            elif num == len(geneInfo['CDS'][i]) - 1: #this is last one
                                sys.stdout.write('%s\t%s%s\n' % (cds[1], partialStop, cds[0]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (cds[1], cds[0]))
                        sys.stdout.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'][i])
                        if geneInfo['db_xref'][i]:
                        	for x in geneInfo['db_xref'][i]:
                        		sys.stdout.write('\t\t\tdb_xref\t%s\n' % x)
                        if geneInfo['note'][i]:
                        	for x in geneInfo['note'][i]:
                        		sys.stdout.write('\t\t\tnote\t%s\n' % x)
                        sys.stdout.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                        sys.stdout.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,i+1))
                        sys.stdout.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,i+1))  
                elif geneInfo['type'] == 'tRNA':
                    if geneInfo['strand'] == '+':
                        for num, exon in enumerate(geneInfo['mRNA'][i]):
                            if num == 0:
                                sys.stdout.write('<%s\t>%s\ttRNA\n' % (exon[0], exon[1]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (exon[0], exon[1]))
                        sys.stdout.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                        if geneInfo['product'] == 'tRNA-Xxx':
                            sys.stdout.write('\t\t\tpseudo\n')
                        if geneInfo['note'] != '':
                            sys.stdout.write('\t\t\tnote\t%s\n' % geneInfo['note'][i])                                    
                    else:
                        sys.stdout.write('<%i\t>%i\tgene\n' % (geneInfo['end'], geneInfo['start']))
                        sys.stdout.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA'][i]):
                            if num == 0:
                                sys.stdout.write('<%s\t>%s\ttRNA\n' % (exon[1], exon[0]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (exon[1], exon[0]))
                        sys.stdout.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                        if geneInfo['product'] == 'tRNA-Xxx':
                            sys.stdout.write('\t\t\tpseudo\n')
                        if geneInfo['note'] != '':
                            sys.stdout.write('\t\t\tnote\t%s\n' % geneInfo['note'][i])
                elif geneInfo['type'] == 'rRNA':
                    if geneInfo['strand'] == '+':
                        sys.stdout.write('<%s\t>%s\t%s\n' % (geneInfo['location'][0],geneInfo['location'][1], geneInfo['type']))
                        sys.stdout.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])   
                    else:
                        sys.stdout.write('<%s\t>%s\t%s\n' % (geneInfo['location'][1],geneInfo['location'][0], geneInfo['type']))
                        sys.stdout.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])

#load into dictionary
Genes = {}
Genes = lib.gff2dict(args.gff3, args.fasta, Genes)
scaffLen = scaffold2Dict(args.fasta)
dict2tbl(Genes, scaffLen)
