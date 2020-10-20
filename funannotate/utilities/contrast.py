#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os
import argparse
import itertools
from Bio import SeqIO
from funannotate.interlap import InterLap
from collections import defaultdict
from collections import OrderedDict
from natsort import natsorted
import numpy as np
import pandas as pd

version = '0.0.1'

###### Functions #########


def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList


def fmtcols(mylist, cols):
    justify = []
    for i in range(0, cols):
        length = max([len(x) for x in mylist[i::cols]])
        length += 2
        ljust = [x.ljust(length) for x in mylist[i::cols]]
        justify.append(ljust)
    justify = flatten(justify)
    num_lines = len(mylist) // cols
    lines = (' '.join(justify[i::num_lines])
             for i in range(0, num_lines))
    return "\n".join(lines)


def translate(cDNA, strand, phase):
    '''
    translate cDNA into protein sequence
    trying to see if I can speed this up over Biopython
    '''
    def _RevComp(s):
        rev_comp_lib = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'M': 'K', 'R': 'Y', 'W': 'W',
                        'S': 'S', 'Y': 'R', 'K': 'M', 'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V', 'X': 'X', 'N': 'N'}
        cseq = ''
        n = len(s)
        s = s.upper()
        for i in range(0, n):
            c = s[n-i-1]
            cseq += rev_comp_lib[c]
        return cseq

    def _split(str, num):
        return [str[start:start+num] for start in range(0, len(str), num)]
    codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
                   'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
                   'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
                   'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
                   'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                   'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
                   'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
                   'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
                   'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                   'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
                   'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
                   'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
                   'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}
    if strand == '-' or strand == -1:
        seq = _RevComp(cDNA)
    else:
        seq = cDNA
    seq = seq[phase:]
    # map seq to proteins
    protSeq = []
    for i in _split(seq, 3):
        if len(i) == 3:
            iSeq = i.upper()
            if iSeq in codon_table:
                aa = codon_table[iSeq]
                protSeq.append(aa)
            else:
                protSeq.append('X')
    return ''.join(protSeq)


def getSeqRegions(SeqRecordDict, header, coordinates):
    # takes SeqRecord dictionary or Index, returns sequence string
    # coordinates is a list of tuples [(1,10), (20,30)]
    result = ''
    sorted_coordinates = sorted(coordinates, key=lambda tup: tup[0])
    for x in sorted_coordinates:
        partial = SeqRecordDict[header][x[0]-1:x[1]]
        result += str(partial.seq)
    return result


def getID(input, type):
    # function to get ID from genbank record.features
    locusTag = None
    ID = None
    Parent = None
    if type == 'gene':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
        else:
            try:
                ID = input.qualifiers['gene'][0]
            except KeyError:
                pass
        return locusTag, ID, locusTag

    elif type == 'mRNA' or type == 'tRNA' or type == 'ncRNA' or type == 'rRNA' or type == 'exon':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
            if locusTag:
                Parent = locusTag
                try:
                    ID = input.qualifiers['transcript_id'][0]
                except KeyError:
                    pass
            else:
                try:
                    locusTag = input.qualifiers['transcript_id'][0]
                    Parent = locusTag
                except KeyError:
                    pass
        else:
            try:
                ID = input.qualifiers['transcript_id'][0]
            except KeyError:
                pass
        if ID:
            if ':' in ID:
                ID = ID.split(':')[-1]
        else:
            try:
                ID = input.qualifiers['standard_name'][0]
            except KeyError:
                pass
        return locusTag, ID, Parent

    elif type == 'CDS':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
            if locusTag:
                Parent = locusTag
                try:
                    ID = input.qualifiers['protein_id'][0]
                except KeyError:
                    pass
            else:
                try:
                    locusTag = input.qualifiers['protein_id'][0]
                    Parent = locusTag
                except KeyError:
                    pass
        else:
            try:
                ID = input.qualifiers['protein_id'][0]
            except KeyError:
                pass
        if ID:
            if ':' in ID:
                ID = ID.split(':')[-1]
        else:
            try:
                ID = input.qualifiers['standard_name'][0]
            except KeyError:
                pass
        return locusTag, ID, Parent


def gff2dict(file, fasta, Genes):
    '''
    general function to take a GFF3 file and return a funannotate standardized dictionary
    locustag: {
    'contig': contigName
    'type': mRNA/rRNA/tRNA/ncRNA
    'location': (start, end) #integer tuple
    'strand': +/-
    'ids': [transcript/protein IDs] #list
    'mRNA':[[(ex1,ex1),(ex2,ex2)]] #list of lists of tuples (start, end)
    'CDS':[[(cds1,cds1),(cds2,cds2)]] #list of lists of tuples (start, end)
    'transcript': [seq1, seq2] #list of mRNA trnascripts
    'cds_transcript': [seq1, seq2] #list of mRNA trnascripts (no UTRs)
    'protein': [protseq1,protseq2] #list of CDS translations
    'codon_start': [1,1] #codon start for translations
    'note': [[first note, second note], [first, second, etc]] #list of lists
    'name': genename
    'product': [hypothetical protein, velvet complex] #list of product definitions
    'go_terms': [[GO:0000001,GO:0000002]] #list of lists
    'db_xref': [[InterPro:IPR0001,PFAM:004384]] #list of lists
    'partialStart': True/False
    'partialStop': True/False
    'source': annotation source
    'phase': [[0,2,1]] list of lists
    '5UTR': [[(),()]] #list of lists of tuples (start, end)
    '3UTR': [[(),()]] #list of lists of tuples (start, end)
    }
    '''
    idParent = {}
    with open(file, 'r') as input:
        for line in input:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split(
                '\t')
            start = int(start)
            end = int(end)
            ID, Parent, Name, Product, GeneFeature = (None,)*5
            Note, DBxref, GO = ([],)*3
            info = attributes.split(';')
            for x in info:
                if x.startswith('ID='):
                    ID = x.replace('ID=', '')
                elif x.startswith('Parent='):
                    Parent = x.replace('Parent=', '')
                elif x.startswith('Name='):
                    Name = x.replace('Name=', '')
                elif x.startswith('Note=') or x.startswith('note='):
                    Note = x.split('ote=')[-1]
                    if ',' in Note:
                        Note = Note.split(',')
                    else:
                        Note = [Note]
                elif x.startswith('DBxref='):
                    DBxref = x.replace('DBxref=', '')
                    if ',' in DBxref:
                        DBxref = DBxref.split(',')
                    else:
                        DBxref = [DBxref]
                elif x.startswith('Ontology_term='):
                    GO = x.replace('Ontology_term=', '')
                    if ',' in GO:
                        GO = GO.split(',')
                    else:
                        GO = [GO]
                elif x.startswith('Product=') or x.startswith('product='):
                    Product = x.split('roduct=')[-1]
                elif x.startswith('description='):
                    Product = x.replace('description=', '')
            if feature == 'gene':
                if not ID in Genes:
                    Genes[ID] = {'name': Name, 'type': None, 'transcript': [], 'cds_transcript': [], 'protein': [], '5UTR': [], '3UTR': [],
                                 'codon_start': [], 'ids': [], 'CDS': [], 'mRNA': [], 'strand': strand,
                                 'location': (start, end), 'contig': contig, 'product': [], 'source': source, 'phase': [],
                                 'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [], 'partialStop': [], 'pseudo': False}
                else:
                    if start < Genes[ID]['location'][0]:
                        Genes[ID]['location'] = (
                            start, Genes[ID]['location'][1])
                    if end > Genes[ID]['location'][1]:
                        Genes[ID]['location'] = (Genes[ID]['location'][0], end)
            else:
                if not ID or not Parent:
                    print("Error, can't find ID or Parent. Malformed GFF file.")
                    print(line)
                    sys.exit(1)
                if feature == 'mRNA' or feature == 'tRNA' or feature == 'rRNA':
                    if not Product:
                        if feature == 'mRNA':
                            Product = 'hypothetical protein'
                    if not Parent in Genes:
                        Genes[Parent] = {'name': Name, 'type': feature, 'transcript': [], 'cds_transcript': [], 'protein': [], '5UTR': [[]], '3UTR': [[]],
                                         'codon_start': [[]], 'ids': [ID], 'CDS': [[]], 'mRNA': [[]], 'strand': strand,
                                         'location': (start, end), 'contig': contig, 'product': [Product], 'source': source, 'phase': [[]],
                                         'db_xref': [DBxref], 'go_terms': [GO], 'note': [Note], 'partialStart': [False], 'partialStop': [False], 'pseudo': False}
                    else:
                        Genes[Parent]['ids'].append(ID)
                        Genes[Parent]['mRNA'].append([])
                        Genes[Parent]['CDS'].append([])
                        Genes[Parent]['phase'].append([])
                        Genes[Parent]['5UTR'].append([])
                        Genes[Parent]['3UTR'].append([])
                        Genes[Parent]['codon_start'].append([])
                        Genes[Parent]['partialStart'].append(False)
                        Genes[Parent]['partialStop'].append(False)
                        Genes[Parent]['product'].append(Product)
                        Genes[Parent]['db_xref'].append(DBxref)
                        Genes[Parent]['go_terms'].append(GO)
                        Genes[Parent]['note'].append(Note)
                        Genes[Parent]['type'] = feature
                        # double check mRNA features are contained in gene coordinates
                        if start < Genes[Parent]['location'][0]:
                            #print('{:} update start: {:} to {:}'.format(Parent, Genes[Parent]['location'][0],start))
                            Genes[Parent]['location'] = (
                                start, Genes[Parent]['location'][1])
                        if end > Genes[Parent]['location'][1]:
                            #print('{:} update stop: {:} to {:}'.format(Parent, Genes[Parent]['location'][1],end))
                            Genes[Parent]['location'] = (
                                Genes[Parent]['location'][0], end)
                    if not ID in idParent:
                        idParent[ID] = Parent
                elif feature == 'exon':
                    if ',' in Parent:
                        parents = Parent.split(',')
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {'name': Name, 'type': None, 'transcript': [], 'cds_transcript': [], 'protein': [], '5UTR': [[]], '3UTR': [[]],
                                                      'codon_start': [[]], 'ids': [p], 'CDS': [], 'mRNA': [[(start, end)]], 'strand': strand,
                                                      'location': None, 'contig': contig, 'product': [], 'source': source, 'phase': [[]],
                                                      'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [False], 'partialStop': [False], 'pseudo': False}
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]['ids'].index(p)
                                Genes[GeneFeature]['mRNA'][i].append(
                                    (start, end))
                elif feature == 'CDS':
                    if ',' in Parent:
                        parents = Parent.split(',')
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {'name': Name, 'type': None, 'transcript': [], 'cds_transcript': [], 'protein': [], '5UTR': [[]], '3UTR': [[]],
                                                      'codon_start': [[]], 'ids': [p], 'CDS': [[(start, end)]], 'mRNA': [], 'strand': strand,
                                                      'location': None, 'contig': contig, 'product': [], 'source': source, 'phase': [[]],
                                                      'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [False], 'partialStop': [False], 'pseudo': False}
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]['ids'].index(p)
                                Genes[GeneFeature]['CDS'][i].append(
                                    (start, end))
                                # add phase
                                Genes[GeneFeature]['phase'][i].append(
                                    int(phase))
                elif feature == 'five_prime_UTR' or feature == 'five_prime_utr':
                    if ',' in Parent:
                        parents = Parent.split(',')
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {'name': Name, 'type': None, 'transcript': [], 'cds_transcript': [], 'protein': [], '5UTR': [[(start, end)]], '3UTR': [[]],
                                                      'codon_start': [[]], 'ids': [p], 'CDS': [], 'mRNA': [[(start, end)]], 'strand': strand,
                                                      'location': None, 'contig': contig, 'product': [], 'source': source, 'phase': [[]],
                                                      'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [False], 'partialStop': [False], 'pseudo': False}
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]['ids'].index(p)
                                Genes[GeneFeature]['5UTR'][i].append(
                                    (start, end))
                elif feature == 'three_prime_UTR' or feature == 'three_prime_utr':
                    if ',' in Parent:
                        parents = Parent.split(',')
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {'name': Name, 'type': None, 'transcript': [], 'cds_transcript': [], 'protein': [], '5UTR': [[]], '3UTR': [[(start, end)]],
                                                      'codon_start': [[]], 'ids': [p], 'CDS': [], 'mRNA': [[(start, end)]], 'strand': strand,
                                                      'location': None, 'contig': contig, 'product': [], 'source': source, 'phase': [[]],
                                                      'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [False], 'partialStop': [False], 'pseudo': False}
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]['ids'].index(p)
                                Genes[GeneFeature]['3UTR'][i].append(
                                    (start, end))
    # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
    SeqRecords = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    for k, v in list(Genes.items()):
        for i in range(0, len(v['ids'])):
            if v['type'] == 'mRNA' or v['type'] == 'tRNA':
                if v['strand'] == '+':
                    sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0])
                else:
                    sortedExons = sorted(
                        v['mRNA'][i], key=lambda tup: tup[0], reverse=True)
                Genes[k]['mRNA'][i] = sortedExons
                mrnaSeq = getSeqRegions(SeqRecords, v['contig'], sortedExons)
                v['transcript'].append(mrnaSeq)
            if v['type'] == 'mRNA':
                if v['strand'] == '+':
                    sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0])
                else:
                    sortedCDS = sorted(
                        v['CDS'][i], key=lambda tup: tup[0], reverse=True)
                # get the codon_start by getting first CDS phase + 1
                indexStart = [x for x, y in enumerate(
                    v['CDS'][i]) if y[0] == sortedCDS[0][0]]
                codon_start = int(v['phase'][i][indexStart[0]]) + 1
                Genes[k]['codon_start'][i] = codon_start
                Genes[k]['CDS'][i] = sortedCDS
                # translate and get protein sequence
                protSeq = None
                cdsSeq = getSeqRegions(SeqRecords, v['contig'], v['CDS'][i])
                v['cds_transcript'].append(cdsSeq)
                protSeq = translate(cdsSeq, v['strand'], v['codon_start'][i]-1)
                v['protein'].append(protSeq)
                if protSeq:
                    if protSeq.endswith('*'):
                        v['partialStop'][i] = False
                    else:
                        v['partialStop'][i] = True
                    if v['codon_start'][i] == 1 and v['protein'][i].startswith('M'):
                        v['partialStart'][i] = False
                    else:
                        v['partialStart'][i] = True
    return Genes


def gb_feature_add2dict(f, record, genes):
    '''
    general function to take a genbank feature from flat file and add to funannotate standardized dictionary
    locustag: {
    'contig': contigName
    'type': mRNA/rRNA/tRNA/ncRNA
    'location': (start, end) #integer tuple
    'strand': +/-
    'ids': [transcript/protein IDs] #list
    'mRNA':[[(ex1,ex1),(ex2,ex2)]] #list of lists of tuples (start, end)
    'CDS':[[(cds1,cds1),(cds2,cds2)]] #list of lists of tuples (start, end)
    'transcript': [seq1, seq2] #list of mRNA trnascripts
    'cds_transcript': [seq1, seq2] list of mRNA (no UTRs)
    'protein': [protseq1,protseq2] #list of CDS translations
    'codon_start': [1,1] #codon start for translations
    'note': [[first note, second note], [first, second, etc]] #list of lists
    'name': genename
    'product': [hypothetical protein, velvet complex] #list of product definitions
    'go_terms': [[GO:0000001,GO:0000002]] #list of lists
    'db_xref': [[InterPro:IPR0001,PFAM:004384]] #list of lists
    'partialStart': True/False
    'partialStop': True/False
    'source': annotation source
    }
    '''
    # get info from features, if there is no locusTag then exit
    if f.type == 'gene' or f.type == 'mRNA' or f.type == 'CDS' or f.type == 'tRNA' or f.type == 'rRNA' or f.type == 'ncRNA':
        locusTag, ID, Parent = getID(f, f.type)
        if not locusTag:
            return genes
    else:
        return genes
    # standard information from every feature
    strand = f.location.strand
    if strand == 1:
        strand = '+'
    elif strand == -1:
        strand = '-'
    start = f.location.nofuzzy_start + 1
    end = f.location.nofuzzy_end
    chr = record.id
    num_parts = len(f.location.parts)
    name, Product = (None,)*2
    Fivepartial, Threepartial = (False,)*2
    DBxref = []
    Note = []
    GO = []
    # parse each type somewhat differently
    if f.type == 'gene':
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': None, 'transcript': [], 'cds_transcript': [], 'protein': [], 'source': 'GenBank',
                               'codon_start': [], 'ids': [], 'CDS': [], 'mRNA': [], 'strand': strand,
                               'location': (int(start), int(end)), 'contig': chr, 'product': [],
                               'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [], 'partialStop': []}
        else:
            genes[locusTag]['location'] = (int(start), int(end))
            genes[locusTag]['strand'] = strand
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    elif f.type == 'tRNA' or f.type == 'rRNA' or f.type == 'ncRNA':
        feature_seq = f.extract(record.seq)
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        try:
            Product = f.qualifiers['product'][0]
            if Product == 'tRNA-OTHER':
                Product = 'tRNA-Xxx'
        except KeyError:
            Product = None
        exonTuples = []
        if num_parts < 2:  # only single exon
            exonTuples.append((int(start), int(end)))
        else:  # more than 1 exon, so loop through
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                exonTuples.append((int(ex_start), int(ex_end)))
        # now we want to sort the positions I think...
        if strand == '+':
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
            if str(f.location.start).startswith('<'):
                Fivepartial = True
            if str(f.location.end).startswith('>'):
                Threepartial = True
        else:
            sortedExons = sorted(
                exonTuples, key=lambda tup: tup[0], reverse=True)
            if str(f.location.start).startswith('<'):
                Threepartial = True
            if str(f.location.end).startswith('>'):
                Fivepartial = True
        # update positions
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': f.type, 'transcript': [feature_seq], 'cds_transcript': [], 'protein': [], 'source': 'GenBank',
                               'codon_start': [], 'ids': [locusTag+'-T1'], 'CDS': [], 'mRNA': [sortedExons], 'strand': strand,
                               'location': (int(start), int(end)), 'contig': chr, 'product': [Product],
                               'db_xref': [DBxref], 'go_terms': [GO], 'note': [Note], 'partialStart': [Fivepartial], 'partialStop': [Threepartial]}
        else:
            genes[locusTag]['mRNA'].append(sortedExons)
            genes[locusTag]['type'] = f.type
            genes[locusTag]['transcript'].append(feature_seq)
            genes[locusTag]['ids'].append(
                locusTag+'-T'+str(len(genes[locusTag]['ids'])+1))
            genes[locusTag]['db_xref'].append(DBxref)
            genes[locusTag]['note'].append(Note)
            genes[locusTag]['go_terms'].append(GO)
            genes[locusTag]['product'].append(Product)
            genes[locusTag]['partialStart'].append(Fivepartial)
            genes[locusTag]['partialStop'].append(Threepartial)
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    elif f.type == 'mRNA':
        feature_seq = f.extract(record.seq)
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        exonTuples = []
        if num_parts < 2:  # only single exon
            exonTuples.append((int(start), int(end)))
        else:  # more than 1 exon, so loop through
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                exonTuples.append((int(ex_start), int(ex_end)))
        # now we want to sort the positions I think...
        if strand == '+':
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
            if str(f.location.start).startswith('<'):
                Fivepartial = True
            if str(f.location.end).startswith('>'):
                Threepartial = True
        else:
            sortedExons = sorted(
                exonTuples, key=lambda tup: tup[0], reverse=True)
            if str(f.location.start).startswith('<'):
                Threepartial = True
            if str(f.location.end).startswith('>'):
                Fivepartial = True
        # update positions
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': f.type, 'transcript': [feature_seq], 'cds_transcript': [], 'protein': [], 'source': 'GenBank',
                               'codon_start': [], 'ids': [], 'CDS': [], 'mRNA': [sortedExons], 'strand': strand,
                               'location': (int(start), int(end)), 'contig': chr, 'product': [],
                               'db_xref': [], 'go_terms': [], 'note': [], 'partialStart': [Fivepartial], 'partialStop': [Threepartial]}
        else:
            genes[locusTag]['mRNA'].append(sortedExons)
            genes[locusTag]['type'] = f.type
            genes[locusTag]['transcript'].append(feature_seq)
            genes[locusTag]['partialStart'].append(Fivepartial)
            genes[locusTag]['partialStop'].append(Threepartial)
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    elif f.type == 'CDS':
        feature_seq = f.extract(record.seq)
        if not ID:
            print(("putative transcript from %s has no ID\n%s" %
                     (locusTag, genes[locusTag])))
            return genes
        try:
            protSeq = f.qualifiers['translation'][0]
        except KeyError:
            print(("%s has no translation" % ID))
            protSeq = ''
        cdsTuples = []
        phase = int(f.qualifiers['codon_start'][0])
        if num_parts < 2:  # only single CDS
            cdsTuples.append((int(start), int(end)))
        else:
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                cdsTuples.append((int(ex_start), int(ex_end)))
        if strand == '+':
            sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0])
        else:
            sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0], reverse=True)
        # check for annotations
        try:
            Product = f.qualifiers['product'][0]
        except KeyError:
            Product = 'hypothetical protein'
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        # note and dbxref are in a dictionary
        for key, value in list(f.qualifiers.items()):
            if key == 'note':
                notes = value[0].split('; ')
                for n in notes:
                    if n.startswith('GO'):
                        GO.append(n)
                    else:
                        Note.append(n)
            elif key == 'db_xref':
                for ref in value:
                    DBxref.append(ref)
        # update dictionary
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': None, 'transcript': [], 'cds_transcript': [feature_seq], 'protein': [], 'source': 'GenBank',
                               'codon_start': [phase], 'ids': [ID], 'CDS': [sortedCDS], 'mRNA': [], 'strand': strand,
                               'location': (int(start), int(end)), 'contig': chr, 'product': [Product],
                               'db_xref': [DBxref], 'go_terms': [GO], 'note': [Note], 'partialStart': [], 'partialStop': []}
        else:
            genes[locusTag]['ids'].append(ID)
            genes[locusTag]['CDS'].append(sortedCDS)
            genes[locusTag]['product'].append(Product)
            genes[locusTag]['protein'].append(protSeq)
            genes[locusTag]['cds_transcript'].append(feature_seq)
            genes[locusTag]['codon_start'].append(phase)
            genes[locusTag]['db_xref'].append(DBxref)
            genes[locusTag]['note'].append(Note)
            genes[locusTag]['go_terms'].append(GO)
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    return genes


def gbk2interlap(input):
    '''
    function to parse GBK file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    '''
    inter = defaultdict(InterLap)
    Genes = {}
    with open(input, 'r') as filein:
        for record in SeqIO.parse(filein, 'genbank'):
            for f in record.features:
                if f.type == 'gene':
                    locusTag, ID, Parent = getID(f, f.type)
                    start = int(f.location.nofuzzy_start)
                    end = int(f.location.nofuzzy_end)
                    inter[record.id].add((start, end, locusTag))
                gb_feature_add2dict(f, record, Genes)
    return inter, Genes


def gff2interlap(input, fasta):
    '''
    function to parse GFF3 file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    '''
    inter = defaultdict(InterLap)
    Genes = {}
    Genes = gff2dict(input, fasta, Genes)
    for k, v in natsorted(list(Genes.items())):
        inter[v['contig']].add((v['location'][0], v['location'][1], k))
    return inter, Genes


def message(loc1, loc2, cdsAED, mrnaAED, protMatches, UTRs, no_change, UTR_added, yardSale, exonChange):
    msg = []
    if not cdsAED or cdsAED == '':
        cds = 0
    else:
        cds = float(cdsAED)
    mrna = float(mrnaAED)
    pos = loc1 == loc2
    # structured message, coordinates, 5prime, 3prime, exon, cds, pident
    if not pos:  # coordinates changed
        msg.append('gene coordinates changed')
    if mrna > 0:
        msg.append('mRNA changed')
    if cds > 0:
        pidentmsg = []
        if protMatches:
            for x in protMatches:
                pidentmsg.append('{0:.0f}%'.format(x))
            msg.append(
                'CDS changed [translation pident: %s]' % ', '.join(pidentmsg))
        msg.append('CDS changed [translation pident: NA]')
    # now work on the counter
    if len(msg) < 1:
        msg = ['no change']
        no_change += 1
    elif any('UTR' in x for x in msg):
        UTR_added += 1
    elif any('mRNA' in x for x in msg):
        exonChange += 1
    else:
        yardSale += 1
    final_message = ';'.join(msg)
    return final_message, no_change, UTR_added, yardSale, exonChange


def pairwiseAlign(query, ref):
    from Bio import pairwise2
    '''
    do global alignment and return pident
    '''
    if query == ref:
        return 100.0, len(ref), len(ref)
    align = pairwise2.align.globalxx(query, ref)
    length = max(len(query), len(ref))
    pident = (align[0][2] / float(length)) * 100
    return pident, align[0][2], len(ref)


def countFeatures(input):
    # given funannotate dictionary, count up some general features
    mRNAg, mRNAt, tRNAg, tRNAt = (0,)*4
    for k, v in natsorted(list(input.items())):
        if v['type'] == 'mRNA':
            mRNAg += 1
            mRNAt += len(v['ids'])
        elif v['type'] == 'tRNA':
            tRNAg += 1
            tRNAt += len(v['ids'])
    return len(input), mRNAg, mRNAt, tRNAg, tRNAt


def compareAnnotations(old, oldformat, new, newformat, fasta, measure_pident, output):
    '''
    function takes two GenBank annotated genomes and compares gene models
    output is a tsv file for each locus and a description of what is different
    can handle multiple transcripts per locus
    '''
    result = {}
    global no_change, identicalCDS, refUnique, queryUnique
    no_change, identicalCDS, refUnique, queryUnique, totalmatches, totallength = (
        0,)*6
    if oldformat == 'gff':
        oldInter, oldGenes = gff2interlap(old, fasta)
    else:
        oldInter, oldGenes = gbk2interlap(old)
    if newformat == 'gff':
        newInter, newGenes = gff2interlap(new, fasta)
    else:
        newInter, newGenes = gbk2interlap(new)
    NumOldLoci, NumOldGenes, NumOldmRNA, NumOldtRNALoci, NumOldtRNA = countFeatures(
        oldGenes)
    NumNewLoci, NumNewGenes, NumNewmRNA, NumNewtRNALoci, NumNewtRNA = countFeatures(
        newGenes)

    # now run some comparisons
    # do the simple stuff first, find models that were deleted
    for contig in oldInter:
        for gene in oldInter[contig]:
            if not gene in newInter[contig]:  # these models are removed
                if not gene[2] in oldGenes:
                    continue
                # populate output dictionary with results
                if not gene[2] in result:
                    refUnique += 1
                    result[gene[2]] = {'contig': oldGenes[gene[2]]['contig'], 'location': oldGenes[gene[2]]['location'], 'ref_type': oldGenes[gene[2]]['type'], 'ref_location': oldGenes[gene[2]]['location'],
                                       'query_location': None, 'query_id': None, 'query_type': None, 'pident': None, 'ref_id': gene[2],
                                       'cdsAED': '1.000', 'exonAED': '1.000', 'ref_transcripts': len(oldGenes[gene[2]]['ids']), 'query_transcripts': 0,
                                       'ref_strand': oldGenes[gene[2]]['strand'], 'query_strand': None}
                    for x in oldGenes[gene[2]]['protein']:
                        totallength += len(x)

    # now go through the updated annotation, comparing to old annot
    for contig in newInter:
        for gene in newInter[contig]:
            # means this is a new model, so add it
            if not gene in oldInter[contig]:
                if not gene[2] in newGenes:
                    continue
                if not gene[2] in result:
                    queryUnique += 1
                    result[gene[2]] = {'contig': newGenes[gene[2]]['contig'], 'location': newGenes[gene[2]]['location'], 'ref_type': None, 'ref_location': None,
                                       'query_location': newGenes[gene[2]]['location'], 'query_id': gene[2], 'query_type': newGenes[gene[2]]['type'], 'pident': None,
                                       'cdsAED': '1.000', 'exonAED': '1.000', 'ref_transcripts': 0, 'query_transcripts': len(newGenes[gene[2]]['ids']), 'ref_id': None,
                                       'ref_strand': None, 'query_strand': newGenes[gene[2]]['strand']}
            else:  # means this is existing model, and need to do some comparisons
                hitList = list(oldInter[contig].find(gene))
                # there might be some overlapping transcripts, so get best hit?
                hit = []
                # get best hit
                for z in hitList:
                    diffs = np.subtract((gene[0], gene[1]), (z[0], z[1]))
                    totaldiffs = abs(diffs[0]) + abs(diffs[1])
                    hit.append((totaldiffs, z[2]))
                besthit = min(hit)

                # get the old annotation
                hitInfo = oldGenes.get(besthit[1])

                # calculate AED
                exonAED = pairwiseAED(
                    newGenes[gene[2]]['mRNA'], hitInfo['mRNA'])
                if newGenes[gene[2]]['type'] == 'mRNA' and hitInfo['type'] == 'mRNA':
                    cdsAED = pairwiseAED(
                        newGenes[gene[2]]['CDS'], hitInfo['CDS'])
                else:
                    cdsAED = '0.000'

                # check translation, to deal with multiple transcripts, lets loop through new
                if measure_pident:
                    protMatches = []
                    if newGenes[gene[2]]['type'] == 'mRNA' and hitInfo['type'] == 'mRNA':
                        for i in range(0, len(newGenes[gene[2]]['ids'])):
                            protMatch = None
                            for y in range(0, len(hitInfo['ids'])):
                                pident, matches, length = pairwiseAlign(
                                    newGenes[gene[2]]['protein'][i], hitInfo['protein'][y])
                                if not protMatch:
                                    protMatch = (pident, matches, length)
                                else:
                                    if pident > protMatch[0]:
                                        protMatch = (pident, matches, length)
                            protMatches.append(protMatch[0])
                            totalmatches += protMatch[1]
                            totallength += protMatch[2]
                else:
                    protMatches = None

                if not besthit[1] in result:
                    result[besthit[1]] = {'contig': newGenes[gene[2]]['contig'], 'location': hitInfo['location'], 'ref_type': hitInfo['type'], 'ref_location': hitInfo['location'],
                                          'query_location': newGenes[gene[2]]['location'], 'query_id': gene[2], 'query_type': newGenes[gene[2]]['type'], 'pident': protMatches,
                                          'cdsAED': cdsAED, 'exonAED': exonAED, 'ref_transcripts': len(hitInfo['ids']), 'query_transcripts': len(newGenes[gene[2]]['ids']),
                                          'ref_strand': hitInfo['strand'], 'query_strand': newGenes[gene[2]]['strand'], 'ref_id': besthit[1]}
                # get some summary stats as you loop through
                if float(exonAED) == 0 and float(cdsAED) == 0:
                    no_change += 1
                elif float(cdsAED) == 0:
                    identicalCDS += 1

    total_cdsAED = []
    total_exonAED = []

    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = sorted(iter(result.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    with open(output, 'w') as out:
        out.write('Reference_Location\tReference_ID\tRef_strand\tRef_Num_Transcripts\tQuery_Location\tQuery_ID\tQuery_strand\tQuery_Num_Transcripts\tmRNA_AED\tCDS_AED\n')
        for k, v in list(sortedGenes.items()):
            Rstart = str(v['location'][0])
            Rend = str(v['location'][1])
            if v['query_id']:
                Qstart = str(v['query_location'][0])
                Qend = str(v['query_location'][1])
            else:
                Qstart = 'None'
                Qend = 'None'
            total_cdsAED.append(float(v['cdsAED']))
            total_exonAED.append(float(v['exonAED']))
            out.write('{:}:{:}-{:}\t{:}\t{:}\t{:}\t{:}:{:}-{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                v['contig'], Rstart, Rend, v['ref_id'], v['ref_strand'], v['ref_transcripts'], v['contig'], Qstart, Qend, v['query_id'], v['query_strand'], v['query_transcripts'], v['exonAED'], v['cdsAED']))
    Avg_cdsAED = sum(total_cdsAED) / float(len(total_cdsAED))
    Avg_exonAED = sum(total_exonAED) / float(len(total_exonAED))
    totalPident = 0.00
    if totalmatches > 0:
        totalPident = (totalmatches / totallength)
    return [NumOldLoci, NumOldGenes, NumOldmRNA, NumOldtRNALoci, NumOldtRNA, refUnique, no_change, identicalCDS, 0.000, 0.000, 1, NumNewLoci, NumNewGenes, NumNewmRNA, NumNewtRNALoci, NumNewtRNA, queryUnique, no_change, identicalCDS, Avg_exonAED, Avg_cdsAED, totalPident]


def findUTRs(cds, mrna, strand):
    '''
    take list of list of CDS coordiantes and compare to list of list of mRNA coordinates to
    determine if 5 prime or 3 prime UTR exist
    '''
    # supporting multiple transcripts, however, they are already matched up and sorted
    UTRs = []
    for i in range(0, len(cds)):
        Fiveprime = False
        Threeprime = False
        refInterlap = InterLap(mrna[i])
        if strand == '+':  # look at first CDS for 5 prime and last CDS for 3 prime
            # means it overlaps with mrNA (which it obviously should)
            if cds[i][0] in refInterlap:
                hit = list(refInterlap.find(cds[i][0]))[0]
                # if first exon, then compare, if not first then there is 5prime UTR
                loc = mrna[i].index(hit)
                if loc == 0:
                    # will return array of exon minus hit at each pos
                    diff = np.subtract(cds[i][0], hit)
                    if diff[0] > 0:
                        Fiveprime = True
                else:
                    Fiveprime = True
            # check for 3 prime UTR
            if cds[i][-1] in refInterlap:
                hit = list(refInterlap.find(cds[i][-1]))[0]
                loc = mrna[i].index(hit)
                if len(mrna[i]) == loc+1:
                    # will return array of exon minus hit at each pos
                    diff = np.subtract(cds[i][-1], hit)
                    if diff[1] < 0:
                        Threeprime = True
                else:
                    Threeprime = True
        else:
            # means it overlaps with mrNA (which it obviously should)
            if cds[i][0] in refInterlap:
                hit = list(refInterlap.find(cds[i][0]))[0]
                # if first exon, then compare, if not first then there is 5prime UTR
                loc = mrna[i].index(hit)
                if loc == 0:
                    # will return array of exon minus hit at each pos
                    diff = np.subtract(cds[i][0], hit)
                    if diff[1] < 0:
                        Fiveprime = True
                else:
                    Fiveprime = True
            # check for 3 prime UTR
            if cds[i][-1] in refInterlap:
                hit = list(refInterlap.find(cds[i][-1]))[0]
                loc = mrna[i].index(hit)
                if len(mrna[i]) == loc+1:
                    # will return array of exon minus hit at each pos
                    diff = np.subtract(cds[i][-1], hit)
                    if diff[0] > 0:
                        Threeprime = True
                else:
                    Threeprime = True
        UTRs.append((Fiveprime, Threeprime))
    return UTRs


def pairwiseAED(query, reference):
    '''
    takes a multiple transcripts and sums AED from lowest pairwise comparison and then calculates
    the average based on number of transcripts in the query
    '''
    AEDsum = []
    pAED = [float(getAED(a, b))
            for a, b in itertools.product(query, reference)]
    # split into parts to get lowest AED
    splitAED = [pAED[i:i+len(query)] for i in range(0, len(pAED), len(query))]
    for pair in splitAED:
        AEDsum.append(min(pair))
    AEDavg = sum(AEDsum) / len(query)
    return '{:.3f}'.format(AEDavg)


def getAED(query, reference):
    '''
    function to calcuate annotation edit distance between two mRNA transcript coordinates
    AED = 1 - (SN + SP / 2)
    SN = fraction of ref predicted
    SP = fraction prediction overlapping the ref
    '''
    def _length(listTup):
        len = 0
        for i in listTup:
            l = abs(i[0] - i[1])
            len += l
        return len
    # check if identical
    if query == reference:
        return '0.000'
    # make sure sorted
    rLen = _length(reference)
    refInterlap = InterLap(reference)
    QueryOverlap = 0
    qLen = 0
    for exon in query:
        qLen += abs(exon[0] - exon[1])
        if exon in refInterlap:  # exon overlaps at least partially with reference
            hit = list(refInterlap.find(exon))
            for h in hit:
                # will return array of exon minus hit at each pos
                diff = np.subtract(exon, h)
                if diff[0] <= 0 and diff[1] >= 0:  # then query exon covers ref exon
                    cov = abs(h[0] - h[1])
                    QueryOverlap += cov
                elif diff[0] <= 0 and diff[1] < 0:  # means query partial covers ref
                    cov = abs(h[0] - exon[1])
                    QueryOverlap += cov
                elif diff[0] > 0 and diff[1] >= 0:  # means query partial covers ref
                    cov = abs(exon[0] - h[1])
                    QueryOverlap += cov
                elif diff[0] > 0 and diff[1] < 1:
                    cov = abs(exon[0] - exon[1])
                    QueryOverlap += cov
    # calculate AED
    if qLen > 0 and rLen > 0:
        SP = QueryOverlap / float(qLen)
        SN = QueryOverlap / float(rLen)
        AED = 1 - ((SN + SP) / 2)
    else:
        AED = 0.000
    return '{:.3f}'.format(AED)


def main(args):
    parser = argparse.ArgumentParser(prog='compare2annotations.py', usage="%(prog)s [options] -q query_annotation -r ref_annotation -o output",
                                     description='''Script is compares query annotation to reference annotation.''',
                                     epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""")
    parser.add_argument('-q', '--query', nargs='+',
                        required=True, help='Genome annotation GBK or GFF3')
    parser.add_argument('-r', '--reference', required=True,
                        help='Genome annotation GBK or GFF3')
    parser.add_argument(
        '-f', '--fasta', help='Genome sequence in FASTA format (if using GFF3)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output comparison file basename')
    parser.add_argument('-c', '--calculate_pident', action='store_true',
                        help='Calculate protein pident at each overlap')
    args = parser.parse_args(args)

    # goal here is to populate annotation into query and reference dictionaries
    # then run comparisons between the query and reference.
    # support multiple queries, so loop through each, build the output text file here
    CombinedResults = {'Stats': ['Total Genes', 'Coding Genes', 'Coding transcripts', 'tRNA Genes', 'tRNA transcripts',
                                 'Unique Genes', 'Identical Genes', 'Identical CDS', 'Avg mRNA AED', 'Avg CDS AED', 'Protein pident']}
    if args.reference.endswith('.gbk') or args.reference.endswith('.gbff') or args.reference.endswith('.gb'):
        Refformat = 'genbank'
    elif args.reference.endswith('.gff3') or args.reference.endswith('.gff'):
        if not args.fasta:
            print('ERROR: Need --fasta if using GFF3 input.')
            sys.exit(1)
        Refformat = 'gff'
    for q in args.query:
        if q.endswith('.gbk') or q.endswith('.gbff') or q.endswith('.gb'):
            Queryformat = 'genbank'
        elif q.endswith('.gff3') or q.endswith('.gff'):
            if not args.fasta:
                print('ERROR: Need --fasta if using GFF3 input.')
                sys.exit(1)
            Queryformat = 'gff'
        # now run comparison and output annotation file
        qbase = os.path.basename(q).rsplit('.', 1)[0]
        if args.output:
            runOutput = args.output + '.' + qbase + '.compare2ref.txt'
        else:
            runOutput = qbase + '.compare2ref.txt'
        result = compareAnnotations(
            args.reference, Refformat, q, Queryformat, args.fasta, args.calculate_pident, runOutput)
        if not 'Reference' in CombinedResults:
            RefResults = result[:11]
            if len(args.query) > 1:
                RefResults[5] = 0
                RefResults[6] = 0
                RefResults[7] = 0
            CombinedResults['Reference'] = RefResults
        CombinedResults[qbase] = result[11:]

    # get together stats from each query:
    df = pd.DataFrame(CombinedResults)
    df.set_index('Stats', inplace=True)
    dfT = df.transpose()
    for x in ['Total Genes', 'Coding Genes', 'Coding transcripts', 'tRNA Genes', 'tRNA transcripts', 'Unique Genes', 'Identical Genes', 'Identical CDS']:
        dfT[x] = pd.Series(["{0:,.0f}".format(val)
                            for val in dfT[x]], index=dfT.index)
    for x in ['Avg mRNA AED', 'Avg CDS AED']:
        dfT[x] = pd.Series(["{0:.3f}".format(val)
                            for val in dfT[x]], index=dfT.index)
    dfT['Protein pident'] = pd.Series(
        ["{0:.1%}".format(val) for val in dfT['Protein pident']], index=dfT.index)
    df2 = dfT.transpose()
    cols = list(df2)
    cols.insert(0, cols.pop(cols.index('Reference')))
    df2 = df2.ix[:, cols]
    df2.to_csv(args.output+'.summary-stats.csv', sep=',')
    print('--------------------------------------------------------------')
    print((df2.to_string(justify='center')))
    print('--------------------------------------------------------------')


if __name__ == "__main__":
    main(sys.argv[1:])
