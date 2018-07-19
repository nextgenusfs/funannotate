#!/usr/bin/env python

import sys, argparse

#setup menu with argparse
parser = argparse.ArgumentParser(prog='stringtie2gff.py', 
    description = '''Script to convert StringTie GTF file to GFF3.''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""")
parser.add_argument('-i', '--input', required=True, help='StringTie GTF file')
args=parser.parse_args()

def gtf2dict(input):
    Genes = {}
    with open(input,'rU') as inFile:
        for line in inFile:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.rstrip()
            #CM002242   StringTie   transcript  4198460 4199001 1000    +   .   gene_id "STRG.18087"; transcript_id "STRG.18087.2"; cov "5.905163"; FPKM "3.279455"; TPM "9.789504";
            #CM002242   StringTie   exon    4198460 4198609 1000    +   .   gene_id "STRG.18087"; transcript_id "STRG.18087.2"; exon_number "1"; cov "6.999466";
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            start = int(start)
            end = int(end)
            ID,transcriptID,exonNum,TPM = (None,)*4
            info = attributes.split(';')
            for x in info:
                x = x.strip()
                x = x.replace('"','')
                if x.startswith('gene_id '):
                    ID = x.replace('gene_id ', '')
                elif x.startswith('transcript_id '):
                    transcriptID = x.replace('transcript_id ', '')
                elif x.startswith('exon_number '):
                    exonNum = x.replace('exon_number ', '')
                elif x.startswith('TPM '):
                	TPM = x.replace('TPM ', '')
            if feature == 'transcript':
                if not ID in Genes:
                    Genes[ID] = {'type': 'mRNA', 'codon_start': [1], 'ids': [transcriptID], 'CDS': [[]], 'mRNA': [[]], 'strand': strand, 
                                'location': (start, end), 'contig': contig, 'source': source, 'tpm': [TPM]}
                else:
                    if start < Genes[ID]['location'][0]:
                        Genes[ID]['location'] = (start,Genes[ID]['location'][1])
                    if end > Genes[ID]['location'][1]:
                        Genes[ID]['location'] = (Genes[ID]['location'][0],end)   
                    Genes[ID]['ids'].append(transcriptID)
                    Genes[ID]['mRNA'].append([])
                    Genes[ID]['CDS'].append([])
                    Genes[ID]['codon_start'].append(1)
                    Genes[ID]['tpm'].append(TPM)
            else:
                if not ID or not transcriptID:
                    print("Error, can't find geneID or transcriptID. Malformed GTF file.")
                    print(line)
                    sys.exit(1)
                if feature == 'exon':
                    if not ID in Genes:
                        Genes[ID] = {'type': 'mRNA', 'codon_start': [1], 'ids': [transcriptID], 'CDS': [[(start,end)]], 'mRNA': [[(start,end)]], 'strand': strand, 
                                'location': (start, end), 'contig': contig, 'source': source, 'tpm': []}
                    else:
                        if transcriptID in Genes[ID]['ids']: #then add exon 
                            i = Genes[ID]['ids'].index(transcriptID)
                            Genes[ID]['mRNA'][i].append((start,end))
                            Genes[ID]['CDS'][i].append((start,end))
    #loop through dictionary and make sure properly sorted exons
    for k,v in Genes.items():
        for i in range(0,len(v['ids'])):
            if v['strand'] == '+':
                sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0])
                sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0])
            else:
                sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0], reverse=True)
                sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0], reverse=True)
            Genes[k]['mRNA'][i] = sortedExons
            Genes[k]['CDS'][i] = sortedCDS
    return Genes

def dict2gff3(input):
    from collections import OrderedDict
    '''
    function to convert funannotate gene dictionary to gff3 output
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    #sort the annotations by contig and start location
    sGenes = sorted(input.iteritems(), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    #then loop through and write GFF3 format
    sys.stdout.write("##gff-version 3\n")
    for k,v in sortedGenes.items():
        sys.stdout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
        for i in range(0,len(v['ids'])):
            #build extra annotations for each transcript if applicable
            extraAnnotations = ''                  
            #now write mRNA feature
            sys.stdout.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};TPM={:}\n".format(v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['tpm'][i]))
            if v['type'] == 'mRNA':
                if '5UTR' in v:
                    #if 5'UTR then write those first
                    num_5utrs = len(v['5UTR'][i])
                    if num_5utrs > 0:
                        for z in range(0,num_5utrs):
                            u_num = z + 1
                            sys.stdout.write("{:}\t{:}\tfive_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr5p{:};Parent={:};\n".format(v['contig'], v['source'], v['5UTR'][i][z][0], v['5UTR'][i][z][1], v['strand'], v['ids'][i], u_num, v['ids'][i]))                          
                #write the exons
                num_exons = len(v['mRNA'][i])
                for x in range(0,num_exons):
                    ex_num = x + 1
                    sys.stdout.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
                #if 3'UTR then write
                if '3UTR' in v:
                    num_3utrs = len(v['3UTR'][i])
                    if num_3utrs > 0:
                        for z in range(0,num_3utrs):
                            u_num = z + 1
                            sys.stdout.write("{:}\t{:}\tthree_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr3p{:};Parent={:};\n".format(v['contig'], v['source'], v['3UTR'][i][z][0], v['3UTR'][i][z][1], v['strand'], v['ids'][i], u_num, v['ids'][i]))                         
            if v['type'] == 'mRNA':
                num_cds = len(v['CDS'][i])
                current_phase = v['codon_start'][i] - 1 #GFF3 phase is 1 less than flat file
                for y in range(0,num_cds):
                    sys.stdout.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(v['contig'], v['source'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                    current_phase = (current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                    if current_phase == 3:
                        current_phase = 0
    
Genes = gtf2dict(args.input)
dict2gff3(Genes)
