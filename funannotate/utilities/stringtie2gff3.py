#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import funannotate.library as lib


def dict2gff3(input):
    from collections import OrderedDict
    '''
    function to convert funannotate gene dictionary to gff3 output
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = sorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # then loop through and write GFF3 format
    sys.stdout.write("##gff-version 3\n")
    for k, v in list(sortedGenes.items()):
        sys.stdout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
            v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
        for i in range(0, len(v['ids'])):
            # build extra annotations for each transcript if applicable
            # now write mRNA feature
            sys.stdout.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};TPM={:}\n".format(
                v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['tpm'][i]))
            if v['type'] == 'mRNA':
                if '5UTR' in v:
                    # if 5'UTR then write those first
                    num_5utrs = len(v['5UTR'][i])
                    if num_5utrs > 0:
                        for z in range(0, num_5utrs):
                            u_num = z + 1
                            sys.stdout.write("{:}\t{:}\tfive_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr5p{:};Parent={:};\n".format(
                                v['contig'], v['source'], v['5UTR'][i][z][0], v['5UTR'][i][z][1], v['strand'], v['ids'][i], u_num, v['ids'][i]))
                # write the exons
                num_exons = len(v['mRNA'][i])
                for x in range(0, num_exons):
                    ex_num = x + 1
                    sys.stdout.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                        v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
                # if 3'UTR then write
                if '3UTR' in v:
                    num_3utrs = len(v['3UTR'][i])
                    if num_3utrs > 0:
                        for z in range(0, num_3utrs):
                            u_num = z + 1
                            sys.stdout.write("{:}\t{:}\tthree_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr3p{:};Parent={:};\n".format(
                                v['contig'], v['source'], v['3UTR'][i][z][0], v['3UTR'][i][z][1], v['strand'], v['ids'][i], u_num, v['ids'][i]))
            if v['type'] == 'mRNA':
                num_cds = len(v['CDS'][i])
                # GFF3 phase is 1 less than flat file
                current_phase = v['codon_start'][i] - 1
                for y in range(0, num_cds):
                    sys.stdout.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(
                        v['contig'], v['source'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                    current_phase = (
                        current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                    if current_phase == 3:
                        current_phase = 0


def main(args):
        # setup menu with argparse
    parser = argparse.ArgumentParser(prog='stringtie2gff.py',
                                     description='''Script to convert StringTie GTF file to GFF3.''',
                                     epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""")
    parser.add_argument('-i', '--input', required=True,
                        help='StringTie GTF file')
    args = parser.parse_args(args)

    Genes = lib.gtf2dict(args.input)
    dict2gff3(Genes)


if __name__ == "__main__":
    main(sys.argv[1:])
