#!/usr/bin/env python

#script to tease apart maker gff into EVM input

import sys

tr = 'transcript_alignments.gff3'
pr = 'protein_alignments.gff3'
gr = 'gene_predictions.gff3'

with open(tr, 'w') as trout:
    with open(pr, 'w') as prout:
        with open(gr, 'w') as grout:
            with open(sys.argv[1], 'rU') as input:
                for line in input:
                    if line.startswith('#'):
                        continue
                    if 'trnascan' in line:
                        continue
                    cols = line.split('\t')
                    if 'maker' in cols[1]:
                        grout.write(line)
                    elif 'protein2genome' in cols[1]:
                        if 'match_part' in cols[2]:
                            cols[2] = 'nucleotide_to_protein_match'
                            cols[5] = '.'
                            prout.write('\t'.join(cols))
                    elif 'est2genome' in cols[1]:
                        if 'match_part' in cols[2]:
                            cols[2] = 'EST_match'
                            cols[5] = '.'
                            trout.write('\t'.join(cols))
                    elif 'cdna2genome' in cols[1]:
                        if 'match_part' in cols[2]:
                            cols[2] = 'EST_match'
                            cols[5] = '.'
                            trout.write('\t'.join(cols))
                    elif 'pred_gff' in cols[1]:
                        if 'match_part' in cols[2]:
                            cols[1] = cols[1].replace('pred_gff:', '')
                            cols[2] = 'EST_match'
                            cols[5] = '100.0'
                            trout.write('\t'.join(cols))