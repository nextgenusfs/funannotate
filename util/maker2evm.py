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
                    if 'ID=trnascan' in line:
                        print line
                    cols = line.split('\t')
                    if 'maker' in cols[1]:
                        grout.write(line)
                    elif 'protein2genome' in cols[1]:
                        if 'match_part' in cols[2]:
                            line = line.replace('match_part', 'nucleotide_to_protein_match')
                            prout.write(line)
                    elif 'est2genome' in cols[1]:
                        if 'match_part' in cols[2]:
                            line = line.replace('match_part', 'EST_match')
                            trout.write(line)
                    elif 'cdna2genome' in cols[1]:
                        if 'match_part' in cols[2]:
                            line = line.replace('match_part', 'EST_match')
                            trout.write(line)