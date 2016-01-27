#!/usr/bin/env python

from Bio import SeqIO
import sys

history = []
with open(sys.argv[3], 'w') as gff:
    with open(sys.argv[2], 'w') as fasta:
        with open(sys.argv[1]) as input:
            SeqRecords = SeqIO.parse(input, 'genbank')
            for record in SeqRecords:
                for f in record.features:
                    if f.type == 'CDS':
                        protID = f.qualifiers['protein_id'][0]
                        locusID = f.qualifiers['locus_tag'][0]
                        start = f.location.nofuzzy_start
                        end = f.location.nofuzzy_end
                        strand = f.location.strand
                        if strand == 1:
                            strand = '+'
                        elif strand == -1:
                            strand = '-'
                        translation = f.qualifiers['translation'][0]
                        product = f.qualifiers['product'][0]
                        chr = record.id
                        if '.' in chr:
                            chr = chr.split('.')[0]
                        if not protID in history:
                            history.append(protID)
                            gff.write("%s\tNCBI\tCDS\t%s\t%s\t.\t%s\t.\tID=%s;Alias=%s;Product=%s;\n" % (chr, start, end, strand, protID, locusID, product))
                            fasta.write(">%s\n%s\n" % (protID, translation))

