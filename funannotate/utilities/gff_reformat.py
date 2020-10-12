#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from natsort import natsorted
from collections import OrderedDict
import funannotate.library as lib


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='gff_reformat.py',
                                     description='''Script to rename gene models GFF3 file.''',
                                     epilog="""Written by Jon Palmer (2020) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-g', '--gff3', required=True,
                        help='Genome annotation GFF3 format')
    parser.add_argument('-f', '--fasta', required=True,
                        help='Genome in FASTA format')
    parser.add_argument('-l', '--locus_tag', default='FUN',
                        help='Basename of gene names')
    parser.add_argument('-n', '--numbering', default=1, type=int,
                        help='Start numbering at')
    parser.add_argument('-o', '--out', required=True, help='Output GFF3')
    args = parser.parse_args(args)

    # load into dictionary
    Genes = {}
    Genes = lib.gff2dict(args.gff3, args.fasta, Genes)
    print('Parsed {:,} gene models from {}'.format(len(Genes), args.gff3))

    # now create ordered dictionary and sort by contig and position
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])

    sGenes = natsorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    counter = args.numbering
    args.locus_tag = args.locus_tag.rstrip('_')
    transcripts = 0
    for k, v in list(sortedGenes.items()):
        locusTag = args.locus_tag+'_'+str(counter).zfill(6)
        renamedGenes[locusTag] = v
        renamedGenes[locusTag]['gene_synonym'].append(k)
        newIds = []
        for i in range(0, len(v['ids'])):
            newIds.append('{}-T{}'.format(locusTag, i+1))
            transcripts += 1
        renamedGenes[locusTag]['ids'] = newIds
        counter += 1

    # write to gff3
    lib.dict2gff3(renamedGenes, args.out)
    print('Sorted and renamed {:,} gene models {:,} transcripts: {}'.format(
        len(renamedGenes), transcripts, args.out))

if __name__ == "__main__":
    main(sys.argv[1:])
