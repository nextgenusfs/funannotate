#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from natsort import natsorted
import funannotate.library as lib
from Bio import SeqIO
from collections import OrderedDict


def scaffold2Dict(input):
    # get scaffold names/lengths
    scaffLen = {}
    with open(input, 'r') as seqin:
        for record in SeqIO.parse(seqin, 'fasta'):
            if not record.id in scaffLen:
                scaffLen[record.id] = len(record.seq)
    return scaffLen

def dicts2tbl(genesDict, scaff2genes, scaffLen, SeqCenter, SeqRefNum,
              annotations=False, external=False, skipList=[]):

    '''
    function to take funannotate annotation dictionaries and convert to NCBI tbl output
    '''
    duplicates = 0
    pseudo = 0
    nocds = 0
    # to parse annotations, will need to have access to GO OBO dictionary
    goDict = {}
    if annotations:
        from goatools import obo_parser
        # location of go.obo
        for item in obo_parser.OBOReader(os.path.join(os.environ["FUNANNOTATE_DB"], 'go.obo')):
            goDict[item.id] = {'name': item.name, 'namespace': item.namespace}

    def _goFormat(id, goDict=goDict):
        # go_function    serine-type endopeptidase activity|0004252||IEA
        # go_process proteolysis|0006508||IEA
        # go_component   nucleus|0005634||IEA
        if id in goDict:
            if goDict[id]['namespace'] == 'biological_process':
                base = 'go_process'
            elif goDict[id]['namespace'] == 'molecular_function':
                base = 'go_function'
            elif goDict[id]['namespace'] == 'cellular_component':
                base = 'go_component'
            reformatted = '\t\t\t{:}\t{:}|{:}||IEA'.format(
                base, goDict[id]['name'], id.replace('GO:', ''))
            return reformatted
        else:
            return False

    for k, v in natsorted(list(scaff2genes.items())):
        sys.stdout.write('>Feature %s\n' % k)
        sys.stdout.write('1\t%s\tREFERENCE\n' % scaffLen.get(k))
        sys.stdout.write('\t\t\t%s\t%s\n' % (SeqCenter, SeqRefNum))
        for genes in v:  # now loop through each gene on the scaffold
            if genes in skipList:
                continue
            # single funannotate standard dictionary
            geneInfo = genesDict.get(genes)
            if 'pseudo' in geneInfo:
                if geneInfo['pseudo']:
                    try:
                        log.debug('{:} is pseudo, skipping'.format(genes))
                    except NameError:
                        print(('{:} is pseudo, skipping'.format(genes)))
                    pseudo += 1
                    continue
            if geneInfo['type'] == 'mRNA' and not geneInfo['CDS']:
                try:
                    log.debug(
                        'Skipping {:} because no CDS found.'.format(genes))
                except NameError:
                    print((
                        'Skipping {:} because no CDS found.'.format(genes)))
                pseudo += 1
                continue
            if geneInfo['type'] == 'mRNA' and not len(geneInfo['ids']) == len(geneInfo['mRNA']) == len(geneInfo['CDS']):
                try:
                    log.debug('Incompatible annotation found: {:}\n{:}'.format(
                        genes, geneInfo))
                except NameError:
                    print(('Incompatible annotation found: {:}\n{:}'.format(
                        genes, geneInfo)))
                duplicates += 1
                continue
            if geneInfo['type'] == 'mRNA' and len(geneInfo['CDS']) == 0:
                nocds += 1
                continue
            if geneInfo['type'] is None:
                continue
            # check for partial models
            if True in geneInfo['partialStart']:
                ps = '<'
            else:
                ps = ''
            if True in geneInfo['partialStop']:
                pss = '>'
            else:
                pss = ''
            # now write gene model
            if geneInfo['strand'] == '+':
                sys.stdout.write('%s%i\t%s%i\tgene\n' % (
                    ps, geneInfo['location'][0], pss, geneInfo['location'][1]))
                if annotations:
                    if geneInfo['name']:
                        sys.stdout.write('\t\t\tgene\t%s\n' % geneInfo['name'])
                    if geneInfo['gene_synonym']:
                        for alias in geneInfo['gene_synonym']:
                            sys.stdout.write('\t\t\tgene_synonym\t%s\n' % alias)
                sys.stdout.write('\t\t\tlocus_tag\t%s\n' % genes)
            else:
                sys.stdout.write('%s%i\t%s%i\tgene\n' % (
                    ps, geneInfo['location'][1], pss, geneInfo['location'][0]))
                if annotations:
                    if geneInfo['name']:
                        sys.stdout.write('\t\t\tgene\t%s\n' % geneInfo['name'])
                    if geneInfo['gene_synonym']:
                        for alias in geneInfo['gene_synonym']:
                            sys.stdout.write('\t\t\tgene_synonym\t%s\n' % alias)
                sys.stdout.write('\t\t\tlocus_tag\t%s\n' % genes)

            # now will output the gene models with -T1, -T2, -T3 annotations based on expression values
            # means need to get the order
            order = []
            # multiple transcripts, so get order of highest TPM
            if len(geneInfo['ids']) > 1:
                tpms = []
                for num, tpm in enumerate(geneInfo['note']):
                    for item in tpm:
                        if item.startswith('TPM:'):
                            value = float(item.split(':')[-1])
                            tpms.append((value, num))
                if len(tpms) > 0:
                    for x in sorted(tpms, reverse=True):
                        order.append(x[1])
                else:
                    order = list(range(0, len(geneInfo['ids'])))
            else:
                order.append(0)
            for num, i in enumerate(order):  # now write mRNA and CDS features
                # if geneInfo['ids'][i].startswith('evm.model'): #if from predict, rename to match locus_tag
                #    protein_id = genes+'-T'+str(num+1)
                # else:
                #    protein_id = geneInfo['ids'][i]
                if external:
                    protein_id = geneInfo['ids'][i]
                else:
                    protein_id = genes+'-T'+str(num+1)
                if geneInfo['type'] == 'mRNA':
                    if geneInfo['partialStart'][i] is False:
                        ps = ''
                    else:
                        ps = '<'
                    if geneInfo['partialStop'][i] is False:
                        pss = ''
                    else:
                        pss = '>'
                    if geneInfo['strand'] == '+':
                        for num, exon in enumerate(geneInfo['mRNA'][i]):
                            # single exon, so slightly differnt method
                            if num == 0 and num == len(geneInfo['mRNA'][i]) - 1:
                                sys.stdout.write('%s%s\t%s%s\tmRNA\n' %
                                          (ps, exon[0], pss, exon[1]))
                            elif num == 0:
                                sys.stdout.write('%s%s\t%s\tmRNA\n' %
                                          (ps, exon[0], exon[1]))
                            # this is last one
                            elif num == len(geneInfo['mRNA'][i]) - 1:
                                sys.stdout.write('%s\t%s%s\n' %
                                          (exon[0], pss, exon[1]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (exon[0], exon[1]))
                        sys.stdout.write('\t\t\tproduct\t%s\n' %
                                  geneInfo['product'][i])
                        sys.stdout.write('\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n' % (protein_id))
                        sys.stdout.write('\t\t\tprotein_id\tgnl|ncbi|%s\n' %
                                  (protein_id))
                        for num, cds in enumerate(geneInfo['CDS'][i]):
                            # single exon, so slightly differnt method
                            if num == 0 and num == len(geneInfo['CDS'][i]) - 1:
                                sys.stdout.write('%s%s\t%s%s\tCDS\n' %
                                          (ps, cds[0], pss, cds[1]))
                            elif num == 0:
                                sys.stdout.write('%s%s\t%s\tCDS\n' %
                                          (ps, cds[0], cds[1]))
                            # this is last one
                            elif num == len(geneInfo['CDS'][i]) - 1:
                                sys.stdout.write('%s\t%s%s\n' %
                                          (cds[0], pss, cds[1]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (cds[0], cds[1]))
                        sys.stdout.write('\t\t\tcodon_start\t%i\n' %
                                  geneInfo['codon_start'][i])
                        if annotations:  # write functional annotation
                            if geneInfo['EC_number'][i]:
                                for EC in geneInfo['EC_number'][i]:
                                    sys.stdout.write('\t\t\tEC_number\t%s\n' % EC)
                            if geneInfo['db_xref'][i]:
                                for xref in geneInfo['db_xref'][i]:
                                    sys.stdout.write('\t\t\tdb_xref\t%s\n' % xref)
                            if geneInfo['go_terms'][i]:
                                for go in geneInfo['go_terms'][i]:
                                    goLine = _goFormat(go)
                                    if goLine:
                                        sys.stdout.write('{:}\n'.format(goLine))
                            if geneInfo['note'][i]:
                                for item in geneInfo['note'][i]:
                                    sys.stdout.write('\t\t\tnote\t%s\n' % item)
                        sys.stdout.write('\t\t\tproduct\t%s\n' %
                                  geneInfo['product'][i])
                        sys.stdout.write('\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n' % (protein_id))
                        sys.stdout.write('\t\t\tprotein_id\tgnl|ncbi|%s\n' %
                                  (protein_id))
                    else:  # means this is on crick strand
                        for num, exon in enumerate(geneInfo['mRNA'][i]):
                            # single exon, so slightly differnt method
                            if num == 0 and num == len(geneInfo['mRNA'][i]) - 1:
                                sys.stdout.write('%s%s\t%s%s\tmRNA\n' %
                                          (ps, exon[1], pss, exon[0]))
                            elif num == 0:
                                sys.stdout.write('%s%s\t%s\tmRNA\n' %
                                          (ps, exon[1], exon[0]))
                            # this is last one
                            elif num == len(geneInfo['mRNA'][i]) - 1:
                                sys.stdout.write('%s\t%s%s\n' %
                                          (exon[1], pss, exon[0]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (exon[1], exon[0]))
                        sys.stdout.write('\t\t\tproduct\t%s\n' %
                                  geneInfo['product'][i])
                        sys.stdout.write('\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n' % (protein_id))
                        sys.stdout.write('\t\t\tprotein_id\tgnl|ncbi|%s\n' %
                                  (protein_id))
                        for num, cds in enumerate(geneInfo['CDS'][i]):
                            # single exon, so slightly differnt method
                            if num == 0 and num == len(geneInfo['CDS'][i]) - 1:
                                sys.stdout.write('%s%s\t%s%s\tCDS\n' %
                                          (ps, cds[1], pss, cds[0]))
                            elif num == 0:
                                sys.stdout.write('%s%s\t%s\tCDS\n' %
                                          (ps, cds[1], cds[0]))
                            # this is last one
                            elif num == (len(geneInfo['CDS'][i]) - 1):
                                sys.stdout.write('%s\t%s%s\n' %
                                          (cds[1], pss, cds[0]))
                            else:
                                sys.stdout.write('%s\t%s\n' % (cds[1], cds[0]))
                        sys.stdout.write('\t\t\tcodon_start\t%i\n' %
                                  geneInfo['codon_start'][i])
                        if annotations:  # write functional annotation
                            if geneInfo['EC_number'][i]:
                                for EC in geneInfo['EC_number'][i]:
                                    sys.stdout.write('\t\t\tEC_number\t%s\n' % EC)
                            if geneInfo['db_xref'][i]:
                                for xref in geneInfo['db_xref'][i]:
                                    sys.stdout.write('\t\t\tdb_xref\t%s\n' % xref)
                            if geneInfo['go_terms'][i]:
                                for go in geneInfo['go_terms'][i]:
                                    goLine = _goFormat(go)
                                    if goLine:
                                        sys.stdout.write('{:}\n'.format(goLine))
                            if geneInfo['note'][i]:
                                for item in geneInfo['note'][i]:
                                    sys.stdout.write('\t\t\tnote\t%s\n' % item)
                        sys.stdout.write('\t\t\tproduct\t%s\n' %
                                  geneInfo['product'][i])
                        sys.stdout.write('\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n' % (protein_id))
                        sys.stdout.write('\t\t\tprotein_id\tgnl|ncbi|%s\n' %
                                  (protein_id))
                elif geneInfo['type'] == 'tRNA':
                    if geneInfo['strand'] == '+':
                        for num, exon in enumerate(geneInfo['mRNA'][i]):
                            if num == 0:
                                sys.stdout.write('%s\t%s\t%s\n' % (
                                    exon[0], exon[1], geneInfo['type']))
                            else:
                                sys.stdout.write('%s\t%s\n' % (exon[0], exon[1]))
                        sys.stdout.write('\t\t\tproduct\t%s\n' %
                                  geneInfo['product'][i])
                        if geneInfo['product'] == 'tRNA-Xxx':
                            sys.stdout.write('\t\t\tpseudo\n')
                    else:
                        for num, exon in enumerate(geneInfo['mRNA'][i]):
                            if num == 0:
                                sys.stdout.write('%s\t%s\t%s\n' % (
                                    exon[1], exon[0], geneInfo['type']))
                            else:
                                sys.stdout.write('%s\t%s\n' % (exon[1], exon[0]))
                        sys.stdout.write('\t\t\tproduct\t%s\n' %
                                  geneInfo['product'][i])
                        if geneInfo['product'] == 'tRNA-Xxx':
                            sys.stdout.write('\t\t\tpseudo\n')
                elif geneInfo['type'] in ['rRNA', 'ncRNA']:
                    if geneInfo['strand'] == '+':
                        sys.stdout.write('%s\t%s\t%s\n' % (
                            geneInfo['location'][0], geneInfo['location'][1], geneInfo['type']))
                        sys.stdout.write('\t\t\tproduct\t%s\n' %
                                  geneInfo['product'][i])
                    else:
                        sys.stdout.write('%s\t%s\t%s\n' % (
                            geneInfo['location'][1], geneInfo['location'][0], geneInfo['type']))
                        sys.stdout.write('\t\t\tproduct\t%s\n' %
                                  geneInfo['product'][i])
    if any(i > 0 for i in [duplicates, pseudo, nocds]):
        try:
            print(('Skipped {:,} annotations: {:,} pseudo genes; {:,} no CDS; {:,} duplicated features'.format(
                sum([pseudo, nocds, duplicates]), pseudo, nocds, duplicates)))
        except NameError:
            print(('Skipped {:,} annotations: {:,} pseudo genes; {:,} no CDS; {:,} duplicated features'.format(
                sum([pseudo, nocds, duplicates]), pseudo, nocds, duplicates)))


def main(args):
        # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='gff2prot.py',
                                     description='''Script to convert GFF3 and FASTA to tbl, proteins, transcripts.''',
                                     epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-g', '--gff3', required=True,
                        help='Genome annotation GFF3 format')
    parser.add_argument('-f', '--fasta', required=True,
                        help='Genome in FASTA format')
    args = parser.parse_args(args)

    # load into dictionary
    Genes = {}
    Genes = lib.gff2dict(args.gff3, args.fasta, Genes)

    # sort the dictionary
    def _sortDict(d):
        return (d[1]['location'][0], d[1]['location'][1])

    # now sort dictionary by contig and location, rename using prefix, translate to protein space to get proper start/stop info
    sGenes = sorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    scaff2genes = {}
    for k, v in list(sortedGenes.items()):
        if not v['contig'] in scaff2genes:
            scaff2genes[v['contig']] = [k]
        else:
            scaff2genes[v['contig']].append(k)

    # get length of scaffolds
    scaffLen = scaffold2Dict(args.fasta)

    # now write table
    dicts2tbl(sortedGenes, scaff2genes, scaffLen, 'CFMR', '12345',
              annotations=True)


if __name__ == "__main__":
    main(sys.argv[1:])
