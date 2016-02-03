#!/usr/bin/env python

import sys
from natsort import natsorted

with open(sys.argv[1], 'rU') as input:
    Genes = {}
    for line in input:
        if not line.startswith(' '):
            continue
        if line.startswith(' Seq name:'):
            scaffold = line.replace(' Seq name: ', '')
            scaffold = scaffold.strip()
        line.strip()
        cols = line.split(' ')
        clean = []
        for i in cols:
            if not i == '':
                clean.append(i)
        if len(clean) > 5:
            if clean[3].startswith('CDS'):
                GeneID = 'gene_' + clean[0]
                strand = clean[1]
                start = clean[4]
                end = clean[6]
                result = (start,end,strand)
                if not GeneID in Genes:
                    Genes[GeneID] = [result]
                else:
                    Genes[GeneID].append(result)
    for k,v in natsorted(Genes.items()):
        genestart = v[0][0]
        strand = v[0][2]
        genestop = v[-1][1]
        sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffold, 'FGENESH', 'gene', genestart, genestop, '.', strand, '.', 'ID='+k+';'))
        sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffold, 'FGENESH', 'mRNA', genestart, genestop, '.', strand, '.', 'ID='+k+'-T1;Parent='+k+';'))
        for i in range(0,len(v)):
            start = v[i][0]
            stop = v[i][1]
            if i == 0: #first CDS, set phase to 0, calculate next phase
                phase = '0'
                diff = int(stop) - int(start) + 1
                next_phase = (int(phase) - diff) % 3
            else:
                phase = next_phase
                diff = int(stop) - int(start) + 1
                next_phase = (int(phase) - diff) % 3

            sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffold, 'FGENESH', 'exon', start, stop, '.', strand, '.', 'ID='+k+':exon'+str(i+1)+';Parent='+k+'-T1;'))
            sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffold, 'FGENESH', 'CDS', start, stop, '.', strand, phase, 'ID=cds.'+k+';Parent='+k+'-T1;'))
