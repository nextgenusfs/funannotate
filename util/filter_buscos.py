#!/usr/bin/env python

#script to reformat Augustus BUSCO results
import sys, os, itertools

if len(sys.argv) < 4:
    print("Usage: filter_buscos.py busco.evm.gff3 full_table_species busco.final.gff3")
    sys.exit(1)

def group_separator(line):
    return line=='\n'

#parse the busco table into dictionary format
busco_complete = {}
with open(sys.argv[2], 'rU') as buscoinput:
    for line in buscoinput:
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        if cols[1] == 'Complete':
            ID = cols[2].replace('evm.model.', '')
            if not ID in busco_complete:
                busco_complete[ID] = (cols[0], cols[3])
            else:
                score = busco_complete.get(ID)[1]
                if float(cols[3]) > float(score):
                    busco_complete[ID] = (cols[0], cols[3])
                    print ID, 'updated dictionary'
                else:
                    print ID, 'is repeated and score is less'

#now parse the evm busco file, group them
results = []
with open(sys.argv[1]) as f:
    for key, group in itertools.groupby(f, group_separator):
        if not key:
            results.append(list(group))

#loop through each gene model, lookup the BUSCO name, and then replace the name with counter based and busco model name
'''
scaffold_1	EVM	gene	18407	18947	.	-	.	ID=evm.TU.scaffold_1.1;Name=EVM%20prediction%20scaffold_1.1
scaffold_1	EVM	mRNA	18407	18947	.	-	.	ID=evm.model.scaffold_1.1;Parent=evm.TU.scaffold_1.1;Name=EVM%20prediction%20scaffold_1.1
scaffold_1	EVM	exon	18772	18947	.	-	.	ID=evm.model.scaffold_1.1.exon1;Parent=evm.model.scaffold_1.1
scaffold_1	EVM	CDS	18772	18947	.	-	0	ID=cds.evm.model.scaffold_1.1;Parent=evm.model.scaffold_1.1
scaffold_1	EVM	exon	18407	18615	.	-	.	ID=evm.model.scaffold_1.1.exon2;Parent=evm.model.scaffold_1.1
scaffold_1	EVM	CDS	18407	18615	.	-	1	ID=cds.evm.model.scaffold_1.1;Parent=evm.model.scaffold_1.1
'''
counter = 0        
with open(sys.argv[3], 'w') as output:
    for i in results:
        counter += 1
        cols = i[0].split('\t')
        ID = cols[8].split(';')[0]
        ID = ID.replace('ID=', '')
        lookup = ID.replace('evm.TU.', '')
        if lookup in busco_complete:
            name = busco_complete.get(lookup)[0]
            geneID = 'gene'+str(counter)
            mrnaID = 'mrna'+str(counter)
            newblock = ''.join(i)
            newblock = newblock.replace('EVM%20prediction%20'+lookup, name)
            newblock = newblock.replace(ID, geneID)
            newblock = newblock.replace('evm.model.'+lookup, mrnaID)
            output.write(newblock+'\n')
