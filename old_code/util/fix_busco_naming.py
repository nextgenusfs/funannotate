#!/usr/bin/env python

#script to reformat Augustus BUSCO results
import sys, os, itertools

if len(sys.argv) < 4:
    print("Usage: fix_busco_naming.py busco_augustus.gff3 full_table_species busco_augustus.fixed.gff3")
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
            if not cols[0] in busco_complete:
                busco_complete[cols[0]] = cols[2]+':'+cols[3]+'-'+cols[4]

#now parse the augustus input file where gene numbers are likely repeated.
results = []
with open(sys.argv[1]) as f:
    for key, group in itertools.groupby(f, group_separator):
        if not key:
            results.append(list(group))

#loop through each gene model, lookup the BUSCO name, and then replace the name with counter based and busco model name
counter = 0        
inverse_busco = {v: k for k, v in busco_complete.items()}
with open(sys.argv[3], 'w') as output:
    for i in results:
        counter += 1
        cols = i[0].split('\t')
        lookup = cols[0]+':'+cols[3]+'-'+cols[4]
        if lookup in inverse_busco:
            name = inverse_busco.get(lookup)
        else:
            name = 'unknown_model'
        ID = cols[8].split(';')[0]
        ID = ID.replace('ID=', '')
        newID = 'gene'+str(counter)
        newblock = ''.join(i)
        newblock = newblock.replace('Augustus%20prediction', name)
        newblock = newblock.replace(ID, newID)
        output.write(newblock+'\n')
