#!/usr/bin/env python
import pdb
import sys
import os
import xml.sax

'''make a table of interpro term counts from multiple sets of species.

The table is useful for comparing the InterPro term content of multiple
species or data sets.

First, run RunIprScan on each set of proteins, specifying a different
output directory for each set of proteins. Then, run this script like
the example shown below. The output is printed on the terminal but you
can redirect it to a text file by using the '>' symbol as shown below.

example:
python ipr_table.py species1/ species2/ ... >output.txt
'''

class IprHandler (xml.sax.ContentHandler):
    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.interpro = set()
        self.seqid = ""
        self.in_protein = False

    def startElement(self, name, attrs):
        if name == 'protein':
            self.in_protein = True
        if name == 'xref' and self.in_protein:
            self.seqid = attrs['id']
            #self.interpro = set()
        if name == 'entry':
            #if not self.interpro:
                #self.interpro = set()
            self.interpro.add(attrs['ac'] + ' ' + attrs['desc'])

    def endElement(self, name):
        if name == 'protein':
            self.in_protein = False



ipr_set = set()
species_data = {}
dir_list = sys.argv[1:]

for dirname in dir_list:
    file_list = os.listdir(dirname)
    ipr_data = {}
    for xmlfile in file_list:
        if not xmlfile.endswith('xml'):
            continue

        parser = xml.sax.make_parser()
        handler = IprHandler()
        parser.setContentHandler(handler)
        parser.parse(open(dirname+'/'+xmlfile))
        for ipr in handler.interpro:
            ipr_set.add(ipr)
        ipr_data[handler.seqid] = handler.interpro

    species_data[dirname] = ipr_data

species_list = species_data.keys()

print 'InterPro Term\t' + '\t'.join(species_list)
for ipr in ipr_set:

    counts = []
    for species in species_list:
#        pdb.set_trace()
        count = 0
        for prot in species_data[species]:
            if ipr in species_data[species][prot]:
                count += 1
        counts.append(str(count))

    print ipr + '\t' + '\t'.join(counts)
