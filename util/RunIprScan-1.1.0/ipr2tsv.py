#!/usr/bin/env python
import pdb
import sys
import os
import xml.sax

if len(sys.argv) < 2:
    print("Usage: ipr2tsv.py runiprscan_results/")
    sys.exit(1)

'''Extract interpro terms from a directory of interpro xml files.

example: ipr2tsv.py runiprcan_results/
'''

#script modified by Jon Palmer for funannotate

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
        #pdb.set_trace()

        for ipr in handler.interpro:
            if not handler.seqid.endswith('-T1'):
                ID = handler.seqid + '-T1'
            else:
                ID = handler.seqid
            interproID = ipr.split(" ")[0]
            interproDesc = ipr.split(" ")[-1]
            print '%s\tdb_xref\tInterPro:%s' % (ID, interproID)
            #print '%s\tnote\tInterPro Description: %s' % (ID, interproDesc)

