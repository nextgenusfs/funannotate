#!/usr/bin/env python
import pdb
import sys
import os
import xml.sax
from goatools import obo_parser

if len(sys.argv) < 2:
    print("Usage: ipr2go.py go.obo.txt runiprscan_results/")
    sys.exit(1)

#load gene ontology into dictionary using goatools
dictGO = obo_parser.GODag(sys.argv[1])
    
'''Extract interpro terms from a directory of interpro xml files.

example: ipr2go.py runiprcan_results/
'''
class IprHandler (xml.sax.ContentHandler):
    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.interpro = set()
        self.go = set()
        self.seqid = ""
        self.in_protein = False

    def startElement(self, name, attrs):
        if name == 'protein':
            self.in_protein = True

        if name == 'xref' and self.in_protein:
            self.seqid = attrs['id']

        if name == 'entry':
            self.interpro.add(attrs['ac'] + ' ' + attrs['desc'])

        if name == 'go-xref':
            self.go.add(attrs['id'] + '\t' + attrs['name'])

    def endElement(self, name):
        if name == 'protein':
            self.in_protein = False

ipr_set = set()
species_data = {}
dir_list = sys.argv[2:]

for dirname in dir_list:
    file_list = os.listdir(dirname)
    ipr_data = {}
    go_data = {}
    for xmlfile in file_list:
        if not xmlfile.endswith('xml'):
            continue

        parser = xml.sax.make_parser()
        handler = IprHandler()
        parser.setContentHandler(handler)
        parser.parse(open(dirname+'/'+xmlfile))
        #pdb.set_trace()

        for go in handler.go:
            if not handler.seqid.endswith('-T1'):
                ID = handler.seqid + '-T1'
            else:
                ID = handler.seqid
            theGO = go.split('\t')[0]
            namespace = dictGO[theGO].namespace
            if namespace == 'biological_process':
                attribute = 'go_process'
            elif namespace == 'molecular_function':
                attribute = 'go_function'
            elif namespace == 'cellular_component':
                attribute = 'go_component'
            description = dictGO[theGO].name
            strippedgo = theGO.replace('GO:', '')
            print '%s\t%s\t%s|%s||IEA' % (ID, attribute, description, strippedgo)
            

