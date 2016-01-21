#!/usr/bin/env python
import pdb
import sys
import os
import xml.sax


'''extract sequences for the given interpro term

usage:
   python extract_ipr.py <dir> <interpro term>

'''
class IprHandler (xml.sax.ContentHandler):
    def __init__(self):
        xml.sax.ContentHandler.__init__(self)

        self.interpro = set()
        self.seqid = ""
        self.seq = ""

        self.in_protein = False
        self.in_sequence = False
        self._buffer = ""

    def startElement(self, name, attrs):
        #print 'starting: ' + name
        if name == 'protein':
            self.in_protein = True
        if name == 'sequence':
            self.in_sequence = True
        if name == 'xref' and self.in_protein:
            self.seqid = attrs['id']
            #self.interpro = set()
        if name == 'entry':
            #if not self.interpro:
                #self.interpro = set()
            self.interpro.add(attrs['ac'])

    def characters(self, content):
            if self.in_sequence:
                self._buffer += content

    def endElement(self, name):
        #print 'ending ' + name
        if name == 'protein':
            self.in_protein = False
        if name == 'sequence':
            self.seq = self._buffer
            #pdb.set_trace()
            self._buffer = ''



file_list = os.listdir(sys.argv[1])

for xmlfile in file_list:
    if not xmlfile.endswith('xml'):
        continue

    parser = xml.sax.make_parser()
    handler = IprHandler()
    parser.setContentHandler(handler)
    parser.parse(open(sys.argv[1]+'/'+xmlfile))

    if sys.argv[2] in handler.interpro:
        #seq = handler.seq
        print '>' + handler.seqid
        print handler.seq
        #pdb.set_trace()

