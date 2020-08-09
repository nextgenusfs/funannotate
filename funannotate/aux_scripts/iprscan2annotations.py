#!/usr/bin/env python

# script written for funannotate by Jon Palmer (2017)
# it will parse an interproscan5 xml file and generate genome annotation file for
# GO terms and IPR terms suitable for GAG

import sys
import xml.etree.cElementTree as etree


if len(sys.argv) < 2:
    print("Usage: iprscan2annotations.py IPRSCAN.xml OUTPUT.annotations.txt")
    sys.exit(1)


def convertGOattribute(namespace):
    if namespace == 'BIOLOGICAL_PROCESS':
        attribute = 'go_process'
    elif namespace == 'MOLECULAR_FUNCTION':
        attribute = 'go_function'
    elif namespace == 'CELLULAR_COMPONENT':
        attribute = 'go_component'
    else:
        print('Error parsing XML GO terms: %s is not a valid term' % namespace)
        sys.exit(1)
    return attribute


with open(sys.argv[2], 'w') as output:
    with open(sys.argv[1]) as xml_file:
        tree = etree.iterparse(xml_file)
        for _, elem in tree:
            if '}' in elem.tag:
                elem.tag = elem.tag.split('}', 1)[1]
            for at in list(elem.attrib.keys()):
                if '}' in at:
                    newat = at.split('}', 1)[1]
                    elem.attrib[newat] = elem.attrib[at]
                    del elem.attrib[at]
        root = tree.root
        # iterate through each of the protein hits
        for hits in root:
            IDs = []
            iprs = []
            gos = []
            signalp = []
            for lv1 in hits:
                if lv1.tag == 'xref':
                    name = lv1.get('id')
                    IDs.append(name)
                if lv1.tag == 'matches':
                    for e in lv1.findall('.//entry'):
                        if not e.get('ac') in iprs:
                            iprs.append(e.get('ac'))
                    for g in lv1.findall('.//go-xref'):
                        cat, goID, desc = g.get(
                            'category'), g.get('id'), g.get('name')
                        cat = convertGOattribute(cat)
                        goHit = (cat, desc, goID)
                        if not goHit in gos:
                            gos.append(goHit)
                    for s in lv1.findall('.//signalp-match'):
                        for lib in s.findall('.//signature-library-release'):
                            if lib.get('library') == "SIGNALP_EUK":
                                for loc in s.findall('.//signalp-location'):
                                    signalp.append(
                                        (loc.get('start'), loc.get('end')))
            # print out annotation file if IPR domains
            if len(iprs) > 0:
                for i in IDs:
                    for x in iprs:
                        output.write('%s\tdb_xref\tInterPro:%s\n' % (i, x))
            if len(gos) > 0:
                for i in IDs:
                    for x in gos:
                        output.write('%s\t%s\t%s|%s||IEA\n' %
                                     (i, x[0], x[1], x[2].replace('GO:', '')))
            # if len(signalp) > 0:
            #    for i in IDs:
            #        for x in signalp:
            #            output.write('%s\tnote\tSECRETED:SignalP(%s-%s)\n' % (i, x[0], x[1]))
