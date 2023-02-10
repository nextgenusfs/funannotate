#!/usr/bin/env python

# script written for funannotate by Jon Palmer (2017)
# it will parse an interproscan5 xml file and generate
# genome annotation file for GO terms and IPR terms

import sys
import os
import xml.etree.cElementTree as etree
from goatools import obo_parser


def convertGOattribute(namespacein):
    namespace = namespacein.upper()
    if namespace == 'BIOLOGICAL_PROCESS':
        attribute = 'go_process'
    elif namespace == 'MOLECULAR_FUNCTION':
        attribute = 'go_function'
    elif namespace == 'CELLULAR_COMPONENT':
        attribute = 'go_component'
    else:
        # print(f'Error parsing XML GO terms: {namespace} is not a valid term')
        attribute = "go_unknown"
        # sys.exit(1)
    return attribute


def main():
    '''Main step of intepro annotations to tab delimited script.'''

    if len(sys.argv) < 2:
        print("Usage: iprscan2annotations.py IPRSCAN.xml OUTPUT.annotations.txt")
        sys.exit(1)

    goDict = {}
    for item in obo_parser.OBOReader(os.path.join(os.environ["FUNANNOTATE_DB"],
                                                  'go.obo')):
        namespace = convertGOattribute(item.namespace)
        goDict[item.id] = {'name': item.name,
                           'namespace': namespace}
        for nm in item.alt_ids:  # also index by alt_id since that may be reported
            goDict[nm] = {'name': item.name,
                          'namespace': namespace}
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
                gos = {}
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
                            if cat is None or desc is None:
                                if goID in goDict:
                                    cat = goDict[goID]['namespace']
                                    desc = goDict[goID]['name']
                                else:
                                    continue
                                    # cat = ""
                                    # desc = ""
                                    # print(f"No GO term {goID} in obo DB")
                            else:
                                cat = convertGOattribute(cat)
                            goHit = (cat, desc, goID)
                            if goID not in gos:
                                gos[goID] = goHit
                        # signalp is processed elsewhere
                        # do we just skip this parsing even?
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
                            output.write(f'{i}\tdb_xref\tInterPro:{x}\n')
                if len(gos) > 0:
                    for i in IDs:
                        for goid in gos:
                            x = gos[goid]
                            GOID = x[2].replace('GO:', '')
                            output.write(f'{i}\t{x[0]}\t{x[1]}|{GOID}||IEA\n')


if __name__ == "__main__":
    main()
