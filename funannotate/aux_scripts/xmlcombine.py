#!/usr/bin/env python
import sys
import os
import os.path
import fnmatch
from xml.etree import cElementTree

cElementTree.register_namespace(
    '', "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5")


def run(xml_files):
    first = None
    for filename in xml_files:
        data = cElementTree.parse(filename).getroot()
        if first is None:
            first = data
        else:
            first.extend(data)
    if first is not None:
        print(cElementTree.tostring(first))


if __name__ == "__main__":
    xml_files = [os.path.join(dirpath, f)
                 for dirpath, dirnames, files in os.walk(sys.argv[1])
                 for f in fnmatch.filter(files, '*.xml')]
    run(xml_files)
