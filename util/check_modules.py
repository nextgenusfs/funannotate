#!/usr/bin/env python

import sys, pkg_resources
import re

def mycmp(version1, version2):
    def normalize(v):
        return [int(x) for x in re.sub(r'(\.0+)*$','', v).split(".")]
    return cmp(normalize(version1), normalize(version2))

funannotate = ['numpy', 'pandas', 'matplotlib', 'scipy', 'scikit-learn', 'psutil', 'natsort', 'goatools', 'seaborn', 'biopython']
min_versions = {'numpy': '1.10.0', 'pandas': '0.16.1', 'matplotlib': '1.5.0', 'scipy': '0.17.0', 'scikit-learn': '0.17.0', 'psutil': '4.0.0', 'natsort': '4.0.0', 'goatools': '0.6.4', 'seaborn': '0.7.0', 'biopython': '1.65'}

vers = sys.version
print "-----------------------"
print "Python Version:"
print "-----------------------"
print vers
print "-----------------------"
print "Modules installed:"
print "-----------------------"
matches = []
for x in funannotate:
    try:
        vers = pkg_resources.get_distribution(x).version
        min = min_versions.get(x)
        vers_test = mycmp(vers, min)
        if vers_test < 0:
            print '%s =! ERROR: v%s is installed, need at least v%s' % (x, vers, min)
        else:
            print x,'=>',vers
    except pkg_resources.DistributionNotFound:
        print x,'=! ERROR. Not installed'
print "-----------------------"