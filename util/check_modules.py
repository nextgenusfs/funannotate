#!/usr/bin/env python

import sys, os, re, pkg_resources, subprocess

def mycmp(version1, version2):
    def normalize(v):
        return [int(x) for x in re.sub(r'(\.0+)*$','', v).split(".")]
    return cmp(normalize(version1), normalize(version2))
    
def check_version1(name):
    try:
        vers = subprocess.Popen([name, '-version'], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
        vers = vers.replace(': ', ' v')
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return (name+': not installed')
    return (vers)
    
def check_version2(name):
    try:
        if name == 'augustus':
            vers = subprocess.Popen([name, '--version'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0].rstrip()
        elif name == 'bamtools':
            vers = subprocess.Popen([name, '--version'], stdout=subprocess.PIPE).communicate()[0].split('\n')[1]
            vers = vers.replace('bamtools ', 'bamtools v')
        elif name == 'gmap':
            vers = subprocess.Popen([name, '--version'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0].split('\n')
            for i in vers:
                if i.startswith('GMAP'):
                    vers = i
                    break
            vers = vers.split(' called')[0].replace('version', 'v')
        else:
            vers = subprocess.Popen([name, '--version'], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
        if 'exonerate' in vers:
            vers = vers.replace('exonerate from ', '')
        if 'AUGUSTUS' in vers:
            vers = vers.split(' is ')[0].replace('(','v').replace(')', '')
        vers = vers.replace('version ', 'v')
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return (name+': not installed')
    return (vers)
    
def check_version3(name):
    try:
        vers = subprocess.Popen([name, '-v'], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
        vers = vers.replace('version open-', 'v')
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return (name+': not installed')
    return (vers)

def check_version4(name):
    try:
        vers = subprocess.Popen([name, 'version'], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
        vers = vers.replace('version ', 'v')
        if name == 'ete3':
            vers = 'ete3 v'+vers
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return (name+': not installed')
    return (vers)

def check_version5(name):
    try:
        if name == 'gmes_petap.pl':
            vers = subprocess.Popen([name], stdout=subprocess.PIPE).communicate()[0].split('\n')
            for i in vers:
                if i.startswith('GeneMark-ES'):
                    vers = i
            vers = vers.replace('version ', 'v')
        elif name == 'blat':
            vers = subprocess.Popen([name], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
            vers = vers.split(' fast')[0]
            vers = vers.split('Standalone ')[-1].replace('v. ', 'v')
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return (name+': not installed')
    return (vers)

def check_version6(name):
    try:
        vers = subprocess.Popen([name, '-h'], stdout=subprocess.PIPE).communicate()[0].split('\n')
        for i in vers:
            if i.startswith('# HMMER'):
                vers = i
                break
        vers = vers.split(';')[0].replace('# ', '')
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return (name+': not installed')
    return (vers)

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
print "External Dependencies:"
print "-----------------------"

programs1 = ['tblastn', 'blastp', 'blastn', 'makeblastdb', 'rmblastn'] #-version
programs2 = ['exonerate', 'bedtools', 'bamtools', 'augustus', 'braker.pl', 'gmap'] #--version
programs3 = ['RepeatModeler', 'RepeatMasker'] #-v
programs4 = ['diamond', 'ete3'] #version
programs5 = ['gmes_petap.pl', 'blat'] #no version option at all, a$$holes
for i in programs1:
    print check_version1(i)
for i in programs2:
    print check_version2(i)
for i in programs3:
    print check_version3(i)
for i in programs4:
    print check_version4(i)
for i in programs5:
    print check_version5(i)
print check_version6('hmmsearch')

