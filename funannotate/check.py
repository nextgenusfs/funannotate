#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import pkg_resources
import subprocess
import errno
from natsort import natsorted
import funannotate.library as lib


def perlVersion():
    proc = subprocess.Popen(['perl', '-e', 'print $];'],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vers = proc.communicate()
    return vers[0].rstrip()


def checkPerlModule(mod):
    proc = subprocess.Popen(['perl', '-M'+mod, '-e', 'print '+mod +
                             '->VERSION . "\n"'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    vers = proc.communicate()
    # stdout in vers[0] and stderr in vers[1]
    if vers[0] != '':
        return vers[0].rstrip()
    elif vers[1].startswith("Can't locate"):
        return False
    else:
        return False
        sys.exit(1)


def checkPyModule(mod):
    try:
        vers = pkg_resources.get_distribution(mod).version
    except pkg_resources.DistributionNotFound:
        vers = False
    return vers


def mycmp(version1, version2):
    def normalize(v):
        return [int(x) for x in re.sub(r'(\.0+)*$', '', v).split(".")]
    return lib.cmp(normalize(version1), normalize(version2))


def check_version1(name):
    vers = False
    try:
        if name == 'java':
            vers = subprocess.Popen(
                [name, '-version'], stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True).communicate()
            vers = vers[1].split('\"')[1]
        else:
            vers = subprocess.Popen(
                [name, '-version'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].split('\n')[0]
            vers = vers.replace(': ', ' ')
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    except BaseException:
        print("\tERROR: %(name)s found but error running %(name)s" % (
            {'name': name}))
        return False
    return (vers)


def check_version2(name):
    vers = False
    try:
        if name == 'augustus':
            vers = subprocess.Popen([name, '--version'], stderr=subprocess.STDOUT,
                                    stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].rstrip()
        elif name == 'bamtools':
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].split('\n')[1]
        elif name == 'gmap':
            vers = subprocess.Popen([name, '--version'], stderr=subprocess.STDOUT,
                                    stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].split('\n')
            for i in vers:
                if i.startswith('GMAP'):
                    vers = i
                    break
            vers = vers.split(' called')[0].replace('version', '')
            vers = vers.replace('GMAP', '').strip()
        elif name == 'hisat2':
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].split('\n')[0]
            vers = vers.split(' ')[-1]
        elif name == 'Trinity':
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].split('\n')[0]
            vers = vers.split('Trinity-v')[-1]
        elif name == 'nucmer':
            verarray = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()
            for res in verarray:
                if re.search("version", res):
                    vers = res.split('version')[-1].strip()
                    break
                elif re.search("^\d+\.", res):
                    vers = res.strip()
                    break
        elif name == 'tbl2asn':
            vers = 'no way to determine, likely 25.X'
        elif name == 'trimal':
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()[0]
            vers = vers.strip()
        elif name == 'mafft':
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()[1]
            vers = vers.strip()
        elif 'Launch_PASA' in name:
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()[1]
            vers = vers.strip()
            vers = vers.split(': ')[-1]
        elif name == 'emapper.py':
            vers = subprocess.Popen([name, '--version'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0]
            vers.strip()
            m = re.match('emapper-(\S+)',vers)
            if m:
                vers = m.group(1)
        else:
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].split('\n')[0]
        if 'exonerate' in vers:
            vers = vers.replace('exonerate from ', '')
        if 'AUGUSTUS' in vers:
            vers = vers.split(' is ')[0].replace('(', '').replace(')', '')
            vers = vers.replace('AUGUSTUS', '').strip()
        vers = vers.replace('version ', '')
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    except BaseException:
        print("\tERROR: %(name)s found but error running %(name)s" % (
            {'name': name}))
        return False
    return (vers)


def check_version3(name):
    vers = False
    try:
        vers = subprocess.Popen(
            [name, '-v'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].split('\n')[0]
        vers = vers.replace('version open-', '')
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    except BaseException:
        print("\tERROR: %(name)s found but error running %(name)s" % (
            {'name': name}))
        return False
    return (vers)


def check_version4(name):
    vers = False
    try:
        if name == 'diamond':
            vers = subprocess.Popen(
                [name, 'version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()
            if vers[1] == '':  # then this is older version and parse the stdout
                vers = vers[0].split('version ')[-1].rstrip()
            else:
                vers = vers[1].split()[1].replace('v', '')
        else:
            vers = subprocess.Popen([name, 'version'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[
                0].split('\n')[0]
            vers = vers.replace('version ', '')
            if name == 'ete3':
                vers = vers.split(' ')[0]
            elif name == 'kallisto':
                vers = vers.split(' ')[-1]
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    except BaseException:
        print("\tERROR: %(name)s found but error running %(name)s" % (
            {'name': name}))
        return False
    return (vers)


def check_version5(name):
    vers = False
    try:
        if name == 'gmes_petap.pl':
            vers = subprocess.Popen([name], stdout=subprocess.PIPE, universal_newlines=True).communicate()[
                0].split('\n')
            for i in vers:
                if i.startswith('GeneMark-ES'):
                    vers = i
            vers = vers.replace('version ', '')
            vers = vers.split(' ')[-1]
        elif name == 'blat':
            vers = subprocess.Popen([name], stdout=subprocess.PIPE, universal_newlines=True).communicate()[
                0].split('\n')[0]
            vers = vers.split(' fast')[0]
            vers = vers.split('Standalone ')[-1].replace('v. ', 'v')
        elif name == 'pslCDnaFilter':
            vers = subprocess.Popen(
                [name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()
            vers = 'no way to determine'
        elif name == 'fasta':
            vers = subprocess.Popen(
                [name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()
            vers = 'no way to determine'
        elif name == 'CodingQuarry':
            vers = subprocess.Popen(
                [name], stdout=subprocess.PIPE, universal_newlines=True).communicate()
            v = vers[0].split('\n')
            for i in v:
                if 'CodingQuarry v.' in i:
                    vers = i.split('v. ')[-1]
        elif name == 'snap':
            vtmp = subprocess.Popen(['snap'], stderr=subprocess.STDOUT,
                                    stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].rstrip().split('\n')
            for i in vtmp:
                if i.startswith('SNAP'):
                    vers = i.split('(version ')[-1].rstrip(')')
        elif name == 'glimmerhmm':
            vtmp = subprocess.Popen(['glimmerhmm'], stderr=subprocess.STDOUT,
                                    stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].rstrip().split('\n')
            vers = '3.0.4'
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    except BaseException:
        print("\tERROR: %(name)s found but error running %(name)s" % (
            {'name': name}))
        return False
    return (vers)


def check_version6(name):
    vers = False
    try:
        if name == 'tRNAscan-SE':
            vers = subprocess.Popen([name, '-h'], stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE, universal_newlines=True).communicate()[1].split('\n')

            for i in vers:
                if i.startswith('tRNAscan-SE'):
                    vers = i
                    break
            vers = vers.split('-SE ')[-1].strip()
        else:
            vers = subprocess.Popen(
                [name, '-h'], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].split('\n')
            for i in vers:
                if i.startswith('# HMMER'):
                    vers = i
                    break
            vers = vers.split(';')[0].replace('# ', '')
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    except BaseException:
        print("\tERROR: %(name)s found but error running %(name)s" % (
            {'name': name}))
        return False
    return (vers)


def check_version7(name):
    vers = False
    try:
        vers = subprocess.Popen([name, '-V'], stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()
        vers = vers.split(' ')[-1]
        if not vers:
            vers = subprocess.Popen([name, '-version'], stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()
            vers = vers.split(' ')[2]
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    except BaseException:
        print("\tERROR: %(name)s found but error running %(name)s" % (
            {'name': name}))
        return False
    return (vers)


def main(args):
    funannotate_perl = ['Getopt::Long', 'Pod::Usage', 'File::Basename', 'threads', 'threads::shared',
                        'Thread::Queue', 'Carp', 'Data::Dumper', 'YAML', 'Hash::Merge', 'Logger::Simple', 'Parallel::ForkManager',
                        'DBI', 'Text::Soundex', 'Scalar::Util::Numeric', 'Tie::File', 'POSIX', 'Storable', 'Clone', 'Bio::Perl',
                        'DBD::mysql', 'JSON', 'LWP::UserAgent', 'DB_File', 'URI::Escape', 'File::Which', 'DBD::SQLite']

    funannotate_python = ['numpy', 'pandas', 'matplotlib', 'scipy', 'scikit-learn',
                          'psutil', 'natsort', 'goatools', 'seaborn', 'biopython', 'requests']

    programs1 = ['tblastn', 'makeblastdb', 'java', 'trimmomatic']  # -version
    programs2 = ['exonerate', 'bedtools', 'bamtools', 'augustus',
                 'samtools', 'gmap', 'hisat2', 'Trinity',
                             'tbl2asn', 'emapper.py', 'minimap2', 'mafft',
                             'trimal', 'stringtie', 'salmon', 'proteinortho', 'tantan']  # --version
    programs3 = []  # -v
    programs4 = ['diamond', 'ete3', 'kallisto']  # version
    programs5 = ['gmes_petap.pl', 'blat', 'pslCDnaFilter', 'fasta',
                 'CodingQuarry', 'snap', 'glimmerhmm']  # no version option at all, a$$holes
    programs6 = ['hmmsearch', 'hmmscan', 'tRNAscan-SE']  # -h
    programs7 = ['signalp']  # -V

    PyVers = sys.version.split(' ')[0]
    PerlVers = perlVersion()
    PyDeps = {}
    PerlDeps = {}
    ExtDeps = {}

    # loop through lists and build dictionary of results so you can print out later
    print("-------------------------------------------------------")
    print(("Checking dependencies for %s" % lib.get_version()))
    print("-------------------------------------------------------")
    global show
    show = False
    if '--show-versions' in sys.argv:
        show = True
    else:
        print("To print all dependencies and versions: funannotate check --show-versions\n")

    print(('You are running Python v %s. Now checking python packages...' % PyVers))
    for mod in funannotate_python:
        if not mod in PyDeps:
            PyDeps[mod] = checkPyModule(mod)
    missing = []
    for k, v in natsorted(list(PyDeps.items())):
        if not v:
            missing.append(k)
        elif show:
            print((k+': '+v))
    if len(missing) > 0:
        for x in missing:
            print((
                '   ERROR: %s not installed, pip install %s or conda install %s' % (x, x, x)))
    else:
        print(("All %i python packages installed" % len(funannotate_python)))
    print("\n")

    for mod in funannotate_perl:
        if not mod in PerlDeps:
            PerlDeps[mod] = checkPerlModule(mod)

    missing = []
    print(('You are running Perl v %s. Now checking perl modules...' % PerlVers))
    for k, v in natsorted(list(PerlDeps.items())):
        if not v:
            missing.append(k)
        elif show:
            print((k+': '+v))
    if len(missing) > 0:
        for x in missing:
            print(('   ERROR: %s not installed, install with cpanm %s ' % (x, x)))
    else:
        print(("All %i Perl modules installed" % len(funannotate_perl)))
    print("\n")

    # check ENV variables
    variables = ['FUNANNOTATE_DB', 'PASAHOME', 'TRINITYHOME',
                 'EVM_HOME', 'AUGUSTUS_CONFIG_PATH', 'GENEMARK_PATH']
    print('Checking Environmental Variables...')
    missing = []
    for var in variables:
        try:
            VARI = os.environ[var]
            if show:
                print(('$%s=%s' % (var, VARI)))
        except KeyError:
            if var == 'TRINITYHOME':
                try:
                    VARI = os.environ['TRINITY_HOME']
                    if show:
                        print(('$%s=%s' % ('TRINITY_HOME', VARI)))
                except KeyError:
                    missing.append(var)
            else:
                missing.append(var)
            pass
    if len(missing) > 0:
        for x in missing:
            print(('\tERROR: %s not set. export %s=/path/to/dir' % (x, x)))
    else:
        print(("All %i environmental variables are set" % (len(variables))))
    print("-------------------------------------------------------")

    if not 'PASAHOME' in missing:
        LAUNCHPASA = os.path.join(os.environ['PASAHOME'], 'Launch_PASA_pipeline.pl')
        programs2.append(LAUNCHPASA)
    print('Checking external dependencies...')
    for prog in programs1:
        if not prog in ExtDeps:
            ExtDeps[prog] = check_version1(prog)
    for prog in programs2:
        if not prog in ExtDeps:
            ExtDeps[prog] = check_version2(prog)
    for prog in programs3:
        if not prog in ExtDeps:
            ExtDeps[prog] = check_version3(prog)
    for prog in programs4:
        if not prog in ExtDeps:
            ExtDeps[prog] = check_version4(prog)
    for prog in programs5:
        if not prog in ExtDeps:
            ExtDeps[prog] = check_version5(prog)
    for prog in programs6:
        if not prog in ExtDeps:
            ExtDeps[prog] = check_version6(prog)
    for prog in programs7:
        if not prog in ExtDeps:
            ExtDeps[prog] = check_version7(prog)

    missing = []
    for k, v in natsorted(list(ExtDeps.items())):
        if not v or v.startswith('dyld:'):
            missing.append(k)
        elif show:
            if 'Launch_PASA_pipeline.pl' in k:
                k = 'PASA'
            print((k+': '+v))
    if len(missing) > 0:
        for x in missing:
            print(('\tERROR: %s not installed' % (x)))
    else:
        print(("All %i external dependencies are installed\n" % (len(ExtDeps))))


if __name__ == "__main__":
    main(sys.argv[1:])
