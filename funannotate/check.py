#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import importlib.metadata
import subprocess
import errno
import shutil
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
        version = importlib.metadata.version(mod)
    except importlib.metadata.PackageNotFoundError:
        # Handle the case where the package metadata is not found
        version = False
    return version


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
            verarray = subprocess.Popen([name, '--version'],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True).communicate()
            for res in verarray:
                if re.search("version", res):
                    vers = res.split('version')[-1].strip()
                    break
                elif re.search(r"^\d+\.", res):
                    vers = res.strip()
                    break
        elif name == 'tbl2asn' or name == 'table2asn':
            (so, se) = subprocess.Popen([name, '-'], stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True).communicate()
            m = re.search(name + r'\s+(\S+)\s+', so+se)
            if m:
                vers = m.group(1)
            else:
                vers = 'no way to determine'
        elif name == 'trimal':
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()[0]
            vers = vers.strip()
        elif name == 'mafft':
            vers = subprocess.Popen(
                [name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()[1]
            vers = vers.strip()
        elif 'Launch_PASA' in name:
            vers = subprocess.Popen([name, '--version'],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True).communicate()[1]
            vers = vers.strip()
            vers = vers.split(': ')[-1]
        elif name == 'emapper.py':
            vers = subprocess.Popen([name, '--version'],
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True).communicate()[0]
            vers.strip()
            m = re.match(r'emapper-(\S+)', vers)
            if m:
                vers = m.group(1)
        elif name == 'pigz':
            (so, se) = subprocess.Popen([name, '--version'],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT,
                                        universal_newlines=True).communicate()
            for str in (so, se):
                m = re.search(r'pigz\s+(\S+)', str)
                if m:
                    vers = m.group(1)
                    break
        elif name == 'signalp':
            vers = subprocess.Popen(['signalp6', '--version'],
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True).communicate()[0]
        else:
            vers = subprocess.Popen([name, '--version'],
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True).communicate()[0].split('\n')[0]
        if 'exonerate' in vers:
            vers = vers.replace('exonerate from ', '')
        if 'SignalP' in vers:
            vers = vers.split(' ')[1]
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
            (so, se) = subprocess.Popen(
                [name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()
            m = re.search(r'pslCDnaFilter \[options\]', so+se)
            if m:
                vers = 'no way to determine'
            else:
                print(f"\tERROR: pslDnaFiler found but error running: {so}{se}")
        elif name == 'fasta':
            (se, so) = subprocess.Popen([name, '-h'],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True).communicate()
            m = re.search(r'version:\s+(\S+)', so+se)
            if m:
                vers = m.group(1)
            else:
                vers = 'no way to determine'
        elif name == 'CodingQuarry':
            vers = subprocess.Popen([name],
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True).communicate()
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


# Rust-accelerated Trinity utilities. These drop-in replace the slower Perl
# equivalents in Trinity's read-partitioning / coverage steps. Trinity resolves
# them by bare name on PATH and silently falls back to the Perl versions when
# they are absent, so their presence is what determines whether the Rust speedup
# is active. Installed/symlinked into the env bin by
# install_scripts/pixi_install_trinity.sh; toggled by toggle_trinity_rust.sh.
TRINITY_RUST_TOOLS = [
    'sam_to_read_coords',
    'extract_reads_per_partition',
    'fragment_coverage_writer',
    'define_coverage_partitions',
]


def check_trinity_rust(tools=None):
    """Report which Rust-accelerated Trinity utilities are available on PATH.

    Returns a dict mapping each tool name to its resolved executable path, or
    None when the tool is not found on PATH.
    """
    if tools is None:
        tools = TRINITY_RUST_TOOLS
    return {tool: lib.which_path(tool) for tool in tools}


def check_evm_rust():
    """Return the path to the Rust EVidenceModeler binary, or None.

    Presence of `evidence_modeler` on PATH is what lets funannotate run the Rust
    EVM engine (FUNANNOTATE_EVM_ENGINE=rust) instead of the Perl EVM pipeline.
    Installed by install_scripts/pixi_install_evm.sh.
    """
    return lib.which_path('evidence_modeler')


def check_pasa_rust():
    """Return the path to the Rust PASA binary (pasa_rust), or None.

    PASA's installer places pasa_rust under $PASAHOME/bin (not necessarily on
    PATH), so look there first and fall back to PATH. Installed by
    install_scripts/pixi_install_pasa.sh.
    """
    pasahome = os.environ.get('PASAHOME')
    if pasahome:
        candidate = os.path.join(pasahome, 'bin', 'pasa_rust')
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return lib.which_path('pasa_rust')


def get_version(program):
    # wrapper function to return version of programs
    if program in ['tblastn', 'makeblastdb', 'java', 'trimmomatic']:
        checker = check_version1
    elif program in ['exonerate', 'bedtools', 'bamtools', 'augustus',
                     'samtools', 'gmap', 'hisat2', 'Trinity',
                     'tbl2asn', 'table2asn', 'emapper.py', 'minimap2', 'mafft',
                     'trimal', 'stringtie', 'salmon', 'proteinortho', 'tantan',
                     'pigz']:
        checker = check_version2
    elif program in ['diamond', 'ete3', 'kallisto']:
        checker = check_version4
    elif program in ['gmes_petap.pl', 'blat', 'pslCDnaFilter', 'fasta',
                     'CodingQuarry', 'snap', 'glimmerhmm']:
        checker = check_version5
    elif program in ['hmmsearch', 'hmmscan', 'tRNAscan-SE']:
        checker = check_version6
    elif program in ['signalp']:
        if shutil.which('signalp6') is not None:
            checker = check_version2
        else:
            checker = check_version7
    else:
        return 'NA'
    return checker(program)


def main(args):
    funannotate_perl = ['Getopt::Long', 'Pod::Usage', 'File::Basename', 'threads', 'threads::shared',
                        'Thread::Queue', 'Carp', 'Data::Dumper', 'YAML', 'Hash::Merge', 'Logger::Simple', 'Parallel::ForkManager',
                        'DBI', 'Text::Soundex', 'Scalar::Util::Numeric', 'Tie::File', 'POSIX', 'Storable', 'Clone',
                        'DBD::mysql', 'JSON', 'LWP::UserAgent', 'DB_File', 'URI::Escape',
                        'File::Which', 'DBD::SQLite', 'local::lib']

    funannotate_python = ['numpy', 'pandas', 'matplotlib', 'scipy', 'scikit-learn',
                          'psutil', 'natsort', 'goatools', 'seaborn', 'biopython', 'requests']

    programs1 = ['tblastn', 'makeblastdb', 'java', 'trimmomatic']  # -version
    programs2 = ['exonerate', 'bedtools', 'bamtools', 'augustus',
                 'samtools', 'gmap', 'hisat2', 'Trinity',
                 'tbl2asn', 'emapper.py', 'minimap2', 'mafft',
                 'trimal', 'stringtie', 'salmon', 'proteinortho', 'tantan',
                 'pigz']  # --version
    programs3 = []  # -v
    programs4 = ['diamond', 'ete3', 'kallisto']  # version
    programs5 = ['gmes_petap.pl', 'blat', 'pslCDnaFilter', 'fasta',
                 'CodingQuarry', 'snap', 'glimmerhmm']  # no version option at all, a$$holes
    programs6 = ['hmmsearch', 'hmmscan', 'tRNAscan-SE']  # -h
    programs7 = []
    if shutil.which('signalp6') is not None:
        programs2.append('signalp')
    else:
        programs7.append('signalp')  # -V5 or lower
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
        if mod not in PyDeps:
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
        if mod not in PerlDeps:
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
                 'AUGUSTUS_CONFIG_PATH', 'GENEMARK_PATH']
    # EVM_HOME is only required if Rust EVM is not available
    has_rust_evm = check_evm_rust() is not None
    if not has_rust_evm:
        variables.append('EVM_HOME')

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

    # EVM_HOME can be set yet point at an incompatible EvidenceModeler. funannotate
    # drives the EVM 1.x Perl pipeline (EvmUtils/partition_EVM_inputs.pl); EVM 2.x
    # dropped those scripts and is not compatible. Without this check the mismatch
    # only surfaces as a cryptic failure deep inside `funannotate predict`
    # (see issues #1071, #1072, #1073, #1078). Flag it up front instead.
    if 'EVM_HOME' in os.environ:
        evm_script = os.path.join(
            os.environ['EVM_HOME'], 'EvmUtils', 'partition_EVM_inputs.pl')
        if not os.path.isfile(evm_script):
            print('\tERROR: $EVM_HOME is set but the EvidenceModeler 1.x scripts were '
                  'not found\n\t       (EvmUtils/partition_EVM_inputs.pl). funannotate '
                  'requires EvidenceModeler 1.x;\n\t       EVM 2.x is not compatible.')
        print("-------------------------------------------------------")

    if 'PASAHOME' not in missing:
        LAUNCHPASA = os.path.join(os.environ['PASAHOME'], 'Launch_PASA_pipeline.pl')
        programs2.append(LAUNCHPASA)
    print('Checking external dependencies...')
    for prog in programs1:
        if prog not in ExtDeps:
            ExtDeps[prog] = check_version1(prog)
    for prog in programs2:
        if prog not in ExtDeps:
            ExtDeps[prog] = check_version2(prog)
    for prog in programs3:
        if prog not in ExtDeps:
            ExtDeps[prog] = check_version3(prog)
    for prog in programs4:
        if prog not in ExtDeps:
            ExtDeps[prog] = check_version4(prog)
    for prog in programs5:
        if prog not in ExtDeps:
            ExtDeps[prog] = check_version5(prog)
    for prog in programs6:
        if prog not in ExtDeps:
            ExtDeps[prog] = check_version6(prog)
    for prog in programs7:
        if prog not in ExtDeps:
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

    # Report status of the custom Rust-optimized engines (EVM, PASA, Trinity).
    # These are optional accelerations with Perl fallbacks, so absence is a NOTE,
    # not an ERROR -- except when FUNANNOTATE_EVM_ENGINE=rust explicitly requests
    # the Rust EVM but its binary is missing.
    print("-------------------------------------------------------")
    print('Checking Rust-optimized engines (EVM, PASA, Trinity)...')

    # EVidenceModeler (Rust)
    evm_rust = check_evm_rust()
    evm_engine = os.environ.get('FUNANNOTATE_EVM_ENGINE', '').lower()
    if evm_rust:
        print('\tEVidenceModeler (Rust): ENABLED')
        if show:
            print('\t\tevidence_modeler: %s' % evm_rust)
    elif evm_engine == 'rust':
        print('\tERROR: FUNANNOTATE_EVM_ENGINE=rust but the Rust `evidence_modeler` '
              'binary was not found\n\t       on PATH. Run `pixi run install-evm`.')
    else:
        print('\tNOTE: Rust EVidenceModeler not found; will use the Perl EVM '
              '(requires $EVM_HOME).')

    # PASA (Rust pasa_rust)
    pasa_rust = check_pasa_rust()
    if pasa_rust:
        print('\tPASA (Rust): ENABLED')
        if show:
            print('\t\tpasa_rust: %s' % pasa_rust)
    else:
        print('\tNOTE: Rust PASA (pasa_rust) not found; PASA will use its '
              'standard pipeline.')

    # Trinity Rust utilities (only relevant when Trinity itself is installed)
    if ExtDeps.get('Trinity'):
        rust_status = check_trinity_rust()
        rust_found = [t for t, p in rust_status.items() if p]
        rust_missing = [t for t, p in rust_status.items() if not p]
        if not rust_missing:
            print('\tTrinity Rust utilities: ENABLED (%i/%i on PATH)'
                  % (len(rust_found), len(rust_status)))
            if show:
                for t in natsorted(rust_found):
                    print('\t\t%s: %s' % (t, rust_status[t]))
        elif rust_found:
            print('\tNOTE: Trinity Rust utilities PARTIAL (%i/%i on PATH); '
                  'Perl fallbacks for the rest'
                  % (len(rust_found), len(rust_status)))
            for t in natsorted(rust_missing):
                print('\t\tmissing: %s' % t)
        else:
            print('\tNOTE: Trinity Rust utilities not installed; '
                  'using the slower Perl equivalents.')
            if show:
                for t in natsorted(rust_missing):
                    print('\t\tmissing: %s' % t)
    print("-------------------------------------------------------")


if __name__ == "__main__":
    main(sys.argv[1:])
