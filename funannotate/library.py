# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from Bio import BiopythonWarning
import os
import uuid
import io
import subprocess
import logging
import sys
import csv
import time
import re
import shutil
import platform
import distro
import multiprocessing
import itertools
import hashlib
import math
import gzip
import operator
import textwrap
import errno
import datetime
from natsort import natsorted
import funannotate.resources as resources
from funannotate.interlap import InterLap
from collections import defaultdict
try:
    from urllib import unquote
except ImportError:
    from urllib.parse import unquote
try:
    from itertools import izip as zip
except ImportError:
    pass
import warnings
from Bio import SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO
warnings.simplefilter('ignore', BiopythonWarning)


# get the working directory, so you can move back into DB folder to find the files you need
global parentdir
parentdir = os.path.join(os.path.dirname(__file__))
GeneMark2GFF = os.path.join(parentdir, 'aux_scripts', 'genemark_gtf2gff3.pl')


class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'


class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).

    '''

    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])


class gzopen(object):
    """Generic opener that decompresses gzipped files
    if needed. Encapsulates an open file or a GzipFile.
    Use the same way you would use 'open()'.
    """

    def __init__(self, fname):
        f = open(fname)
        # Read magic number (the first 2 bytes) and rewind.
        magic_number = f.read(2)
        f.seek(0)
        # Encapsulated 'self.f' is a file or a GzipFile.
        if magic_number == b'\x1f\x8b':
            self.f = gzip.GzipFile(fileobj=f)
        else:
            self.f = f

    # Define '__enter__' and '__exit__' to use in
    # 'with' blocks. Always close the file and the
    # GzipFile if applicable.
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        try:
            self.f.fileobj.close()
        except AttributeError:
            pass
        finally:
            self.f.close()

    # Reproduce the interface of an open file
    # by encapsulation.
    def __getattr__(self, name):
        return getattr(self.f, name)

    def __iter__(self):
        return iter(self.f)

    def __next__(self):
        return next(self.f)


def createdir(name):
    try:
        os.makedirs(name)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass


def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)


def len_without_format(text):
    try:
        return len(remove_formatting(text))
    except TypeError:
        return len(str(text))


def remove_formatting(text):
    return re.sub('\033.*?m', '', text)


def colour(text, text_colour):
    bold_text = 'bold' in text_colour
    text_colour = text_colour.replace('bold', '')
    underline_text = 'underline' in text_colour
    text_colour = text_colour.replace('underline', '')
    text_colour = text_colour.replace('_', '')
    text_colour = text_colour.replace(' ', '')
    text_colour = text_colour.lower()
    if 'red' in text_colour:
        coloured_text = RED
    elif 'green' in text_colour:
        coloured_text = GREEN
    elif 'yellow' in text_colour:
        coloured_text = YELLOW
    elif 'dim' in text_colour:
        coloured_text = DIM
    else:
        coloured_text = ''
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def green(text):
    return GREEN + text + END_FORMATTING


def bold_green(text):
    return GREEN + BOLD + text + END_FORMATTING


def red(text):
    return RED + text + END_FORMATTING


def magenta(text):
    return MAGENTA + text + END_FORMATTING


def bold_red(text):
    return RED + BOLD + text + END_FORMATTING


def bold(text):
    return BOLD + text + END_FORMATTING


def bold_underline(text):
    return BOLD + UNDERLINE + text + END_FORMATTING


def underline(text):
    return UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def dim_underline(text):
    return DIM + UNDERLINE + text + END_FORMATTING


def bold_yellow(text):
    return YELLOW + BOLD + text + END_FORMATTING


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def bold_red_underline(text):
    return RED + BOLD + UNDERLINE + text + END_FORMATTING


def print_table(table, alignments='', max_col_width=30, col_separation=3, indent=2,
                row_colour=None, sub_colour=None, row_extra_text=None, leading_newline=False,
                subsequent_indent='', return_str=False, header_format='underline',
                hide_header=False, fixed_col_widths=None, left_align_header=True,
                bottom_align_header=True, verbosity=1):
    """
    Args:
        table: a list of lists of strings (one row is one list, all rows should be the same length)
        alignments: a string of L and R, indicating the alignment for each row
        max_col_width: values longer than this will be wrapped
        col_separation: the number of spaces between columns
        indent: the number of spaces between the table and the left side of the terminal
        row_colour: a dictionary of row indices and their colour names
        sub_colour: a dictionary of values to colour names for which the text colour will be set
        row_extra_text: a dictionary of row indices and extra text to display after the row
        leading_newline: if True, the function will print a blank line above the table
        subsequent_indent: this string will be added to the start of wrapped text lines
        return_str: if True, this function will return a string of the table instead of printing it
        header_format: the formatting (colour, underline, etc) of the header line
        hide_header: if True, the header is not printed
        fixed_col_widths: a list to specify exact column widths (automatic if not used)
        left_align_header: if False, the header will follow the column alignments
        bottom_align_header: if False, the header will align to the top, like other rows
        verbosity: the table will only be logged if the logger verbosity is >= this value
    """
    # this function is written by Ryan Wick in Unicycler code
    # modified to not support colors
    column_count = len(table[0])
    table = [x[:column_count] for x in table]
    table = [x + [''] * (column_count - len(x)) for x in table]
    if row_colour is None:
        row_colour = {}
    if sub_colour is None:
        sub_colour = {}
    if row_extra_text is None:
        row_extra_text = {}
    if leading_newline:
        print('')

    # Ensure the alignments string is the same length as the column count
    alignments += 'L' * (column_count - len(alignments))
    alignments = alignments[:column_count]

    if fixed_col_widths is not None:
        col_widths = fixed_col_widths
    else:
        col_widths = [0] * column_count
        for row in table:
            col_widths = [min(max(col_widths[i], len_without_format(x)), max_col_width)
                          for i, x in enumerate(row)]
    separator = ' ' * col_separation
    indenter = ' ' * indent
    full_table_str = ''
    for i, row in enumerate(table):
        row = [str(x) for x in row]
        if hide_header and i == 0:
            continue

        if fixed_col_widths is not None:
            wrapped_row = []
            for col, fixed_width in zip(row, fixed_col_widths):
                wrapper = textwrap.TextWrapper(subsequent_indent=subsequent_indent,
                                               width=fixed_width)
                wrapped_row.append(wrapper.wrap(col))
        else:
            wrapper = textwrap.TextWrapper(
                subsequent_indent=subsequent_indent, width=max_col_width)
            wrapped_row = [wrapper.wrap(x) for x in row]

        row_rows = max(len(x) for x in wrapped_row)
        if i == 0 and bottom_align_header:
            wrapped_row = [[''] * (row_rows - len(x)) + x for x in wrapped_row]

        for j in range(row_rows):
            row_line = [x[j] if j < len(x) else '' for x in wrapped_row]
            aligned_row = []
            for value, col_width, alignment in zip(row_line, col_widths, alignments):
                if alignment == 'L' or (i == 0 and left_align_header):
                    aligned_row.append(value.ljust(col_width))
                elif alignment == 'C':
                    aligned_row.append(value.center(col_width))
                else:
                    aligned_row.append(value.rjust(col_width))
            row_str = separator.join(aligned_row)
            if i in row_extra_text:
                row_str += row_extra_text[i]
            if i == 0 and header_format:
                row_str = colour(row_str, header_format)
            if i in row_colour:
                row_str = colour(row_str, row_colour[i])
            for text, colour_name in list(sub_colour.items()):
                row_str = row_str.replace(text, colour(text, colour_name))
            if j < row_rows - 1 and UNDERLINE in row_str:
                row_str = re.sub('\033\[4m', '', row_str)
            if return_str:
                full_table_str += indenter + row_str + '\n'
            else:
                print((indenter + row_str))
    if return_str:
        return full_table_str


def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=open(
            os.devnull, 'w'), cwd=parentdir).communicate()[0]
        return out
    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', '--short', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = False
    return GIT_REVISION


def Funzip(input, output, cpus):
    '''
    function to unzip as fast as it can, pigz -> bgzip -> gzip
    '''
    if which('pigz'):
        cmd = ['pigz', '--decompress', '-c', '-p', str(cpus), input]
    else:
        cmd = ['gzip', '--decompress', '-c', input]
    try:
        runSubprocess2(cmd, '.', log, output)
    except NameError:
        with open(output, 'w') as outfile:
            subprocess.call(cmd, stdout=outfile)


def Fzip(input, output, cpus):
    '''
    function to zip as fast as it can, pigz -> bgzip -> gzip
    '''
    if which('pigz'):
        cmd = ['pigz', '-c', '-p', str(cpus), input]
    else:
        cmd = ['gzip', '-c', input]
    try:
        runSubprocess2(cmd, '.', log, output)
    except NameError:
        with open(output, 'w') as outfile:
            subprocess.call(cmd, stdout=outfile)


def Fzip_inplace(input, cpus):
    '''
    function to zip as fast as it can, pigz -> bgzip -> gzip
    '''
    if which('pigz'):
        cmd = ['pigz', '-f', '-p', str(cpus), input]
    else:
        cmd = ['gzip', '-f', input]
    try:
        runSubprocess(cmd, '.', log)
    except NameError:
        subprocess.call(cmd)

# RNA seq mediated modules


def concatenateReads(input, output):
    '''
    Since I can't seem to get the comma separated lists to work with subprocess modules, just
    concatenate FASTQ files in order and use a single file, input should be a list of FASTQ files
    using system cat here so that gzipped files are concatenated correctly
    '''
    cmd = ['cat']
    cmd = cmd + input
    runSubprocess2(cmd, '.', log, output)


def which2(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def open_pipe(command, mode='r', buff=1024*1024):
    import subprocess
    import signal
    if 'r' in mode:
        return subprocess.Popen(command, shell=True, bufsize=buff,
                                stdout=subprocess.PIPE, universal_newlines=True,
                                preexec_fn=lambda: signal.signal(
                                    signal.SIGPIPE, signal.SIG_DFL)
                                ).stdout
    elif 'w' in mode:
        return subprocess.Popen(command, shell=True, bufsize=buff, universal_newlines=True,
                                stdin=subprocess.PIPE).stdin
    return None


NORMAL = 0
PROCESS = 1
PARALLEL = 2

WHICH_BZIP2 = which2("bzip2")
WHICH_PBZIP2 = which2("pbzip2")


def open_bz2(filename, mode='r', buff=1024*1024, external=PARALLEL):
    if external is None or external == NORMAL:
        import bz2
        return bz2.BZ2File(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_BZIP2:
            return open_bz2(filename, mode, buff, NORMAL)
        if 'r' in mode:
            return open_pipe("bzip2 -dc " + filename, mode, buff)
        elif 'w' in mode:
            return open_pipe("bzip2 >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PBZIP2:
            return open_bz2(filename, mode, buff, PROCESS)
        if 'r' in mode:
            return open_pipe("pbzip2 -dc " + filename, mode, buff)
        elif 'w' in mode:
            return open_pipe("pbzip2 >" + filename, mode, buff)
    return None


WHICH_GZIP = which2("gzip")
WHICH_PIGZ = which2("pigz")


def open_gz(filename, mode='r', buff=1024*1024, external=PARALLEL):
    if external is None or external == NORMAL:
        import gzip
        return gzip.GzipFile(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_GZIP:
            return open_gz(filename, mode, buff, NORMAL)
        if 'r' in mode:
            return open_pipe("gzip -dc " + filename, mode, buff)
        elif 'w' in mode:
            return open_pipe("gzip >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PIGZ:
            return open_gz(filename, mode, buff, PROCESS)
        if 'r' in mode:
            return open_pipe("pigz -dc " + filename, mode, buff)
        elif 'w' in mode:
            return open_pipe("pigz >" + filename, mode, buff)
    return None


WHICH_XZ = which2("xz")


def open_xz(filename, mode='r', buff=1024*1024, external=PARALLEL):
    if WHICH_XZ:
        if 'r' in mode:
            return open_pipe("xz -dc " + filename, mode, buff)
        elif 'w' in mode:
            return open_pipe("xz >" + filename, mode, buff)
    return None


def zopen(filename, mode='r', buff=1024*1024, external=PARALLEL):
    """
    Open pipe, zipped, or unzipped file automagically

    # external == 0: normal zip libraries
    # external == 1: (zcat, gzip) or (bzcat, bzip2)
    # external == 2: (pigz -dc, pigz) or (pbzip2 -dc, pbzip2)
    """
    if 'r' in mode and 'w' in mode:
        return None
    if filename.startswith('!'):
        return open_pipe(filename[1:], mode, buff)
    elif filename.endswith('.bz2'):
        return open_bz2(filename, mode, buff, external)
    elif filename.endswith('.gz'):
        return open_gz(filename, mode, buff, external)
    elif filename.endswith('.xz'):
        return open_xz(filename, mode, buff, external)
    else:
        return open(filename, mode, buff)
    return None


def execute(cmd):
    DEVNULL = open(os.devnull, 'w')
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             universal_newlines=True, stderr=DEVNULL)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def getDiamondVersion():
    vers = subprocess.Popen(['diamond', 'version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()
    if vers[1] == '':  # then this is older version and parse the stdout
        vers = vers[0].split('version ')[-1].rstrip()
    else:
        vers = vers[1].split()[1].replace('v', '')
    return vers


def CheckDiamondDB(database):
    diamond_version = getDiamondVersion()
    DBvers = None
    for line in execute(['diamond', 'dbinfo', '-d', database]):
        if line.startswith('Database format version ='):
            DBvers = int(line.strip().split('= ')[-1])
    if not DBvers:
        log.error('Could not determine diamond database version')
        return False
    runVers = None
    if diamond_version < '0.9.10':
        return False
    elif diamond_version < '0.9.25':
        runVers = 2
    else:
        runVers = 3
    if runVers >= DBvers:
        return True
    else:
        return False


def CheckFASTQandFix(forward, reverse, cpus=2):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    # open and check first header, if okay exit, if not fix
    file1 = FastqGeneralIterator(zopen(forward, 'rt'))
    file2 = FastqGeneralIterator(zopen(reverse, 'rt'))
    check = True
    for read1, read2 in zip(file1, file2):
        if ' ' in read1[0] and ' ' in read2[0]:
            # std illumina, exit
            if read1[0].split(' ', 1)[1].startswith('1') and read2[0].split(' ', 1)[1].startswith('2'):
                break
            else:
                log.debug("R1 header: {} and R2 header: {} are not 1 and 2 as expected".format(read1[0],read2[0]))
                check = False
                break
        elif read1[0].endswith('/1') and read2[0].endswith('/2'):  # also acceptable
            break
        else:  # it is not okay missing paired information
            log.debug("R1 header: {} and R2 header: {} are missing pairing as expected".format(read1[0],read2[0]))
            check = False
            break
    file1.close()
    file2.close()
    if not check:
        log.error('ERROR: FASTQ headers are not properly paired, see logfile and reformat your FASTQ headers')
        sys.exit(1)
        '''
        # now need to fix these reads
        log.info(
            "PE reads do not conform to Trinity naming convention (need either /1 /2 or std illumina), fixing...")
        # work on forward reads first
        if forward.endswith('.gz'):
            Funzip(forward, forward+'.bak', cpus)
            SafeRemove(forward)
        else:
            os.rename(forward, forward+'.bak')
        # now add ending to reads
        with open(forward+'.fix', 'w') as forwardfix:
            for title, seq, qual in FastqGeneralIterator(open(forward+'.bak')):
                title = title+'/1'
                forwardfix.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        Fzip(forward+'.fix', forward, cpus)
        SafeRemove(forward+'.bak')
        SafeRemove(forward+'.fix')
        # now work on reverse reads
        if reverse.endswith('.gz'):
            Funzip(reverse, reverse+'.bak', cpus)
        else:
            os.rename(reverse, reverse+'.bak')
        with open(reverse+'.fix', 'w') as reversefix:
            for title, seq, qual in FastqGeneralIterator(open(reverse+'.bak')):
                title = title+'/2'
                reversefix.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        # zip back up to original file
        Fzip(reverse+'.fix', reverse, cpus)
        SafeRemove(reverse+'.bak')
        SafeRemove(reverse+'.fix')
        '''
    else:
        log.debug('FASTQ headers seem compatible with Trinity')
    return 0


def SafeRemove(input):
    if os.path.isdir(input):
        shutil.rmtree(input)
    elif os.path.isfile(input):
        os.remove(input)
    else:
        return


def runSubprocess(cmd, dir, logfile):
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        logfile.error('CMD ERROR: {}'.format(' '.join(cmd)))
        if stdout:
            logfile.error(stdout.decode("utf-8"))
        if stderr:
            logfile.error(stderr.decode("utf-8"))
        sys.exit(1)
    else:
        if stdout:
            logfile.debug(stdout.decode("utf-8"))
        if stderr:
            logfile.debug(stderr.decode("utf-8"))


def runSubprocess2(cmd, dir, logfile, output):
    # function where output of cmd is STDOUT, capture STDERR in logfile
    logfile.debug(' '.join(cmd))
    with open(output, 'w') as out:
        proc = subprocess.Popen(cmd, cwd=dir, stdout=out,
                                stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if proc.returncode != 0:
        logfile.error('CMD ERROR: {}'.format(' '.join(cmd)))
        if stderr:
            logfile.error(stderr)
        sys.exit(1)
    else:
        if stderr:
            if stderr[0] is not None:
                logfile.debug(stderr)


def runSubprocess3(cmd, dir, logfile):
    # function where STDOUT pipes to FNULL, capture STDERR in logfile
    FNULL = open(os.devnull, 'w')
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=FNULL, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        logfile.debug(stderr)


def runSubprocess4(cmd, dir, logfile):
    # function where STDOUT and STDERR pipes to FNULL
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        logfile.error('CMD ERROR: {}'.format(' '.join(cmd)))
        if stdout:
            print(stdout)
        if stderr:
            print(stderr)
        sys.exit(1)


def runSubprocess5(cmd, dir, logfile, input, output):
    # function where STDOUT to file, STDIN as input, STDERR pipes to logfile
    logfile.debug(' '.join(cmd))
    with open(input) as infile:
        with open(output, 'w') as out:
            proc = subprocess.Popen(cmd, cwd=dir, stdin=infile, stdout=out,
                                    stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if proc.returncode != 0:
        logfile.error('CMD ERROR: {}'.format(' '.join(cmd)))
        if stderr:
            logfile.error(stderr)
        sys.exit(1)
    else:
        if stderr:
            if stderr[0] is not None:
                logfile.debug(stderr)


def runSubprocess6(cmd, dir, logfile, logfile2):
    # function where cmd captured in logfile, but both stdout and stdin piped to additional logfile
    logfile.debug(' '.join(cmd))
    with open(logfile2, 'w') as logout:
        proc = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            logfile.error('CMD ERROR: {}'.format(' '.join(cmd)))
            if stdout:
                logfile.error(stdout)
            if stderr:
                logfile.error(stderr)
            sys.exit(1)
        else:
            if stdout:
                logout.write(stdout)
            if stderr:
                logout.write(stderr)


def runSubprocess7(cmd, dir, logfile, output):
    # function where output of cmd is STDOUT, capture STDERR in logfile
    logfile.debug(' '.join(cmd))
    with open(output, 'a') as out:
        proc = subprocess.Popen(cmd, cwd=dir, stdout=out,
                                stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if proc.returncode != 0:
        logfile.error('CMD ERROR: {}'.format(' '.join(cmd)))
        if stderr:
            logfile.error(stderr)
        sys.exit(1)
    else:
        if stderr:
            if stderr[0] is not None:
                logfile.debug(stderr)


def runSubprocess8(cmd, dir, logfile, output):
    # function where output of cmd is STDOUT, capture STDERR in FNULL
    logfile.debug(' '.join(cmd))
    with open(output, 'w') as out:
        proc = subprocess.Popen(cmd, cwd=dir, stdout=out,
                                stderr=subprocess.PIPE)
        stderr = proc.communicate()
        if proc.returncode != 0:
            logfile.error('CMD ERROR: {}'.format(' '.join(cmd)))
            if stderr:
                logfile.error(stderr)
            sys.exit(1)


def evmGFFvalidate(input, evmpath, logfile):
    Validator = os.path.join(evmpath, 'EvmUtils', 'gff3_gene_prediction_file_validator.pl')
    cmd = ['perl', Validator, os.path.realpath(input)]
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = proc.communicate()
    if not stderr:
        return True
    else:
        logfile.error(stderr.rstrip())
        return False


def hashfile(afile, hasher, blocksize=65536):
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.digest()


def sha256_check(file1, file2):
    files = [file1, file2]
    output = [(fname, hashfile(open(fname, 'rb'), hashlib.sha256()))
              for fname in files]
    if output[0][1] == output[1][1]:
        return True
    else:
        return False


def readBlocks(source, pattern):
    buffer = []
    for line in source:
        try:
            line = line.decode('utf-8')
        except AttributeError:
            line = line
        if line.startswith(pattern):
            if buffer:
                yield buffer
            buffer = [line]
        else:
            buffer.append(line)
    yield buffer


def readBlocks2(source, startpattern, endpattern):
    buffer = []
    for line in source:
        try:
            line = line.decode('utf-8')
        except AttributeError:
            line = line
        if line.startswith(startpattern) or line.endswith(endpattern):
            if buffer:
                yield buffer
            buffer = [line]
        else:
            buffer.append(line)
    yield buffer


def empty_line_sep(line):
    return line == '\n'


def get_parent_dir(directory):
    return os.path.dirname(directory)


def getSize(filename):
    st = os.stat(filename)
    return st.st_size


def checkinputs(filename):
    if not os.path.isfile(filename):
        log.error("%s is not a valid file, exiting" % filename)
        sys.exit(1)
    size = getSize(filename)
    if size < 2:  # this is 1 character...
        log.error("%s appears to be empty, exiting" % filename)
        sys.exit(1)


def make_tarfile(output_filename, source_dir):
    import tarfile
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def multipleReplace(text, wordDict):
    for key in wordDict:
        text = text.replace(key, wordDict[key])
    return text


def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None


def which(name):
    try:
        with open(os.devnull) as devnull:
            diff = ['tbl2asn', 'dustmasker', 'mafft', 'signalp',
                    'proteinortho', 'ete3', 'phyml', 'phobius.pl', 'tantan']
            if not any(name in x for x in diff):
                subprocess.Popen([name], stdout=devnull,
                                 stderr=devnull, universal_newlines=True).communicate()
            else:
                if name == 'signalp':
                    subprocess.Popen([name, '-V'], stdout=devnull,
                                     stderr=devnull, universal_newlines=True).communicate()
                elif name == 'dustmasker':
                    subprocess.Popen(
                        [name, '-version-full'], stdout=devnull, stderr=devnull, universal_newlines=True).communicate()
                elif name == 'tbl2asn':
                    subprocess.Popen(
                        [name, '--help'], stdout=devnull, stderr=devnull, universal_newlines=True).communicate()
                elif name == 'raxmlHPC-PTHREADS':
                    subprocess.Popen(
                        [name, '-version'], stdout=devnull, stderr=devnull, universal_newlines=True).communicate()
                elif name == 'ete3':
                    subprocess.Popen(
                        [name, 'version'], stdout=devnull, stderr=devnull, universal_newlines=True).communicate()
                elif name == 'phobius.pl':
                    subprocess.Popen([name, '-h'], stdout=devnull,
                                     stderr=devnull, universal_newlines=True).communicate()
                else:
                    subprocess.Popen(
                        [name, '--version'], stdout=devnull, stderr=devnull, universal_newlines=True).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


def vers_tblastn():
    p1 = subprocess.Popen(['tblastn', '-version'],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    vers = p1.communicate()[0].split('+')[0]
    vers = vers.split(' ')[-1]
    return vers


def CheckDependencies(input):
    missing = []
    for p in input:
        if which(p) is False:
            missing.append(p)
    if missing != []:
        error = ", ".join(missing)
        try:
            log.error(
                "Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        except NameError:
            print("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        sys.exit(1)


def checkannotations(input):
    if input and os.path.isfile(input):
        filesize = getSize(input)
        if int(filesize) < 1:
            return False
        else:
            return True
    elif input and os.path.islink(input):
        return True
    else:
        return False


def line_count(fname):
    with open(fname) as f:
        i = -1
        for i, l in enumerate(f):
            pass
    return i + 1


def countfasta(input):
    count = 0
    with open(input, 'r') as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def getGeneBasename(fastafile):
    bases = []
    with open(fastafile, 'r') as input:
        for line in input:
            line = line.replace('\n', '')
            if line.startswith('>'):
                line = line.replace('>', '')
                transcript, gene = line.split(' ')
                if '_' in gene:
                    Base = line.split('_')[0]+'_'
                elif '-' in gene:
                    Base = line.split('-')[0]
                else:
                    Base = gene
                if not Base in bases:
                    bases.append(Base)
    return bases


def get_version():
    from pkg_resources import get_distribution
    __version__ = get_distribution('funannotate').version
    return __version__


def ver_tuple(z):
    return tuple([int(x) for x in z.split('.') if x.isdigit()])


def cmp(a, b):
    return (a > b) - (a < b)


def ver_cmp(a, b):
    return cmp(ver_tuple(a), ver_tuple(b))


def versionCheck(a, b):
    if ver_cmp(a, b) == -1:
        return False
    else:
        return True


def checkAugustusFunc():
    '''
    function to try to test Augustus installation is working, note segmentation fault still results in a pass
    '''
    functional = False
    p1 = subprocess.Popen(['augustus', '--version'], stderr=subprocess.STDOUT,
                          stdout=subprocess.PIPE, universal_newlines=True).communicate()
    stdout, stderr = p1
    if isinstance(stdout, str):
        try:
            stdout = stdout.decode('ascii', 'ignore').encode('ascii')
        except AttributeError:
            pass
    version = stdout.split(' is ')[0]
    model = os.path.join(parentdir, 'config', 'EOG092C0B3U.prfl')
    if not os.path.isfile(model):
        log.error("Testing Augustus Error: installation seems wrong, can't find prfl model")
        sys.exit(1)
    profile = '--proteinprofile='+model
    proc = subprocess.Popen(['augustus', '--species=anidulans', profile, os.path.join(parentdir, 'config', 'busco_test.fa')],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = proc.communicate()
    stderr = stderr.strip()
    if isinstance(stdout, str):
        try:
            stdout = stdout.decode('ascii', 'ignore').encode('ascii')
        except AttributeError:
            pass
    stdout = stdout.strip().split('\n')
    if stderr.startswith('augustus: ERROR'):
        print(stderr)
        return version, functional
    else:
        for line in stdout:
            line = line.strip()
            if line.startswith('# start gene g1'):
                functional = True
    return version, functional


def maker2evm(inputfile, outputdir):
    tr = os.path.join(outputdir, 'transcript_alignments.gff3')
    pr = os.path.join(outputdir, 'protein_alignments.gff3')
    gr = os.path.join(outputdir, 'gene_predictions.gff3')
    with open(tr, 'w') as trout:
        with open(pr, 'w') as prout:
            with open(gr, 'w') as grout:
                with open(inputfile, 'r') as input:
                    for line in input:
                        if line.startswith('#'):
                            continue
                        if 'trnascan' in line:
                            continue
                        cols = line.split('\t')
                        if 'maker' in cols[1]:
                            grout.write(line)
                        elif 'protein2genome' in cols[1]:
                            if 'match_part' in cols[2]:
                                cols[2] = 'nucleotide_to_protein_match'
                                cols[5] = '.'
                                prout.write('\t'.join(cols))
                        elif 'est2genome' in cols[1]:
                            if 'match_part' in cols[2]:
                                cols[2] = 'EST_match'
                                cols[5] = '.'
                                trout.write('\t'.join(cols))
                        elif 'cdna2genome' in cols[1]:
                            if 'match_part' in cols[2]:
                                cols[2] = 'EST_match'
                                cols[5] = '.'
                                trout.write('\t'.join(cols))
                        elif 'pred_gff' in cols[1]:
                            if 'match_part' in cols[2]:
                                cols[1] = cols[1].replace('pred_gff:', '')
                                cols[2] = 'EST_match'
                                cols[5] = '100.0'
                                trout.write('\t'.join(cols))


def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList


def fmtcols(mylist, cols):
    justify = []
    for i in range(0, cols):
        length = max([len(x) for x in mylist[i::cols]])
        length += 2
        ljust = [x.ljust(length) for x in mylist[i::cols]]
        justify.append(ljust)
    justify = flatten(justify)
    num_lines = len(mylist) / cols
    lines = (' '.join(justify[i::num_lines])
             for i in range(0, num_lines))
    return "\n".join(lines)


def list_columns(obj, cols=4, columnwise=True, gap=4):
    """
    Print the given list in evenly-spaced columns.

    Parameters
    ----------
    obj : list
        The list to be printed.
    cols : int
        The number of columns in which the list should be printed.
    columnwise : bool, default=True
        If True, the items in the list will be printed column-wise.
        If False the items in the list will be printed row-wise.
    gap : int
        The number of spaces that should separate the longest column
        item/s from the next column. This is the effective spacing
        between columns based on the maximum len() of the list items.
    """

    sobj = [str(item) for item in obj]
    if cols > len(sobj):
        cols = len(sobj)
    max_len = max([len(item) for item in sobj])
    if columnwise:
        cols = int(math.ceil(float(len(sobj)) / float(cols)))
    plist = [sobj[i: i+cols] for i in range(0, len(sobj), cols)]
    if columnwise:
        if not len(plist[-1]) == cols:
            plist[-1].extend(['']*(len(sobj) - len(plist[-1])))
        plist = list(zip(*plist))
    printer = '\n'.join([
        ''.join([c.ljust(max_len + gap) for c in p])
        for p in plist])
    return printer


def roundup(x):
    return x if x % 100 == 0 else x + 100 - x % 100


def maxabs(a, axis=None):
    import numpy as np
    """Return slice of a, keeping only those values that are furthest away
    from 0 along axis"""
    maxa = a.max(axis=axis)
    mina = a.min(axis=axis)
    p = abs(maxa) > abs(mina)  # bool, or indices where +ve values win
    n = abs(mina) > abs(maxa)  # bool, or indices where -ve values win
    if axis is None:
        if p:
            return maxa
        else:
            return mina
    shape = list(a.shape)
    shape.pop(axis)
    out = np.zeros(shape, dtype=a.dtype)
    out[p] = maxa[p]
    out[n] = mina[n]
    return out


def setupLogging(LOGNAME):
    global log
    if 'darwin' in sys.platform:
        stdoutformat = logging.Formatter(
            colr.GRN+'%(asctime)s'+colr.END+': %(message)s', datefmt='[%b %d %I:%M %p]')
    else:
        stdoutformat = logging.Formatter(
            '%(asctime)s: %(message)s', datefmt='[%b %d %I:%M %p]')
    fileformat = logging.Formatter(
        '%(asctime)s: %(message)s', datefmt='[%x %H:%M:%S]')
    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(stdoutformat)
    log.addHandler(sth)
    fhnd = logging.FileHandler(LOGNAME)
    fhnd.setLevel(logging.DEBUG)
    fhnd.setFormatter(fileformat)
    log.addHandler(fhnd)


def renameGFF(input, newname, output):
    contigs = set()
    with open(output, 'w') as outfile:
        with open(input, 'r') as infile:
            for line in infile:
                if line.startswith('>'):  # remove any fasta sequences
                    continue
                if line.startswith('#'):
                    outfile.write(line)
                else:
                    cols = line.split('\t')
                    # make sure it has correct columns to be GFF
                    if len(cols) == 9:
                        contigs.add(cols[0])
                        outfile.write('{}\t{}\t{}'.format(cols[0], newname,
                                                          '\t'.join(cols[2:])))
    return contigs


def countGFFgenes(input):
    count = 0
    if os.path.exists(input):
        with open(input, 'r') as f:
            for line in f:
                if "\tgene\t" in line:
                    count += 1
    return count


def countEVMpredictions(input):
    Counts = {'total': 0}
    with open(input, 'r') as f:
        for line in f:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.strip()
            contig, source, feature, start, end, blank, strand, score, info = line.split(
                '\t')
            if feature == 'gene':
                Counts['total'] += 1
                if not source in Counts:
                    Counts[source] = 1
                else:
                    Counts[source] += 1
    return Counts


def countGMAPtranscripts(input):
    count = 0
    with open(input, 'r') as f:
        for line in f:
            if line.startswith('###'):
                count += 1
    return count


def runMultiProgress(function, inputList, cpus, progress=True):
    # setup pool
    p = multiprocessing.Pool(cpus)
    # setup results and split over cpus
    tasks = len(inputList)
    results = []
    for i in inputList:
        results.append(p.apply_async(function, [i]))
    # refresh pbar every 5 seconds
    if progress:
        while True:
            incomplete_count = sum(1 for x in results if not x.ready())
            if incomplete_count == 0:
                break
            sys.stdout.write("     Progress: %.2f%% \r" %
                            (float(tasks - incomplete_count) / tasks * 100))
            sys.stdout.flush()
            time.sleep(1)
    p.close()
    p.join()


def runMultiNoProgress(function, inputList, cpus):
    # setup pool
    p = multiprocessing.Pool(cpus)
    # setup results and split over cpus
    results = []
    for i in inputList:
        results.append(p.apply_async(function, [i]))
    p.close()
    p.join()


def cleanProteins(inputList, output):
    # expecting a list of protein fasta files for combining/cleaning headers
    # make sure you aren't duplicated sequences names
    # dropping proteins less than 50 amino acids
    seen = set()
    with open(output, 'w') as out:
        for x in inputList:
            with open(x, 'r') as input:
                for rec in SeqIO.parse(input, 'fasta'):
                    if len(rec.seq) < 50:
                        continue
                    # explicitly check for swissprot and jgi
                    if rec.id.startswith('sp|') or rec.id.startswith('jgi|'):
                        ID = rec.id.split('|')[-1]
                    else:
                        ID = rec.id
                    # now clean up the shit
                    badshit = [':', ';', '/', '\\', '.', ',', '%']
                    for i in badshit:
                        if i in ID:
                            ID = ID.replace(i, '_')
                    if not ID in seen:
                        seen.add(ID)
                    else:
                        # means that ID has already been used, so add a number to it, auto increment
                        counter = 1
                        while ID in seen:
                            oldnum = counter-1
                            ID = ID.replace('_'+str(oldnum),
                                            '') + '_'+str(counter)
                            counter += 1
                        seen.add(ID)
                    out.write('>%s\n%s\n' % (ID, rec.seq))


def genemark2busco(genemark, bedfile, output):
    #function to load coords into Interlap from bedfile, then pull out
    #genemark EVM gff3 format
    counter = 0
    inter = bed2interlap(bedfile)
    with open(output, 'w') as outfile:
        with open(genemark, 'r') as infile:
            for gene_model in readBlocks(infile, '\n'):
                if len(gene_model) < 2:
                    continue
                if gene_model[0] == '\n':
                    cols = gene_model[1].split('\t')
                else:
                    cols = gene_model[0].split('\t')
                coords = [int(cols[3]), int(cols[4])]
                chr = cols[0]
                if interlapIntersect(coords, chr, inter):
                    counter += 1
                    outfile.write('{}'.format(''.join(gene_model)))
    return counter

def evidence2busco(evidence, bedfile, output):
    counter = 0
    inter = bed2interlap(bedfile)
    with open(output, 'w') as outfile:
        with open(evidence, 'r') as infile:
            for hit in readBlocks(infile, '\n'):
                hit = [x for x in hit if x != '\n']
                if len(hit) == 1:
                    start = int(hit[0].split('\t')[3])
                    end = int(hit[0].split('\t')[4])
                    coords = [start, end]
                    chr = hit[0].split('\t')[0]
                elif len(hit) > 1:
                    start = int(hit[0].split('\t')[3])
                    end = int(hit[-1].split('\t')[4])
                    chr = hit[0].split('\t')[0]
                    if start < end:
                        coords = [start, end]
                    else:
                        coords = [end, start]
                else:
                    continue
                if interlapIntersect(coords, chr, inter):
                    counter += 1
                    outfile.write('{}\n'.format(''.join(hit)))
    return counter


def fix_busco_naming(busco_infile, genome, augustus, gffout, ploidy=1,
                     proteins=False):
    def group_separator(line):
        return line == '\n'
    # parse the busco table into dictionary format
    busco_complete = {}
    passing = ['Complete']
    if ploidy > 1:
        passing.append('Duplicated')
    with open(busco_infile, 'r') as buscoinput:
        for line in buscoinput:
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            if cols[1] in passing:
                if not cols[0] in busco_complete:
                    busco_complete[cols[0]] = cols[2]+':'+cols[3]+'-'+cols[4]
    # now parse the augustus input file where gene numbers are likely repeated.
    results = []
    with open(augustus) as f:
        for key, group in itertools.groupby(f, group_separator):
            if not key:
                results.append(list(group))
    # loop through each gene model, lookup the BUSCO name, and then replace the name with counter based and busco model name
    tmpOut = augustus+'.intermediate'
    counter = 0
    inverse_busco = {v: k for k, v in list(busco_complete.items())}
    with open(tmpOut, 'w') as output:
        for i in results:
            counter += 1
            cols = i[0].split('\t')
            lookup = cols[0]+':'+cols[3]+'-'+cols[4]
            if lookup in inverse_busco:
                name = inverse_busco.get(lookup)
            else:
                name = 'unknown_model'
            ID = cols[8].split(';')[0]
            ID = ID.replace('ID=', '')
            newID = 'gene'+str(counter)
            newblock = ''.join(i)
            newblock = newblock.replace('Augustus%20prediction', name)
            newblock = newblock.replace(ID, newID)
            output.write(newblock+'\n')
    #write to GFF3 properly and fix CDS
    Genes = {}
    Genes = gff2dict(tmpOut, genome, Genes)
    dict2gff3(Genes, gffout)
    if proteins:
        dict2proteins(Genes, proteins)


def gb2output(input, output1, output2, output3):
    with open(output1, 'w') as proteins:
        with open(output2, 'w') as transcripts:
            with open(output3, 'w') as scaffolds:
                with open(input, 'r') as gbk:
                    SeqRecords = SeqIO.parse(gbk, 'genbank')
                    for record in SeqRecords:
                        scaffolds.write(">%s\n%s\n" % (record.id, record.seq))
                        for f in record.features:
                            if f.type == "CDS":
                                proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], softwrap(
                                    f.qualifiers['translation'][0].rstrip('*'))))
                            if f.type == "mRNA":
                                feature_seq = f.extract(record.seq)
                                transcripts.write(">%s\n%s\n" % (
                                    f.qualifiers['locus_tag'][0], softwrap(feature_seq)))


def sortGFF(input, output, order):
    cmd = ['bedtools', 'sort', '-header', '-faidx', order, '-i', input]
    with open(output, 'w') as out:
        proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        if stderr[0] is None:
            if stderr[1] != '':
                log.error(
                    "Sort GFF failed, unreferenced scaffold present in gene predictions, check logfile")
                sys.exit(1)


def sortBedproper(input, output):
    # sort BED file same as GFF3 files
    data = []
    with open(input, 'r') as infile:
        for line in infile:
            if line.startswith('\n'):
                continue
            line = line.rstrip()
            cols = line.split('\t')
            data.append(cols)
    # we can now sort
    sort_data = natsorted(data, key=lambda x: (x[0], int(x[1])))
    # now we can write back out to file
    with open(output, 'w') as outfile:
        for x in sort_data:
            outfile.write('{}\n'.format('\t'.join(x)))


def sortGFFproper(input, output):
    # function to sort GFF3 file but maintain gene, mrna, exon, cds order
    data = []
    features = set()
    comments = []
    with open(input, 'r') as infile:
        for line in infile:
            if line.startswith('\n'):
                continue
            if line.startswith('#'):
                comments.append(line)
                continue
            line = line.rstrip()
            cols = line.split('\t')
            data.append(cols)
            features.add(cols[2])
    # build sort order dictionary for features
    order_map = {'gene': 0, 'mRNA': 1, 'transcript': 2, 'tRNA': 3, 'ncRNA': 4,
                 'rRNA': 5, 'pseudogene': 6, 'five_prime_utr': 7,
                 'five_prime_UTR': 8, 'exon': 9, 'CDS': 10,
                 'three_prime_utr': 11, 'three_prime_UTR': 12}
    idx = len(order_map)
    for x in features:
        if x not in order_map:
            order_map[x] = idx
            idx += 1
    # we can now sort
    sort_data = natsorted(data, key=lambda x: (x[0], int(x[3]), order_map[x[2]]))
    # now we can write back out to file
    with open(output, 'w') as outfile:
        for y in comments:
            outfile.write(y)
        for x in sort_data:
            outfile.write('{}\n'.format('\t'.join(x)))


def checkGenBank(input):
    count = 0
    with open(input, 'r') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                if f.type == 'CDS':
                    count += 1
    if count == 0:
        return False
    else:
        return True


def countGenBank(input):
    cds = 0
    trna = 0
    dnas = 0
    with open(input, 'r') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            dnas += 1
            for f in record.features:
                if f.type == 'CDS':
                    cds += 1
                elif f.type == 'tRNA':
                    trna += 1
    return dnas, cds, trna


def checkFastaHeaders(input, limit):
    length = 0
    names = []
    with open(input, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                line = line.replace('\n', '')
                ID = line.replace('>', '').strip()
                names.append(ID)
                # subtract one character for fasta carrot
                headlen = len(line) - 1
                if headlen > length:
                    length = headlen
    if length > int(limit):
        return (False, names)
    else:
        return (True, names)


def analyzeAssembly(input, header_max=16):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    bad_names = []
    IUPAC = {'A', 'C', 'G', 'T', 'R', 'Y',
             'S', 'W', 'K', 'M', 'B', 'D',
             'H', 'V', 'N'}
    nuc_errors = {}
    suspect = {}
    with open(input, 'r') as infile:
        for title, seq in SimpleFastaParser(infile):
            if len(title) > header_max:
                bad_names.append(title)
            # get number of each contig
            characters = {}
            for nuc in seq:
                nuc = nuc.upper()
                if not nuc in characters:
                    characters[nuc] = 1
                else:
                    characters[nuc] += 1
            # check for non IUPAC characters
            errors = []
            for k, v in characters.items():
                if k not in IUPAC:
                    errors.append((k, v))
            if len(errors) > 0:
                nuc_errors[title] = errors
            # if there are less than 4 characters in scaffolds, its suspect
            if len(characters) < 4:
                suspect[title] = characters
    return bad_names, nuc_errors, suspect


def BamHeaderTest(genome, mapping):
    # get list of fasta headers from genome
    genome_headers = []
    with open(genome, 'r') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            if rec.id not in genome_headers:
                genome_headers.append(rec.id)
    # get list of fasta headers from BAM
    bam_headers = []
    cmd = ['samtools', 'idxstats', os.path.realpath(mapping)]
    for line in execute(cmd):
        line = line.rstrip()
        chr, length, mapped, unmapped = line.split('\t')[:4]
        if chr != '*':
            bam_headers.append(chr)
    # now compare lists, basically if BAM headers not in genome headers, then output bad names to logfile and return FALSE
    genome_headers = set(genome_headers)
    diffs = [x for x in bam_headers if x not in genome_headers]
    if len(diffs) > 0:
        log.debug(
            "ERROR: These BAM headers not found in genome FASTA headers\n%s" % ','.join(diffs))
        return False
    else:
        return True


def mapCount(input, location_dict, output):
    Counts = {}
    for aln in execute(['samtools', 'view', os.path.realpath(input)]):
        cols = aln.split('\t')
        if not cols[2] in Counts:
            Counts[cols[2]] = 1
        else:
            Counts[cols[2]] += 1
    with open(output, 'w') as outfile:
        outfile.write("#mRNA-ID\tgene-ID\tLocation\tTPM\n")
        for k, v in natsorted(list(location_dict.items())):
            if k in Counts:
                tpm = Counts.get(k)
            else:
                tpm = 0
            geneID = v[0]
            location = v[1]
            outfile.write('{:}\t{:}\t{:}\t{:.2f}\n'.format(
                k, geneID, location, float(tpm)))


def tokenizeString(aString, separators):
    # separators is an array of strings that are being used to split the the string.
    # sort separators in order of descending length
    separators.sort(key=len)
    listToReturn = []
    i = 0
    while i < len(aString):
        theSeparator = ""
        for current in separators:
            if current == aString[i:i+len(current)]:
                theSeparator = current
        if theSeparator != "":
            listToReturn += [theSeparator]
            i = i + len(theSeparator)
        else:
            if listToReturn == []:
                listToReturn = [""]
            if(listToReturn[-1] in separators):
                listToReturn += [""]
            listToReturn[-1] += aString[i]
            i += 1
    return listToReturn


def bam2gff3(input, output):
    count = 0
    with open(output, 'w') as gffout:
        gffout.write('##gff-version 3\n')
        for aln in execute(['samtools', 'view', os.path.realpath(input)]):
            cols = aln.split('\t')
            if cols[1] == '0':
                strand = '+'
            elif cols[1] == '16':
                strand = '-'
            else:
                continue
            cs = None
            nm = None
            tags = cols[11:]
            if not tags:
                continue
            for x in tags:
                if x.startswith('cs:'):
                    cs = x.replace('cs:Z:', '')
                if x.startswith('NM:'):
                    nm = int(x.split(':')[-1])
            if nm is None or cs is None:
                continue
            matches = 0
            ProperSplice = True
            splitter = []
            exons = [int(cols[3])]
            position = int(cols[3])
            query = [1]
            querypos = 0
            num_exons = 1
            gaps = 0
            splitter = tokenizeString(cs, [':', '*', '+', '-', '~'])
            for i, x in enumerate(splitter):
                if x == ':':
                    matches += int(splitter[i+1])
                    position += int(splitter[i+1])
                    querypos += int(splitter[i+1])
                elif x == '-':
                    gaps += 1
                elif x == '+':
                    gaps += 1
                    querypos += len(splitter[i+1])
                elif x == '~':
                    if cols[1] == '0':
                        if splitter[i+1].startswith('gt') and splitter[i+1].endswith('ag'):
                            ProperSplice = True
                        elif splitter[i+1].startswith('at') and splitter[i+1].endswith('ac'):
                            ProperSplice = True
                        else:
                            ProperSplice = False
                    elif cols[1] == '16':
                        if splitter[i+1].startswith('ct') and splitter[i+1].endswith('ac'):
                            ProperSplice = True
                        elif splitter[i+1].startswith('gt') and splitter[i+1].endswith('at'):
                            ProperSplice = True
                        else:
                            ProperSplice = False
                    num_exons += 1
                    exons.append(position)
                    query.append(querypos)
                    query.append(querypos+1)
                    intronLen = int(splitter[i+1][2:-2])
                    position += intronLen
                    exons.append(position)
            # add last Position
            exons.append(position)
            query.append(len(cols[9]))
            # convert exon list into list of exon tuples
            exons = list(zip(exons[0::2], exons[1::2]))
            queries = list(zip(query[0::2], query[1::2]))
            if ProperSplice:
                mismatches = nm - gaps
                pident = 100 * (matches / (matches + mismatches))
                if pident < 80:
                    continue
                count += 1
                for i, exon in enumerate(exons):
                    start = exon[0]
                    end = exon[1]-1
                    if strand == '+':
                        qstart = queries[i][0]
                        qend = queries[i][1]
                    else:
                        qstart = len(cols[9]) - queries[i][1] + 1
                        qend = len(cols[9]) - queries[i][0] + 1
                    gffout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:.2f}\t{:}\t{:}\tID={:};Target={:} {:} {:}\n'.format(
                        cols[2], 'genome', 'cDNA_match', start, end, pident, strand, '.', cols[0], cols[0], qstart, qend))
    return count


def bam2ExonsHints(input, gff3, hints):
    count = 0
    with open(gff3, 'w') as gffout:
        gffout.write('##gff-version 3\n')
        with open(hints, 'w') as hintsout:
            num = -1
            for aln in execute(['samtools', 'view', os.path.realpath(input)]):
                num += 1
                cols = aln.split('\t')
                if cols[1] == '0':
                    strand = '+'
                elif cols[1] == '16':
                    strand = '-'
                else:
                    continue
                cs = None
                nm = None
                tags = cols[11:]
                for x in tags:
                    if x.startswith('cs:'):
                        cs = x.replace('cs:Z:', '')
                    if x.startswith('NM:'):
                        nm = int(x.split(':')[-1])
                if nm is None or cs is None:
                    continue
                matches = 0
                ProperSplice = True
                splitter = []
                exons = [int(cols[3])]
                position = int(cols[3])
                query = [1]
                querypos = 0
                num_exons = 1
                gaps = 0
                splitter = tokenizeString(cs, [':', '*', '+', '-', '~'])
                for i, x in enumerate(splitter):
                    if x == ':':
                        matches += int(splitter[i+1])
                        position += int(splitter[i+1])
                        querypos += int(splitter[i+1])
                    elif x == '-':
                        gaps += 1
                    elif x == '+':
                        gaps += 1
                        querypos += len(splitter[i+1])
                    elif x == '~':
                        if cols[1] == 0:
                            if splitter[i+1].startswith('gt') and splitter[i+1].endswith('ag'):
                                ProperSplice = True
                            elif splitter[i+1].startswith('at') and splitter[i+1].endswith('ac'):
                                ProperSplice = True
                            else:
                                ProperSplice = False
                                break
                        elif cols[1] == 16:
                            if splitter[i+1].startswith('ct') and splitter[i+1].endswith('ac'):
                                ProperSplice = True
                            elif splitter[i+1].startswith('gt') and splitter[i+1].endswith('at'):
                                ProperSplice = True
                            else:
                                ProperSplice = False
                                break
                        num_exons += 1
                        exons.append(position)
                        query.append(querypos)
                        query.append(querypos+1)
                        intronLen = int(splitter[i+1][2:-2])
                        position += intronLen
                        exons.append(position)
                # add last Position
                exons.append(position)
                query.append(len(cols[9]))

                # convert exon list into list of exon tuples
                exons = list(zip(exons[0::2], exons[1::2]))
                queries = list(zip(query[0::2], query[1::2]))
                introns = []
                if len(exons) > 1:
                    for x, y in enumerate(exons):
                        try:
                            introns.append((y[1], exons[x+1][0]-1))
                        except IndexError:
                            pass
                if ProperSplice:
                    mismatches = nm - gaps
                    pident = 100 * (matches / (matches + mismatches))
                    if pident < 80:
                        continue
                    feature = 'EST_match'
                    if pident > 95:
                        feature = 'cDNA_match'
                    count += 1
                    for i, exon in enumerate(exons):
                        start = exon[0]
                        end = exon[1]-1
                        qstart = queries[i][0]
                        qend = queries[i][1]
                        if i == 0 or i == len(exons)-1:
                            gffout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:.2f}\t{:}\t{:}\tID=minimap2_{:};Target={:} {:} {:} {:}\n'.format(
                                cols[2], 'genome', feature, start, end, pident, strand, '.', num+1, cols[0], qstart, qend, strand))
                            hintsout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tgrp=minimap2_{:};pri=4;src=E\n'.format(
                                cols[2], 'b2h', 'ep', start, end, 0, strand, '.', num+1, cols[0]))
                        else:
                            gffout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:.2f}\t{:}\t{:}\tID=minimap2_{:};Target={:} {:} {:} {:}\n'.format(
                                cols[2], 'genome', feature, start, end, pident, strand, '.', num+1, cols[0], qstart, qend, strand))
                            hintsout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tgrp=minimap2_{:};pri=4;src=E\n'.format(
                                cols[2], 'b2h', 'exon', start, end, 0, strand, '.', num+1, cols[0]))
                    if len(introns) > 0:
                        for z in introns:
                            hintsout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tgrp=minimap2_{:};pri=4;src=E\n'.format(
                                cols[2], 'b2h', 'intron', z[0], z[1], 1, strand, '.', num+1, cols[0]))
    return count


def combineTranscripts(minimap, gmap, output):
    '''
    function to combine minimap GFF3 and gmap GFF3 files
    need to rename GMAP as you loop through and GFF3 from gmap is kind of messed up.
    '''
    with open(output, 'w') as out:
        if minimap:
            with open(minimap, 'r') as mini:
                for line in mini:
                    out.write(line)
        else:
            out.write('##gff-version 3\n')
        with open(gmap, 'r') as gmap_in:
            for i, aln in enumerate(readBlocks(gmap_in, '###')):
                for x in aln:
                    if not x.startswith('#'):
                        contig, source, feature, start, end, score, strand, phase, attributes = x.split(
                            '\t')
                        info = attributes.split(';')
                        for y in info:
                            if y.startswith('Target='):
                                Target = y
                        out.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID=gmap_{:};{:}\n'.format(
                            contig, source, feature, start, end, score, strand, phase, i+1, Target))


def RevComp(s):
    rev_comp_lib = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'M': 'K', 'R': 'Y', 'W': 'W',
                    'S': 'S', 'Y': 'R', 'K': 'M', 'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V', 'X': 'X', 'N': 'N'}
    cseq = ''
    n = len(s)
    s = s.upper()
    for i in range(0, n):
        c = s[n-i-1]
        cseq += rev_comp_lib[c]
    return cseq


def translate(cDNA, strand, phase):
    '''
    translate cDNA into protein sequence
    trying to see if I can speed this up over Biopython
    '''
    def _split(str, num):
        return [str[start:start+num] for start in range(0, len(str), num)]
    codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
                   'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
                   'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
                   'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
                   'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                   'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
                   'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
                   'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
                   'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                   'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
                   'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
                   'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
                   'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}
    if strand == '-' or strand == -1:
        seq = RevComp(cDNA)
    else:
        seq = cDNA
    seq = seq[phase:]
    # map seq to proteins
    protSeq = []
    for i in _split(seq, 3):
        if len(i) == 3:
            iSeq = i.upper()
            if iSeq in codon_table:
                aa = codon_table[iSeq]
                protSeq.append(aa)
            else:
                protSeq.append('X')
    return ''.join(protSeq)


def extend2stop(seqDict, header, coordinates, strand, phase, protLen):
    '''
    try to extend a CDS lacking a stop to find a stop codon
    it will extend a CDS up to 20 codons (60 bp) from the current
    frame to find a stop codon, if none is found it will return
    the original coordinates
    '''
    sorted_coordinates = sorted(coordinates, key=lambda tup: tup[0])
    if strand == '+':
        newStop = sorted_coordinates[-1][1]+60
        if newStop > len(seqDict[header]):
            newStop = len(seqDict[header])
        lastTup = (sorted_coordinates[-1][0], newStop)
        if len(sorted_coordinates) > 1:
            newCoords = sorted_coordinates[:-1]
            newCoords.append(lastTup)
        else:
            newCoords = [lastTup]
        updateCDS = getSeqRegions(seqDict, header, newCoords)
        updateProt = translate(updateCDS, strand, phase)
        if '*' in updateProt:
            num = (updateProt.find('*') - protLen + 1) * 3
            finalTup = (sorted_coordinates[-1][0],
                        sorted_coordinates[-1][1]+num)
            if len(sorted_coordinates) > 1:
                finalCoords = sorted_coordinates[:-1]
                finalCoords.append(finalTup)
            else:
                finalCoords = [finalTup]
            return True, finalCoords
        else:
            return False, coordinates
    else:
        newStop = sorted_coordinates[0][0]-60
        if newStop < 1:
            newStop = 1
        lastTup = (newStop, sorted_coordinates[0][1])
        newCoords = [lastTup]
        if len(sorted_coordinates) > 1:
            newCoords += sorted_coordinates[1:]
        updateCDS = getSeqRegions(seqDict, header, newCoords)
        updateProt = translate(updateCDS, strand, phase)
        if '*' in updateProt:
            num = (updateProt.find('*') - protLen + 1) * 3
            finalTup = (sorted_coordinates[0][0]-num, sorted_coordinates[0][1])
            finalCoords = [finalTup]
            if len(sorted_coordinates) > 1:
                finalCoords += sorted_coordinates[1:]
            finalSort = sorted(
                finalCoords, key=lambda tup: tup[0], reverse=True)
            return True, finalSort
        else:
            return False, coordinates


def getSeqRegions(SeqRecordDict, header, coordinates):
    # takes SeqRecord dictionary or Index, returns sequence string
    # coordinates is a list of tuples [(1,10), (20,30)]
    result = ''
    sorted_coordinates = sorted(coordinates, key=lambda tup: tup[0])
    for x in sorted_coordinates:
        partial = SeqRecordDict[header][x[0]-1:x[1]]
        result += str(partial.seq)
    return result


def convertgff2tbl(gff, prefix, fasta, prots, trans, tblout, external=False):
    from collections import OrderedDict
    '''
    function to convert directly from gff to tbl
    '''

    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])

    # load GFF annotations into funannotate dictionary
    Genes = {}
    Genes = gff2dict(gff, fasta, Genes)
    # get scaffold names/lengths
    scaffLen = {}
    with open(fasta, 'r') as seqin:
        for record in SeqIO.parse(seqin, 'fasta'):
            if not record.id in scaffLen:
                scaffLen[record.id] = len(record.seq)
    # get partialStart/stop info and load scaffold dictionary with coordinates of Genes
    sGenes = sorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    scaff2genes = {}
    counter = 1
    for k, v in list(sortedGenes.items()):
        if not prefix:
            locusTag = k
        else:
            locusTag = prefix+'_'+str(counter).zfill(6)
        if not locusTag in renamedGenes:
            renamedGenes[locusTag] = v
            if not v['contig'] in scaff2genes:
                scaff2genes[v['contig']] = [locusTag]
            else:
                scaff2genes[v['contig']].append(locusTag)
        counter += 1
    if external:
        log.info('Found {:,} gene models from GFF3 annotation'.format(len(sortedGenes)))
    dicts2tbl(renamedGenes, scaff2genes, scaffLen, 'CFMR', '12345', [],
              tblout, external=external)
    # transcript to geneID dictionary
    geneDB = {}
    for k, v in list(renamedGenes.items()):
        for x in v['ids']:
            if not x in geneDB:
                geneDB[x] = k

    # write to protein and transcripts
    with open(prots, 'w') as protout:
        with open(trans, 'w') as tranout:
            for k, v in natsorted(list(Genes.items())):
                if v['pseudo'] and v['pseudo'] is True:
                    continue
                for i, x in enumerate(v['ids']):
                    try:
                        Transcript = str(v['transcript'][i])
                    except IndexError:
                        print((k, v))
                    if v['strand'] == '-':
                        Transcript = RevComp(Transcript)
                    tranout.write('>%s %s\n%s\n' % (x, k, softwrap(Transcript)))
                    if v['type'] == 'mRNA':
                        Prot = v['protein'][i]
                        if Prot.endswith('*'):
                            Prot = Prot.rstrip('*')
                        protout.write('>%s %s\n%s\n' % (x, k, softwrap(Prot)))

    return len(Genes), geneDB


def tblfilter(input, remove, output):
    '''
    function to take an NCBI tbl file and drop gene models present in remove file
    '''
    # get items to remove list
    removeModels = []
    with open(remove, 'r') as file:
        for line in file:
            if line.startswith('#') or line.startswith('\n'):
                continue
            line = line.strip()
            if not line in removeModels:
                removeModels.append(line)
    # now loop through tbl file and get line positions of gene models
    found = []
    with open(output, 'w') as outfile:
        with open(input, 'r') as infile:
            for gene in readBlocks2(infile, '>Feature', '\tgene\n'):
                if gene[0].startswith('>Feature'):
                    outfile.write(''.join(gene))
                else:
                    locusTag = None
                    for x in gene:
                        if x.startswith('\t\t\tlocus_tag\t'):
                            locusTag = x.split('\t')[-1].rstrip()
                    if locusTag and not locusTag in removeModels:
                        outfile.write(''.join(gene))
                    else:
                        if not locusTag:
                            log.debug(
                                'LocusTag not found parsing NCBI Tbl file (this should not happen)')
                            print(gene)
                        else:
                            found.append(locusTag)
    log.debug("Removed %i out of %i gene models from annotation" %
              (len(found), len(removeModels)))
    s = set(found)
    diff = [x for x in removeModels if x not in s]
    if len(diff) > 0:
        log.debug('Could not find %i gene models:\n %s' %
                  (len(diff), ','.join(diff)))


def annotations2dict(input, geneDB={}, custom=False):
    Annotations = {}
    with open(input, 'r') as all_annots:
        for line in all_annots:
            line = line.replace('\n', '')
            ID, refDB, description = line.split('\t')
            if description == '':  # there is nothing here, so skip
                continue
            if refDB == 'name' or refDB == 'product':
                if len(geneDB) == 0:
                    if '-T' in ID:
                        geneID = ID.split('-T')[0]
                    else:
                        geneID = ID
                else:
                    if ID in geneDB:
                        geneID = geneDB[ID]
                    else:
                        geneID = ID
            else:
                geneID = ID
            if not geneID in Annotations:
                Annotations[geneID] = {refDB: [description]}
            else:
                if not refDB in Annotations[geneID]:
                    Annotations[geneID][refDB] = [description]
                else:
                    Annotations[geneID][refDB].append(description)
    if custom:
        log.info("Parsing custom annotations from {:}".format(custom))
        with open(custom, 'r') as custom_annots:
            for line in custom_annots:
                line = line.rstrip()
                try:
                    if line.count('\t') != 2:
                        continue
                except UnicodeDecodeError:
                    log.error('Error parsing the custom annotations:')
                    print(line)
                    sys.exit(1)
                ID, refDB, description = line.split('\t')
                if description == '':
                    continue
                if refDB in ['name', 'product', 'gene_synonym']:
                    if len(geneDB) == 0:
                        if '-T' in ID:
                            geneID = ID.split('-T')[0]
                        else:
                            geneID = ID
                    else:
                        if ID in geneDB:
                            geneID = geneDB[ID]
                        else:
                            geneID = ID
                else:
                    geneID = ID
                if not geneID in Annotations:
                    Annotations[geneID] = {refDB: [description]}
                else:
                    if not refDB in Annotations[geneID]:
                        Annotations[geneID][refDB] = [description]
                    elif refDB == 'name':
                        previousNames = Annotations[geneID][refDB]
                        if not 'gene_synonym' in Annotations[geneID]:
                            Annotations[geneID]['gene_synonym'] = previousNames
                        else:
                            Annotations[geneID]['gene_synonym'] += previousNames
                        Annotations[geneID][refDB] = [description]
                    elif refDB == 'product':
                        Annotations[geneID][refDB] = [description]
                    else:
                        Annotations[geneID][refDB].append(description)
    # make sure no synonyms are repeated
    for k, v in natsorted(Annotations.items()):
        if 'gene_synonym' in v and 'name' in v:
            synonym_set = set(v['gene_synonym'])
            cleaned = [x for x in synonym_set if x not in v['name']]
            Annotations[k]['gene_synonym'] = cleaned
        elif 'gene_synonm' in v:
            synonym_set = set(v['gene_synonym'])
            Annotations[k]['gene_synonym'] = list(synonym_set)
    return Annotations


def updateTBL(input, annotDict, output, prefix=False, newtag=False):
    '''
    general function to parse ncbi tbl format and add functional annotation
    '''
    log.debug('Parsing tbl file: {:}'.format(os.path.abspath(input)))
    tmpoutput = output+'.tmp'
    with open(input, 'r') as infile:
        with open(tmpoutput, 'w') as outfile:
            for gene in readBlocks2(infile, '>Feature', '\tgene\n'):
                transcriptsSeen = []
                # transcriptNum = 0
                if gene[0].startswith('>Feature'):
                    outfile.write(''.join(gene))
                else:
                    locusTag, locusTagIndex, LocusType, geneAnnot, transcriptAnnot = (
                        None,)*5
                    for i, x in enumerate(gene):
                        if x.startswith('\t\t\tlocus_tag\t'):
                            locusTag = x.split('\t')[-1].rstrip()
                            locusTagIndex = i
                    if not locusTagIndex:
                        outfile.write(''.join(gene))
                        continue
                    try:
                        locusType = gene[locusTagIndex+1].split('\t')[-1].rstrip()
                    except IndexError:
                        print(gene)
                    except TypeError:
                        print(gene)
                    if locusType in ['tRNA', 'ncRNA', 'rRNA']:
                        outfile.write(''.join(gene))
                    elif locusType == 'mRNA':
                        if locusTag in annotDict:
                            geneAnnot = annotDict.get(locusTag)
                        else:
                            geneAnnot = {}
                        for line in gene:
                            if line.startswith('\t\t\tlocus_tag\t'):
                                if 'name' in geneAnnot:
                                    outfile.write('\t\t\tgene\t%s\n' %
                                                  geneAnnot['name'][0])
                                if 'gene_synonym' in geneAnnot:
                                    for z in set(geneAnnot['gene_synonym']):
                                        outfile.write('\t\t\tgene_synonym\t%s\n' % z)
                                outfile.write(line)
                            elif line.startswith('\t\t\tproduct\t'):
                                if not 'product' in geneAnnot:
                                    outfile.write(line)
                            elif line.startswith('\t\t\ttranscript_id\t'):
                                ID = line.split('|')[-1]
                                ID = ID.split('_mrna')[0]
                                if not ID in transcriptsSeen:
                                    transcriptsSeen.append(ID)
                                transcriptNum = len(transcriptsSeen)
                                if ID in annotDict:
                                    transcriptAnnot = annotDict.get(ID)
                                if 'product' in geneAnnot:
                                    Description = geneAnnot['product'][0]
                                    if transcriptNum > 1:
                                        Description = Description + ', variant {:}'.format(transcriptNum)
                                    outfile.write('\t\t\tproduct\t%s\n' % Description)
                                outfile.write(line)
                            elif line.startswith('\t\t\tcodon_start\t'):
                                outfile.write(line)
                                if transcriptAnnot:
                                    for item in transcriptAnnot:
                                        if item in ['name', 'product', 'gene_synonym']:
                                            continue
                                        for x in set(transcriptAnnot[item]):
                                            outfile.write('\t\t\t%s\t%s\n' % (item, x))
                            else:
                                outfile.write(line)
    if newtag:
        with open(output, 'w') as outfile:
            with open(tmpoutput, 'r') as infile:
                for line in infile:
                    if line.startswith('\t\t\tlocus_tag\t'):
                        line = line.replace('\t'+prefix, '\t'+newtag)
                    elif line.startswith('\t\t\ttranscript_id\t') or line.startswith('\t\t\tprotein_id\t'):
                        line = line.replace('|'+prefix, '|'+newtag)
                    outfile.write(line)
        os.remove(tmpoutput)
    else:
        os.rename(tmpoutput, output)


def bed2gff3(input, output):
    '''
    convert repeats bed file into GFF3 format
    Contig245   36  69  Repeat_1
    Contig245   265 288 Repeat_2
    Contig245   477 493 Repeat_3
    Contig245   780 797 Repeat_4
    Contig245   997 1016    Repeat_5
    '''
    with open(output, 'w') as outfile:
        outfile.write("##gff-version 3\n")
        with open(input, 'r') as bedfile:
            for line in bedfile:
                line = line.strip()
                if line.startswith('\n'):
                    continue
                contig, start, end, name = line.split('\t')
                start = int(start) + 1  # bed is 0-based, gff 1-based
                outfile.write(
                    '{:}\tRepeatMasker\tdispersed_repeat\t{:}\t{:}\t.\t+\t.\tID={:}\n'.format(contig, start, end, name))


def findUTRs(cds, mrna, strand):
    import numpy
    FiveUTR = []
    ThreeUTR = []
    if cds != mrna:
        inter = InterLap()
        inter.add(cds)
        for i, x in enumerate(mrna):
            if not x in inter:
                loc = (list(inter)[0][0], list(inter)[-1][1])
                diff = numpy.subtract(x, loc)
                if diff[0] < 0 and diff[1] < 0:
                    if strand == '+':
                        FiveUTR.append(x)
                    else:
                        ThreeUTR.append(x)
                elif diff[0] > 0 and diff[1] > 0:
                    if strand == '+':
                        ThreeUTR.append(x)
                    else:
                        FiveUTR.append(x)
            else:
                hit = list(inter.find(x))
                if x == hit[0]:
                    continue
                else:
                    diff = numpy.subtract(x, hit[0])
                    if strand == '+':
                        if int(diff[0]) < 1 and int(diff[1]) == 0:
                            FiveUTR.append((x[0], hit[0][0]-1))
                        elif int(diff[1]) > 1 and int(diff[0]) == 0:
                            ThreeUTR.append((hit[0][1]+1, x[1]))
                        elif int(diff[0]) < 1 and int(diff[1]) > 1:
                            FiveUTR.append((x[0], hit[0][0]-1))
                            ThreeUTR.append((hit[0][1]+1, x[1]))
                    else:
                        if diff[0] == 0 and diff[1] > 0:
                            FiveUTR.append((hit[0][1]+1, x[1]))
                        elif diff[0] < 0 and diff[1] == 0:
                            ThreeUTR.append((x[0], hit[0][0]-1))
                        elif diff[0] < 0 and diff[1] > 0:
                            FiveUTR.append((hit[0][1]+1, x[1]))
                            ThreeUTR.append((x[0], hit[0][0]-1))
    return FiveUTR, ThreeUTR


def dict2nucleotides2(input, prots, trans, cdstrans):
    '''
    function to generate protein and transcripts from dictionary
    '''
    # write to protein and transcripts
    with open(prots, 'w') as protout:
        with open(trans, 'w') as tranout:
            with open(cdstrans, 'w') as cdsout:
                for k, v in natsorted(list(input.items())):
                    if 'pseudo' in v:
                        if v['pseudo']:
                            continue
                    if v['type'] == 'mRNA' and not v['CDS']:
                        continue
                    if v['type'] == 'mRNA' and not len(v['ids']) == len(v['mRNA']) == len(v['CDS']):
                        continue
                    for i, x in enumerate(v['ids']):
                        try:
                            Transcript = str(v['transcript'][i])
                            if v['strand'] == '-':
                                Transcript = RevComp(Transcript)
                            tranout.write('>{:} {:}\n{:}\n'.format(
                                x, k, softwrap(Transcript)))
                        except IndexError:
                            pass
                        try:
                            CDStranscript = str(v['cds_transcript'][i])
                            if v['strand'] == '-':
                                CDStranscript = RevComp(CDStranscript)
                            cdsout.write('>{:} {:}\n{:}\n'.format(
                                x, k, softwrap(CDStranscript)))
                        except IndexError:
                            pass
                        if v['type'] == 'mRNA':
                            try:
                                Prot = v['protein'][i]
                            except IndexError:
                                print(('ERROR', k, v))
                                sys.exit(1)
                            if Prot.endswith('*'):
                                Prot = Prot[:-1]
                            protout.write('>{:} {:}\n{:}\n'.format(
                                x, k, softwrap(Prot)))

def simpleFastaStats(fasta):
    from Bio.SeqUtils import GC
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    contigs = []
    with open(fasta, 'r') as infile:
        for header, seq in SimpleFastaParser(infile):
            contigs.append(seq)
    contigs = sorted(contigs, key=lambda x: len(x), reverse=True)
    lengths = [len(x) for x in contigs]
    pctGC = round(GC(''.join(contigs)), 2)
    #N50
    totalLen = sum(lengths)
    n50_len = totalLen*.50
    n90_len = totalLen*.90
    n50 = None
    n90 = None
    runningSum = 0
    for y, c in enumerate(lengths):
        runningSum += c
        if not n50 and runningSum >= n50_len:
            n50 = c
        if not n90 and runningSum >= n90_len:
            n90 = c
    l50 = lengths.index(n50) + 1
    l90 = lengths.index(n90) + 1
    numContigs = len(contigs)
    avg_length = '{:.2f}'.format(totalLen / float(len(contigs)))
    return numContigs, totalLen, pctGC, float(avg_length), n50, l50, n90, l90


def databases2json(FUNDB):
    resources = {}
    dbfile = os.path.join(FUNDB, 'funannotate-db-info.txt')
    if not os.path.isfile(dbfile):
        return resources
    else:
        with open(dbfile, 'r') as infile:
            for line in infile:
                line = line.rstrip()
                cols = line.split('\t')
                resources[cols[0]] = {
                    'type': cols[1],
                    'version': cols[3],
                    'date': cols[4],
                    'num-records': cols[5]
                }
        return resources


def annotation_summary(fasta, output, gff=False, tbl=False, pct=0.90,
                       transcripts=False, proteins=False, previous=False,
                       database='.', command='', organism=''):
    '''
    function to output annotation stats from GFF3 or TBL files
    '''
    import json
    stats = {
        'format': 'annotation',
        'command': command,
        'organism': organism,
        'software':{
            'name': 'funannotate',
            'version': get_version(),
            'date': datetime.datetime.today().strftime('%Y-%m-%d'),
            'resources': databases2json(database)
        },
        'assembly': {
            'num_contigs': 0,
            'length': 0,
            'mean_length': 0,
            'N50': 0,
            'L50': 0,
            'N90': 0,
            'L90': 0,
            'GC_content': 0,
        },
        'annotation': {
            'genes': 0,
            'common_name': 0,
            'mRNA': 0,
            'tRNA': 0,
            'ncRNA': 0,
            'rRNA': 0,
            'avg_gene_length': 0,
            'transcript-level': {
                'CDS_transcripts': 0,
                'CDS_five_utr': 0,
                'CDS_three_utr': 0,
                'CDS_no_utr': 0,
                'CDS_five_three_utr': 0,
                'CDS_complete': 0,
                'CDS_no-start': 0,
                'CDS_no-stop': 0,
                'CDS_no-start_no-stop': 0,
                'total_exons': 0,
                'total_cds_exons': 0,
                'multiple_exon_transcript': 0,
                'single_exon_transcript': 0,
                'avg_exon_length': 0,
                'avg_protein_length': 0,
                'functional': {
                    'go_terms': 0,
                    'interproscan': 0,
                    'eggnog': 0,
                    'pfam': 0,
                    'cazyme': 0,
                    'merops': 0,
                    'busco': 0,
                    'secretion': 0
                }
            }
        }

    }
    if previous:  # load some stats that cant calculate from annotation
        with open(previous, 'r') as infile:
            previousStats = json.load(infile)
        try:
            stats['annotation']['transcript-level']['pct_exon_overlap_protein_evidence'] = previousStats['annotation']['transcript-level']['pct_exon_overlap_protein_evidence']
        except KeyError:
            pass
        try:
            stats['annotation']['transcript-level']['pct_exon_overlap_transcript_evidence'] = previousStats['annotation']['transcript-level']['pct_exon_overlap_transcript_evidence']
        except KeyError:
            pass
    num, tot, gc, avg, n50, l50, n90, l90 = simpleFastaStats(fasta)
    stats['assembly']['num_contigs'] = num
    stats['assembly']['length'] = tot
    stats['assembly']['GC_content'] = gc
    stats['assembly']['mean_length'] = avg
    stats['assembly']['N50'] = n50
    stats['assembly']['L50'] = l50
    stats['assembly']['N90'] = n90
    stats['assembly']['L90'] = l90
    Genes = {}
    if tbl:
        Genes = tbl2dict(tbl, fasta, Genes)
    elif gff:
        Genes = gff2dict(gff, fasta, Genes)
    if len(Genes) > 0:
        protLengths = []
        geneLengths = []
        exonLengths = []
        for k, v in Genes.items():
            stats['annotation']['genes'] += 1
            gLength = v['location'][1] - v['location'][0]
            geneLengths.append(gLength)
            if v['type'] == 'tRNA':
                stats['annotation']['tRNA'] += 1
            elif v['type'] == 'rRNA':
                stats['annotation']['rRNA'] += 1
            elif v['type'] == 'ncRNA':
                stats['annotation']['ncRNA'] += 1
            if v['name']:
                stats['annotation']['common_name'] += 1
            for i in range(0, len(v['ids'])):
                if v['type'] == 'mRNA':
                    stats['annotation']['mRNA'] += 1
                    stats['annotation']['transcript-level']['CDS_transcripts'] += 1
                    pLen = len(v['protein'][i])
                    if v['protein'][i].endswith('*'):
                        pLen -= 1
                    protLengths.append(pLen)
                    if len(v['mRNA'][i]) > 1:
                        stats['annotation']['transcript-level']['multiple_exon_transcript'] += 1
                        for y in v['mRNA'][i]:
                            exon_length = y[1] - y[0]
                            exonLengths.append(exon_length)
                    else:
                        stats['annotation']['transcript-level']['single_exon_transcript'] += 1
                    stats['annotation']['transcript-level']['total_exons'] += len(v['mRNA'][i])
                    stats['annotation']['transcript-level']['total_exons'] += len(v['5UTR'][i])
                    stats['annotation']['transcript-level']['total_exons'] += len(v['3UTR'][i])
                    stats['annotation']['transcript-level']['total_cds_exons'] += len(v['CDS'][i])
                    if v['partialStart'][i] and v['partialStop'][i]:
                        stats['annotation']['transcript-level']['CDS_no-start_no-stop'] += 1
                    elif v['partialStart'][i]:
                        stats['annotation']['transcript-level']['CDS_no-start'] += 1
                    elif v['partialStop'][i]:
                        stats['annotation']['transcript-level']['CDS_no-stop'] += 1
                    else:
                        stats['annotation']['transcript-level']['CDS_complete'] += 1
                    if len(v['5UTR'][i]) > 0 and len(v['3UTR'][i]) > 0:
                        stats['annotation']['transcript-level']['CDS_five_three_utr'] += 1
                    elif len(v['3UTR'][i]) > 0:
                        stats['annotation']['transcript-level']['CDS_three_utr'] += 1
                    elif len(v['5UTR'][i]) > 0:
                        stats['annotation']['transcript-level']['CDS_three_utr'] += 1
                    else:
                        stats['annotation']['transcript-level']['CDS_no_utr'] += 1
                    if v['go_terms'][i]:
                        stats['annotation']['transcript-level']['functional']['go_terms'] += 1
                    if any(s.startswith('PFAM:') for s in v['db_xref'][i]):
                        stats['annotation']['transcript-level']['functional']['pfam'] += 1
                    if any(s.startswith('InterPro:') for s in v['db_xref'][i]):
                        stats['annotation']['transcript-level']['functional']['interproscan'] += 1
                    if any(s.startswith('EggNog:') for s in v['note'][i]):
                        stats['annotation']['transcript-level']['functional']['eggnog'] += 1
                    if any(s.startswith('CAZy:') for s in v['note'][i]):
                        stats['annotation']['transcript-level']['functional']['cazyme'] += 1
                    if any(s.startswith('MEROPS:') for s in v['note'][i]):
                        stats['annotation']['transcript-level']['functional']['merops'] += 1
                    if any(s.startswith('BUSCO:') for s in v['note'][i]):
                        stats['annotation']['transcript-level']['functional']['busco'] += 1
                    if any(s.startswith('SECRETED:') for s in v['note'][i]):
                        stats['annotation']['transcript-level']['functional']['secretion'] += 1
        stats['annotation']['avg_gene_length'] = round(sum(geneLengths) / float(len(geneLengths)), 2)
        stats['annotation']['transcript-level']['avg_protein_length'] = round(sum(protLengths) / float(len(protLengths)), 2)
        stats['annotation']['transcript-level']['avg_exon_length'] = round(sum(exonLengths) / float(len(exonLengths)), 2)
        exonBED = 'tmp.exon.{}.bed'.format(os.getpid())
        if transcripts or proteins:
            exonCount = 0
            bedtools_cmd = ['bedtools', 'intersect', '-a', exonBED,
                            '-u', '-f', str(pct), '-s', '-b']
            with open(exonBED, 'w') as outfile:
                for k, v in Genes.items():
                    for i in range(0, len(v['ids'])):
                        for z, x in enumerate(v['mRNA'][i]):
                            exonCount += 1
                            outfile.write('{}\t{}\t{}\t{}.exon{}\t.\t{}\n'.format(
                                v['contig'], x[0]-1, x[1], v['ids'][i], z+1, v['strand']))
        if transcripts:  # calculate exons covered by transcripts
            cmd = bedtools_cmd + [transcripts]
            overlapCount = 0
            for line in execute(cmd):
                overlapCount += 1
            pctOverlap = '{:.2f}'.format(overlapCount/exonCount*100)
            stats['annotation']['transcript-level']['pct_exon_overlap_transcript_evidence'] = float(pctOverlap)
        if proteins:  # calculate exons covered by proteins
            cmd = bedtools_cmd + [proteins]
            overlapCount = 0
            for line in execute(cmd):
                overlapCount += 1
            pctOverlap = '{:.2f}'.format(overlapCount/exonCount*100)
            stats['annotation']['transcript-level']['pct_exon_overlap_protein_evidence'] = float(pctOverlap)
        if os.path.isfile(exonBED):
            os.remove(exonBED)
    # write to json format
    with open(output, 'w') as outfile:
        json.dump(stats, outfile, indent=4)


def tbl2allout(input, fasta, GFF, Proteins, Transcripts, cdsTranscripts, DNA):
    '''
    function to convert NCBI tbl format directly to other formats; this will be a replacement
    for Genbank derived output files and correctly parse/print the transcript/proteins
    '''
    Genes = {}
    Genes = tbl2dict(input, fasta, Genes)
    # write GFF
    dict2gff3(Genes, GFF)
    # write to protein and transcripts
    dict2nucleotides2(Genes, Proteins, Transcripts, cdsTranscripts)
    # copy over DNA fasta file
    shutil.copyfile(fasta, DNA)


def tbl2dict(input, fasta, Genes):
    '''
    need a method to convert directly from NCBI tbl format to several output formats
    to avoid conversion problems with GBK files that have mutliple transcripts
    if can load funannotate dictionary directly from tbl format, then can write the other
    formats directly
    '''
    with open(input, 'r') as infile:
        contig = ''
        for item in readBlocks2(infile, '>Feature', '\tgene\n'):
            if item[0].startswith('>Feature'):  # this will be contig header block
                contig = item[0].rstrip().split(' ')[-1]
            else:  # these are all gene model blocks
                geneID, Name, type, start, end, fivepartial, threepartial, strand, location = (
                    None,)*9
                codon_start = []
                transcriptID = []
                proteinID = []
                synonyms = []
                product = []
                first, firstpartial, second, secondpartial = (False,)*4
                position = None
                # check number of transcripts
                tNum = 0
                for z in item:
                    if z.startswith('\t\t\ttranscript_id'):
                        tNum += 1
                if tNum > 0:
                    tNum = int(tNum / 2)
                else:
                    tNum = 1
                # setup lists for transcripts
                mRNA = [[] for y in range(tNum)]
                CDS = [[] for y in range(tNum)]
                note = [[] for y in range(tNum)]
                dbxref = [[] for y in range(tNum)]
                ECnum = [[] for y in range(tNum)]
                go_terms = [[] for y in range(tNum)]
                fivepartial = [False, ]*tNum
                threepartial = [False, ]*tNum
                currentNum = 0
                for x in item:
                    exonF, exonR, cdsF, cdsR, cols = (None,)*5
                    if x.endswith('\tgene\n') and not position:
                        cols = x.strip().split('\t')
                        position = 'gene'
                        if cols[0].startswith('<'):
                            first = int(cols[0].split('<')[-1])
                        else:
                            first = int(cols[0])
                        if cols[1].startswith('>'):
                            second = int(cols[1].split('>')[-1])
                        else:
                            second = int(cols[1])
                        if first < second:
                            start = first
                            end = second
                            strand = '+'
                        else:
                            start = second
                            end = first
                            strand = '-'
                        location = (start, end)
                    elif x.startswith('\t\t\tgene\t'):
                        Name = x.strip().split('\t')[-1]
                    elif x.startswith('\t\t\tlocus_tag\t'):
                        geneID = x.strip().split('\t')[-1]
                    elif x.endswith('\ttRNA\n') and x.count('\t') == 2 and position == 'gene':
                        type = 'tRNA'
                        position = 'tRNA'
                        cols = x.strip().split('\t')
                        exonF = int(cols[0].replace('<', ''))
                        exonR = int(cols[1].replace('>', ''))
                        if strand == '+':
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif x.endswith('\tncRNA\n') and x.count('\t') == 2 and position == 'gene':
                        type = 'ncRNA'
                        position = 'ncRNA'
                        cols = x.strip().split('\t')
                        exonF = int(cols[0].replace('<', ''))
                        exonR = int(cols[1].replace('>', ''))
                        if strand == '+':
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif x.endswith('\trRNA\n') and x.count('\t') == 2 and position == 'gene':
                        type = 'rRNA'
                        position = 'rRNA'
                        cols = x.strip().split('\t')
                        exonF = int(cols[0].replace('<', ''))
                        exonR = int(cols[1].replace('>', ''))
                        if strand == '+':
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif x.endswith('\tmRNA\n') and x.count('\t') == 2:
                        if position == 'CDS':
                            currentNum += 1
                        elif position == 'gene':
                            type = 'mRNA'
                        position = 'mRNA'
                        cols = x.strip().split('\t')
                        exonF = int(cols[0].replace('<', ''))
                        exonR = int(cols[1].replace('>', ''))
                        if strand == '+':
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif x.endswith('\tCDS\n') and x.count('\t') == 2:
                        position = 'CDS'
                        cols = x.strip().split('\t')
                        cdsF = int(cols[0].replace('<', ''))
                        cdsR = int(cols[1].replace('>', ''))
                        if strand == '+':
                            CDS[currentNum].append((cdsF, cdsR))
                        else:
                            CDS[currentNum].append((cdsR, cdsF))
                    elif x.startswith('\t\t\tcodon_start\t'):
                        cNum = int(x.strip().split('\t')[-1])
                        codon_start.append(cNum)
                    elif x.startswith('\t\t\tproduct\t') and position != 'mRNA':
                        product.append(x.strip().split('\t')[-1])
                    elif x.startswith('\t\t\ttranscript_id\t'):
                        tID = x.strip().split('|')[-1]
                        if '_mrna' in tID:
                            tID = tID.replace('_mrna', '')
                        if not tID in transcriptID:
                            transcriptID.append(tID)
                    elif x.startswith('\t\t\tprotein_id\t'):
                        pID = x.strip().split('|')[-1]
                        if not pID in proteinID:
                            proteinID.append(pID)
                    elif x.startswith('\t\t\tgene_synonym\t'):
                        synonyms.append(x.strip().split('\t')[-1])
                    elif x.startswith('\t\t\tgo_'):  # go terms
                        go_terms[currentNum].append(
                            'GO:{:}'.format(x.strip().split('|')[1]))
                    elif x.startswith('\t\t\tnote\t'):
                        note[currentNum].append(x.strip().split('\t')[-1])
                    elif x.startswith('\t\t\tdb_xref\t'):
                        dbxref[currentNum].append(x.strip().split('\t')[-1])
                    elif x.startswith('\t\t\tEC_number\t'):
                        ECnum[currentNum].append(x.strip().split('\t')[-1])
                    elif position == 'mRNA' and x.count('\t') == 1:
                        cols = x.strip().split('\t')
                        exonF = int(cols[0].replace('<', ''))
                        exonR = int(cols[1].replace('>', ''))
                        if strand == '+':
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif position in ['tRNA', 'ncRNA', 'rRNA'] and x.count('\t') == 1:
                        cols = x.strip().split('\t')
                        exonF = int(cols[0].replace('<', ''))
                        exonR = int(cols[1].replace('>', ''))
                        if strand == '+':
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif position == 'CDS' and x.count('\t') == 1:
                        cols = x.strip().split('\t')
                        cdsF = int(cols[0].replace('<', ''))
                        cdsR = int(cols[1].replace('>', ''))
                        if strand == '+':
                            CDS[currentNum].append((cdsF, cdsR))
                        else:
                            CDS[currentNum].append((cdsR, cdsF))
                if not geneID in Genes:
                    if type in ['tRNA', 'ncRNA', 'rRNA']:
                        Genes[geneID] = {'name': Name, 'type': type,
                                         'transcript': [],
                                         'cds_transcript': [],
                                         'protein': [], '5UTR': [[]],
                                         '3UTR': [[]],
                                         'codon_start': codon_start,
                                         'ids': [geneID+'-T1'], 'CDS': CDS,
                                         'mRNA': mRNA, 'strand': strand,
                                         'gene_synonym': synonyms,
                                         'location': location,
                                         'contig': contig,
                                         'product': product,
                                         'source': 'funannotate', 'phase': [],
                                         'db_xref': dbxref,
                                         'go_terms': go_terms,
                                         'EC_number': ECnum, 'note': note,
                                         'partialStart': [True],
                                         'partialStop': [True],
                                         'pseudo': False
                                         }
                    else:
                        Genes[geneID] = {'name': Name, 'type': type,
                                         'transcript': [], 'cds_transcript': [],
                                         'protein': [], '5UTR': [], '3UTR': [],
                                         'codon_start': codon_start,
                                         'ids': proteinID, 'CDS': CDS,
                                         'mRNA': mRNA, 'strand': strand,
                                         'gene_synonym': synonyms,
                                         'location': location,
                                         'contig': contig, 'product': product,
                                         'source': 'funannotate', 'phase': [],
                                         'db_xref': dbxref,
                                         'go_terms': go_terms,
                                         'EC_number': ECnum, 'note': note,
                                         'partialStart': fivepartial,
                                         'partialStop': threepartial,
                                         'pseudo': False
                                         }
    # now we need to sort coordinates, get protein/transcript sequences and capture UTRs
    SeqRecords = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    for k, v in list(Genes.items()):
        for i in range(0, len(v['ids'])):
            if v['type'] in ['mRNA', 'tRNA', 'ncRNA']:
                if v['strand'] == '+':
                    sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0])
                else:
                    sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0],
                                         reverse=True)
                Genes[k]['mRNA'][i] = sortedExons
                mrnaSeq = getSeqRegions(SeqRecords, v['contig'], sortedExons)
                Genes[k]['transcript'].append(mrnaSeq)
            if v['type'] == 'mRNA':
                if v['strand'] == '+':
                    sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0])
                else:
                    sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0],
                                       reverse=True)
                # get the codon_start by getting first CDS phase + 1
                cdsSeq = getSeqRegions(SeqRecords, v['contig'], sortedCDS)
                Genes[k]['cds_transcript'].append(cdsSeq)
                Genes[k]['CDS'][i] = sortedCDS
                protSeq, codon_start = (None,)*2
                protSeq = translate(cdsSeq, v['strand'], v['codon_start'][i]-1)
                if protSeq:
                    Genes[k]['protein'].append(protSeq)
                    if protSeq.endswith('*'):
                        Genes[k]['partialStop'][i] = False
                    else:
                        Genes[k]['partialStop'][i] = True
                    if v['codon_start'][i] == 1 and protSeq.startswith('M'):
                        Genes[k]['partialStart'][i] = False
                    else:
                        Genes[k]['partialStart'][i] = True
                # get UTRs
                try:
                    FiveUTR, ThreeUTR = findUTRs(sortedCDS, sortedExons,
                                                 v['strand'])
                    Genes[k]['5UTR'].append(FiveUTR)
                    Genes[k]['3UTR'].append(ThreeUTR)
                except ValueError:
                    print(('ERROR', k, v))
    return Genes


def dicts2tbl(genesDict, scaff2genes, scaffLen, SeqCenter, SeqRefNum, skipList,
              output, annotations=False, external=False):
    '''
    function to take funannotate annotation dictionaries and convert to NCBI tbl output
    '''
    duplicates = 0
    pseudo = 0
    nocds = 0
    # to parse annotations, will need to have access to GO OBO dictionary
    goDict = {}
    if annotations:
        from goatools import obo_parser
        # location of go.obo
        for item in obo_parser.OBOReader(os.path.join(os.environ["FUNANNOTATE_DB"], 'go.obo')):
            goDict[item.id] = {'name': item.name, 'namespace': item.namespace}

    def _goFormat(id, goDict=goDict):
        # go_function    serine-type endopeptidase activity|0004252||IEA
        # go_process proteolysis|0006508||IEA
        # go_component   nucleus|0005634||IEA
        if id in goDict:
            if goDict[id]['namespace'] == 'biological_process':
                base = 'go_process'
            elif goDict[id]['namespace'] == 'molecular_function':
                base = 'go_function'
            elif goDict[id]['namespace'] == 'cellular_component':
                base = 'go_component'
            reformatted = '\t\t\t{:}\t{:}|{:}||IEA'.format(
                base, goDict[id]['name'], id.replace('GO:', ''))
            return reformatted
        else:
            return False

    with open(output, 'w') as tbl:
        for k, v in natsorted(list(scaff2genes.items())):
            tbl.write('>Feature %s\n' % k)
            tbl.write('1\t%s\tREFERENCE\n' % scaffLen.get(k))
            tbl.write('\t\t\t%s\t%s\n' % (SeqCenter, SeqRefNum))
            for genes in v:  # now loop through each gene on the scaffold
                if genes in skipList:
                    continue
                # single funannotate standard dictionary
                geneInfo = genesDict.get(genes)
                if 'pseudo' in geneInfo:
                    if geneInfo['pseudo']:
                        try:
                            log.debug('{:} is pseudo, skipping'.format(genes))
                        except NameError:
                            print(('{:} is pseudo, skipping'.format(genes)))
                        pseudo += 1
                        continue
                if geneInfo['type'] == 'mRNA' and not geneInfo['CDS']:
                    try:
                        log.debug(
                            'Skipping {:} because no CDS found.'.format(genes))
                    except NameError:
                        print((
                            'Skipping {:} because no CDS found.'.format(genes)))
                    pseudo += 1
                    continue
                if geneInfo['type'] == 'mRNA' and not len(geneInfo['ids']) == len(geneInfo['mRNA']) == len(geneInfo['CDS']):
                    try:
                        log.debug('Incompatible annotation found: {:}\n{:}'.format(
                            genes, geneInfo))
                    except NameError:
                        print(('Incompatible annotation found: {:}\n{:}'.format(
                            genes, geneInfo)))
                    duplicates += 1
                    continue
                if geneInfo['type'] == 'mRNA' and len(geneInfo['CDS']) == 0:
                    nocds += 1
                    continue
                if geneInfo['type'] is None:
                    continue
                # check for partial models
                if True in geneInfo['partialStart']:
                    ps = '<'
                else:
                    ps = ''
                if True in geneInfo['partialStop']:
                    pss = '>'
                else:
                    pss = ''
                # now write gene model
                if geneInfo['strand'] == '+':
                    tbl.write('%s%i\t%s%i\tgene\n' % (
                        ps, geneInfo['location'][0], pss, geneInfo['location'][1]))
                    if annotations:
                        if geneInfo['name']:
                            tbl.write('\t\t\tgene\t%s\n' % geneInfo['name'])
                        if geneInfo['gene_synonym']:
                            for alias in geneInfo['gene_synonym']:
                                tbl.write('\t\t\tgene_synonym\t%s\n' % alias)
                    tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                else:
                    tbl.write('%s%i\t%s%i\tgene\n' % (
                        ps, geneInfo['location'][1], pss, geneInfo['location'][0]))
                    if annotations:
                        if geneInfo['name']:
                            tbl.write('\t\t\tgene\t%s\n' % geneInfo['name'])
                        if geneInfo['gene_synonym']:
                            for alias in geneInfo['gene_synonym']:
                                tbl.write('\t\t\tgene_synonym\t%s\n' % alias)
                    tbl.write('\t\t\tlocus_tag\t%s\n' % genes)

                # now will output the gene models with -T1, -T2, -T3 annotations based on expression values
                # means need to get the order
                order = []
                # multiple transcripts, so get order of highest TPM
                if len(geneInfo['ids']) > 1:
                    tpms = []
                    for num, tpm in enumerate(geneInfo['note']):
                        for item in tpm:
                            if item.startswith('TPM:'):
                                value = float(item.split(':')[-1])
                                tpms.append((value, num))
                    if len(tpms) > 0:
                        for x in sorted(tpms, reverse=True):
                            order.append(x[1])
                    else:
                        order = list(range(0, len(geneInfo['ids'])))
                else:
                    order.append(0)
                for num, i in enumerate(order):  # now write mRNA and CDS features
                    # if geneInfo['ids'][i].startswith('evm.model'): #if from predict, rename to match locus_tag
                    #    protein_id = genes+'-T'+str(num+1)
                    # else:
                    #    protein_id = geneInfo['ids'][i]
                    if external:
                        protein_id = geneInfo['ids'][i]
                    else:
                        protein_id = genes+'-T'+str(num+1)
                    if geneInfo['type'] == 'mRNA':
                        if geneInfo['partialStart'][i] is False:
                            ps = ''
                        else:
                            ps = '<'
                        if geneInfo['partialStop'][i] is False:
                            pss = ''
                        else:
                            pss = '>'
                        if geneInfo['strand'] == '+':
                            for num, exon in enumerate(geneInfo['mRNA'][i]):
                                # single exon, so slightly differnt method
                                if num == 0 and num == len(geneInfo['mRNA'][i]) - 1:
                                    tbl.write('%s%s\t%s%s\tmRNA\n' %
                                              (ps, exon[0], pss, exon[1]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tmRNA\n' %
                                              (ps, exon[0], exon[1]))
                                # this is last one
                                elif num == len(geneInfo['mRNA'][i]) - 1:
                                    tbl.write('%s\t%s%s\n' %
                                              (exon[0], pss, exon[1]))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                            tbl.write('\t\t\tproduct\t%s\n' %
                                      geneInfo['product'][i])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n' % (protein_id))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s\n' %
                                      (protein_id))
                            for num, cds in enumerate(geneInfo['CDS'][i]):
                                # single exon, so slightly differnt method
                                if num == 0 and num == len(geneInfo['CDS'][i]) - 1:
                                    tbl.write('%s%s\t%s%s\tCDS\n' %
                                              (ps, cds[0], pss, cds[1]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tCDS\n' %
                                              (ps, cds[0], cds[1]))
                                # this is last one
                                elif num == len(geneInfo['CDS'][i]) - 1:
                                    tbl.write('%s\t%s%s\n' %
                                              (cds[0], pss, cds[1]))
                                else:
                                    tbl.write('%s\t%s\n' % (cds[0], cds[1]))
                            tbl.write('\t\t\tcodon_start\t%i\n' %
                                      geneInfo['codon_start'][i])
                            if annotations:  # write functional annotation
                                if geneInfo['EC_number'][i]:
                                    for EC in geneInfo['EC_number'][i]:
                                        tbl.write('\t\t\tEC_number\t%s\n' % EC)
                                if geneInfo['db_xref'][i]:
                                    for xref in geneInfo['db_xref'][i]:
                                        tbl.write('\t\t\tdb_xref\t%s\n' % xref)
                                if geneInfo['go_terms'][i]:
                                    for go in geneInfo['go_terms'][i]:
                                        goLine = _goFormat(go)
                                        if goLine:
                                            tbl.write('{:}\n'.format(goLine))
                                if geneInfo['note'][i]:
                                    for item in geneInfo['note'][i]:
                                        tbl.write('\t\t\tnote\t%s\n' % item)
                            tbl.write('\t\t\tproduct\t%s\n' %
                                      geneInfo['product'][i])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n' % (protein_id))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s\n' %
                                      (protein_id))
                        else:  # means this is on crick strand
                            for num, exon in enumerate(geneInfo['mRNA'][i]):
                                # single exon, so slightly differnt method
                                if num == 0 and num == len(geneInfo['mRNA'][i]) - 1:
                                    tbl.write('%s%s\t%s%s\tmRNA\n' %
                                              (ps, exon[1], pss, exon[0]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tmRNA\n' %
                                              (ps, exon[1], exon[0]))
                                # this is last one
                                elif num == len(geneInfo['mRNA'][i]) - 1:
                                    tbl.write('%s\t%s%s\n' %
                                              (exon[1], pss, exon[0]))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[1], exon[0]))
                            tbl.write('\t\t\tproduct\t%s\n' %
                                      geneInfo['product'][i])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n' % (protein_id))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s\n' %
                                      (protein_id))
                            for num, cds in enumerate(geneInfo['CDS'][i]):
                                # single exon, so slightly differnt method
                                if num == 0 and num == len(geneInfo['CDS'][i]) - 1:
                                    tbl.write('%s%s\t%s%s\tCDS\n' %
                                              (ps, cds[1], pss, cds[0]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tCDS\n' %
                                              (ps, cds[1], cds[0]))
                                # this is last one
                                elif num == (len(geneInfo['CDS'][i]) - 1):
                                    tbl.write('%s\t%s%s\n' %
                                              (cds[1], pss, cds[0]))
                                else:
                                    tbl.write('%s\t%s\n' % (cds[1], cds[0]))
                            tbl.write('\t\t\tcodon_start\t%i\n' %
                                      geneInfo['codon_start'][i])
                            if annotations:  # write functional annotation
                                if geneInfo['EC_number'][i]:
                                    for EC in geneInfo['EC_number'][i]:
                                        tbl.write('\t\t\tEC_number\t%s\n' % EC)
                                if geneInfo['db_xref'][i]:
                                    for xref in geneInfo['db_xref'][i]:
                                        tbl.write('\t\t\tdb_xref\t%s\n' % xref)
                                if geneInfo['go_terms'][i]:
                                    for go in geneInfo['go_terms'][i]:
                                        goLine = _goFormat(go)
                                        if goLine:
                                            tbl.write('{:}\n'.format(goLine))
                                if geneInfo['note'][i]:
                                    for item in geneInfo['note'][i]:
                                        tbl.write('\t\t\tnote\t%s\n' % item)
                            tbl.write('\t\t\tproduct\t%s\n' %
                                      geneInfo['product'][i])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n' % (protein_id))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s\n' %
                                      (protein_id))
                    elif geneInfo['type'] == 'tRNA':
                        if geneInfo['strand'] == '+':
                            for num, exon in enumerate(geneInfo['mRNA'][i]):
                                if num == 0:
                                    tbl.write('%s\t%s\t%s\n' % (
                                        exon[0], exon[1], geneInfo['type']))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                            tbl.write('\t\t\tproduct\t%s\n' %
                                      geneInfo['product'][i])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tbl.write('\t\t\tpseudo\n')
                        else:
                            for num, exon in enumerate(geneInfo['mRNA'][i]):
                                if num == 0:
                                    tbl.write('%s\t%s\t%s\n' % (
                                        exon[1], exon[0], geneInfo['type']))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[1], exon[0]))
                            tbl.write('\t\t\tproduct\t%s\n' %
                                      geneInfo['product'][i])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tbl.write('\t\t\tpseudo\n')
                    elif geneInfo['type'] in ['rRNA', 'ncRNA']:
                        if geneInfo['strand'] == '+':
                            tbl.write('%s\t%s\t%s\n' % (
                                geneInfo['location'][0], geneInfo['location'][1], geneInfo['type']))
                            tbl.write('\t\t\tproduct\t%s\n' %
                                      geneInfo['product'][i])
                        else:
                            tbl.write('%s\t%s\t%s\n' % (
                                geneInfo['location'][1], geneInfo['location'][0], geneInfo['type']))
                            tbl.write('\t\t\tproduct\t%s\n' %
                                      geneInfo['product'][i])
    if any(i > 0 for i in [duplicates, pseudo, nocds]):
        try:
            print(('Skipped {:,} annotations: {:,} pseudo genes; {:,} no CDS; {:,} duplicated features'.format(
                sum([pseudo, nocds, duplicates]), pseudo, nocds, duplicates)))
        except NameError:
            print(('Skipped {:,} annotations: {:,} pseudo genes; {:,} no CDS; {:,} duplicated features'.format(
                sum([pseudo, nocds, duplicates]), pseudo, nocds, duplicates)))


def GFF2tbl(evm, trnascan, fasta, scaffLen, prefix, Numbering, SeqCenter,
            SeqRefNum, tblout):
    from collections import OrderedDict
    '''
    function to take EVM protein models and tRNA scan GFF to produce a GBK tbl file as well
    as a new GFF3 file. The function will also rename locus_id if passed.
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][1])

    # load GFF into dictionary
    Genes = {}
    Genes = gff2dict(evm, fasta, Genes)
    Genes = gff2dict(trnascan, fasta, Genes)

    # now sort dictionary by contig and location, rename using prefix, translate to protein space to get proper start/stop info
    sGenes = natsorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    scaff2genes = {}
    count = Numbering
    for k, v in list(sortedGenes.items()):
        if prefix:
            locusTag = prefix+'_'+str(count).zfill(6)
        else:
            locusTag = k
        renamedGenes[locusTag] = v
        if not v['contig'] in scaff2genes:
            scaff2genes[v['contig']] = [locusTag]
        else:
            scaff2genes[v['contig']].append(locusTag)
        count += 1

    # write tbl outputfile
    dicts2tbl(renamedGenes, scaff2genes, scaffLen,
              SeqCenter, SeqRefNum, [], tblout)


def checkRefSeq(input):
    refseq = False
    with open(input, 'r') as infile:
        for record in SeqIO.parse(infile, 'genbank'):
            if 'RefSeq' in record.annotations['keywords']:
                refseq = True
            break
    return refseq


def getGBKinfo(input):
    accession = None
    organism = None
    strain = None
    isolate = None
    gb_gi = None
    WGS_accession = None
    version = None
    with open(input, 'r') as infile:
        for record in SeqIO.parse(infile, 'genbank'):
            try:
                WGS_accession = 'WGS:' + \
                    record.annotations['contig'].split(
                        ':')[0].replace('join(', '')[:4]
            except KeyError:
                pass
            try:
                accession = record.annotations['accessions'][0]
            except KeyError:
                pass
            try:
                organism = record.annotations['organism'].replace(
                    'Unclassified.', '').rstrip()
            except KeyError:
                pass
            try:
                gb_gi = record.annotations['gi']
            except KeyError:
                pass
            try:
                version = record.annotations['sequence_version']
            except KeyError:
                pass
            for f in record.features:
                if f.type == "source":
                    isolate = f.qualifiers.get("isolate", [None])[0]
                    strain = f.qualifiers.get("strain", [None])[0]
            break
    return organism, strain, isolate, accession, WGS_accession, gb_gi, version


def getGBKLocusTag(input):
    LocusTags = []
    with open(input, 'r') as infile:
        for record in SeqIO.parse(infile, 'genbank'):
            for f in record.features:
                if f.type == 'gene':
                    ID = f.qualifiers['locus_tag'][0]
                    if not ID in LocusTags:
                        LocusTags.append(ID)
    lastTag = natsorted(LocusTags)[-1]
    if not '_' in lastTag:
        print('ERROR: underscore "_" not found in locus_tag, exiting.')
        sys.exit(1)
    tag, count = lastTag.rsplit('_', 1)
    justify = len(count)
    return tag, count, justify


def gb2dna(input, output):
    with open(output, 'w') as outfile:
        with open(input, 'r') as infile:
            for record in SeqIO.parse(infile, 'genbank'):
                outfile.write(">%s\n%s\n" %
                              (record.id, softwrap(str(record.seq))))


def getID(input, type):
    # function to get ID from genbank record.features
    locusTag = None
    ID = None
    Parent = None
    if type == 'gene':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
        else:
            try:
                ID = input.qualifiers['gene'][0]
            except KeyError:
                pass
        return locusTag, ID, locusTag

    elif type in ['mRNA', 'tRNA', 'ncRNA', 'rRNA', 'misc_RNA', 'exon']:
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
            if locusTag:
                Parent = locusTag
                try:
                    ID = input.qualifiers['transcript_id'][0]
                except KeyError:
                    pass
            else:
                try:
                    locusTag = input.qualifiers['transcript_id'][0]
                    Parent = locusTag
                except KeyError:
                    pass
        else:
            try:
                ID = input.qualifiers['transcript_id'][0]
            except KeyError:
                pass
        if ID:
            if ':' in ID:
                ID = ID.split(':')[-1]
        else:
            try:
                ID = input.qualifiers['standard_name'][0]
            except KeyError:
                pass
        return locusTag, ID, Parent

    elif type == 'CDS':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
            if locusTag:
                Parent = locusTag
                try:
                    ID = input.qualifiers['protein_id'][0]
                except KeyError:
                    pass
            else:
                try:
                    locusTag = input.qualifiers['protein_id'][0]
                    Parent = locusTag
                except KeyError:
                    pass
        else:
            try:
                ID = input.qualifiers['protein_id'][0]
            except KeyError:
                ID = locusTag
        if ID:
            if ':' in ID:
                ID = ID.split(':')[-1]
        else:
            try:
                ID = input.qualifiers['standard_name'][0]
            except KeyError:
                pass
        return locusTag, ID, Parent


def gb2nucleotides(input, prots, trans, dna):
    '''
    function to generate protein, transcripts, and contigs from genbank file
    '''
    genes = {}
    with open(dna, 'w') as dnaout:
        with open(input, 'r') as filein:
            for record in SeqIO.parse(filein, 'genbank'):
                dnaout.write(">%s\n%s\n" %
                             (record.id, softwrap(str(record.seq))))
                for f in record.features:
                    gb_feature_add2dict(f, record, genes)
    # write to protein and transcripts
    dict2nucleotides(genes, prots, trans)
    return len(genes)


def dict2proteins(input, prots):
    with open(prots, 'w') as protout:
        for k, v in natsorted(list(input.items())):
            if 'pseudo' in v:
                if v['pseudo']:
                    continue
            if v['type'] == 'mRNA' and not v['CDS']:
                continue
            if v['type'] == 'mRNA' and not len(v['ids']) == len(v['mRNA']) == len(v['CDS']):
                continue
            for i, x in enumerate(v['ids']):
                if v['type'] == 'mRNA':
                    Prot = v['protein'][i]
                    protout.write('>{:} {:}\n{:}\n'.format(
                        x, k, softwrap(Prot)))


def dict2nucleotides(input, prots, trans):
    '''
    function to generate protein and transcripts from dictionary
    '''
    # write to protein and transcripts
    with open(prots, 'w') as protout:
        with open(trans, 'w') as tranout:
            for k, v in natsorted(list(input.items())):
                if 'pseudo' in v:
                    if v['pseudo']:
                        continue
                if v['type'] == 'mRNA' and not v['CDS']:
                    continue
                if v['type'] == 'mRNA' and not len(v['ids']) == len(v['mRNA']) == len(v['CDS']):
                    continue
                for i, x in enumerate(v['ids']):
                    try:
                        Transcript = str(v['transcript'][i])
                        tranout.write('>{:} {:}\n{:}\n'.format(
                            x, k, softwrap(Transcript)))
                    except IndexError:
                        pass
                    if v['type'] == 'mRNA':
                        Prot = v['protein'][i]
                        protout.write('>{:} {:}\n{:}\n'.format(
                            x, k, softwrap(Prot)))


def gb2gffnuc(input, gff, prots, trans, dna):
    '''
    function to generate protein, transcripts, and contigs from genbank file
    '''
    genes = {}
    with open(dna, 'w') as dnaout:
        with open(input, 'r') as filein:
            for record in SeqIO.parse(filein, 'genbank'):
                dnaout.write(">{:}\n{:}\n".format(
                    record.id, softwrap(str(record.seq))))
                for f in record.features:
                    gb_feature_add2dict(f, record, genes)
    # write gff3 output
    dict2gff3(genes, gff)
    # write to protein and transcripts
    dict2nucleotides(genes, prots, trans)
    return len(genes)


def gb2parts(input, tbl, gff, prots, trans, dna):
    '''
    function returns a dictionary of all gene models from a genbank file this function
    can handle multiple transcripts per locus/gene
    '''
    genes = {}
    scaff2genes = {}
    scaffLen = {}
    with open(dna, 'w') as dnaout:
        with open(input, 'r') as filein:
            for record in SeqIO.parse(filein, 'genbank'):
                dnaout.write(">{:}\n{:}\n".format(
                    record.id, softwrap(str(record.seq))))
                Contig = record.id
                if not Contig in scaffLen:
                    scaffLen[Contig] = len(record.seq)
                for f in record.features:
                    if f.type == 'gene':
                        locusTag, ID, Parent = getID(f, f.type)
                        if not Contig in scaff2genes:
                            scaff2genes[Contig] = [locusTag]
                        else:
                            scaff2genes[Contig].append(locusTag)
                    gb_feature_add2dict(f, record, genes)

    # write tbl output
    dicts2tbl(genes, scaff2genes, scaffLen, 'CFMR', '12345', [], tbl)
    # write gff3 output
    dict2gff3_old(genes, gff)
    # write to protein and transcripts
    dict2nucleotides(genes, prots, trans)
    return len(genes)


def gb_feature_add2dict(f, record, genes):
    '''
    general function to take a genbank feature from flat file and add to funannotate standardized dictionary
    locustag: {
    'contig': contigName
    'type': mRNA/rRNA/tRNA/ncRNA
    'location': (start, end) #integer tuple
    'strand': +/-
    'ids': [transcript/protein IDs] #list
    'mRNA':[[(ex1,ex1),(ex2,ex2)]] #list of lists of tuples (start, end)
    'CDS':[[(cds1,cds1),(cds2,cds2)]] #list of lists of tuples (start, end)
    'transcript': [seq1, seq2] #list of mRNA trnascripts
    'cds_transcript': [seq1, seq2] list of mRNA (no UTRs)
    'protein': [protseq1,protseq2] #list of CDS translations
    'protein_id': [id,id] #from NCBI
    'codon_start': [1,1] #codon start for translations
    'note': [[first note, second note], [first, second, etc]] #list of lists
    'name': genename
    'product': [hypothetical protein, velvet complex] #list of product definitions
    'go_terms': [[GO:0000001,GO:0000002]] #list of lists
    'db_xref': [[InterPro:IPR0001,PFAM:004384]] #list of lists
    'partialStart': True/False
    'partialStop': True/False
    'source': annotation source
    'pseudo': True/False
    }
    '''
    # get info from features, if there is no locusTag then exit
    if f.type and f.type in ['gene', 'mRNA', 'CDS', 'tRNA', 'rRNA', 'ncRNA', 'exon', 'misc_RNA']:
        try:
            locusTag, ID, Parent = getID(f, f.type)
        except TypeError:
            print('ERROR parsing GBK record')
            print(f)
            sys.exit(1)
        if not locusTag:
            return genes
    else:
        return genes
    # check for mismatching funannotate ID locus tag basename
    if ID and '-T' in ID:  # then this is from funannotate, okay to modify - this is to capture apparent tbl2asn local error
        # there is a problem, update locusTag with basename of ID
        if ID.split('-T')[0] != locusTag:
            locusTag = ID.split('-T')[0]
    # standard information from every feature
    strand = f.location.strand
    if strand == 1:
        strand = '+'
    elif strand == -1:
        strand = '-'
    start = f.location.nofuzzy_start + 1
    end = f.location.nofuzzy_end
    chr = record.id
    num_parts = len(f.location.parts)
    name, Product = (None,)*2
    Fivepartial, Threepartial = (False,)*2
    DBxref = []
    Note = []
    GO = []
    EC = []
    synonyms = []
    pseudo = False
    if 'pseudo' in f.qualifiers:
        pseudo = True
    # parse each type somewhat differently
    if f.type == 'gene':
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        if 'gene_synonym' in f.qualifiers:
            for z in f.qualifiers['gene_synonym']:
                synonyms.append(z)
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': None, 'transcript': [],
                               'cds_transcript': [], 'protein': [],
                               'source': 'GenBank', '5UTR': [[]], '3UTR': [[]],
                               'codon_start': [], 'ids': [], 'CDS': [],
                               'mRNA': [], 'strand': strand,
                               'location': (int(start), int(end)),
                               'contig': chr, 'product': [],
                               'gene_synonym': synonyms, 'EC_number': [],
                               'db_xref': [], 'go_terms': [], 'note': [],
                               'partialStart': [], 'partialStop': [],
                               'protein_id': [], 'pseudo': pseudo}
        else:
            genes[locusTag]['location'] = (int(start), int(end))
            genes[locusTag]['strand'] = strand
            genes[locusTag]['gene_synonym'] = synonyms
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    elif f.type in ['tRNA', 'rRNA', 'ncRNA', 'misc_RNA']:
        feature_seq = f.extract(record.seq)
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        try:
            Product = f.qualifiers['product'][0]
            if Product == 'tRNA-OTHER':
                Product = 'tRNA-Xxx'
        except KeyError:
            Product = None
        exonTuples = []
        if num_parts < 2:  # only single exon
            exonTuples.append((int(start), int(end)))
        else:  # more than 1 exon, so loop through
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                exonTuples.append((int(ex_start), int(ex_end)))
        # now we want to sort the positions I think...
        if strand == '+':
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
            if str(f.location.start).startswith('<'):
                Fivepartial = True
            if str(f.location.end).startswith('>'):
                Threepartial = True
        else:
            sortedExons = sorted(
                exonTuples, key=lambda tup: tup[0], reverse=True)
            if str(f.location.start).startswith('<'):
                Threepartial = True
            if str(f.location.end).startswith('>'):
                Fivepartial = True
        # update positions
        if f.type == 'misc_RNA':
            feature = 'ncRNA'
        else:
            feature = f.type
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': feature,
                               'transcript': [feature_seq],
                               'cds_transcript': [], 'protein': [],
                               'source': 'GenBank', '5UTR': [[]], '3UTR': [[]],
                               'codon_start': [], 'ids': [locusTag+'-T1'],
                               'CDS': [], 'mRNA': [sortedExons],
                               'strand': strand,
                               'location': (int(start), int(end)),
                               'contig': chr, 'product': [Product],
                               'protein_id': [], 'pseudo': pseudo,
                               'gene_synonym': synonyms, 'EC_number': [EC],
                               'db_xref': [DBxref], 'go_terms': [GO],
                               'note': [Note], 'partialStart': [Fivepartial],
                               'partialStop': [Threepartial]}
        else:
            genes[locusTag]['mRNA'].append(sortedExons)
            genes[locusTag]['type'] = feature
            genes[locusTag]['transcript'].append(feature_seq)
            genes[locusTag]['cds_transcript'].append(None)
            genes[locusTag]['protein'].append(None)
            genes[locusTag]['ids'].append(
                locusTag+'-T'+str(len(genes[locusTag]['ids'])+1))
            genes[locusTag]['db_xref'].append(DBxref)
            genes[locusTag]['note'].append(Note)
            genes[locusTag]['go_terms'].append(GO)
            genes[locusTag]['EC_number'].append(EC)
            genes[locusTag]['product'].append(Product)
            genes[locusTag]['partialStart'].append(Fivepartial)
            genes[locusTag]['partialStop'].append(Threepartial)
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    elif f.type == 'mRNA':
        feature_seq = f.extract(record.seq)
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        exonTuples = []
        if num_parts < 2:  # only single exon
            exonTuples.append((int(start), int(end)))
        else:  # more than 1 exon, so loop through
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                exonTuples.append((int(ex_start), int(ex_end)))
        # now we want to sort the positions I think...
        if strand == '+':
            sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
            if str(f.location.start).startswith('<'):
                Fivepartial = True
            if str(f.location.end).startswith('>'):
                Threepartial = True
        else:
            sortedExons = sorted(
                exonTuples, key=lambda tup: tup[0], reverse=True)
            if str(f.location.start).startswith('<'):
                Threepartial = True
            if str(f.location.end).startswith('>'):
                Fivepartial = True
        # update positions
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': f.type,
                               'transcript': [feature_seq],
                               'cds_transcript': [], 'protein': [],
                               'source': 'GenBank', '5UTR': [[]], '3UTR': [[]],
                               'codon_start': [], 'ids': [], 'CDS': [],
                               'mRNA': [sortedExons], 'strand': strand,
                               'location': (int(start), int(end)),
                               'contig': chr, 'product': [], 'protein_id': [],
                               'pseudo': pseudo, 'gene_synonym': synonyms,
                               'EC_number': [],
                               'db_xref': [], 'go_terms': [],
                               'note': [], 'partialStart': [Fivepartial],
                               'partialStop': [Threepartial]}
        else:
            genes[locusTag]['mRNA'].append(sortedExons)
            genes[locusTag]['type'] = f.type
            genes[locusTag]['transcript'].append(feature_seq)
            genes[locusTag]['partialStart'].append(Fivepartial)
            genes[locusTag]['partialStop'].append(Threepartial)
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    elif f.type == 'exon':  # assuming need to overwrite mRNA feature then?
        if len(genes[locusTag]['mRNA']) == 0:
            genes[locusTag]['mRNA'] = []
            genes[locusTag]['transcript'] = []
            feature_seq = f.extract(record.seq)
            try:
                name = f.qualifiers['gene'][0]
            except KeyError:
                pass
            exonTuples = []
            if num_parts < 2:  # only single exon
                exonTuples.append((int(start), int(end)))
            else:  # more than 1 exon, so loop through
                for i in range(0, num_parts):
                    ex_start = f.location.parts[i].nofuzzy_start + 1
                    ex_end = f.location.parts[i].nofuzzy_end
                    exonTuples.append((int(ex_start), int(ex_end)))
            # now we want to sort the positions I think...
            if strand == '+':
                sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
                if str(f.location.start).startswith('<'):
                    Fivepartial = True
                if str(f.location.end).startswith('>'):
                    Threepartial = True
            else:
                sortedExons = sorted(
                    exonTuples, key=lambda tup: tup[0], reverse=True)
                if str(f.location.start).startswith('<'):
                    Threepartial = True
                if str(f.location.end).startswith('>'):
                    Fivepartial = True
            # update positions
            if not locusTag in genes:
                genes[locusTag] = {'name': name, 'type': f.type,
                                'transcript': [feature_seq],
                                'cds_transcript': [], 'protein': [],
                                'source': 'GenBank', '5UTR': [[]], '3UTR': [[]],
                                'codon_start': [], 'ids': [], 'CDS': [],
                                'mRNA': [sortedExons], 'strand': strand,
                                'location': (int(start), int(end)),
                                'contig': chr, 'product': [], 'protein_id': [],
                                'db_xref': [], 'go_terms': [], 'note': [],
                                'gene_synonym': synonyms, 'EC_number': [],
                                'partialStart': [Fivepartial],
                                'partialStop': [Threepartial], 'pseudo': pseudo}
            else:
                genes[locusTag]['mRNA'].append(sortedExons)
                genes[locusTag]['transcript'].append(feature_seq)
                genes[locusTag]['partialStart'].append(Fivepartial)
                genes[locusTag]['partialStop'].append(Threepartial)
    elif f.type == 'CDS' and 'codon_start' in f.qualifiers:
        feature_seq = f.extract(record.seq)
        if not ID:
            try:
                log.info("putative transcript from %s has no ID\n(%s %s %s)" % (
                    locusTag, locusTag, ID, Parent))
            except NameError:
                print(("putative transcript from %s has no ID\n(%s %s %s)" %
                      (locusTag, locusTag, ID, Parent)))
            return genes
        try:
            protSeq = f.qualifiers['translation'][0]
        except KeyError:
            try:
                log.debug("%s has no translation" % ID)
            except NameError:
                print(("%s has no translation" % ID))
            protSeq = ''
        cdsTuples = []
        phase = int(f.qualifiers['codon_start'][0])
        if num_parts < 2:  # only single CDS
            cdsTuples.append((int(start), int(end)))
        else:
            for i in range(0, num_parts):
                ex_start = f.location.parts[i].nofuzzy_start + 1
                ex_end = f.location.parts[i].nofuzzy_end
                cdsTuples.append((int(ex_start), int(ex_end)))
        if strand == '+':
            sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0])
        else:
            sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0], reverse=True)
        # check for annotations
        try:
            Product = f.qualifiers['product'][0]
        except KeyError:
            Product = 'hypothetical protein'
        try:
            name = f.qualifiers['gene'][0]
        except KeyError:
            pass
        # note and dbxref are in a dictionary
        for key, value in list(f.qualifiers.items()):
            if key == 'note':
                notes = value[0].split('; ')
                for n in notes:
                    if n.startswith('GO'):
                        GO.append(n)
                    else:
                        Note.append(n)
            elif key == 'db_xref':
                for ref in value:
                    DBxref.append(ref)
            elif key == 'EC_number':
                for x in value:
                    EC.append(x)
        # update dictionary
        if not locusTag in genes:
            genes[locusTag] = {'name': name, 'type': 'mRNA',
                               'transcript': [], '5UTR': [[]], '3UTR': [[]],
                               'cds_transcript': [feature_seq],
                               'protein': [], 'source': 'GenBank',
                               'codon_start': [phase], 'ids': [locusTag+'-T1'],
                               'CDS': [sortedCDS], 'mRNA': [],
                               'strand': strand,
                               'location': (int(start), int(end)),
                               'contig': chr, 'product': [Product],
                               'gene_synonym': synonyms, 'EC_number': [EC],
                               'protein_id': [ID],
                               'db_xref': [DBxref], 'go_terms': [GO],
                               'note': [Note], 'partialStart': [],
                               'partialStop': [], 'pseudo': pseudo}
        else:
            genes[locusTag]['protein_id'].append(ID)
            genes[locusTag]['ids'].append(
                locusTag+'-T'+str(len(genes[locusTag]['ids'])+1))
            genes[locusTag]['CDS'].append(sortedCDS)
            genes[locusTag]['5UTR'].append([])
            genes[locusTag]['3UTR'].append([])
            genes[locusTag]['product'].append(Product)
            genes[locusTag]['protein'].append(protSeq)
            genes[locusTag]['cds_transcript'].append(feature_seq)
            genes[locusTag]['codon_start'].append(phase)
            genes[locusTag]['db_xref'].append(DBxref)
            genes[locusTag]['note'].append(Note)
            genes[locusTag]['go_terms'].append(GO)
            genes[locusTag]['EC_number'].append(EC)
            if not genes[locusTag]['type']:
                genes[locusTag]['type'] = 'mRNA'
            if not genes[locusTag]['name']:
                genes[locusTag]['name'] = name
    return genes


def bed2interlapNames(bedfile):
    # load interlap object from a bed file
    inter = defaultdict(InterLap)
    with open(bedfile, 'r') as infile:
        for line in infile:
            line = line.strip()
            chr, start, end, name = line.split('\t')[:4]
            inter[chr].add((int(start), int(end), name))
    return inter


def bed2interlap(bedfile):
    # load interlap object from a bed file
    inter = defaultdict(InterLap)
    with open(bedfile, 'r') as infile:
        for line in infile:
            line = line.strip()
            chr, start, end = line.split('\t')[:3]
            inter[chr].add((int(start), int(end)))
    return inter


def interlapIntersect(coords, contig, interObj):
    # return interlap coords of an intersection
    if coords in interObj[contig]:
        return True
    else:
        return False


def gff2interlap(input, fasta):
    '''
    function to parse GFF3 file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    '''
    inter = defaultdict(InterLap)
    Genes = {}
    Genes = gff2dict(input, fasta, Genes)
    for k, v in natsorted(list(Genes.items())):
        inter[v['contig']].add((v['location'][0], v['location'][1], k))
    return inter, Genes


def gff2interlapDict(input, fasta, inter, Dict):
    '''
    function to parse GFF3 file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    '''
    Genes = {}
    Genes = gff2dict(input, fasta, Genes, gap_filter=True)
    for k, v in natsorted(list(Genes.items())):
        inter[v['contig']].add(
            (v['location'][0], v['location'][1], v['strand'], k))
    # merge dictionary and return
    Dict = merge_dicts(Dict, Genes)
    return inter, Dict


def merge_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def exonerate2hints(file, outfile):
    # mimic exonerate2hints from GFF3 exonerate file
    # CDSpart +/- 15 bp to each match
    # intron as is
    '''
    #gff3 via EVM
    scaffold_20 exonerate   nucleotide_to_protein_match 225035  225823  82.13   +   .   ID=match.11677.2;Target=VC83_07547 1 96
    scaffold_20 exonerate   nucleotide_to_protein_match 53957   54342   92.93   +   .   ID=match.11677.3;Target=VC83_02595 1 129
    scaffold_20 exonerate   nucleotide_to_protein_match 54397   54904   92.93   +   .   ID=match.11677.3;Target=VC83_02595 130 299
    scaffold_107    exonerate   nucleotide_to_protein_match 77634   78119   89.95   -   .   ID=match.11677.5;Target=VC83_08471 1 163
    scaffold_107    exonerate   nucleotide_to_protein_match 77501   77546   89.95   -   .   ID=match.11677.5;Target=VC83_08471 163 178
    scaffold_107    exonerate   nucleotide_to_protein_match 77385   77422   89.95   -   .   ID=match.11677.5;Target=VC83_08471 179 191

    #corresponding exonerate2hints
    scaffold_20 xnt2h   CDSpart 225050  225808  .   +   .   src=XNT;grp=VC83_07547;pri=4
    scaffold_20 xnt2h   CDSpart 53972   54327   .   +   .   src=XNT;grp=VC83_02595;pri=4
    scaffold_20 xnt2h   intron  54343   54396   .   +   .   src=XNT;grp=VC83_02595;pri=4
    scaffold_20 xnt2h   CDSpart 54412   54889   .   +   .   src=XNT;grp=VC83_02595;pri=4
    scaffold_107    xnt2h   CDSpart 77649   78104   .   -   .   src=XNT;grp=VC83_08471;pri=4
    scaffold_107    xnt2h   intron  77547   77633   .   -   .   src=XNT;grp=VC83_08471;pri=4
    scaffold_107    xnt2h   CDSpart 77516   77531   .   -   .   src=XNT;grp=VC83_08471;pri=4
    scaffold_107    xnt2h   intron  77423   77500   .   -   .   src=XNT;grp=VC83_08471;pri=4
    scaffold_107    xnt2h   CDSpart 77400   77407   .   -   .   src=XNT;grp=VC83_08471;pri=4

    '''
    Genes = {}
    with open(file, 'r') as input:
        for line in input:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split(
                '\t')
            start = int(start)
            end = int(end)
            ID, Target = (None,)*2
            info = attributes.split(';')
            for x in info:
                if x.startswith('ID='):
                    ID = x.replace('ID=', '')
                elif x.startswith('Target='):
                    Target = x.replace('Target=', '').split(' ')[0]
            if not ID in Genes:
                Genes[ID] = {'id': ID, 'target': Target, 'loc': [
                    (start, end)], 'strand': strand, 'contig': contig}
            else:
                Genes[ID]['loc'].append((start, end))
    # now lets sort through and write hints file
    with open(outfile, 'w') as output:
        for k, v in natsorted(list(Genes.items())):
            if v['strand'] == '+':
                sortedCDS = sorted(v['loc'], key=lambda tup: tup[0])
                for i, x in enumerate(sortedCDS):  # loop through tuples
                    output.write('{:}\txnt2h\tCDSpart\t{:}\t{:}\t.\t{:}\t.\tsrc=XNT;grp={:};pri=4\n'.format(
                        v['contig'], x[0]-15, x[1]+15, v['strand'], v['target']))
                    if len(sortedCDS) > 1:
                        try:
                            output.write('{:}\txnt2h\tintron\t{:}\t{:}\t.\t{:}\t.\tsrc=XNT;grp={:};pri=4\n'.format(
                                v['contig'], x[1]+1, sortedCDS[i+1][0]-1, v['strand'], v['target']))
                        except IndexError:
                            pass
            else:
                sortedCDS = sorted(
                    v['loc'], key=lambda tup: tup[0], reverse=True)
                for i, x in enumerate(sortedCDS):  # loop through tuples
                    output.write('{:}\txnt2h\tCDSpart\t{:}\t{:}\t.\t{:}\t.\tsrc=XNT;grp={:};pri=4\n'.format(
                        v['contig'], x[0]+15, x[1]-15, v['strand'], v['target']))
                    if len(sortedCDS) > 1:
                        try:
                            output.write('{:}\txnt2h\tintron\t{:}\t{:}\t.\t{:}\t.\tsrc=XNT;grp={:};pri=4\n'.format(
                                v['contig'], sortedCDS[i+1][1]+1, x[0]-1, v['strand'], v['target']))
                        except IndexError:
                            pass


def alignments2dict(input, Genes):
    '''
    function to take a transcript_alignments file and create dictionary
    structure for each alignment
    '''
    with open(input, 'r') as infile:
        for line in infile:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split(
                '\t')
            start = int(start)
            end = int(end)
            ID, Target, Extra = (None,)*3
            for x in attributes.split(';'):
                if x.startswith('ID='):
                    ID = x.replace('ID=', '')
                elif x.startswith('Target='):
                    Target, Extra = x.split(' ', 1)
                    Target = Target.replace('Target=', '')
            if not ID:
                continue
            if not ID in Genes:
                Genes[ID] = {'mRNA': [(start, end)], 'strand': strand, 'pident': [score],
                             'location': (start, end), 'contig': contig, 'extra': [Extra]}
            else:
                if contig != Genes[ID]['contig']:
                    log.debug('ERROR: {:} mapped to multiple contigs: {:} and {:}'.format(ID, contig, Genes[ID]['contig']))
                    continue
                elif strand != Genes[ID]['strand']:
                    log.debug('ERROR: {:} mapped has different strands'.format(ID))
                    continue
                else:
                    Genes[ID]['mRNA'].append((start, end))
                    Genes[ID]['pident'].append(score)
                    Genes[ID]['extra'].append(Extra)
                    # double check mRNA features are contained in gene coordinates
                    if start < Genes[ID]['location'][0]:
                        Genes[ID]['location'] = (
                            start, Genes[ID]['location'][1])
                    if end > Genes[ID]['location'][1]:
                        Genes[ID]['location'] = (
                            Genes[ID]['location'][0], end)
    return Genes

def introns_from_exons(input):
    introns = []
    if len(input) > 1:
        for x, y in enumerate(input):
            try:
                introns.append((y[1]+1, input[x+1][0]-1))
            except IndexError:
                pass
    return introns

def dict2hints(input, hints):
    from collections import OrderedDict
    '''
    function to take simple alignments dictionary and ouput augustus hints file
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = natsorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    with open(hints, 'w') as hintsout:
        for k, v in list(sortedGenes.items()):
            sortedExons = sorted(v['mRNA'], key=lambda tup: tup[0])
            introns = introns_from_exons(sortedExons)
            for i, exon in enumerate(sortedExons):
                if i == 0 or i == len(sortedExons)-1:
                    hintsout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tgrp={:};pri=4;src=E\n'.format(
                            v['contig'], 'b2h', 'ep', exon[0], exon[1], 0, v['strand'], '.', k))
                else:
                    hintsout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tgrp={:};pri=4;src=E\n'.format(
                            v['contig'], 'b2h', 'exon', exon[0], exon[1], 0, v['strand'], '.', k))
            if len(introns) > 0:
                for z in introns:
                    hintsout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tgrp={:};pri=4;src=E\n'.format(
                            v['contig'], 'b2h', 'intron', z[0], z[1], 1, v['strand'], '.', k))


def dict2transcriptgff3(input, output):
    from collections import OrderedDict
    '''
    function to take simple alignments dictionary and ouput GFF3 transcripts file
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = natsorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    with open(output, 'w') as outfile:
        outfile.write('##gff-version 3\n')
        for k, v in list(sortedGenes.items()):
            for i, exon in enumerate(v['mRNA']):
                outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:};Target={:} {:}\n'.format(
                    v['contig'], 'genome', 'cDNA_match', exon[0], exon[1], v['pident'][i], v['strand'], '.',
                    k, k, v['extra'][i]))


def harmonize_transcripts(genome, alignments, gfffile, hintsfile, evidence=None, tmpdir='.', cpus=1, maxintron=3000):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    '''
    function to check if evidence transcripts are missing from existing alignments and/or
    write the augustus hints file
    '''
    Genes = {}
    Genes = alignments2dict(alignments, Genes)
    log.info('Parsed {:,} transcript alignments from: {:}'.format(len(Genes), alignments))
    if evidence:  # if nothing here then just move on
        uniqueTranscripts = os.path.join(tmpdir, 'transcript_evidence_unique.fasta')
        seqcount = 0
        with open(uniqueTranscripts, 'w') as fasta_outfile:
            for file in evidence:
                with open(file, 'r') as fasta_infile:
                    for title, seq in SimpleFastaParser(fasta_infile):
                        if ' ' in title:
                            id = title.split(' ')[0]
                        else:
                            id = title
                        if not id in Genes:
                            fasta_outfile.write('>{:}\n{:}\n'.format(title, softwrap(seq)))
                            seqcount += 1
        if seqcount > 0:
            log.info('Aligning {:,} unique transcripts [not found in exising alignments] with minimap2'.format(seqcount))
            minimapBAM = os.path.join(tmpdir, 'transcript_evidence_unique.bam')
            minimapGFF = os.path.join(tmpdir, 'transcript_evidence_unique.gff3')
            minimap2Align(uniqueTranscripts, genome, cpus, maxintron, minimapBAM)
            mappedReads = bam2gff3(str(minimapBAM), minimapGFF)
            if mappedReads > 0:
                log.info('Mapped {:,} of these transcripts to the genome, adding to alignments'.format(mappedReads))
                Genes = alignments2dict(minimapGFF, Genes)
            else:
                log.info('Mapped 0 of these transcripts to the genome')
    log.info('Creating transcript EVM alignments and Augustus transcripts hintsfile')
    dict2transcriptgff3(Genes, gfffile)
    dict2hints(Genes, hintsfile)


def gff2dict(file, fasta, Genes, debug=False, gap_filter=False):
    '''
    general function to take a GFF3 file and return a funannotate standardized dictionary
    locustag: {
    'contig': contigName
    'type': mRNA/rRNA/tRNA/ncRNA
    'location': (start, end) #integer tuple
    'strand': +/-
    'ids': [transcript/protein IDs] #list
    'mRNA':[[(ex1,ex1),(ex2,ex2)]] #list of lists of tuples (start, end)
    'CDS':[[(cds1,cds1),(cds2,cds2)]] #list of lists of tuples (start, end)
    'transcript': [seq1, seq2] #list of mRNA trnascripts
    'cds_transcript': [seq1, seq2] #list of mRNA trnascripts (no UTRs)
    'protein': [protseq1,protseq2] #list of CDS translations
    'codon_start': [1,1] #codon start for translations
    'note': [[first note, second note], [first, second, etc]] #list of lists
    'name': genename
    'product': [hypothetical protein, velvet complex] #list of product definitions
    'gene_synonym': Aliases
    'EC_number': [[ec number]]
    'go_terms': [[GO:0000001,GO:0000002]] #list of lists
    'db_xref': [[InterPro:IPR0001,PFAM:004384]] #list of lists
    'partialStart': True/False
    'partialStop': True/False
    'source': annotation source
    'phase': [[0,2,1]] list of lists
    '5UTR': [[(),()]] #list of lists of tuples (start, end)
    '3UTR': [[(),()]] #list of lists of tuples (start, end)
    }
    '''
    idParent = {}
    SeqRecords = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    with open(file, 'r') as input:
        for line in input:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            if feature not in ['gene', 'mRNA', 'exon', 'CDS', 'tRNA',
                               'ncRNA', 'rRNA', 'pseudogene', 'five_prime_UTR',
                               'five_prime_utr', 'three_prime_UTR',
                               'three_prime_utr', 'transcript']:
                continue
            if not contig in SeqRecords:
                continue
            start = int(start)
            end = int(end)
            ID, Parent, Name, Product, GeneFeature, gbkey = (None,)*6
            Note, DBxref, GO, synonyms, ECnum = ([],)*5
            info = attributes.split(';')
            for x in info:
                if x.startswith('ID='):
                    ID = x.replace('ID=', '')
                elif x.startswith('Parent='):
                    Parent = x.replace('Parent=', '')
                elif x.startswith('Name='):
                    Name = x.replace('Name=', '')
                elif x.startswith('Note=') or x.startswith('note='):
                    Note = x.split('ote=')[-1]
                    if ',' in Note:
                        Note = Note.split(',')
                    else:
                        Note = [Note]
                elif x.startswith('Dbxref='):
                    DBxref = x.replace('Dbxref=', '')
                    if ',' in DBxref:
                        DBxref = DBxref.split(',')
                    else:
                        DBxref = [DBxref]
                elif x.startswith('Ontology_term='):
                    GO = x.replace('Ontology_term=', '')
                    if ',' in GO:
                        GO = GO.split(',')
                    else:
                        GO = [GO]
                elif x.startswith('EC_number='):
                    ECnum = x.split('=',1)[-1]
                    if ',' in ECnum:
                        ECnum = ECnum.split(',')
                    else:
                        ECnum = [ECnum]
                elif x.startswith('Product=') or x.startswith('product='):
                    Product = unquote(x.split('roduct=')[-1])
                elif x.startswith('description='):
                    Product = unquote(x.replace('description=', ''))
                elif x.startswith('Alias='):
                    synonyms = x.replace('Alias=', '')
                    synonyms = synonyms.split(',')
                elif x.startswith('gbkey='):  # genbank uses
                    gbkey = x.split('=', 1)[-1]
            if feature == 'gene' or feature == 'pseudogene':
                if not ID in Genes:
                    if feature == 'pseudogene':
                        pseudoFlag = True
                    else:
                        pseudoFlag = False
                    Genes[ID] = {'name': Name, 'type': None, 'transcript': [],
                                 'cds_transcript': [], 'protein': [], '5UTR': [],
                                 '3UTR': [], 'gene_synonym': synonyms,
                                 'codon_start': [], 'ids': [], 'CDS': [],
                                 'mRNA': [], 'strand': strand,
                                 'EC_number': [],
                                 'location': (start, end), 'contig': contig,
                                 'product': [], 'source': source, 'phase': [],
                                 'db_xref': [], 'go_terms': [], 'note': [],
                                 'partialStart': [], 'partialStop': [],
                                 'pseudo': pseudoFlag}
                else:
                    if start < Genes[ID]['location'][0]:
                        Genes[ID]['location'] = (
                            start, Genes[ID]['location'][1])
                    if end > Genes[ID]['location'][1]:
                        Genes[ID]['location'] = (Genes[ID]['location'][0], end)
            else:
                if not ID or not Parent:
                    sys.stderr.write("Error, can't find ID or Parent. Malformed GFF file.\n")
                    sys.stderr.write(line)
                    sys.exit(1)
                if feature in ['mRNA', 'transcript', 'tRNA', 'ncRNA', 'rRNA']:
                    if gbkey and gbkey == 'misc_RNA':
                        feature = 'ncRNA'
                    if not Product:
                        if feature in ['mRNA', 'transcript']:
                            Product = 'hypothetical protein'
                    if not Parent in Genes:
                        Genes[Parent] = {'name': Name, 'type': feature,
                                         'transcript': [], 'cds_transcript': [],
                                         'protein': [], '5UTR': [[]], '3UTR': [[]],
                                         'codon_start': [[]], 'ids': [ID],
                                         'CDS': [[]], 'mRNA': [[]], 'strand': strand,
                                         'location': (start, end), 'contig': contig,
                                         'product': [Product], 'source': source,
                                         'phase': [[]], 'gene_synonym': synonyms,
                                         'db_xref': [DBxref], 'go_terms': [GO],
                                         'EC_number': [ECnum],
                                         'note': [Note], 'partialStart': [False],
                                         'partialStop': [False], 'pseudo': False}
                    else:
                        Genes[Parent]['ids'].append(ID)
                        Genes[Parent]['mRNA'].append([])
                        Genes[Parent]['CDS'].append([])
                        Genes[Parent]['phase'].append([])
                        Genes[Parent]['5UTR'].append([])
                        Genes[Parent]['3UTR'].append([])
                        Genes[Parent]['codon_start'].append([])
                        Genes[Parent]['partialStart'].append(False)
                        Genes[Parent]['partialStop'].append(False)
                        Genes[Parent]['product'].append(Product)
                        Genes[Parent]['db_xref'].append(DBxref)
                        Genes[Parent]['EC_number'].append(ECnum)
                        Genes[Parent]['gene_synonym'] += synonyms
                        Genes[Parent]['go_terms'].append(GO)
                        Genes[Parent]['note'].append(Note)
                        Genes[Parent]['type'] = feature
                        # double check mRNA features are contained in gene coordinates
                        if start < Genes[Parent]['location'][0]:
                            # print('{:} update start: {:} to {:}'.format(Parent, Genes[Parent]['location'][0],start))
                            Genes[Parent]['location'] = (
                                start, Genes[Parent]['location'][1])
                        if end > Genes[Parent]['location'][1]:
                            # print('{:} update stop: {:} to {:}'.format(Parent, Genes[Parent]['location'][1],end))
                            Genes[Parent]['location'] = (
                                Genes[Parent]['location'][0], end)
                    if not ID in idParent:
                        idParent[ID] = Parent
                elif feature == 'exon':
                    if ',' in Parent:
                        parents = Parent.split(',')
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {'name': Name, 'type': None,
                                                      'transcript': [], 'cds_transcript': [],
                                                      'protein': [], '5UTR': [[]], '3UTR': [[]],
                                                      'codon_start': [[]], 'ids': [p],
                                                      'CDS': [], 'mRNA': [[(start, end)]], 'strand': strand,
                                                      'location': None, 'contig': contig,
                                                      'product': [], 'source': source, 'phase': [[]],
                                                      'db_xref': [], 'go_terms': [],
                                                      'EC_number': [],
                                                      'note': [], 'partialStart': [False],
                                                      'partialStop': [False], 'pseudo': False,
                                                      'gene_synonym': synonyms}
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]['ids'].index(p)
                                Genes[GeneFeature]['mRNA'][i].append(
                                    (start, end))
                elif feature == 'CDS':
                    if ',' in Parent:
                        parents = Parent.split(',')
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {'name': Name, 'type': None,
                                                      'transcript': [], 'cds_transcript': [],
                                                      'protein': [], '5UTR': [[]], '3UTR': [[]],
                                                      'codon_start': [[]], 'ids': [p],
                                                      'CDS': [[(start, end)]], 'mRNA': [], 'strand': strand,
                                                      'location': None, 'contig': contig,
                                                      'product': [], 'source': source, 'phase': [[]],
                                                      'db_xref': [], 'go_terms': [],
                                                      'EC_number': [],
                                                      'note': [], 'partialStart': [False],
                                                      'partialStop': [False], 'pseudo': False,
                                                      'gene_synonym': synonyms}
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]['ids'].index(p)
                                Genes[GeneFeature]['CDS'][i].append(
                                    (start, end))
                                # add phase
                                try:
                                    Genes[GeneFeature]['phase'][i].append(int(phase))
                                except ValueError:
                                    Genes[GeneFeature]['phase'][i].append('?')
                elif feature == 'five_prime_UTR' or feature == 'five_prime_utr':
                    if ',' in Parent:
                        parents = Parent.split(',')
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {'name': Name, 'type': None,
                                                      'transcript': [], 'cds_transcript': [],
                                                      'protein': [], '5UTR': [[(start, end)]], '3UTR': [[]],
                                                      'codon_start': [[]], 'ids': [p],
                                                      'CDS': [], 'mRNA': [[(start, end)]], 'strand': strand,
                                                      'location': None, 'contig': contig,
                                                      'product': [], 'source': source, 'phase': [[]],
                                                      'db_xref': [], 'go_terms': [],
                                                      'EC_number': [],
                                                      'note': [], 'partialStart': [False],
                                                      'partialStop': [False], 'pseudo': False,
                                                      'gene_synonym': synonyms,}
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]['ids'].index(p)
                                Genes[GeneFeature]['5UTR'][i].append(
                                    (start, end))
                elif feature == 'three_prime_UTR' or feature == 'three_prime_utr':
                    if ',' in Parent:
                        parents = Parent.split(',')
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {'name': Name, 'type': None,
                                                      'transcript': [], 'cds_transcript': [],
                                                      'protein': [], '5UTR': [[]], '3UTR': [[(start, end)]],
                                                      'codon_start': [[]], 'ids': [p],
                                                      'CDS': [], 'mRNA': [[(start, end)]], 'strand': strand,
                                                      'location': None, 'contig': contig,
                                                      'product': [], 'source': source, 'phase': [[]],
                                                      'db_xref': [], 'go_terms': [],
                                                      'EC_number': [],
                                                      'note': [], 'partialStart': [False],
                                                      'partialStop': [False], 'pseudo': False,
                                                      'gene_synonym': synonyms}
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]['ids'].index(p)
                                Genes[GeneFeature]['3UTR'][i].append(
                                    (start, end))
    # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
    for k, v in list(Genes.items()):
        for i in range(0, len(v['ids'])):
            if v['type'] in ['mRNA', 'tRNA', 'ncRNA', 'rRNA']:
                if v['strand'] == '+':
                    sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0])
                else:
                    sortedExons = sorted(
                        v['mRNA'][i], key=lambda tup: tup[0], reverse=True)
                Genes[k]['mRNA'][i] = sortedExons
                mrnaSeq = getSeqRegions(SeqRecords, v['contig'], sortedExons)
                if gap_filter:
                    mrnaSeq, Genes[k]['mRNA'][i] = start_end_gap(mrnaSeq, Genes[k]['mRNA'][i])
                v['transcript'].append(mrnaSeq)
            if v['type'] == 'mRNA':
                if not v['CDS'][i]:
                    sys.stderr.write('ERROR: ID={:} has no CDS features, removing gene model\n'.format(k))
                    del Genes[k]
                    continue
                if v['strand'] == '+':
                    sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0])
                else:
                    sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0], reverse=True)
                #get the codon_start by getting first CDS phase + 1
                indexStart = [x for x, y in enumerate(v['CDS'][i]) if y[0] == sortedCDS[0][0]]
                cdsSeq = getSeqRegions(SeqRecords, v['contig'], sortedCDS)
                if gap_filter:
                    cdsSeq, v['CDS'][i] = start_end_gap(cdsSeq, v['CDS'][i])
                Genes[k]['cds_transcript'].append(cdsSeq)
                Genes[k]['CDS'][i] = sortedCDS
                protSeq, codon_start = (None,)*2
                try:
                    currentphase = v['phase'][i]
                except IndexError:
                    pass
                if '?' in v['phase'][i]: #dont know the phase -- malformed GFF3, try to find best CDS
                    translateResults = []
                    for y in [1,2,3]:
                        protSeq = translate(cdsSeq, v['strand'], y-1)
                        if not protSeq:
                            log.debug('Translation of {:} using {:} phase failed'.format(v['ids'][i], y-1))
                            continue
                        numStops = protSeq.count('*')
                        if protSeq[-1] == '*':
                            numStops -= 1
                        translateResults.append((y, numStops, protSeq))
                    sortedResults = sorted(translateResults, key=lambda tup: tup[1])
                    codon_start = sortedResults[0][0]
                    protSeq = sortedResults[0][2]
                else:
                    try:
                        codon_start = int(v['phase'][i][indexStart[0]]) + 1
                    except IndexError:
                        pass
                    #translate and get protein sequence
                    protSeq = translate(cdsSeq, v['strand'], codon_start-1)
                Genes[k]['codon_start'][i] = codon_start
                v['protein'].append(protSeq)
                if protSeq:
                    if protSeq.endswith('*'):
                        v['partialStop'][i] = False
                    else:
                        v['partialStop'][i] = True
                    if v['codon_start'][i] == 1 and v['protein'][i].startswith('M'):
                        v['partialStart'][i] = False
                    else:
                        v['partialStart'][i] = True
        # since its possible updated the mRNA/CDS fields, double check that gene coordinates are ok
        if k not in Genes:
            continue
        all_mRNA_coords = [item for sublist in v['mRNA'] for item in sublist]
        try:
            Genes[k]['location'] = (min(all_mRNA_coords, key=lambda item: item[0])[0], max(all_mRNA_coords, key=lambda item: item[1])[1])
        except ValueError:
            continue
        # clean up any repeated synonym
        if len(v['gene_synonym']) > 1:
            uniqueSynonyms = set(v['gene_synonym'])
            Genes[k]['gene_synonym'] = list(uniqueSynonyms)
    return Genes


def start_end_gap(seq, coords):
    if seq.startswith('N'):
        oldLen = len(seq)
        seq = seq.lstrip('N')
        numLeftStripped = oldLen - len(seq)
        coords[0] = (coords[0][0]+numLeftStripped, coords[0][1])
    if seq.endswith('N'):
        oldLen = len(seq)
        seq = seq.rstrip('N')
        numRightStripped = oldLen - len(seq)
        coords[-1] = (coords[-1][0], coords[-1][1]-numRightStripped)
    return seq, coords


def simplifyGO(inputList):
    simple = []
    for x in inputList:
        if x.startswith('GO:'):
            simple.append(x.strip())
        elif ' ' in x:
            simple.append(x.split(' ')[1])
    return simple


def dict2gff3(input, output, debug=False):
    from collections import OrderedDict
    '''
    function to convert funannotate gene dictionary to gff3 output
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = natsorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # then loop through and write GFF3 format
    with open(output, 'w') as gffout:
        gffout.write("##gff-version 3\n")
        for k, v in list(sortedGenes.items()):
            if 'pseudo' in v:
                if v['pseudo']:
                    continue
            if v['type'] == 'mRNA' and not v['CDS']:
                continue
            if v['type'] == 'mRNA' and not len(v['ids']) == len(v['mRNA']) == len(v['CDS']):
                continue
            if v['type'] == 'mRNA' and len(v['CDS']) == 0:
                continue
            if v['type'] is None:
                continue
            if v['name']:
                if 'gene_synonym' in v and len(v['gene_synonym']) > 0:
                    gffout.write(
                        "{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};Alias={:};\n".format(
                            v['contig'], v['source'],v['location'][0],
                            v['location'][1], v['strand'], k, v['name'],
                            ','.join(v['gene_synonym'])))
                else:
                    gffout.write(
                        "{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(
                            v['contig'], v['source'], v['location'][0],
                            v['location'][1], v['strand'], k, v['name']))
            else:
                if 'gene_synonym' in v and len(v['gene_synonym']) > 0:
                    gffout.write(
                        "{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Alias={:};\n".format(
                            v['contig'], v['source'], v['location'][0],
                            v['location'][1], v['strand'], k,
                            ','.join(v['gene_synonym'])))
                else:
                    gffout.write(
                        "{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                            v['contig'], v['source'], v['location'][0],
                            v['location'][1], v['strand'], k))

            for i in range(0, len(v['ids'])):
                # make sure coordinates are sorted
                if v['strand'] == '+':
                    sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0])
                    sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0])
                    if '5UTR' in v and v['5UTR'][i]:
                        sortedFive = sorted(
                            v['5UTR'][i], key=lambda tup: tup[0])
                    if '3UTR' in v and v['3UTR'][i]:
                        sortedThree = sorted(
                            v['3UTR'][i], key=lambda tup: tup[0])
                else:
                    sortedExons = sorted(
                        v['mRNA'][i], key=lambda tup: tup[0], reverse=True)
                    sortedCDS = sorted(
                        v['CDS'][i], key=lambda tup: tup[0], reverse=True)
                    if '5UTR' in v and v['5UTR'][i]:
                        sortedFive = sorted(
                            v['5UTR'][i], key=lambda tup: tup[0], reverse=True)
                    if '3UTR' in v and v['3UTR'][i]:
                        sortedThree = sorted(
                            v['3UTR'][i], key=lambda tup: tup[0], reverse=True)
                # build extra annotations for each transcript if applicable
                extraAnnotations = ''
                if 'gene_synonym' in v and len(v['gene_synonym']) > 0:
                    extraAnnotations = extraAnnotations + \
                        'Alias={:};'.format(','.join(v['gene_synonym']))
                if len(v['go_terms'][i]) > 0:
                    go_annotations = simplifyGO(v['go_terms'][i])
                    extraAnnotations = extraAnnotations + \
                        'Ontology_term={:};'.format(','.join(go_annotations))
                if len(v['db_xref'][i]) > 0:
                    extraAnnotations = extraAnnotations + \
                        'Dbxref={:};'.format(','.join(v['db_xref'][i]))
                if 'EC_number' in v and len(v['EC_number'][i]) > 0:
                    extraAnnotations = extraAnnotations + \
                        'EC_number={:};'.format(','.join(v['EC_number'][i]))
                if len(v['note'][i]) > 0:
                    CleanedNote = []  # need to make sure no commas or semi-colons in these data else will cause problems in parsing GFF3 output downstream
                    for x in v['note'][i]:
                        if ';' in x:
                            x = x.replace(';', '.')
                        if ':' in x:
                            base, values = x.split(':', 1)
                            if not ',' in values:
                                CleanedNote.append(base+':'+values)
                            else:
                                for y in values.split(','):
                                    CleanedNote.append(base+':'+y)
                        else:
                            CleanedNote.append(x.replace(',', ''))

                    extraAnnotations = extraAnnotations + \
                        'note={:};'.format(','.join(CleanedNote))
                # now write mRNA feature
                gffout.write(
                    "{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};{:}\n".format(
                        v['contig'], v['source'], v['type'], v['location'][0],
                        v['location'][1], v['strand'], v['ids'][i], k,
                        v['product'][i], extraAnnotations))
                if v['type'] in ['mRNA', 'tRNA', 'ncRNA']:
                    if '5UTR' in v and v['5UTR'][i]:
                        # if 5'UTR then write those first
                        num_5utrs = len(v['5UTR'][i])
                        if num_5utrs > 0:
                            for z in range(0, num_5utrs):
                                u_num = z + 1
                                gffout.write("{:}\t{:}\tfive_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr5p{:};Parent={:};\n".format(
                                    v['contig'], v['source'], sortedFive[z][0], sortedFive[z][1], v['strand'], v['ids'][i],
                                    u_num, v['ids'][i]))
                    # write the exons
                    num_exons = len(v['mRNA'][i])
                    for x in range(0, num_exons):
                        ex_num = x + 1
                        gffout.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                                     v['contig'], v['source'], sortedExons[x][0], sortedExons[x][1], v['strand'],
                                     v['ids'][i], ex_num, v['ids'][i]))
                    # if 3'UTR then write
                    if '3UTR' in v and v['3UTR'][i]:
                        num_3utrs = len(v['3UTR'][i])
                        if num_3utrs > 0:
                            for z in range(0, num_3utrs):
                                u_num = z + 1
                                gffout.write("{:}\t{:}\tthree_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr3p{:};Parent={:};\n".format(
                                             v['contig'], v['source'], sortedThree[z][0], sortedThree[z][1], v['strand'],
                                             v['ids'][i], u_num, v['ids'][i]))
                if v['type'] == 'mRNA':
                    num_cds = len(v['CDS'][i])
                    # GFF3 phase is 1 less than flat file
                    current_phase = v['codon_start'][i] - 1
                    for y in range(0, num_cds):
                        gffout.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(
                                     v['contig'], v['source'], sortedCDS[y][0], sortedCDS[y][1], v['strand'],
                                     current_phase, v['ids'][i], v['ids'][i]))
                        current_phase = (
                            current_phase - (int(sortedCDS[y][1]) - int(sortedCDS[y][0]) + 1)) % 3
                        if current_phase == 3:
                            current_phase = 0


def dict2gff3_old(input, output):
    from collections import OrderedDict
    '''
    function to convert funannotate gene dictionary to gff3 output
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = sorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # then loop through and write GFF3 format
    with open(output, 'w') as gffout:
        gffout.write("##gff-version 3\n")
        for k, v in list(sortedGenes.items()):
            if 'pseudo' in v:
                if v['pseudo']:
                    continue
            if v['type'] == 'mRNA' and not v['CDS']:
                continue
            if v['type'] == 'mRNA' and not len(v['ids']) == len(v['mRNA']) == len(v['CDS']):
                continue
            if v['type'] == 'mRNA' and len(v['CDS']) == 0:
                continue
            if v['type'] is None:
                continue
            if v['name']:
                gffout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(
                    v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k, v['name']))
            else:
                gffout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                    v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
            for i in range(0, len(v['ids'])):
                # build extra annotations for each transcript if applicable
                extraAnnotations = ''
                if len(v['go_terms'][i]) > 0:
                    extraAnnotations = extraAnnotations + \
                        'Ontology_term={:};'.format(','.join(v['go_terms'][i]))
                if len(v['db_xref'][i]) > 0:
                    extraAnnotations = extraAnnotations + \
                        'Dbxref={:};'.format(','.join(v['db_xref'][i]))
                if len(v['note'][i]) > 0:
                    extraAnnotations = extraAnnotations + \
                        'note={:};'.format(','.join(v['note'][i]))
                # now write mRNA feature
                gffout.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};{:}\n".format(
                    v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i], extraAnnotations))
                if v['type'] == 'mRNA' or v['type'] == 'tRNA':
                    if '5UTR' in v:
                        # if 5'UTR then write those first
                        num_5utrs = len(v['5UTR'][i])
                        if num_5utrs > 0:
                            for z in range(0, num_5utrs):
                                u_num = z + 1
                                gffout.write("{:}\t{:}\tfive_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr5p{:};Parent={:};\n".format(
                                    v['contig'], v['source'], v['5UTR'][i][z][0], v['5UTR'][i][z][1], v['strand'], v['ids'][i], u_num, v['ids'][i]))
                    # write the exons
                    num_exons = len(v['mRNA'][i])
                    for x in range(0, num_exons):
                        ex_num = x + 1
                        gffout.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                            v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
                    # if 3'UTR then write
                    if '3UTR' in v:
                        num_3utrs = len(v['3UTR'][i])
                        if num_3utrs > 0:
                            for z in range(0, num_3utrs):
                                u_num = z + 1
                                gffout.write("{:}\t{:}\tthree_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr3p{:};Parent={:};\n".format(
                                    v['contig'], v['source'], v['3UTR'][i][z][0], v['3UTR'][i][z][1], v['strand'], v['ids'][i], u_num, v['ids'][i]))
                if v['type'] == 'mRNA':
                    num_cds = len(v['CDS'][i])
                    # GFF3 phase is 1 less than flat file
                    current_phase = v['codon_start'][i] - 1
                    for y in range(0, num_cds):
                        gffout.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(
                            v['contig'], v['source'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                        current_phase = (
                            current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                        if current_phase == 3:
                            current_phase = 0


def dict2gff3noUTRs(input, output):
    from collections import OrderedDict
    '''
    function to convert funannotate gene dictionary to gff3 output, no UTRs!
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = sorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # then loop through and write GFF3 format
    with open(output, 'w') as gffout:
        gffout.write("##gff-version 3\n")
        for k, v in list(sortedGenes.items()):
            if v['name']:
                gffout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(
                    v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k, v['name']))
            else:
                gffout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                    v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
            for i in range(0, len(v['ids'])):
                # build extra annotations for each transcript if applicable
                extraAnnotations = ''
                if len(v['go_terms'][i]) > 0:
                    extraAnnotations = extraAnnotations + \
                        'Ontology_term={:};'.format(','.join(v['go_terms'][i]))
                if len(v['db_xref'][i]) > 0:
                    extraAnnotations = extraAnnotations + \
                        'Dbxref={:};'.format(','.join(v['db_xref'][i]))
                if len(v['note'][i]) > 0:
                    extraAnnotations = extraAnnotations + \
                        'note={:};'.format(','.join(v['note'][i]))
                # now write mRNA feature
                gffout.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};{:}\n".format(
                    v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i], extraAnnotations))
                if v['type'] == 'tRNA':
                    # write the exons and CDS features
                    num_exons = len(v['mRNA'][i])
                    for x in range(0, num_exons):
                        ex_num = x + 1
                        gffout.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                            v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
                elif v['type'] == 'mRNA':
                    num_cds = len(v['CDS'][i])
                    # GFF3 phase is 1 less than flat file
                    current_phase = v['codon_start'][i] - 1
                    for y in range(0, num_cds):
                        ex_num = y + 1
                        gffout.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                            v['contig'], v['source'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
                        gffout.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(
                            v['contig'], v['source'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                        current_phase = (
                            current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                        if current_phase == 3:
                            current_phase = 0


def gtf2dict(input):
    Genes = {}
    with open(input, 'r') as inFile:
        for line in inFile:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.rstrip()
            # CM002242   StringTie   transcript  4198460 4199001 1000    +   .   gene_id "STRG.18087"; transcript_id "STRG.18087.2"; cov "5.905163"; FPKM "3.279455"; TPM "9.789504";
            # CM002242   StringTie   exon    4198460 4198609 1000    +   .   gene_id "STRG.18087"; transcript_id "STRG.18087.2"; exon_number "1"; cov "6.999466";
            contig, source, feature, start, end, score, strand, phase, attributes = line.split(
                '\t')
            start = int(start)
            end = int(end)
            ID, transcriptID, TPM = (None,)*3
            info = attributes.split(';')
            for x in info:
                x = x.strip()
                x = x.replace('"', '')
                if x.startswith('gene_id '):
                    ID = x.replace('gene_id ', '')
                elif x.startswith('transcript_id '):
                    transcriptID = x.replace('transcript_id ', '')
                elif x.startswith('TPM '):
                    TPM = x.replace('TPM ', '')

            if feature == 'transcript':
                if not ID in Genes:
                    Genes[ID] = {'type': 'mRNA', 'codon_start': [1], 'ids': [transcriptID], 'CDS': [[]], 'mRNA': [[]], 'strand': strand,
                                 'location': (start, end), 'contig': contig, 'source': source, 'tpm': [TPM]}
                else:
                    if start < Genes[ID]['location'][0]:
                        Genes[ID]['location'] = (
                            start, Genes[ID]['location'][1])
                    if end > Genes[ID]['location'][1]:
                        Genes[ID]['location'] = (Genes[ID]['location'][0], end)
                    Genes[ID]['ids'].append(transcriptID)
                    Genes[ID]['mRNA'].append([])
                    Genes[ID]['CDS'].append([])
                    Genes[ID]['codon_start'].append(1)
                    Genes[ID]['tpm'].append(TPM)
            else:
                if not ID or not transcriptID:
                    print(
                        "Error, can't find geneID or transcriptID. Malformed GTF file.")
                    print(line)
                    sys.exit(1)
                if feature == 'exon':
                    if not ID in Genes:
                        Genes[ID] = {'type': 'mRNA', 'codon_start': [1], 'ids': [transcriptID], 'CDS': [[(start, end)]], 'mRNA': [[(start, end)]], 'strand': strand,
                                     'location': (start, end), 'contig': contig, 'source': source, 'tpm': []}
                    else:
                        if transcriptID in Genes[ID]['ids']:  # then add exon
                            i = Genes[ID]['ids'].index(transcriptID)
                            Genes[ID]['mRNA'][i].append((start, end))
                            Genes[ID]['CDS'][i].append((start, end))
    # loop through dictionary and make sure properly sorted exons
    for k, v in list(Genes.items()):
        for i in range(0, len(v['ids'])):
            if v['strand'] == '+':
                sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0])
                sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0])
            else:
                sortedExons = sorted(
                    v['mRNA'][i], key=lambda tup: tup[0], reverse=True)
                sortedCDS = sorted(
                    v['CDS'][i], key=lambda tup: tup[0], reverse=True)
            Genes[k]['mRNA'][i] = sortedExons
            Genes[k]['CDS'][i] = sortedCDS
    return Genes


def Stringtie_dict2gff3(input, output):
    from collections import OrderedDict
    '''
    function to convert funannotate gene dictionary to gff3 output
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = sorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # then loop through and write GFF3 format
    with open(output, 'w') as outfile:
        outfile.write("##gff-version 3\n")
        for k, v in list(sortedGenes.items()):
            outfile.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
            for i in range(0, len(v['ids'])):
                # build extra annotations for each transcript if applicable
                # now write mRNA feature
                outfile.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};TPM={:}\n".format(
                    v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['tpm'][i]))
                if v['type'] == 'mRNA':
                    if '5UTR' in v:
                        # if 5'UTR then write those first
                        num_5utrs = len(v['5UTR'][i])
                        if num_5utrs > 0:
                            for z in range(0, num_5utrs):
                                u_num = z + 1
                                outfile.write("{:}\t{:}\tfive_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr5p{:};Parent={:};\n".format(
                                    v['contig'], v['source'], v['5UTR'][i][z][0], v['5UTR'][i][z][1], v['strand'], v['ids'][i], u_num, v['ids'][i]))
                    # write the exons
                    num_exons = len(v['mRNA'][i])
                    for x in range(0, num_exons):
                        ex_num = x + 1
                        outfile.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                            v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
                    # if 3'UTR then write
                    if '3UTR' in v:
                        num_3utrs = len(v['3UTR'][i])
                        if num_3utrs > 0:
                            for z in range(0, num_3utrs):
                                u_num = z + 1
                                outfile.write("{:}\t{:}\tthree_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr3p{:};Parent={:};\n".format(
                                    v['contig'], v['source'], v['3UTR'][i][z][0], v['3UTR'][i][z][1], v['strand'], v['ids'][i], u_num, v['ids'][i]))
                if v['type'] == 'mRNA':
                    num_cds = len(v['CDS'][i])
                    # GFF3 phase is 1 less than flat file
                    current_phase = v['codon_start'][i] - 1
                    for y in range(0, num_cds):
                        outfile.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(
                            v['contig'], v['source'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                        current_phase = (
                            current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                        if current_phase == 3:
                            current_phase = 0


def Quarry2GFF3(input, output):
    with open(output, 'w') as outfile:
        outfile.write(("##gff-version 3\n"))
        exonCounts = {}
        GeneCount = 1
        with open(input, 'r') as infile:
            for line in infile:
                line = line.strip()
                contig, source, feature, start, end, score, strand, phase, attributes = line.split(
                    '\t')
                source = 'CodingQuarry'
                ID, Parent, Name = (None,)*3
                info = attributes.split(';')
                for x in info:
                    if x.startswith('ID='):
                        ID = x.replace('ID=', '')
                    elif x.startswith('Parent='):
                        Parent = x.replace('Parent=', '')
                if ID and ' ' in ID:
                    ID = ID.split(' ')[0]
                if Parent and ' ' in Parent:
                    Parent = Parent.split(' ')[0]
                if feature == 'gene':
                    geneID = 'gene_'+str(GeneCount)
                    transID = 'transcript_'+str(GeneCount)+'-T1'
                    outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:};Name={:};Alias={:};\n'.format(
                        contig, source, feature, start, end, score, strand, phase, geneID, geneID, ID))
                    outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:};Parent={:};Alias={:};\n'.format(
                        contig, source, 'mRNA', start, end, '.', strand, '.', transID, geneID, ID))
                    GeneCount += 1
                elif feature == 'CDS':
                    if not transID in exonCounts:
                        exonCounts[transID] = 1
                    else:
                        exonCounts[transID] += 1
                    num = exonCounts.get(transID)
                    outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:}.exon{:};Parent={:};\n'.format(
                        contig, source, 'exon', start, end, '.', strand, '.', transID, num, transID))
                    outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:}.cds;Parent={:};\n'.format(
                        contig, source, feature, start, end, score, strand, phase, transID, transID))


def runStringtie(bamfile, cpus, output):
    '''
    Function to run stringtie from bamfile
    Note that when given the bamfile, no way to determine strandeness so will run unstranded
    '''
    cmd = ['stringtie', '-p', str(cpus), os.path.realpath(bamfile)]
    runSubprocess2(cmd, '.', log, os.path.abspath(output))


def runCodingQuarry(genome, stringtie, cpus, output):
    '''
    run CodingQuarry from stringtie GTF input file
    '''
    # first get basename directory as need to create tmp CodingQuarry dir
    basedir = os.path.dirname(genome)
    tmpdir = os.path.join(basedir, 'CodingQuarry')
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # convert GTF to GFF3 file
    stringtieGFF3 = os.path.join(basedir, 'stringtie.gff3')
    Genes = gtf2dict(stringtie)
    Stringtie_dict2gff3(Genes, stringtieGFF3)
    # now setup command and run from tmpdir folder
    cmd = ['CodingQuarry', '-p',
           str(cpus), '-f', os.path.realpath(genome), '-t', os.path.realpath(stringtieGFF3)]
    runSubprocess(cmd, tmpdir, log)
    # capture results and reformat to proper GFF3
    result = os.path.join(tmpdir, 'out', 'PredictedPass.gff3')
    if not checkannotations(result):
        log.error('CodingQuarry failed, moving on without result, check logfile')
        return False
    else:
        Quarry2GFF3(result, output)
        return True


def runCodingQuarryTrained(genome, species, tmpdir, cpus, output):
    # now setup command and run from tmpdir folder
    log.info(
        'CodingQuarry prediction is running using {:} paremeters'.format(species))
    cmd = ['CodingQuarry', '-p',
           str(cpus), '-f', os.path.realpath(genome), '-s', species]
    log.debug(' '.join(cmd))
    myENV = os.environ
    if 'QUARRY_PATH' in myENV:
        del myENV['QUARRY_PATH']
    FNULL = open(os.devnull, 'w')
    p1 = subprocess.Popen(cmd, stdout=FNULL, stderr=FNULL,
                          cwd=tmpdir, env=dict(myENV))
    p1.communicate()

    # capture results and reformat to proper GFF3
    result = os.path.join(tmpdir, 'out', 'PredictedPass.gff3')
    if not checkannotations(result):
        log.error('CodingQuarry failed, moving on without result, check logfile')
        return False
    else:
        Quarry2GFF3(result, output)
        return True


def dict2gtf(input, output):
    from collections import OrderedDict

    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # sort the annotations by contig and start location
    sGenes = sorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    with open(output, 'w') as gtfout:
        for k, v in list(sortedGenes.items()):
            if v['type'] != 'mRNA':
                continue
            if 'pseudo' in v:
                if v['pseudo']:
                    continue
            if v['type'] == 'mRNA' and not v['CDS']:
                continue
            if v['type'] == 'mRNA' and not len(v['ids']) == len(v['mRNA']) == len(v['CDS']):
                continue
            for i in range(0, len(v['ids'])):
                # create attributes string
                attributes = 'gene_id "{:}"; transcript_id "{:}";'.format(
                    k, v['ids'][i])
                # if v['name']:
                #    attributes = attributes + ' Name "{:}";'.format(v['name'])
                if len(v['5UTR'][i]) > 0:
                    for utr in v['5UTR'][i]:
                        gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                            v['contig'], v['source'], '5UTR', utr[0], utr[1], 0, v['strand'], 0, attributes))
                if not v['partialStart'][i]:
                    if v['strand'] == '+':
                        startCodon = (v['CDS'][i][0][0], v['CDS'][i][0][0]+2)
                    else:
                        startCodon = (v['CDS'][i][0][1]-2, v['CDS'][i][0][1])
                    gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                        v['contig'], v['source'], 'start_codon', startCodon[0], startCodon[1], 0, v['strand'], 0, attributes))
                for x, cds in enumerate(v['CDS'][i]):
                    if v['partialStop'][i]:  # then just write the whole CDS as no reason to move codon back
                        gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                            v['contig'], v['source'], 'CDS', cds[0], cds[1], 0, v['strand'], v['phase'][i][x], attributes))
                    else:
                        if v['strand'] == '+':
                            if x == len(v['CDS'][i])-1:  # this is last one
                                gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                                    v['contig'], v['source'], 'CDS', cds[0], cds[1]-3, 0, v['strand'], v['phase'][i][x], attributes))
                                gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                                    v['contig'], v['source'], 'stop_codon', cds[1]-2, cds[1], 0, v['strand'], 0, attributes))
                            else:
                                gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                                    v['contig'], v['source'], 'CDS', cds[0], cds[1], 0, v['strand'], v['phase'][i][x], attributes))
                        else:
                            if x == len(v['CDS'][i])-1:  # this is last one
                                gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                                    v['contig'], v['source'], 'CDS', cds[0]+3, cds[1], 0, v['strand'], v['phase'][i][x], attributes))
                                gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                                    v['contig'], v['source'], 'stop_codon', cds[0], cds[0]+2, 0, v['strand'], 0, attributes))
                            else:
                                gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                                    v['contig'], v['source'], 'CDS', cds[0], cds[1], 0, v['strand'], v['phase'][i][x], attributes))
                if len(v['3UTR'][i]) > 0:
                    for utr in v['3UTR'][i]:
                        gtfout.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(
                            v['contig'], v['source'], '3UTR', utr[0], utr[1], 0, v['strand'], 0, attributes))
            gtfout.write('\n')


def gff3_to_gtf(input, genome, output):
    Genes = {}
    Genes = gff2dict(input, genome, Genes)
    dict2gtf(Genes, output)


def gb2allout(input, GFF, Proteins, Transcripts, DNA):
    '''
    function to split GBK file into parts, need to be able to deal with multiple transcripts and get naming correct
    assumption is that the mRNA and CDS features from multiple transcripts are in order, i.e. the first mRNA feature
    you see corresponds to first CDS feature, etc. **hopefully this is an okay assumption**
    '''
    # idea is to populate the dictionary first, then write GFF, proteins, transcripts, can write DNA on first pass
    genes = {}
    with open(DNA, 'w') as scaffolds:
        with open(input, 'r') as gbk:
            for record in SeqIO.parse(gbk, 'genbank'):
                scaffolds.write(">{:}\n{:}\n".format(
                    record.id, softwrap(str(record.seq))))
                for f in record.features:
                    gb_feature_add2dict(f, record, genes)
    # write GFF
    dict2gff3_old(genes, GFF)
    # write to protein and transcripts
    dict2nucleotides(genes, Proteins, Transcripts)


def minimap2Align(transcripts, genome, cpus, intron, output):
    '''
    function to align transcripts to genome using minimap2
    huge speed increase over gmap + blat
    '''
    FNULL = open(os.devnull, 'w')
    bamthreads = int(round(int(cpus) / 2))
    if bamthreads > 4:
        bamthreads = 4
    minimap2_cmd = ['minimap2', '-ax', 'splice', '-t',
                    str(cpus), '--cs', '-u', 'b', '-G', str(intron), genome,
                    transcripts]
    samtools_cmd = ['samtools', 'sort', '--reference', genome,
                    '-@', str(bamthreads), '-o', output, '-']
    log.debug('{} | {}'.format(' '.join(minimap2_cmd), ' '. join(samtools_cmd)))
    p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    p2 = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=FNULL, stdin=p1.stdout)
    p1.stdout.close()
    p2.communicate()


def iso_seq_minimap2(transcripts, genome, cpus, intron, output):
    '''
    function to align PB iso-seq data
    '''
    FNULL = open(os.devnull, 'w')
    bamthreads = int(round(int(cpus) / 2))
    if bamthreads > 4:
        bamthreads = 4
    minimap2_cmd = ['minimap2', '-ax', 'splice', '-t',
                    str(cpus), '--cs', '-uf', '-C5', '-G', str(intron), genome,
                    transcripts]
    samtools_cmd = ['samtools', 'sort', '--reference', genome,
                    '-@', str(bamthreads), '-o', output, '-']
    log.debug('{} | {}'.format(' '.join(minimap2_cmd), ' '. join(samtools_cmd)))
    p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    p2 = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=FNULL, stdin=p1.stdout)
    p1.stdout.close()
    p2.communicate()


def nanopore_cDNA_minimap2(transcripts, genome, cpus, intron, output):
    '''
    function to nanopore 2d cDNA
    '''
    FNULL = open(os.devnull, 'w')
    bamthreads = int(round(int(cpus) / 2))
    if bamthreads > 4:
        bamthreads = 4
    minimap2_cmd = ['minimap2', '-ax', 'splice', '-t',
                    str(cpus), '--cs', '-G', str(intron), genome, transcripts]
    samtools_cmd = ['samtools', 'sort', '--reference', genome,
                    '-@', str(bamthreads), '-o', output, '-']
    log.debug('{} | {}'.format(' '.join(minimap2_cmd), ' '. join(samtools_cmd)))
    p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    p2 = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=FNULL, stdin=p1.stdout)
    p1.stdout.close()
    p2.communicate()


def nanopore_mRNA_minimap2(transcripts, genome, cpus, intron, output):
    '''
    function to nanopore direct mRNA reads
    '''
    FNULL = open(os.devnull, 'w')
    bamthreads = int(round(int(cpus) / 2))
    if bamthreads > 4:
        bamthreads = 4
    minimap2_cmd = ['minimap2', '-ax', 'splice', '-t',
                    str(cpus), '--cs', '-uf', '-k14', '-G', str(intron),
                    genome, transcripts]
    samtools_cmd = ['samtools', 'sort', '--reference', genome,
                    '-@', str(bamthreads), '-o', output, '-']
    log.debug('{} | {}'.format(' '.join(minimap2_cmd), ' '. join(samtools_cmd)))
    p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
    p2 = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=FNULL, stdin=p1.stdout)
    p1.stdout.close()
    p2.communicate()


def mergeBAMs(*args, **kwargs):
    cmd = ['samtools', 'merge', '-@', str(kwargs['cpus']), kwargs['output']]
    cmd = cmd + list(args)
    runSubprocess(cmd, '.', log)


def catFiles(*args, **kwargs):
    cmd = ['cat']
    cmd = cmd + list(args)
    runSubprocess2(cmd, '.', log, kwargs['output'])


def runGMAP(transcripts, genome, cpus, intron, tmpdir, output):
    # first build genome database
    build_log = os.path.join(tmpdir, 'gmap-build.log')
    with open(build_log, 'w') as logfile:
        subprocess.call(['gmap_build', '-D', tmpdir, '-d', 'genome',
                         '-k', '13', genome], stdout=logfile, stderr=logfile)
    # now map transcripts
    map_log = os.path.join(tmpdir, 'gmap-map.log')
    with open(map_log, 'w') as logfile:
        with open(output, 'w') as out:
            subprocess.call(['gmap', '--cross-species', '-f', '3', '-K', str(intron), '-n', '1', '-t', str(
                cpus), '-B', '5', '-D', tmpdir, '-d', 'genome', transcripts], stdout=out, stderr=logfile)


def runBUSCO(input, Database, cpus, tmpdir, output):
    # run busco in protein mapping mode
    if (sys.version_info > (3, 0)):
        BUSCO = os.path.join(parentdir,
                             'aux_scripts', 'funannotate-BUSCO2.py')
    else:
        BUSCO = os.path.join(parentdir,
                             'aux_scripts', 'funannotate-BUSCO2-py2.py')
    cmd = [BUSCO, '-i', input, '-m', 'proteins', '-l',
           Database, '-o', 'busco', '-c', str(cpus), '-f']
    runSubprocess(cmd, tmpdir, log)
    # now parse output and write to annotation file
    with open(output, 'w') as out:
        with open(os.path.join(tmpdir, 'run_busco', 'full_table_busco.tsv'), 'r') as busco:
            for line in busco:
                if line.startswith('#'):
                    continue
                col = line.split('\t')
                # if diploid these should show up, but problematic for drawing trees....
                if col[1] == 'Complete' or col[1] == 'Duplicated':
                    out.write("%s\tnote\tBUSCO:%s\n" % (col[2], col[0]))


def dupBUSCO2gff(ID, base_folder, locationID):
    hmmerfolder = os.path.join(base_folder, 'hmmer_output')
    geneID = ''
    AugFile = ''
    GFFfile = os.path.join(base_folder, 'augustus_output', 'gffs', ID+'.gff')
    if geneID == '':
        for file in os.listdir(hmmerfolder):
            if file.startswith(ID):
                with open(os.path.join(hmmerfolder, file), 'r') as hmmer:
                    for line in hmmer:
                        if not line.startswith('#'):
                            longID = line.split()[0]
                            longID = longID.replace(']', '')
                            partsID = longID.split('[')
                            if locationID == partsID[1]:
                                geneID = partsID[0]
                                AugFile = os.path.join(
                                    base_folder, 'augustus_output', 'predicted_genes', file)
                                break
    # so now should have gene name, get the GFF from augustus
    with open(GFFfile, 'w') as gffout:
        with open(AugFile, 'r') as augustus:
            for pred in readBlocks(augustus, '# start gene'):
                if pred[0].startswith('# This output'):
                    continue
                if pred[0].startswith('##gff-version 3'):
                    continue
                if pred[0].startswith('# Please cite'):
                    continue
                if geneID in pred[0]:
                    for x in pred:
                        if not x.startswith('#'):
                            gffout.write(x)


def getCompleteBuscos(input, ploidy=1):
    busco_complete = {}
    with open(input, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('#'):
                continue
            passing = ['Complete']
            if ploidy > 1:
                passing.append('Duplicated')
            cols = line.split('\t')
            if cols[1] in passing:
                busco, status, gene, score, length = cols
                if gene not in busco_complete:
                    busco_complete[gene] = busco
    return busco_complete


def filterGFF3(keepDict, genome, gff3, output):
    #load into Dictionary
    Genes = {}
    Genes = gff2dict(gff3, genome, Genes)
    filtered = {}
    for k,v in Genes.items():
        if v['ids'][0] in keepDict:
            filtered[k] = v
    dict2gff3(filtered, output)


def parseBUSCO2genome(input, ploidy, ContigSizes, output):
    # input is BUSCO output, ploidy is integer, ContigSizes is dictionary, output is a bedfile, function returns dictionary
    busco_complete = {}
    hits = {}
    with open(output, 'w') as bedfile:
        with open(input, 'r') as buscoinput:
            for line in buscoinput:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                if cols[1] == 'Complete' or cols[1] == 'Duplicated':
                    contig = cols[2]
                    start = cols[3]
                    end = cols[4]
                    score = cols[5]
                    length = cols[6]
                    ID = contig+':'+start+'-'+end
                    if cols[1] == 'Complete':
                        if not cols[0] in hits:
                            hits[cols[0]] = (
                                ID, score, contig, start, end, length)
                    if ploidy > 1:
                        if cols[1] == 'Duplicated':
                            if not cols[0] in hits:
                                hits[cols[0]] = (
                                    ID, score, contig, start, end, length)
                                dupBUSCO2gff(
                                    cols[0], os.path.dirname(input), ID)
                            else:
                                oldscore = float(hits.get(cols[0])[1])
                                if float(score) > oldscore:
                                    hits[cols[0]] = (
                                        ID, score, contig, start, end, length)
                                    dupBUSCO2gff(
                                        cols[0], os.path.dirname(input), ID)
            for k, v in natsorted(list(hits.items())):
                # validate locations for bedfile, move 100 bp in each direction for bedfile
                start = int(v[3]) - 100
                if start < 1:  # negative no good
                    start = 1
                end = int(v[4]) + 100
                # check it doesn't go past contig length
                if end > ContigSizes.get(contig):
                    end = ContigSizes.get(contig)
                bedfile.write('%s\t%i\t%i\t%s\n' % (contig, start, end, k))
                busco_complete[k] = v[0]
    return busco_complete


def RepeatBlast(input, cpus, evalue, DataBase, tmpdir, output, diamond=True):
    # run blastp against repeats
    blast_tmp = os.path.join(tmpdir, 'repeats.xml')
    if diamond:
        blastdb = os.path.join(DataBase, 'repeats.dmnd')
        cmd = ['diamond', 'blastp', '--sensitive', '--query', input, '--threads', str(cpus),
               '--out', blast_tmp, '--db', blastdb, '--evalue', str(evalue), '--max-target-seqs', '1', '--outfmt', '5']
    else:
        blastdb = os.path.join(DataBase, 'REPEATS')
        cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus),
               '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input]
    runSubprocess4(cmd, '.', log)
    # parse results
    with open(output, 'w') as out:
        with open(blast_tmp, 'r') as results:
            for qresult in SearchIO.parse(results, "blast-xml"):
                hits = qresult.hits
                ID = qresult.id
                num_hits = len(hits)
                if num_hits > 0:
                    length = 0
                    for i in range(0, len(hits[0].hsps)):
                        length += hits[0].hsps[i].aln_span
                    pident = hits[0].hsps[0].ident_num / float(length)
                    out.write("%s\t%s\t%f\t%s\n" %
                              (ID, hits[0].id, pident, hits[0].hsps[0].evalue))


def eggnog2dict(annotations):
    # load in annotation dictionary
    EggNog = {}
    with open(annotations, 'r') as input:
        reader = csv.reader(input, delimiter='\t')
        for line in reader:
            EggNog[line[1]] = line[5]
    return EggNog


def number_present(s):
    return any(i.isdigit() for i in s)


def capfirst(x):
    return x[0].upper() + x[1:]


def item2index(inputList, item):
    # return the index of an item in the input list
    item_index = None
    for x in inputList:
        if item in x:
            item_index = inputList.index(x)
    return item_index


def getEggNogHeaders(input):
    IDi, DBi, OGi, Genei, COGi, Desci = (None,)*6
    with open(input, 'r') as infile:
        for line in infile:
            line = line.replace('\n', '')
            if line.startswith('#query_name'):  # this is HEADER
                headerCols = line.split('\t')
                IDi = item2index(headerCols, 'query_name')
                Genei = item2index(headerCols, 'predicted_gene_name')
                DBi = item2index(headerCols, 'Annotation_tax_scope')
                OGi = item2index(headerCols, 'OGs')
                COGi = item2index(headerCols, 'COG cat')
                Desci = item2index(headerCols, 'eggNOG annot')
                break
    return IDi, DBi, OGi, Genei, COGi, Desci


def parseEggNoggMapper(input, output):
    Definitions = {}
    # indexes from header file
    IDi, DBi, OGi, Genei, COGi, Desci = getEggNogHeaders(input)
    # take annotations file from eggnog-mapper and create annotations
    with open(output, 'w') as out:
        with open(input, 'r') as infile:
            for line in infile:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                cols = line.split('\t')
                ID = cols[IDi]
                DB = cols[DBi].split('[')[0]
                OGs = cols[OGi].split(',')
                NOG = ''
                for x in OGs:
                    if DB in x:
                        NOG = 'ENOG41' + x.split('@')[0]
                Gene = ''
                if cols[Genei] != '':
                    if not '_' in cols[Genei] and not '.' in cols[Genei] and number_present(cols[Genei]):
                        Gene = cols[Genei]
                Description = cols[Desci]
                if NOG == '':
                    continue
                if not NOG in Definitions:
                    Definitions[NOG] = Description
                out.write("%s\tnote\tEggNog:%s\n" % (ID, NOG))
                if cols[COGi] != '':
                    out.write("%s\tnote\tCOG:%s\n" %
                              (ID, cols[COGi].replace(' ', '')))
                if Gene != '':
                    product = Gene.lower()+'p'
                    product = capfirst(product)
                    out.write("%s\tname\t%s\n" % (ID.split('-T')[0], Gene))
                    out.write("%s\tproduct\t%s\n" % (ID, product))
                    if Description != '':
                        out.write("%s\tnote\t%s\n" % (ID, Description))
    return Definitions


def batch_iterator(iterator, batch_size):
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def fasta2chunks(input, chunks, tmpdir, output):
    # split the input fasta file into 20 chunks to process
    with open(input, 'r') as seqs:
        SeqCount = countfasta(input)
        SeqRecords = SeqIO.parse(seqs, 'fasta')
        chunks = SeqCount / int(chunks)
        # divide into chunks, store in tmp file
        folder = os.path.join(tmpdir, output)
        if not os.path.exists(folder):
            os.makedirs(folder)
        else:
            shutil.rmtree(folder)
            os.makedirs(folder)
        for i, batch in enumerate(batch_iterator(SeqRecords, chunks)):
            filename = "chunk_%i.fa" % (i+1)
            tmpout = os.path.join(folder, filename)
            handle = open(tmpout, "w")
            SeqIO.write(batch, handle, "fasta")
            handle.close()


def signalP(input, tmpdir, output):
    # split input file into chunks, 20 should mean < 200 proteins per chunk
    from funannotate.check import check_version7
    version = check_version7('signalp')
    if '.' in version:
        version = int(version.split('.')[0])
    if version > 4:
        cmd = ['signalp', '-stdout', '-org', 'euk', '-format', 'short', '-fasta']
    else:
        cmd = ['signalp', '-t', 'euk', '-f', 'short']
    fasta2chunks(input, 40, tmpdir, 'signalp_tmp')
    for file in os.listdir(os.path.join(tmpdir, 'signalp_tmp')):
        if file.startswith('chunk'):
            file = os.path.join(tmpdir, 'signalp_tmp', file)
            tmp_out = file.replace('.fa', '.signalp.out')
            cmd1 = cmd + [file]
            runSubprocess2(cmd1, '.', log, tmp_out)
    # now concatenate all outputs
    if os.path.isfile(output):
        os.remove(output)
    with open(output, 'a') as finalout:
        for file in os.listdir(os.path.join(tmpdir, 'signalp_tmp')):
            if file.endswith('.signalp.out'):
                file = os.path.join(tmpdir, 'signalp_tmp', file)
                with open(file) as infile:
                    finalout.write(infile.read())
    # cleanup tmp directory
    shutil.rmtree(os.path.join(tmpdir, 'signalp_tmp'))


def parseSignalP(sigP, secretome_annot):
    sigpDict = {}
    version = 4
    with open(sigP, 'r') as results:
        for line in results:
            line = line.rstrip()
            if line.startswith('#'):
                if line.startswith('# SignalP-5'):
                    version = 5
                continue
            if version < 5:
                col = line.split(' ')  # not tab delimited
                col = [_f for _f in col if _f]  # clean up empty spaces
                if col[9] == 'Y':  # then there is signal peptide
                    ID = col[0]
                    end = int(col[2]) - 1
                    sigpDict[ID] = [end, '', '']
            else:  # version 5 has different format and tab delimited hooray!
                if '\t' in line:
                    cols = line.split('\t')
                    if cols[1] != 'OTHER':  # then signal peptide
                        ID, prediction, score1, score2, position = cols[:5]
                        components = position.split()
                        pos = components[2].split('-')[0]
                        prob = components[-1]
                        aa = components[3].replace('.', '')
                        sigpDict[ID] = [pos, aa, prob]

    with open(secretome_annot, 'w') as secout:
        for k, v in natsorted(list(sigpDict.items())):
            if v[1] != '':
                secout.write("{:}\tnote\tSECRETED:SignalP(1-{:},cutsite={:},prob={:})\n".format(k, v[0], v[1], v[2]))
            else:
                secout.write("{:}\tnote\tSECRETED:SignalP(1-{:})\n".format(k, v[0]))


def parsePhobiusSignalP(phobius, sigP, membrane_annot, secretome_annot):
    # give directory of annotate_misc, first get phobius results
    '''
    This is what phobius results look like
    ID  TM  SP  Prediction
    VE00_00001  0   0   o
    VE00_00002  2   0   i198-219o283-301i
    VE00_00003  0   0   o
    VE00_00004  0   Y   n8-18c23/24o
    VE00_00005  12  0   i49-69o89-107i119-138o144-167i179-200o212-234i280-299o319-341i348-366o378-398i410-430o442-465i
    '''
    pSecDict = {}
    pTMDict = {}
    sigpDict = {}
    # parsing short format phobius
    with open(phobius, 'r') as input1:
        for line in input1:
            line = line.rstrip()
            if line.startswith('ID') or line.startswith('SEQ'):
                continue
            if '\t' in line:
                cols = line.split('\t')
            else:
                cols = line.split()
            geneID = cols[0]
            if int(cols[1]) > 0:  # then found TM domain
                annot = cols[3]
                if not geneID in pTMDict:
                    pTMDict[geneID] = 'TransMembrane:'+cols[1]+' ('+annot+')'
            if cols[2] == 'Y':  # then sig pep discovered
                location = cols[3].split('/')[0]
                clevage = location.split('c')[-1]
                if not geneID in pSecDict:
                    pSecDict[geneID] = [clevage, '', '']
    if sigP:  # will be passed FALSE if signalP data missing
        # parse signalp output and turn into annotation file
        version = 4
        with open(sigP, 'r') as results:
            for line in results:
                line = line.rstrip()
                if line.startswith('#'):
                    if line.startswith('# SignalP-5'):
                        version = 5
                    continue
                if version < 5:
                    col = line.split(' ')  # not tab delimited
                    col = [_f for _f in col if _f]  # clean up empty spaces
                    if col[9] == 'Y':  # then there is signal peptide
                        ID = col[0]
                        end = int(col[2]) - 1
                        sigpDict[ID] = [end, '', '']
                else:  # version 5 has different format and tab delimited hooray!
                    if '\t' in line:
                        cols = line.split('\t')
                        if cols[1] != 'OTHER':  # then signal peptide
                            ID, prediction, score1, score2, position = cols[:5]
                            components = position.split()
                            pos = components[2].split('-')[0]
                            prob = components[-1]
                            aa = components[3].replace('.', '')
                            sigpDict[ID] = [pos, aa, prob]
    else:
        sigpDict = pSecDict
    # write annotation files
    with open(membrane_annot, 'w') as memout:
        for k, v in natsorted(list(pTMDict.items())):
            memout.write("%s\tnote\t%s\n" % (k, v))
    with open(secretome_annot, 'w') as secout:
        for k, v in natsorted(list(sigpDict.items())):
            if v[1] != '':
                secout.write("{:}\tnote\tSECRETED:SignalP(1-{:},cutsite={:},prob={:})\n".format(k, v[0], v[1], v[2]))
            else:
                secout.write("{:}\tnote\tSECRETED:SignalP(1-{:})\n".format(k, v[0]))


def n_lower_chars(string):
    return sum(1 for c in string if c.islower())


def CheckAugustusSpecies(input):
    # get the possible species from augustus
    augustus_list = []
    for i in os.listdir(os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], 'species')):
        if not i.startswith('.'):
            augustus_list.append(i)
    augustus_list = set(augustus_list)
    if input in augustus_list:
        return True
    else:
        return False


def CheckFunannotateSpecies(input, db):
    # get the possible species from funannotateDB dir -- on install mirrored Augustus
    species_list = []
    for i in os.listdir(os.path.join(db, 'trained_species')):
        if not i.startswith('.'):
            species_list.append(i)
    species_list = set(species_list)
    if input in species_list:
        return True
    else:
        return False


def SortRenameHeaders(input, output):
    # sort records and write temp file
    with open(output, 'w') as out:
        with open(input, 'r') as input:
            records = list(SeqIO.parse(input, 'fasta'))
            records.sort(cmp=lambda x, y: cmp(len(y), len(x)))
            counter = 1
            for rec in records:
                rec.name = ''
                rec.description = ''
                rec.id = 'scaffold_' + str(counter)
                counter += 1
            SeqIO.write(records, out, 'fasta')


def validate_tRNA(input, genes, gaps, output):
    # run bedtools intersect to keep only input that dont intersect with either genes or gaps
    sortedInput = os.path.abspath(input)+'.sorted.gff3'
    #sortGFFproper(input, sortedInput)
    cmd1 = ['bedtools', 'sort', '-i', input]
    with open(sortedInput, 'w') as outfile:
        subprocess.call(cmd1, stdout=outfile)
    sortedGenes = os.path.abspath(genes)+'.sorted.gff3'
    #sortGFFproper(genes, sortedGenes)
    cmd2 = ['bedtools', 'sort', '-i', genes]
    with open(sortedGenes, 'w') as outfile:
        subprocess.call(cmd2, stdout=outfile)
    if gaps:
        sortedGaps = os.path.abspath(gaps)+'.sorted.gff3'
        #sortGFFproper(gaps, sortedGaps)
        cmd3 = ['bedtools', 'sort', '-i', gaps]
        with open(sortedGaps, 'w') as outfile:
            subprocess.call(cmd3, stdout=outfile)
    cmd = ['bedtools', 'intersect', '-sorted', '-v', '-a', sortedInput, '-b', sortedGenes]
    if gaps:
        cmd.append(sortedGaps)
    tmpOut = os.path.abspath(output)+'.tmp'
    runSubprocess2(cmd, '.', log, tmpOut)
    # now sort properly
    sortGFFproper(tmpOut, output)
    os.remove(tmpOut)


# via https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
def list2groups(L):
    if len(L) < 1:
        return
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:  # Part of the group, bump the end
            last = n
        else:  # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last  # Yield the last group


def checkMask(genome, bedfile):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    # load contig names and sizes into dictionary, get masked repeat stats
    GenomeLength = 0
    maskedSize = 0
    masked = {}
    ContigSizes = {}
    with open(genome, 'r') as input:
        for header, Seq in SimpleFastaParser(input):
            if ' ' in header:
                ID = header.split(' ')[0]
            else:
                ID = header
            if not ID in masked:
                masked[ID] = []
            if not ID in ContigSizes:
                ContigSizes[ID] = len(Seq)
            GenomeLength += len(Seq)
            maskedSize += n_lower_chars(Seq)
            for i, c in enumerate(Seq):
                if c.islower():
                    masked[ID].append(i)  # 0 based
    if maskedSize == 0:  # not softmasked, return False
        with open(bedfile, 'w') as bedout:
            bedout.write('')
        return ContigSizes, GenomeLength, maskedSize, 0.0
    else:
        counter = 1
        with open(bedfile, 'w') as bedout:
            for k, v in natsorted(list(masked.items())):
                repeats = list(list2groups(v))
                for item in repeats:
                    if len(item) == 2:
                        bedout.write('{:}\t{:}\t{:}\tRepeat_{:}\n'.format(
                            k, item[0], item[1], counter))
                        counter += 1
    percentMask = maskedSize / float(GenomeLength)
    return ContigSizes, GenomeLength, maskedSize, percentMask


def maskingstats2bed(input, counter, alock):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    masked = []
    gaps = []
    maskedSize = 0
    bedfilename = input.replace('.fasta', '.bed')
    gapfilename = input.replace('.fasta', '.gaps')
    with open(input, 'r') as infile:
        for header, Seq in SimpleFastaParser(infile):
            if ' ' in header:
                ID = header.split(' ')[0]
            else:
                ID = header
            for i, c in enumerate(Seq):
                if c == 'N' or c == 'n':
                    masked.append(i)
                    maskedSize += 1
                    gaps.append(i)
                elif c.islower():
                    masked.append(i)  # 0 based
                    maskedSize += 1

    if maskedSize > 0:  # not softmasked, return False
        with open(bedfilename, 'w') as bedout:
            repeats = list(list2groups(masked))
            for item in repeats:
                if len(item) == 2:
                    bedout.write('{:}\t{:}\t{:}\tRepeat_\n'.format(
                        ID, item[0], item[1]))
    if len(gaps) > 0:
        with open(gapfilename, 'w') as gapout:
            bedGaps = list(list2groups(gaps))
            for item in bedGaps:
                if len(item) == 2:
                    gapout.write(
                        '{:}\t{:}\t{:}\tassembly-gap_\n'.format(ID, item[0], item[1]))
    with alock:
        counter.value += maskedSize


def mask_safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        maskingstats2bed(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def checkMasklowMem(genome, bedfile, gapsfile, cpus):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    # load contig names and sizes into dictionary, get masked repeat stats
    ContigSizes = {}
    tmpdir = os.path.join(os.path.dirname(genome), 'mask_'+str(uuid.uuid4()))
    os.makedirs(tmpdir)
    file_list = []
    with open(genome, 'r') as input:
        for header, Seq in SimpleFastaParser(input):
            if ' ' in header:
                ID = header.split(' ')[0]
            else:
                ID = header
            if not ID in ContigSizes:
                ContigSizes[ID] = len(Seq)
            with open(os.path.join(tmpdir, ID+'.fasta'), 'w') as fastaout:
                fastaout.write('>{:}\n{:}\n'.format(ID, Seq))
            file_list.append(os.path.join(tmpdir, ID+'.fasta'))
    # num = 1
    p = multiprocessing.Pool(processes=cpus)
    TotalMask = multiprocessing.Manager().Value('i', 0)
    lock = multiprocessing.Manager().Lock()
    result = []
    for i in file_list:
        result.append(p.apply_async(mask_safe_run, [i, TotalMask, lock]))
    p.close()
    p.join()
    repeatNum = 1
    gapNum = 1
    with open(bedfile, 'w') as bedout:
        for file in natsorted(os.listdir(tmpdir)):
            if file.endswith('.bed'):
                with open(os.path.join(tmpdir, file), 'r') as infile:
                    for line in infile:
                        line = line.replace(
                            'Repeat_', 'Repeat_'+str(repeatNum))
                        bedout.write(line)
                        repeatNum += 1
    with open(gapsfile, 'w') as gapout:
        for file in natsorted(os.listdir(tmpdir)):
            if file.endswith('.gaps'):
                with open(os.path.join(tmpdir, file), 'r') as infile:
                    for line in infile:
                        line = line.replace(
                            'assembly-gap_', 'assembly-gap_'+str(gapNum))
                        gapout.write(line)
                        gapNum += 1

    SafeRemove(tmpdir)
    GenomeLength = sum(ContigSizes.values())
    percentMask = TotalMask.value / float(GenomeLength)
    return ContigSizes, GenomeLength, TotalMask.value, percentMask


def RunGeneMarkES(command, input, ini, maxintron, softmask, cpus, tmpdir, output, fungus):
    # make directory to run script from
    outdir = os.path.join(tmpdir, 'genemark')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if cpus > 64:
        cpus = 64
    contigs = os.path.abspath(input)
    log.info("Running GeneMark-ES on assembly")
    cmd = [command, '--ES', '--max_intron', str(maxintron), '--soft_mask', str(
        softmask), '--cores', str(cpus), '--sequence', contigs]
    if fungus == 'fungus':
        cmd = cmd + ['--fungus']
    if ini:
        cmd = cmd + ['--ini_mod', os.path.abspath(ini)]
    runSubprocess3(cmd, outdir, log)
    # rename results and grab mod file
    try:
        os.rename(os.path.join(outdir, 'output', 'gmhmm.mod'),
                  os.path.join(tmpdir, 'gmhmm.mod'))
    except OSError:
        log.error("GeneMark-ES failed: {:} file missing, please check logfiles.".format(
            os.path.join(outdir, 'output', 'gmhmm.mod')))
    # convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join(outdir, 'genemark.gtf')
    if checkannotations(gm_gtf):
        # log.info("Converting GeneMark GTF file to GFF3")
        with open(output, 'w') as out:
            subprocess.call([GeneMark2GFF, gm_gtf], stdout=out)


def RunGeneMarkET(command, input, ini, evidence, maxintron, softmask, cpus, tmpdir, output, fungus):
    # make directory to run script from
    outdir = os.path.join(tmpdir, 'genemark')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if cpus > 64:
        cpus = 64
    contigs = os.path.abspath(input)
    # get only intron information from evidence
    hintsfile = os.path.join(tmpdir, 'genemark.intron-hints.gff')
    with open(hintsfile, 'w') as hints:
        with open(evidence, 'r') as evid:
            for line in evid:
                if '\tintron\t' in line and '\tb2h\t' in line:
                    tmprow = line.split("\t")
                    tmprow[5] = "500"  # for intron hint score to be 500
                    hints.write("\t".join(tmprow))
    log.info("Running GeneMark-ET on assembly")
    cmd = [command, '--ET', os.path.abspath(hintsfile), '--max_intron', str(
        maxintron), '--soft_mask', str(softmask), '--cores', str(cpus), '--sequence', contigs]
    if fungus == 'fungus':
        cmd = cmd + ['--fungus']
    if ini:
        cmd = cmd + ['--ini_mod', os.path.abspath(ini)]
    runSubprocess3(cmd, outdir, log)
    # rename results and grab mod file
    try:
        os.rename(os.path.join(outdir, 'output', 'gmhmm.mod'),
                  os.path.join(tmpdir, 'gmhmm.mod'))
    except OSError:
        log.error("GeneMark-ET failed: {:} file missing, please check logfiles.".format(
            os.path.join(outdir, 'output', 'gmhmm.mod')))
    # convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join(outdir, 'genemark.gtf')
    if checkannotations(gm_gtf):
        # log.info("Converting GeneMark GTF file to GFF3")
        with open(output, 'w') as out:
            subprocess.call([GeneMark2GFF, gm_gtf], stdout=out)


def dict2glimmer(input, output):
    # take funannotate dictionary convert to glimmer training format
    with open(output, 'w') as outfile:
        for k, v in list(input.items()):
            for i in range(0, len(v['ids'])):
                for c in v['CDS'][i]:
                    if v['strand'] == '+':
                        outfile.write('{:} {:} {:}\n'.format(
                            v['contig'], c[0], c[1]))
                    else:
                        outfile.write('{:} {:} {:}\n'.format(
                            v['contig'], c[1], c[0]))
            outfile.write('\n')


def glimmer2gff3(input, output):
    '''
    scaffold_39     GlimmerHMM      mRNA    23692   25015   .       +       .       ID=scaffold_39.path1.gene12;Name=scaffold_39.path1.gene12
    scaffold_39     GlimmerHMM      CDS     23692   23886   .       +       0       ID=scaffold_39.cds12.1;Parent=scaffold_39.path1.gene12;Name=scaffold_39.path1.gene12;Note=initial-exon
    scaffold_39     GlimmerHMM      CDS     24282   24624   .       +       0       ID=scaffold_39.cds12.2;Parent=scaffold_39.path1.gene12;Name=scaffold_39.path1.gene12;Note=internal-exon
    scaffold_39     GlimmerHMM      CDS     24711   25015   .       +       2       ID=scaffold_39.cds12.3;Parent=scaffold_39.path1.gene12;Name=scaffold_39.path1.gene12;Note=final-exon
    scaffold_39     GlimmerHMM      mRNA    25874   27899   .       -       .       ID=scaffold_39.path1.gene13;Name=scaffold_39.path1.gene13
    scaffold_39     GlimmerHMM      CDS     25874   26973   .       -       2       ID=scaffold_39.cds13.1;Parent=scaffold_39.path1.gene13;Name=scaffold_39.path1.gene13;Note=final-exon
    scaffold_39     GlimmerHMM      CDS     27257   27899   .       -       0       ID=scaffold_39.cds13.2;Parent=scaffold_39.path1.gene13;Name=scaffold_39.path1.gene13;Note=initial-exon
    '''
    with open(output, 'w') as outfile:
        outfile.write(("##gff-version 3\n"))
        exonCounts = {}
        GeneCount = 1
        skipList = []
        idsSeen = {}
        with open(input, 'r') as infile:
            for i, line in enumerate(infile):
                if line.startswith('##sequence-region'):
                    idsSeen = {}
                if line.startswith('#') or line.startswith('\n'):
                    continue
                line = line.strip()
                if line.count('\t') < 8:
                    print('ERROR parsing GlimmerHMM Raw output in line {}:\n   {}'.format(i+1, line))
                    continue
                contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
                ID, Parent, Name = (None,)*3
                info = attributes.split(';')
                for x in info:
                    if x.startswith('ID='):
                        ID = x.replace('ID=', '')
                    elif x.startswith('Parent='):
                        Parent = x.replace('Parent=', '')
                if Parent and Parent in skipList:
                    continue
                if feature == 'mRNA':
                    genelen = int(end) - int(start)
                    if genelen < 150:
                        if not ID in skipList:
                            skipList.append(ID)
                        continue
                    geneID = 'glimmerG_'+str(GeneCount)
                    transID = 'glimmerT_'+str(GeneCount)+'-T1'
                    idsSeen[ID] = (geneID, transID)
                    outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:};Alias={:};\n'.format(
                        contig, source, 'gene', start, end, score, strand, phase, geneID, ID))
                    outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:};Parent={:};Alias={:};\n'.format(
                        contig, source, 'mRNA', start, end, '.', strand, '.', transID, geneID, ID))
                    GeneCount += 1
                elif feature == 'CDS':
                    if Parent in idsSeen:
                        geneID, transID = idsSeen.get(Parent)
                        if not transID in exonCounts:
                            exonCounts[transID] = 1
                        else:
                            exonCounts[transID] += 1
                        num = exonCounts.get(transID)
                        outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:}.exon{:};Parent={:};\n'.format(
                            contig, source, 'exon', start, end, '.', strand, '.', transID, num, transID))
                        outfile.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\tID={:}.cds;Parent={:};\n'.format(
                            contig, source, feature, start, end, score, strand, phase, transID, transID))
                    else:
                        print('ERROR parsing GlimmerHMM Raw output in line {}:\n   {}'.format(i+1, line))


def runGlimmerHMM(fasta, gff3, dir, output):
    '''
    wrapper to run GlimmerHMM training followed by prediction
    input is GFF3 format high quality models, i.e. from PASA/transdecoder
    output is standard GFF3 format
    '''
    # generate training directory ontop of the dir that is passed
    tmpdir = os.path.join(dir, 'glimmerhmm')
    if os.path.isdir(tmpdir):
        SafeRemove(tmpdir)

    # generate glimmer training input
    # load gff3 into dictionary
    Genes = {}
    Genes = gff2dict(os.path.abspath(gff3), os.path.abspath(fasta), Genes)
    glimmExons = os.path.join(dir, 'glimmer.exons')
    dict2glimmer(Genes, glimmExons)

    # now run trainGlimmerHMM
    cmd = ['trainGlimmerHMM', os.path.abspath(
        fasta), os.path.abspath(glimmExons), '-d', tmpdir]
    runSubprocess4(cmd, '.', log)  # runSubproces4 --> stdout/stderr to devnull

    # now run GlimmerHMM prediction
    # glimmhmm.pl <glimmerhmm_program> <fasta_file> <train_dir> <options>
    glimmerRaw = os.path.abspath(os.path.join(dir, 'glimmerHMM.output.raw'))
    cmd = ['perl', which_path('glimmhmm.pl'), which_path(
        'glimmerhmm'), os.path.abspath(fasta), os.path.abspath(tmpdir), '-g']
    runSubprocess2(cmd, dir, log, glimmerRaw)

    # now convert to proper GFF3 format
    glimmer2gff3(glimmerRaw, output)

    return os.path.abspath(tmpdir)


def runGlimmerHMMTrained(fasta, training, dir, output):
    glimmerRaw = os.path.abspath(os.path.join(dir, 'glimmerHMM.output.raw'))
    cmd = ['perl', which_path('glimmhmm.pl'), which_path(
        'glimmerhmm'), os.path.abspath(fasta), os.path.abspath(training), '-g']
    runSubprocess2(cmd, dir, log, glimmerRaw)
    # now convert to proper GFF3 format
    glimmer2gff3(glimmerRaw, output)


def glimmer_run_check(Result, training, weights):
    if checkannotations(Result):
        log.info('Using existing GlimmerHMM results: {:}'.format(Result))
        return False
    if not checkannotations(training):
        log.info(
            'GlimmerHMM training failed, empty training set: {:}'.format(training))
        return False
    if weights < 1:
        log.info(
            'Skipping GlimmerHMM prediction as weight set to {:}'.format(weights))
        return False
    programs = ['trainGlimmerHMM', 'glimmerhmm', 'glimmhmm.pl']
    for x in programs:
        if not which_path(x):
            log.info(
                'GlimmerHMM failed, dependency not in $PATH: {:}'.format(x))
            return False
    return True


def dict2zff(scaffoldDict, GeneDict, output):
    # take funannotate dictionary convert to zff training format
    with open(output, 'w') as outfile:
        for k, v in natsorted(list(scaffoldDict.items())):
            outfile.write('>{:}\n'.format(k))
            for genes in v:
                gd = GeneDict.get(genes)
                for i in range(0, len(gd['ids'])):
                    for num, c in enumerate(gd['CDS'][i]):
                        if gd['strand'] == '+':
                            start = c[0]
                            end = c[1]
                        else:
                            start = c[1]
                            end = c[0]
                        if num == 0:
                            outfile.write('Einit\t{:}\t{:}\t{:}\n'.format(
                                start, end, gd['ids'][i]))
                        elif num == len(gd['CDS'][i])-1:
                            outfile.write('Eterm\t{:}\t{:}\t{:}\n'.format(
                                start, end, gd['ids'][i]))
                        else:
                            outfile.write('Exon\t{:}\t{:}\t{:}\n'.format(
                                start, end, gd['ids'][i]))


def zff2gff3(input, fasta, output):
    '''
    >scaffold_40
    Einit   7104    7391    -       14.809  0       0       1       scaffold_40-snap.1
    Eterm   6728    7039    -       1.974   0       0       2       scaffold_40-snap.1
    Einit   8935    9070    +       9.578   0       1       0       scaffold_40-snap.2
    Exon    9119    9206    +       10.413  2       2       0       scaffold_40-snap.2
    Exon    9254    9389    +       21.529  1       0       2       scaffold_40-snap.2
    Eterm   9439    10128   +       42.769  0       0       0       scaffold_40-snap.2
    Einit   11784   12139   -       38.847  0       2       2       scaffold_40-snap.3
    Eterm   11185   11761   -       72.324  1       0       0       scaffold_40-snap.3
    Einit   13191   13250   -       7.662   0       0       1       scaffold_40-snap.4
    Eterm   12498   13019   -       63.296  0       0       1       scaffold_40-snap.4
    Einit   16359   16608   +       41.592  0       1       2       scaffold_40-snap.5
    Exon    16628   16712   +       13.780  2       2       0       scaffold_40-snap.5
    Exon    16795   17012   +       26.393  1       1       1       scaffold_40-snap.5
    Eterm   17224   17381   +       8.331   2       0       2       scaffold_40-snap.5
    >scaffold_41
    Exon    65      951     -       169.146 1       1       0       scaffold_41-snap.1
    '''
    # need to load/generate a funannotate dictionary, then output to gff3 format
    Genes = {}
    contig = ''
    with open(input, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('#') or line.startswith('\n'):
                continue
            elif line.startswith('>'):
                contig = line[1:]
            else:
                feature, start, end, strand, score, fiveo, threeo, phase, ID = line.split(
                    '\t')
                start = int(start)
                end = int(end)
                # phase = int(phase)
                phase = '?'  # phase in GFF3 doesn't seem to be accurate, so guess it by translation of all 3 frames
                if not ID in Genes:
                    Genes[ID] = {'name': None, 'type': 'mRNA', 'transcript': [], 'cds_transcript': [], 'protein': [], '5UTR': [[]], '3UTR': [[]],
                                 'codon_start': [[]], 'ids': [ID+'-T1'], 'CDS': [[(start, end)]], 'mRNA': [[(start, end)]], 'strand': strand,
                                 'location': (start, end), 'contig': contig, 'product': [[]], 'source': 'snap', 'phase': [[phase]],
                                 'db_xref': [[]], 'EC_number': [[]], 'gene_synonym': [], 'go_terms': [[]], 'note': [[]], 'partialStart': [[]], 'partialStop': [[]], 'pseudo': False}
                else:
                    Genes[ID]['CDS'][0].append((start, end))
                    Genes[ID]['mRNA'][0].append((start, end))
                    Genes[ID]['phase'][0].append(phase)
                    if start < Genes[ID]['location'][0]:
                        Genes[ID]['location'] = (
                            start, Genes[ID]['location'][1])
                    if end > Genes[ID]['location'][1]:
                        Genes[ID]['location'] = (Genes[ID]['location'][0], end)
    # translate, check partial, etc
    SeqRecords = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    for k, v in list(Genes.items()):
        i = 0
        if v['strand'] == '+':
            sortedExons = sorted(v['mRNA'][i], key=lambda tup: tup[0])
            sortedCDS = sorted(v['CDS'][i], key=lambda tup: tup[0])
        else:
            sortedExons = sorted(
                v['mRNA'][i], key=lambda tup: tup[0], reverse=True)
            sortedCDS = sorted(
                v['CDS'][i], key=lambda tup: tup[0], reverse=True)
        Genes[k]['mRNA'][i] = sortedExons
        mrnaSeq = getSeqRegions(SeqRecords, v['contig'], sortedExons)
        Genes[k]['transcript'].append(mrnaSeq)

        # get the codon_start by getting first CDS phase + 1
        indexStart = [x for x, y in enumerate(
            v['CDS'][i]) if y[0] == sortedCDS[0][0]]
        cdsSeq = getSeqRegions(SeqRecords, v['contig'], sortedCDS)
        Genes[k]['cds_transcript'].append(cdsSeq)
        Genes[k]['CDS'][i] = sortedCDS
        protSeq, codon_start = (None,)*2
        if '?' in v['phase'][i]:  # dont know the phase -- malformed GFF3, try to find best CDS
            translateResults = []
            for y in [1, 2, 3]:
                protSeq = translate(cdsSeq, v['strand'], y-1)
                numStops = protSeq.count('*')
                if protSeq[-1] == '*':
                    numStops -= 1
                translateResults.append((y, numStops, protSeq))
            sortedResults = sorted(translateResults, key=lambda tup: tup[1])
            codon_start = sortedResults[0][0]
            protSeq = sortedResults[0][2]
        else:
            codon_start = int(v['phase'][i][indexStart[0]]) + 1
            # translate and get protein sequence
            protSeq = translate(cdsSeq, v['strand'], codon_start-1)
        Genes[k]['codon_start'][i] = codon_start
        if protSeq:
            Genes[k]['protein'].append(protSeq)
            if protSeq.endswith('*'):
                Genes[k]['partialStop'][i] = False
            else:
                Genes[k]['partialStop'][i] = True
            if codon_start == 1 and protSeq.startswith('M'):
                Genes[k]['partialStart'][i] = False
            else:
                Genes[k]['partialStart'][i] = True

    # now write to GFF3
    dict2gff3(Genes, output)


def cq_run_check(cqResult, bam, stringtie, weight):
    if checkannotations(cqResult):
        log.info('Using existing CodingQuarry results: {:}'.format(cqResult))
        return False
    if weight < 1:
        log.info(
            'Skipping CodingQuarry prediction as weight set to {:}'.format(weight))
        return False
    if not bam and not stringtie:
        log.info('Skipping CodingQuarry as there are no RNA-seq data')
        return False
    # check if dependencies installed
    programs = []
    if stringtie and checkannotations(stringtie):
        programs = ['CodingQuarry']
    elif bam and checkannotations(bam):
        programs = ['stringtie', 'CodingQuarry']
    for x in programs:
        if not which_path(x):
            log.info(
                'CodingQuarry failed, dependency not in $PATH: {:}'.format(x))
            return False
    # if you get here should be good
    return True


def snap_run_check(snapResult, training, weight):
    if checkannotations(snapResult):
        log.info('Using existing snap results: {:}'.format(snapResult))
        return False
    if not checkannotations(training):
        log.info(
            'Snap training failed, empty training set: {:}'.format(training))
        return False
    if weight < 1:
        log.info(
            'Skipping snap prediction as weight set to {:}'.format(weight))
        return False
    programs = ['fathom', 'snap', 'forge', 'hmm-assembler.pl']
    for x in programs:
        if not which_path(x):
            log.info('Snap failed, dependency not in $PATH: {:}'.format(x))
            return False
    return True


def runSnap(fasta, gff3, minintron, maxintron, dir, output):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    '''
    wrapper to run Snap training followed by prediction
    input is GFF3 format high quality models, i.e. from PASA/transdecoder
    output is standard GFF3 format
    '''
    from collections import OrderedDict
    snapHMM = os.path.join(dir, 'snap-trained.hmm')
    snapRaw = os.path.join(dir, 'snap-prediction.zff')
    if not checkannotations(snapRaw):
        # generate training directory ontop of the dir that is passed
        tmpdir = os.path.join(dir, 'snaptrain')
        if os.path.isdir(tmpdir):
            SafeRemove(tmpdir)
        os.makedirs(tmpdir)

        # load gff3 into dictionary
        Genes = {}
        Genes = gff2dict(os.path.abspath(gff3), os.path.abspath(fasta), Genes)
        scaff2genes = {}

        # sort the dictionary
        def _sortDict(d):
            return (d[1]['contig'], d[1]['location'][0])
        sGenes = sorted(iter(Genes.items()), key=_sortDict)
        sortedGenes = OrderedDict(sGenes)
        scaff2genes = {}
        for k, v in list(sortedGenes.items()):
            if not v['contig'] in scaff2genes:
                scaff2genes[v['contig']] = [k]
            else:
                scaff2genes[v['contig']].append(k)

        # get only scaffolds that have gene models for training
        log.debug('{:} gene models to train snap on {:} scaffolds'.format(
            len(sGenes), len(scaff2genes)))
        trainingFasta = os.path.join(dir, 'snap-training.scaffolds.fasta')
        with open(trainingFasta, 'w') as outfile:
            with open(os.path.abspath(fasta), 'r') as infile:
                for title, seq in SimpleFastaParser(infile):
                    if title in list(scaff2genes.keys()):
                        outfile.write('>{:}\n{:}\n'.format(
                            title, softwrap(seq)))

        # convert to ZFF format
        origzff = os.path.join(dir, 'snap.training.zff')
        dict2zff(scaff2genes, Genes, origzff)

        # now run SNAP training
        cmd = ['fathom', os.path.abspath(origzff), os.path.abspath(
            trainingFasta), '-categorize', '1000', '-min-intron', str(minintron), '-max-intron', str(maxintron)]
        runSubprocess(cmd, tmpdir, log)

        cmd = ['fathom', 'uni.ann', 'uni.dna', '-export', '1000', '-plus']
        runSubprocess(cmd, tmpdir, log)

        cmd = ['forge', 'export.ann', 'export.dna']
        runSubprocess(cmd, tmpdir, log)

        cmd = ['perl', which_path('hmm-assembler.pl'), 'snap-trained', tmpdir]
        runSubprocess2(cmd, '.', log, snapHMM)

        # now run SNAP prediction
        cmd = ['snap', os.path.abspath(snapHMM), os.path.abspath(fasta)]
        runSubprocess2(cmd, '.', log, snapRaw)

    # convert zff to proper gff3
    zff2gff3(snapRaw, fasta, output)

    return os.path.abspath(snapHMM)


def runSnapTrained(fasta, hmm, dir, output):
    snapRaw = os.path.join(dir, 'snap-prediction.zff')
    # now run SNAP prediction
    cmd = ['snap', os.path.abspath(hmm), os.path.abspath(fasta)]
    runSubprocess2(cmd, '.', log, snapRaw)
    # convert zff to proper gff3
    zff2gff3(snapRaw, fasta, output)


def MemoryCheck():
    import psutil
    mem = psutil.virtual_memory()
    RAM = int(mem.total)
    return round(RAM / 1024000000)


def systemOS():
    if sys.platform == 'darwin':
        system_os = 'MacOSX ' + platform.mac_ver()[0]
    elif sys.platform == 'linux':
        linux_version = distro.linux_distribution()
        system_os = linux_version[0] + ' '+linux_version[1]
    else:
        system_os = sys.platform
    return system_os


def SystemInfo():
    system_os = systemOS()
    python_vers = str(
        sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])
    log.info("OS: %s, %i cores, ~ %i GB RAM. Python: %s" %
             (system_os, multiprocessing.cpu_count(), MemoryCheck(), python_vers))


def runtRNAscan(input, tmpdir, output, cpus=1, precalc=False):
    tRNAout = os.path.join(tmpdir, 'tRNAscan.out')
    tRNAlenOut = os.path.join(tmpdir, 'tRNAscan.len-filtered.out')
    if not precalc:
        if os.path.isfile(tRNAout):  # tRNAscan can't overwrite file, so check
            os.remove(tRNAout)
        cmd = ['tRNAscan-SE', '-o', tRNAout, '--thread', str(cpus), input]
        log.debug(' '.join(cmd))
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            log.error('CMD ERROR: {}'.format(' '.join(cmd)))
        if stdout:
            log.debug(stdout.decode("utf-8"))
        if stderr:
            log.debug(stderr.decode("utf-8"))
    else:
        shutil.copyfile(precalc, tRNAout)
    if not checkannotations(tRNAout):
        log.info('tRNAscan-SE seems to have failed, check logfile for error. You can pass precalculated results to --trnascan')
        return False
    # enforce NCBI length rules
    with open(tRNAlenOut, 'w') as lenOut:
        with open(tRNAout, 'r') as infile:
            for line in infile:
                if line.startswith('Sequence') or line.startswith('Name') or line.startswith('--------'):
                    lenOut.write('%s' % line)
                else:
                    cols = line.split('\t')
                    start = cols[2]
                    end = cols[3]
                    if int(start) < int(end):
                        length = abs(int(end) - int(start))
                    else:
                        length = abs(int(start) - int(end))
                    if length < 50 or length > 150:
                        continue
                    else:
                        lenOut.write('%s' % line)

    # now convert to GFF3
    trna2gff = os.path.join(parentdir, 'aux_scripts', 'trnascan2gff3.pl')
    with open(output, 'w') as out:
        subprocess.call(['perl', trna2gff, '--input', tRNAlenOut], stdout=out)
    return True


def runtbl2asn(folder, template, discrepency, organism, isolate, strain, parameters, version):
    '''
    function to run NCBI tbl2asn
    '''
    # get funannotate version
    fun_version = get_version()
    # input should be a folder
    if not os.path.isdir(folder):
        log.error("tbl2asn error: %s is not a directory, exiting" % folder)
        sys.exit(1)
    # based on organism, isolate, strain, construct meta info for -j flag
    if not organism:
        log.error("tbl2asn error: organism not specified")
        sys.exit(1)
    meta = "[organism=" + organism + "]"
    if isolate:
        isolate_meta = "[isolate=" + isolate + "]"
        meta = meta + " " + isolate_meta
    if strain:
        strain_meta = "[strain=" + strain + "]"
        meta = meta + " " + strain_meta
    cmd = ['tbl2asn', '-y', '"Annotated using '+fun_version+'"', '-N',
           str(version), '-p', folder, '-t', template, '-M', 'n', '-Z', discrepency, '-j', '"'+meta+'"', '-V', 'b', '-c', 'fx', '-T', '-a', 'r10u']
    # check for custom parameters
    if parameters:
        params = parameters.split(' ')
        cmd = cmd + params
    runSubprocess(cmd, '.', log)
    return ' '.join(cmd)


def gb2smurf(input, prot_out, smurf_out):
    with open(smurf_out, 'w') as smurf:
        with open(prot_out, 'w') as proteins:
            with open(input, 'r') as gbk:
                SeqRecords = SeqIO.parse(gbk, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        name = re.sub('[^0-9]', '', record.name)
                        if f.type == "CDS":
                            proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], softwrap(
                                f.qualifiers['translation'][0].rstrip('*'))))
                            locus_tag = f.qualifiers.get(
                                "locus_tag", ["No ID"])[0]
                            product_name = f.qualifiers.get(
                                "product", ["No Description"])[0]
                            mystart = f.location.start
                            myend = f.location.end
                            strand = f.location.strand
                            if strand == 1:
                                smurf.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip(
                                    "0"), int(mystart), int(myend), product_name))
                            else:
                                smurf.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip(
                                    "0"), int(myend), int(mystart), product_name))


def GAGprotClean(input, output):
    '''
    gag.py v1 had headers like:
    >>evm.model.Contig100.1 protein
    gag.py v2 has headers like:
    >protein|evm.model.scaffold_1.169 ID=evm.model.scaffold_1.169|Parent=evm.TU.scaffold_1.169|Name=EVM%20prediction%20scaffold_1.169
    '''
    with open(output, 'w') as outfile:
        with open(input, 'ru') as infile:
            for rec in SeqIO.parse(infile, 'fasta'):
                if rec.id.startswith('protein|'):
                    ID = rec.id.replace('protein|', '').split(' ')[0]
                else:
                    ID = rec.id.split(' ')[0]
                rec.id = ID
                rec.name = ''
                rec.description = ''
                SeqIO.write(rec, outfile, 'fasta')


def OldRemoveBadModels(proteins, gff, length, repeats, BlastResults, tmpdir, output):
    # first run bedtools to intersect models where 90% of gene overlaps with repeatmasker region
    repeat_temp = os.path.join(tmpdir, 'genome.repeats.to.remove.gff')
    cmd = ['bedtools', 'intersect', '-f', '0.9', '-a', gff, '-b', repeats]
    runSubprocess2(cmd, '.', log, repeat_temp)
    # now remove those proteins that do not have valid starts, less then certain length, and have internal stops
    remove = []
    reason = {}
    # parse the results from bedtools and add to remove list
    with open(repeat_temp, 'r') as input:
        for line in input:
            if "\tgene\t" in line:
                ninth = line.split('ID=')[-1]
                ID = ninth.split(";")[0]
                remove.append(ID)
                if not ID in reason:
                    reason[ID] = 'remove_reason=repeat_overlap;'
    # parse the results from BlastP search of transposons
    with open(BlastResults, 'r') as input:
        for line in input:
            col = line.split('\t')
            remove.append(col[0])
            if not col[0] in reason:
                ID = col[0].replace('evm.model.', 'evm.TU.')
                reason[ID] = 'remove_reason=repeat_match;'
            else:
                ID = col[0].replace('evm.model.', 'evm.TU.')
                reason[ID] = 'remove_reason=repeat_overalap|repeat_match;'

    # I'm only seeing these models with GAG protein translations, so maybe that is a problem? skip enforcing start with M
    with open(proteins, 'r') as input:
        SeqRecords = SeqIO.parse(input, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)[:-1]
            ID = rec.id.replace('evm.model.', 'evm.TU.')
            if len(Seq) < int(length):
                remove.append(ID)
                if not ID in reason:
                    reason[ID] = 'remove_reason=seq_too_short;'
            if 'XX' in Seq:
                remove.append(ID)
                if not rec.id in reason:
                    reason[ID] = 'remove_reason=model_span_gap;'
    remove = [w.replace('evm.TU.', '') for w in remove]
    remove = [w.replace('evm.model.', '') for w in remove]
    remove = set(remove)
    if len(remove) > 0:
        remove_match = re.compile(r'\b\evm.(.*?:%s)[\.;]\b' % '|'.join(remove))
        with open(output, 'w') as out:
            with open(os.path.join(tmpdir, 'bad_models.gff'), 'w') as out2:
                with open(gff, 'r') as GFF:
                    for line in GFF:
                        if '\tstart_codon\t' in line:
                            continue
                        if '\tstop_codon\t' in line:
                            continue
                        matchLine = remove_match.search(line)
                        if not matchLine:
                            # remove the Name attribute as it sticks around in GBK file
                            line = re.sub(';Name=.*$', ';', line)
                            out.write(line)
                        else:
                            # print matchLine.group()
                            # print line
                            if "\tgene\t" in line:
                                bad_ninth = line.split('ID=')[-1]
                                bad_ID = bad_ninth.split(";")[0]
                                bad_reason = reason.get(bad_ID)
                                if bad_reason:
                                    line = line.replace(
                                        '\n', ';'+bad_reason+'\n')
                                    # print bad_reason
                                else:
                                    log.debug(
                                        "%s was removed in removeBadModels function for unknown reason, please check manually" % bad_ID)
                                    line = line.replace(
                                        '\n', ';remove_reason=unknown;\n')
                                    # print 'uknown'
                            out2.write(line)
    else:  # if nothing to remove, just print out GFF
        with open(output, 'w') as out:
            with open(gff, 'r') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue
                    # remove the Name attribute as it sticks around in GBK file
                    line = re.sub(';Name=.*$', ';', line)
                    out.write(line)


def RemoveBadModels(proteins, gff, length, repeats, BlastResults, tmpdir, methods, output):
    reason = {}
    tooShort = 0
    repeat = 0
    gapspan = 0
    if 'overlap' in methods:
        # first run bedtools to intersect models where 90% of gene overlaps with repeatmasker region
        repeat_temp = os.path.join(tmpdir, 'genome.repeats.to.remove.gff')
        gffSorted = os.path.abspath(gff)+'.sorted.gff'
        bedSorted = os.path.abspath(repeats)+'.sorted.bed'
        #sortBedproper(repeats, bedSorted)
        cmd1 = ['bedtools', 'sort', '-i', repeats]
        with open(bedSorted, 'w') as bedout:
            subprocess.call(cmd1, stdout=bedout)
        #sortGFFproper(gff, gffSorted)
        cmd2 = ['bedtools', 'sort', '-i', gff]
        with open(gffSorted, 'w') as gffout:
            subprocess.call(cmd2, stdout=gffout)
        cmd = ['bedtools', 'intersect', '-sorted', '-f', '0.9', '-a', gffSorted, '-b', bedSorted]
        runSubprocess2(cmd, '.', log, repeat_temp)
        # parse the results from bedtools and add to remove list
        with open(repeat_temp, 'r') as input:
            for line in input:
                if "\tgene\t" in line:
                    ninth = line.split('ID=')[-1]
                    ID = ninth.split(";")[0]
                    if not ID in reason:
                        reason[ID] = 'remove_reason=repeat_overlap;'
    if 'blast' in methods:
        # parse the results from BlastP search of transposons
        with open(BlastResults, 'r') as input:
            for line in input:
                col = line.split('\t')
                if not col[0] in reason:
                    ID = col[0].replace('evm.model.', 'evm.TU.')
                    reason[ID] = 'remove_reason=repeat_match;'
                else:
                    ID = col[0].replace('evm.model.', 'evm.TU.')
                    reason[ID] = 'remove_reason=repeat_overlap|repeat_match;'
    # always do these checks
    # Look for models that are too short
    with open(proteins, 'r') as input:
        SeqRecords = SeqIO.parse(input, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)[:-1]
            ID = rec.id.replace('evm.model.', 'evm.TU.')
            if len(Seq) < int(length):
                if not ID in reason:
                    reason[ID] = 'remove_reason=seq_too_short;'
            if 'XX' in Seq:
                if not rec.id in reason:
                    reason[ID] = 'remove_reason=model_span_gap;'
    # now read the EVM gene models in Blocks so you can parse gene ID
    numTotal = len(reason)
    for k, v in reason.items():
        if 'model_span_gap' in v:
            gapspan += 1
        elif 'seq_too_short' in v:
            tooShort += 1
        else:
            repeat += 1
    if numTotal > 0:
        log.info("Found {:,} gene models to remove: {:,} too short; {:,} span gaps; {:,} transposable elements".format(
            numTotal, tooShort, gapspan, repeat))
        with open(output, 'w') as out:
            with open(os.path.join(tmpdir, 'bad_models.gff'), 'w') as out2:
                with open(gff, 'r') as GFF:
                    for gene_model in readBlocks(GFF, '\n'):
                        if len(gene_model) > 1:
                            if gene_model[0].startswith('\n'):
                                ID = gene_model[1].split(
                                    'ID=')[-1].split(';')[0]
                            else:
                                ID = gene_model[0].split(
                                    'ID=')[-1].split(';')[0]
                            if ID in reason:
                                out2.write('#%s removed; %s\n' %
                                           (ID, reason.get(ID)))
                                for line in gene_model:
                                    if not line.startswith('\n'):
                                        out2.write('%s' % (line))
                            else:
                                for line in gene_model:
                                    # remove the Name attribute as it sticks around in GBK file
                                    line = re.sub(';Name=.*$', ';', line)
                                    out.write('%s' % (line))


def CleantRNAtbl(GFF, TBL, output):
    # clean up genbank tbl file from gag output
    # try to read through GFF file, make dictionary of tRNA genes and products
    TRNA = {}
    matches = []
    with open(GFF, 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            line = line.replace('\n', '')
            scaffold, source, feature, start, end, score, orientation, phase, info = line.split(
                '\t')
            if feature == 'tRNA':
                ID = info.split(';')[0].replace('ID=', '')
                ID = ID.replace('-T1', '')
                product = info.split('product=')[-1]
                TRNA[ID] = product
                matches.append(product)
    matches = set(matches)
    tRNAmatch = re.compile(r'\t\t\tproduct\t%s\n' % '|'.join(matches))
    with open(output, 'w') as out:
        with open(TBL, 'r') as input:
            for line in input:
                if line.startswith('\t\t\tlocus_tag\t'):
                    out.write(line)
                    geneID = line.split('locus_tag\t')[-1].replace('\n', '')
                    if geneID in TRNA:
                        CurrentProduct = TRNA.get(geneID)
                        if 'tRNA-Xxx' == CurrentProduct:
                            out.write("\t\t\tpseudo\n")
                elif line.startswith("\t\t\tproduct\ttRNA-Xxx"):
                    out.write(line)
                    out.write("\t\t\tpseudo\n")
                    next(input)
                    next(input)
                elif tRNAmatch.search(line):
                    out.write(line)
                    next(input)
                    next(input)
                else:  # otherwise just write line
                    out.write(line)


def getFailedProductNames(input, GeneDict):
    # input is NCBI tbl2asn discrepency report, parse to get suspect product names
    failed = {}
    with open(input, 'r') as discrep:
        for block in readBlocks(discrep, 'DiscRep_'):
            if 'DiscRep_SUB:SUSPECT_PRODUCT_NAMES::' in block[0]:
                reason = []
                for item in block:
                    if item.startswith('DiscRep_SUB:'):
                        bad = item.split('::')[-1].rstrip()
                        if 'features' in bad.lower():
                            bad = bad.split('features ')[-1]
                        reason.append(bad)
                    elif item.startswith('genome:'):
                        gene = item.split('\t')[-1].strip()
                        if gene.startswith('DiscRep'):
                            continue
                        if gene in GeneDict:
                            hit = GeneDict.get(gene)
                            if not hit[0] in failed:
                                failed[hit[0]] = (hit[1], gene, reason)
    return failed


def ParseErrorReport(input, Errsummary, val, Discrep, output, keep_stops):
    errors = []
    gapErrors = []
    remove = []
    with open(Errsummary) as summary:
        for line in summary:
            if 'ERROR' in line:
                # there are probably other errors you are unaware of....
                if 'SEQ_DESCR.OrganismIsUndefinedSpecies' in line or 'SEQ_DESCR.BadOrgMod' in line or 'SEQ_FEAT.MissingTrnaAA' in line or 'SEQ_INST.TerminalNs' in line:
                    pass
                elif 'SEQ_FEAT.NoStop' in line:
                    if keep_stops:
                        pass
                    else:
                        err = line.split(" ")[-1].rstrip()
                        errors.append(err)
                elif 'SEQ_FEAT.FeatureBeginsOrEndsInGap' in line:
                    err = line.split(" ")[-1].rstrip()
                    gapErrors.append(err)
                else:
                    err = line.split(" ")[-1].rstrip()
                    errors.append(err)
    # parse the discrepency report and look for overlapping genes, so far, all have been tRNA's in introns, so just get those for now.
    with open(Discrep, 'r') as discrep:
        # process discrepency report into blocks, then look for block headers where overlapping genes are, remove only tRNA models right now
        for block in readBlocks(discrep, 'DiscRep_'):
            if 'DiscRep_ALL:OVERLAPPING_GENES::' in block[0] or 'DiscRep_SUB:RNA_CDS_OVERLAP::' in block[0]:
                for item in block:
                    if item.startswith('genome:tRNA'):
                        gene = item.split('\t')[-1].replace('\n', '')
                        if gene.startswith('DiscRep'):
                            continue
                        tRNA = gene + '_tRNA'
                        exon = gene + '_exon'
                        remove.append(gene)
                        remove.append(tRNA)
                        remove.append(exon)
            if 'DiscRep_ALL:FIND_OVERLAPPED_GENES::' in block[0]:
                for item in block:
                    gene = item.split('\t')[-1].replace('\n', '')
                    if gene.startswith('DiscRep'):
                        continue
                    tRNA = gene + '_tRNA'
                    exon = gene + '_exon'
                    remove.append(gene)
                    remove.append(tRNA)
                    remove.append(exon)

    # there are no errors, then just remove stop/start codons and move on
    if len(errors) < 1 and len(remove) < 1:
        with open(output, 'w') as out:
            with open(input, 'r') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue
                    out.write(line)
    else:
        with open(val) as validate:
            for line in validate:
                if any(x in line for x in errors):
                    mRNA = line.split("ncbi|")[-1].replace(']', '').rstrip()
                    gene = mRNA.replace('evm.model', 'evm.TU')
                    exon = mRNA + '.exon'
                    mRNA = mRNA + ';'
                    remove.append(mRNA)
                    remove.append(gene)
                    remove.append(exon)
                # this is only picking up tRNAs right now, which "probably" is all that it needs to.....but u never know
                if any(x in line for x in gapErrors):
                    cols = line.split(' ')
                    if 'Gene:' in cols:
                        gene = line.split('Gene: ')[-1]
                        gene = gene.split(' ')[0]
                        tRNA = gene + '_tRNA'
                        exon = gene + '_exon'
                        remove.append(gene)
                        remove.append(tRNA)
                        remove.append(exon)
        # make sure no empty strings
        remove = list([_f for _f in remove if _f])
        remove = set(remove)
        remove_match = re.compile(r'\b(?:%s)+\b' % '|'.join(remove))
        with open(output, 'w') as out:
            with open(input, 'r') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue
                    if not remove_match.search(line):
                        if '\tgene\t' in line:
                            line = line.replace('Name=;', '')
                        out.write(line)


def antismash_version(input):
    # choose v4, v5 or v6 parser
    version = 4
    with open(input, 'r') as infile:
        for rec in SeqIO.parse(infile, 'genbank'):
            if 'structured_comment' in rec.annotations:
                if 'antiSMASH-Data' in rec.annotations['structured_comment']:
                    version = int(
                        rec.annotations['structured_comment']['antiSMASH-Data']['Version'].split('.')[0])
            break
    return version


def ParseAntiSmash(input, tmpdir, output, annotations):
    smash_version = antismash_version(input)
    log.info("Now parsing antiSMASH v{:} results, finding SM clusters".format(smash_version))
    BackBone = {}
    SMCOGs = {}
    bbSubType = {}
    bbDomains = {}
    smProducts = {}
    backboneCount = 0
    clusterCount = 0
    cogCount = 0
    # parse antismash genbank to get clusters in bed format and slice the record for each cluster prediction
    with open(output, 'w') as antibed:
        with open(input, 'r') as input:
            SeqRecords = SeqIO.parse(input, 'genbank')
            for rec_num,record in enumerate(SeqRecords):
                for f in record.features:
                    locusTag, ID, Parent = (None,)*3
                    if smash_version < 6:
                        baseName = 'Cluster'
                        if '_' in record.id:
                            try:
                                numericalContig = '{}_{}'.format(baseName, int(record.id.rsplit('_', 1)[-1]))
                            except ValueError:
                                if '.' in record.id:
                                    numericalContig = '{}_{}'.format(baseName, int(record.id.rsplit('.', 1)[0].rsplit('_', 1)[-1]))
                        else:  # just get the numbers
                            numericalContig = '{}_{}'.format(baseName, int(''.join(filter(str.isdigit, record.id))))
                    else:
                        numericalContig = 'Cluster'
                    # parse v4 differently than version 5
                    if smash_version == 4:
                        if f.type == "cluster":
                            clusterCount += 1
                            chr = record.id
                            start = f.location.nofuzzy_start
                            end = f.location.nofuzzy_end
                            clusternum = f.qualifiers.get(
                                "note")[0].replace("Cluster number: ", "")
                            antibed.write("%s\t%s\t%s\tCluster_%s\t0\t+\n" %
                                          (chr, start, end, clusternum))
                        Domains = []
                        if f.type == "CDS":
                            locusTag, ID, Parent = getID(f, f.type)
                            if not ID:
                                continue
                            ID = ID.replace('ncbi_', '')
                            if f.qualifiers.get('sec_met'):
                                for k, v in list(f.qualifiers.items()):
                                    if k == 'sec_met':
                                        for i in v:
                                            if i.startswith('Type:'):
                                                type = i.replace('Type: ', '')
                                                backboneCount += 1
                                                BackBone[ID] = type
                                            if i.startswith('NRPS/PKS subtype:'):
                                                subtype = i.replace(
                                                    'NRPS/PKS subtype: ', '')
                                                bbSubType[ID] = subtype
                                            if i.startswith('NRPS/PKS Domain:'):
                                                doms = i.replace(
                                                    'NRPS/PKS Domain: ', '')
                                                doms = doms.split('. ')[0]
                                                Domains.append(doms)
                                    bbDomains[ID] = Domains
                            for k, v in list(f.qualifiers.items()):
                                if k == 'note':
                                    for i in v:
                                        if i.startswith('smCOG:'):
                                            COG = i.replace('smCOG: ', '')
                                            COG = COG.split(' (')[0]
                                            SMCOGs[ID] = COG
                                            cogCount += 1
                                        elif not i.startswith('smCOG tree'):
                                            notes = i
                                            smProducts[ID] = notes
                    elif smash_version >= 5:
                        if f.type == "protocluster":
                            clusterCount += 1
                            chr = record.id
                            start = f.location.nofuzzy_start
                            # if '<' in start:
                            #    start = start.replace('<', '')
                            end = f.location.nofuzzy_end
                            # if '>' in end:
                            #    end = end.replace('>', '')
                            clusternum = int(f.qualifiers.get(
                                "protocluster_number")[0])
                            if smash_version >= 6:
                                antibed.write("{:}\t{:}\t{:}\t{:}_{:}\t0\t+\n".format(
                                    chr, start, end, numericalContig, clusternum))
                            else:
                                antibed.write("{:}\t{:}\t{:}\t{:}.{:}\t0\t+\n".format(
                                    chr, start, end, numericalContig, clusternum))
                        Domains = []
                        if f.type == "CDS":
                            locusTag, ID, Parent = getID(f, f.type)
                            if not ID:
                                continue
                            ID = ID.replace('ncbi_', '')
                            if f.qualifiers.get('NRPS_PKS'):
                                for k, v in list(f.qualifiers.items()):
                                    if k == 'NRPS_PKS':
                                        for i in v:
                                            if i.startswith('type:'):
                                                type = i.replace('type: ', '')
                                                backboneCount += 1
                                                BackBone[ID] = type
                                            if i.startswith('NRPS_PKS subtype:'):
                                                subtype = i.replace(
                                                    'NRPS_PKS subtype: ', '')
                                                bbSubType[ID] = subtype
                                            if i.startswith('Domain:'):
                                                doms = i.replace(
                                                    'Domain: ', '')
                                                doms = doms.split('. ')[0]
                                                Domains.append(doms)
                                    bbDomains[ID] = Domains
                            for k, v in list(f.qualifiers.items()):
                                if k == 'gene_functions':
                                    for i in v:
                                        if '(smcogs)' in i:
                                            COG = i.split(
                                                '(smcogs)')[-1].strip()
                                            COG = COG.split(' (')[0]
                                            SMCOGs[ID] = COG
                                            cogCount += 1
                                elif k == 'gene_kind':
                                    if 'biosynthetic' in v:
                                        backboneCount += 1

    # if smash_version == 4:
    log.info("Found %i clusters, %i biosynthetic enyzmes, and %i smCOGs predicted by antiSMASH" % (
        clusterCount, backboneCount, cogCount))

    # now generate the annotations to add to genome
    with open(annotations, 'w') as out:
        # add product annotations - use bbSubType --> BackBone
        for k, v in natsorted(list(BackBone.items())):
            ID = k
            if k in bbSubType:
                hit = bbSubType.get(k)
                if hit == 'NRPS':
                    hit = 'Nonribosomal Peptide Synthase (NRPS)'
                if hit == 'Type I Iterative PKS':
                    hit = 'Type I Iterative Polyketide synthase (PKS)'
            else:
                hit = v
            if hit == 'terpene':
                hit = 'terpene cyclase'
            elif hit == 'other':
                hit = 'putative secondary metabolism biosynthetic enzyme'
            elif hit == 'indole':
                hit = 'aromatic prenyltransferase (DMATS family)'
            elif hit == 'alkaloid' or hit == 'lignan' or hit == 'saccharide' or hit == 'polyketide':
                hit = 'putative ' + hit + ' biosynthetic cluster'
            elif hit == 'putative':
                hit = 'putative uncategorized biosynthetic cluster'
            elif '-' in hit:
                hit = 'putative ' + hit + ' biosynthetic cluster'
            if hit != 'none':
                out.write("%s\tproduct\t%s\n" % (ID, hit))
        # add annots from smProducts
        for k, v in list(smProducts.items()):
            ID = k
            if v != 'none' and not 'BLAST' in v:
                sys.stdout.write("%s\tproduct\t%s\n" % (ID, v))
        # add smCOGs into note section
        for k, v in list(SMCOGs.items()):
            ID = k
            if v != 'none':
                out.write("%s\tnote\t%s\n" % (ID, v))

    return bbDomains, bbSubType, BackBone


def GetClusterGenes(input, GFF, genome, annotations):
    # load clusters into InterLap
    interClust = bed2interlapNames(input)

    # load GFF3 into Dictionary
    Genes = {}
    Genes = gff2dict(GFF, genome, Genes)

    # loop through genes and check if in Clusters
    dictClusters = {}
    for k, v in natsorted(Genes.items()):
        if v['type'] == 'mRNA':
            if v['location'] in interClust[v['contig']]:
                best_hit = list(interClust[v['contig']].find(v['location']))[0]
                clusterName = best_hit[2]
                if not clusterName in dictClusters:
                    dictClusters[clusterName] = v['ids']
                else:
                    dictClusters[clusterName] += v['ids']
    # write the output file
    with open(annotations, 'w') as annotout:
        for k, v in list(dictClusters.items()):
            for i in v:
                annotout.write("%s\tnote\tantiSMASH:%s\n" % (i, k))

    return dictClusters


def splitFASTA(input, outputdir):
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    with open(input, 'r') as InputFasta:
        SeqRecords = SeqIO.parse(InputFasta, 'fasta')
        for record in SeqRecords:
            name = str(record.id)
            outputfile = os.path.join(outputdir, name+'.fa')
            with open(outputfile, 'w') as output:
                SeqIO.write(record, output, 'fasta')


def genomeStats(input):
    from Bio.SeqUtils import GC
    lengths = []
    GeeCee = []
    Genes = 0
    tRNA = 0
    Prots = 0
    locus_tag = ''
    organism = None
    isolate = None
    strain = None
    uniqueIso = None
    with open(input, 'r') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            lengths.append(len(record.seq))
            GeeCee.append(str(record.seq))
            organism = record.annotations['organism'].replace(
                ' Unclassified.', '')
            for f in record.features:
                if f.type == "source":
                    isolate = f.qualifiers.get("isolate", [None])[0]
                    strain = f.qualifiers.get("strain", [None])[0]
                if f.type == "CDS":
                    Prots += 1
                if f.type == "gene":
                    Genes += 1
                    if Genes == 1:
                        locus_tag = f.qualifiers.get("locus_tag")[
                            0].split('_')[0]
                if f.type == "tRNA":
                    tRNA += 1
    if strain:
        log.info("working on %s %s" % (organism, strain))
        uniqueIso = strain.replace(' ', '')
    elif isolate:
        log.info("working on %s %s" % (organism, isolate))
        uniqueIso = isolate.replace(' ', '')
    else:
        log.info("working on %s" % organism)
    GenomeSize = sum(lengths)
    LargestContig = max(lengths)
    ContigNum = len(lengths)
    AvgContig = int(round(GenomeSize / ContigNum))
    pctGC = round(GC("".join(GeeCee)), 2)

    # now get N50
    lengths.sort()
    nlist = []
    for x in lengths:
        nlist += [x]*x
    if len(nlist) % 2 == 0:
        medianpos = int(len(nlist) / 2)
        N50 = int((nlist[medianpos] + nlist[medianpos-1]) / 2)
    else:
        medianpos = int(len(nlist) / 2)
        N50 = int(nlist[medianpos])
    # return values in a list
    return [organism, uniqueIso, locus_tag, "{0:,}".format(GenomeSize)+' bp', "{0:,}".format(LargestContig)+' bp', "{0:,}".format(AvgContig)+' bp', "{0:,}".format(ContigNum), "{0:,}".format(N50)+' bp', "{:.2f}".format(pctGC)+'%', "{0:,}".format(Genes), "{0:,}".format(Prots), "{0:,}".format(tRNA)]


def MEROPS2dict(input):
    dict = {}
    with open(input, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                cols = line.split(' ')
                ID = cols[0].replace('>', '')
                family = cols[1].replace('\n', '')
                dict[ID] = family
    return dict


def getEggNogfromNote(input):
    dict = {}
    with open(input, 'r') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    try:
                        ID = f.qualifiers['locus_tag'][0]
                    except KeyError:
                        log.debug("%s has no locus_tag, skipping")
                        continue
                    for k, v in list(f.qualifiers.items()):
                        if k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith('EggNog:'):
                                    hit = i.replace('EggNog:', '')
                                    if not ID in dict:
                                        dict[ID] = hit
    return dict


def getStatsfromNote(input, word, Database):
    dict = {}
    meropsDict = MEROPS2dict(os.path.join(Database, 'merops.formatted.fa'))
    with open(input, 'r') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    try:
                        ID = f.qualifiers['locus_tag'][0]
                    except KeyError:
                        log.debug("%s has no locus_tag, skipping")
                        continue
                    for k, v in list(f.qualifiers.items()):
                        if k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith(word+':'):
                                    hit = i.replace(word+':', '')
                                    if hit.startswith('MER'):  # change to family name
                                        hit = meropsDict.get(hit)
                                    if not hit in dict:
                                        dict[hit] = [ID]
                                    else:
                                        dict[hit].append(ID)
    return dict


def getSMBackbones(input):
    dict = {'NRPS': 0, 'PKS': 0, 'Hybrid': 0}
    with open(input, 'r') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                if f.type == 'CDS':
                    product = f.qualifiers['product'][0]
                    if not product == 'hypothetical protein':
                        if product == "Hybrid PKS-NRPS":
                            dict['Hybrid'] += 1
                        if product == "Nonribosomal Peptide Synthase (NRPS)":
                            dict['NRPS'] += 1
                        if 'Polyketide synthase (PKS)' in product:
                            dict['PKS'] += 1
    return dict


def parseGOterms(input, folder, genome):
    with open(os.path.join(folder, 'associations.txt'), 'a') as assoc:
        with open(os.path.join(folder, genome+'.txt'), 'w') as terms:
            with open(input, 'r') as gbk:
                SeqRecords = SeqIO.parse(gbk, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        if f.type == 'CDS':
                            try:
                                ID = f.qualifiers['locus_tag'][0]
                            except KeyError:
                                log.debug("%s has no locus_tag, skipping")
                                continue
                            GOS = []
                            for k, v in list(f.qualifiers.items()):
                                if k == 'note':
                                    notes = v[0].split('; ')
                                    for i in notes:
                                        if i.startswith('GO'):
                                            go_term = i.split(' ')[1]
                                            GOS.append(go_term)
                            if GOS:
                                assoc.write("%s\t%s\n" % (ID, ";".join(GOS)))
                                terms.write("%s\n" % ID)


def getStatsfromDbxref(input, word):
    dict = {}
    with open(input, 'r') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    try:
                        ID = f.qualifiers['locus_tag'][0]
                    except KeyError:
                        log.debug("%s has no locus_tag, skipping")
                        continue
                    for k, v in list(f.qualifiers.items()):
                        if k == 'db_xref':
                            for i in v:
                                if i.startswith(word+':'):
                                    hit = i.replace(word+':', '')
                                    if not hit in dict:
                                        dict[hit] = [ID]
                                    else:
                                        dict[hit].append(ID)
    return dict


def getGBKannotation(input, Database):
    '''
    Function will loop through GBK file pulling out funannotate functional annotation
    and returning a list of dictionaries for each annotation class
    '''
    # convert merops on the fly, need database
    meropsDict = MEROPS2dict(os.path.join(Database, 'merops.formatted.fa'))
    SMs = {'NRPS': 0, 'PKS': 0, 'Hybrid': 0}
    pfams = {}
    iprs = {}
    nogs = {}
    cogs = {}
    merops = {}
    cazys = {}
    secreted = {}
    membrane = {}
    buscos = {}
    secmet = {}
    with open(input, 'r') as infile:
        for record in SeqIO.parse(infile, 'genbank'):
            for f in record.features:
                locusTag, ID, Parent = (None,)*3
                if f.type == 'CDS':
                    locusTag, ID, Parent = getID(f, f.type)
                    if not ID:
                        continue
                    product = f.qualifiers['product'][0]
                    if product == "Hybrid PKS-NRPS":
                        SMs['Hybrid'] += 1
                    if product == "Nonribosomal Peptide Synthase (NRPS)":
                        SMs['NRPS'] += 1
                    if 'Polyketide synthase (PKS)' in product:
                        SMs['PKS'] += 1
                    for k, v in list(f.qualifiers.items()):
                        if k == 'db_xref':
                            for i in v:
                                if i.startswith('PFAM:'):
                                    hit = i.replace('PFAM:', '')
                                    if not hit in pfams:
                                        pfams[hit] = [ID]
                                    else:
                                        pfams[hit].append(ID)
                                elif i.startswith('InterPro:'):
                                    hit = i.replace('InterPro:', '')
                                    if not hit in iprs:
                                        iprs[hit] = [ID]
                                    else:
                                        iprs[hit].append(ID)
                        if k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith('EggNog:'):
                                    hit = i.replace('EggNog:', '')
                                    if not ID in nogs:
                                        nogs[ID] = hit
                                elif i.startswith('BUSCO:'):
                                    hit = i.replace('BUSCO:', '')
                                    if not hit in buscos:
                                        buscos[hit] = [ID]
                                    else:
                                        buscos[hit].append(ID)
                                elif i.startswith('MEROPS:'):  # change to family name
                                    hit = i.replace('MEROPS:', '')
                                    hit = meropsDict.get(hit)
                                    if not hit in merops:
                                        merops[hit] = [ID]
                                    else:
                                        merops[hit].append(ID)
                                elif i.startswith('CAZy:'):
                                    hit = i.replace('CAZy:', '')
                                    if not hit in cazys:
                                        cazys[hit] = [ID]
                                    else:
                                        cazys[hit].append(ID)
                                elif i.startswith('COG:'):
                                    hit = i.replace('COG:', '')
                                    hits = hit.split(',')
                                    for x in hits:
                                        if not x in cogs:
                                            cogs[x] = [ID]
                                        else:
                                            cogs[x].append(ID)
                                elif i.startswith('SECRETED:'):
                                    hit = i.replace('SECRETED:', '')
                                    if not hit in secreted:
                                        secreted[hit] = [ID]
                                    else:
                                        secreted[hit].append(ID)
                                elif i.startswith('TransMembrane:'):
                                    hit = i.replace('TransMembrane:', '')
                                    if not hit in membrane:
                                        membrane[hit] = [ID]
                                    else:
                                        membrane[hit].append(ID)
                                elif i.startswith('antiSMASH:'):
                                    hit = i.replace('antiSMASH:', '')
                                    if not hit in secmet:
                                        secmet[hit] = [ID]
                                    else:
                                        secmet[hit].append(ID)
    return [pfams, iprs, nogs, buscos, merops, cazys, cogs, secreted, membrane, secmet, SMs]


def annotationtable(input, Database, HeaderNames, InterProDict, output):
    from collections import OrderedDict
    '''
    Function will create a tsv annotation table from GenBank file
    trying to capture all annotation in a parsable tsv file or
    something that could be imported into excel
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # convert merops on the fly, need database
    meropsDict = MEROPS2dict(os.path.join(Database, 'merops.formatted.fa'))
    # get note new/unique note names
    uniqueNotes = OrderedDict()
    for x in HeaderNames:
        if not x in ['BUSCO', 'CAZy', 'COG', 'EggNog', 'SECRETED', 'GO', 'MEROPS', 'TransMembrane']:
            uniqueNotes[x] = []
    # load genbank into funannotate dictionary (required as we need transcript/cds/etc)
    Genes = {}
    with open(input, 'r') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                gb_feature_add2dict(f, record, Genes)
    SeqRecords = SeqIO.to_dict(SeqIO.parse(input, 'genbank'))
    sGenes = natsorted(Genes.items(), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # input should be fully annotation GBK file from funannotate
    with open(output, 'w') as outfile:
        header = ['GeneID', 'TranscriptID', 'Feature', 'Contig', 'Start',
            'Stop', 'Strand', 'Name', 'Product', 'Alias/Synonyms', 'EC_number',
            'BUSCO', 'PFAM', 'InterPro', 'EggNog', 'COG', 'GO Terms',
            'Secreted', 'Membrane', 'Protease', 'CAZyme']
        header += uniqueNotes.keys()
        header += ['Notes', 'gDNA', 'mRNA', 'CDS-transcript', 'Translation']
        outfile.write('%s\n' % '\t'.join(header))
        for k,v in sortedGenes.items():
            for i in range(0,len(v['ids'])):
                # for each new feature, start with empty lists
                pfams = []
                iprs = []
                GOS = v['go_terms'][i]
                nogs = []
                cogs = []
                merops = []
                cazys = []
                secreted = []
                membrane = []
                therest = []
                buscos = []
                ecnum = []
                alias = []
                for key,value in uniqueNotes.items():
                    uniqueNotes[key] = []
                # now grab the data
                for y in v['db_xref'][i]:
                    if y.startswith('PFAM:'):
                        hit = y.replace('PFAM:', '')
                        pfams.append(hit)
                    elif y.startswith('InterPro:'):
                        hit = y.replace('InterPro:', '')
                        # look up description in dictionary
                        desc = InterProDict.get(hit)
                        iprs.append('{:} {:}'.format(hit, desc))
                for y in v['gene_synonym']:
                    alias.append(y)
                for y in v['EC_number'][i]:
                    ecnum.append(y)
                for y in v['note'][i]:
                    if y.startswith('EggNog:'):
                        hit = y.replace('EggNog:', '')
                        nogs.append(hit)
                    elif y.startswith('BUSCO:'):
                        hit = y.replace('BUSCO:', '')
                        buscos.append(hit)
                    elif y.startswith('MEROPS:'):  # change to family name
                        hit = y.replace('MEROPS:', '')
                        if hit in meropsDict:
                            hit = meropsDict.get(hit)
                            merops.append(hit)
                        else:
                            log.error("MEROPS database inconsistency: %s not found" % hit)
                    elif y.startswith('CAZy:'):
                        hit = y.replace('CAZy:', '')
                        cazys.append(hit)
                    elif y.startswith('COG:'):
                        hit = y.replace('COG:', '')
                        hits = hit.split(',')
                        for x in hits:
                            desc = x + ':'+ resources.COGS.get(x)
                            cogs.append(desc)
                    elif y.startswith('SECRETED:'):
                        hit = y.replace('SECRETED:', '')
                        secreted.append(hit)
                    elif y.startswith('TransMembrane:'):
                        hit = y.replace('TransMembrane:', '')
                        membrane.append(hit)
                    elif y.startswith(tuple(uniqueNotes.keys())):
                        try:
                            n = y.split(':')[0]
                            hit = y.split(':', 1)[1]
                            uniqueNotes[n].append(hit)
                        except IndexError:
                            hit = y
                            therest.append(hit)
                    else:  # capture everything else
                        hit = y
                        therest.append(hit)

                # bring together output
                result = [k, v['ids'][i], v['type'], v['contig'],
                          str(v['location'][0]), str(v['location'][1]),
                          v['strand'], v['name'],
                          v['product'][i],';'.join(alias),
                          ';'.join(ecnum),';'.join(buscos),
                          ';'.join(pfams),';'.join(iprs),
                          ';'.join(nogs),';'.join(cogs),
                          ';'.join(GOS),
                          ';'.join(secreted),
                          ';'.join(membrane),
                          ';'.join(merops),
                          ';'.join(cazys)
                          ]
                for key,value in uniqueNotes.items():
                    result.append(';'.join(value))
                gDNA = getSeqRegions(SeqRecords, v['contig'], [v['location']])
                try:
                    Transcript = str(v['transcript'][i])
                except IndexError:
                    if v['cds_transcript'][i]:
                        Transcript = str(v['cds_transcript'][i])
                    else:
                        print('{:} has no mrna or cds transcript'.format(k))
                        pass
                if v['type'] == 'mRNA':
                    CDSTranscript = str(v['cds_transcript'][i])
                    Protein = v['protein'][i]
                else:
                    CDSTranscript = ''
                    Protein = ''
                if v['strand'] == '-':
                    gDNA = RevComp(gDNA)
                    Transcript = RevComp(Transcript)
                    CDSTranscript = RevComp(CDSTranscript)
                result += [';'.join(therest), gDNA, Transcript,
                           CDSTranscript, Protein]
                # convert any None's to empty string
                result = ['' if x is None else x for x in result]
                # write to file
                outfile.write('%s\n' % '\t'.join(result))


def annotationtableOld(input, Database, output):
    '''
    Function will create a tsv annotation table from GenBank file
    trying to capture all annotation in a parsable tsv file or
    something that could be imported into excel
    '''
    # convert merops on the fly, need database
    meropsDict = MEROPS2dict(os.path.join(Database, 'merops.formatted.fa'))
    # input should be fully annotation GBK file from funannotate
    with open(output, 'w') as outfile:
        header = ['GeneID', 'Feature', 'Contig', 'Start', 'Stop', 'Strand', 'Name', 'Product', 'BUSCO', 'PFAM',
                  'InterPro', 'EggNog', 'COG', 'GO Terms', 'Secreted', 'Membrane', 'Protease', 'CAZyme', 'Notes', 'Translation']
        outfile.write('%s\n' % '\t'.join(header))
        for record in SeqIO.parse(input, 'genbank'):
            Contig = record.id
            for f in record.features:
                if f.type in ['tRNA', 'ncRNA', 'rRNA']:
                    ID = f.qualifiers['locus_tag'][0]
                    Start = f.location.nofuzzy_start
                    End = f.location.nofuzzy_end
                    strand = f.location.strand
                    if strand == 1:
                        Strand = '+'
                    elif strand == -1:
                        Strand = '-'
                    try:
                        Product = f.qualifiers['product'][0]
                    except KeyError:
                        Product = "None"
                    result = [ID, f.type, Contig, str(Start), str(
                        End), Strand, '', Product, '', '', '', '', '', '', '', '', '', '', '', '']
                    outfile.write('%s\n' % '\t'.join(result))
                if f.type == 'CDS':
                    ID = f.qualifiers['locus_tag'][0]
                    Start = f.location.nofuzzy_start
                    End = f.location.nofuzzy_end
                    strand = f.location.strand
                    if strand == 1:
                        Strand = '+'
                    elif strand == -1:
                        Strand = '-'
                    try:
                        Product = f.qualifiers['product'][0]
                    except KeyError:
                        Product = 'hypothetical protein'
                    try:
                        Name = f.qualifiers['gene'][0]
                    except KeyError:
                        Name = ''
                    try:
                        Translation = f.qualifiers['translation'][0]
                    except KeyError:
                        Translation = ''
                    pfams = []
                    iprs = []
                    GOS = []
                    nogs = []
                    cogs = []
                    merops = []
                    cazys = []
                    secreted = []
                    membrane = []
                    therest = []
                    buscos = []
                    for k, v in list(f.qualifiers.items()):
                        if k == 'db_xref':
                            for i in v:
                                if i.startswith('PFAM:'):
                                    hit = i.replace('PFAM:', '')
                                    pfams.append(hit)
                                elif i.startswith('InterPro:'):
                                    hit = i.replace('InterPro:', '')
                                    iprs.append(hit)
                        elif k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith('GO'):
                                    go_term = i.split(' ')[1]
                                    GOS.append(go_term)
                                elif i.startswith('EggNog:'):
                                    hit = i.replace('EggNog:', '')
                                    nogs.append(hit)
                                elif i.startswith('BUSCO:'):
                                    hit = i.replace('BUSCO:', '')
                                    buscos.append(hit)
                                elif i.startswith('MEROPS:'):  # change to family name
                                    hit = i.replace('MEROPS:', '')
                                    if hit in meropsDict:
                                        hit = meropsDict.get(hit)
                                        merops.append(hit)
                                    else:
                                        log.error(
                                            "MEROPS database inconsistency: %s not found" % hit)
                                elif i.startswith('CAZy:'):
                                    hit = i.replace('CAZy:', '')
                                    cazys.append(hit)
                                elif i.startswith('COG:'):
                                    hit = i.replace('COG:', '')
                                    hits = hit.split(',')
                                    for x in hits:
                                        desc = x + ':' + resources.COGS.get(x)
                                        cogs.append(desc)
                                elif i.startswith('SECRETED:'):
                                    hit = i.replace('SECRETED:', '')
                                    secreted.append(hit)
                                elif i.startswith('TransMembrane:'):
                                    hit = i.replace('TransMembrane:', '')
                                    membrane.append(hit)
                                else:  # capture everything else
                                    hit = i
                                    therest.append(hit)
                    result = [ID, 'CDS', Contig, str(Start), str(End), Strand, Name, Product, ';'.join(buscos), ';'.join(pfams), ';'.join(iprs), ';'.join(
                        nogs), ';'.join(cogs), ';'.join(GOS), ';'.join(secreted), ';'.join(membrane), ';'.join(merops), ';'.join(cazys), ';'.join(therest), Translation]
                    outfile.write('%s\n' % '\t'.join(result))


def ncbiCheckErrors(error, validation, genename, fixOut):
    ncbi_error = 0
    actual_error = 0
    with open(error, 'r') as errors:
        for line in errors:
            line = line.strip()
            if 'ERROR' in line:
                num = line.split(' ')[0]
                ncbi_error += int(num)
    # if errors in summary, then parse validation report, only get errors with gene names
    if ncbi_error > 0:
        # see if we can get the gene models that need to be fixed
        needFixing = {}
        with open(validation, 'r') as validationFile:
            for line in validationFile:
                line = line.strip()
                if line.startswith('ERROR') and genename in line:
                    actual_error += 1
                    parts = line.split(' ')
                    for x in parts:
                        if genename in x:
                            ID = x.split('|')[-1]
                    if '-' in ID:
                        ID = ID.split('-')[0]
                    reason = line.split(' FEATURE:')[0]
                    reason = reason.split('] ')[-1]
                    if not ID in needFixing:
                        needFixing[ID] = reason
        if actual_error > 0:
            log.info("There are %i gene models that need to be fixed." %
                     actual_error)
            print('-------------------------------------------------------')
            with open(fixOut, 'w') as fix:
                fix.write('#GeneID\tError Message\n')
                for k, v in natsorted(list(needFixing.items())):
                    fix.write('%s\t%s\n' % (k, v))
                    print(('%s\t%s' % (k, v)))
    return actual_error


def convert2counts(input):
    import pandas as pd
    Counts = []
    for i in range(0, len(input)):
        dict = {}
        for k, v in list(input[i].items()):
            dict[k] = len(v)
        Counts.append(dict)
    df = pd.DataFrame(Counts)
    df.fillna(0, inplace=True)  # fill in zeros for missing data
    return df


def gb2proteinortho(input, folder, name):
    gffOut = os.path.join(folder, name+'.gff')
    FastaOut = os.path.join(folder, name+'.faa')
    Transcripts = os.path.join(folder, name+'.transcripts.fa')
    genes = {}
    with open(input, 'r') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                gb_feature_add2dict(f, record, genes)
    # now output the files you need
    with open(gffOut, 'w') as gff:
        with open(FastaOut, 'w') as fasta:
            with open(Transcripts, 'w') as transcripts:
                for k, v in natsorted(list(genes.items())):
                    if v['type'] == 'mRNA':
                        for i, item in enumerate(v['ids']):
                            transcripts.write(">{:} {:} codon_start={:} strand={:}\n{:}\n".format(
                                item, k, v['codon_start'][i], v['strand'], v['cds_transcript'][i]))
                            fasta.write(">%s %s\n%s\n" %
                                        (item, k, v['protein'][i]))
                            gff.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(
                                v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], item, k, v['product'][i]))


def drawStackedBar(panda, type, labels, ymax, output, colors=False):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    import seaborn as sns
    import numpy as np
    from funannotate.stackedBarGraph import StackedBarGrapher as StackedBarGrapher
    # stackedbargraph from summary data
    SBG = StackedBarGrapher()
    # labels
    d_labels = panda.index.values
    # y-ticks
    ticks = np.linspace(0, ymax, 6)
    ticks = list(ticks)
    nums = [int(x) for x in ticks]
    vals = [str(x) for x in nums]
    yticks = [nums, vals]
    # colors
    if not colors:
        color_palette = sns.hls_palette(
            len(panda.columns), l=.4, s=.8).as_hex()
        color_palette = [str(x).upper() for x in color_palette]
    else:
        color_palette = colors
    # set up plot
    sns.set_style('darkgrid')
    sns.set_context('paper')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    YLabel = "Number of "+type
    SBG.stackedBarPlot(ax, panda, color_palette, xLabels=panda.index.values,
                       endGaps=True, gap=0.25, xlabel="Genomes", ylabel=YLabel, yTicks=yticks)
    plt.title(type+" summary")
    # get the legend
    legends = []
    i = 0
    for column in panda.columns:
        legends.append(mpatches.Patch(
            color=color_palette[i], label=panda.columns.values[i] + ": " + labels.get(panda.columns.values[i])))
        i += 1
    lgd = ax.legend(handles=legends, fontsize=6, loc='upper left',
                    bbox_to_anchor=(1.02, 1), borderaxespad=0)
    plt.ylim([0, ymax])
    # set the font size - i wish I knew how to do this proportionately.....but setting to something reasonable.
    for item in ax.get_xticklabels():
        item.set_fontsize(8)
    # setup the plot
    fig.subplots_adjust(bottom=0.4)
    fig.savefig(output, format='pdf', bbox_extra_artists=(
        lgd,), bbox_inches='tight')
    plt.close(fig)


def drawHeatmap(df, color, output, labelsize, annotate):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    # get size of table
    width = len(df.columns) / 2
    height = len(df.index) / 4
    fig, ax = plt.subplots(figsize=(width, height))
    cbar_ax = fig.add_axes(shrink=0.4)
    if annotate:
        sns.heatmap(df, linewidths=0.5, cmap=color, ax=ax,
                    fmt="d", annot_kws={"size": 4}, annot=True)
    else:
        sns.heatmap(df, linewidths=0.5, cmap=color, ax=ax, annot=False)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    for item in ax.get_xticklabels():
        item.set_fontsize(8)
    for item in ax.get_yticklabels():
        item.set_fontsize(int(labelsize))
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)


def donutplot(df, LongName, output, colors=False):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import seaborn as sns
    # create data
    longnames = []
    for x in df.columns.tolist():
        if x in LongName:
            longnames.append(LongName.get(x))
        else:
            longnames.append(x)
    names = df.columns.tolist()
    data = df.values.tolist()
    species = df.index.values
    # get size of table
    categories = len(df.columns)
    total = len(df.index)
    Rows = total // 2
    Rows += total % 2
    Position = list(range(1, total+1))
    # get colors figured out
    if not colors:
        color_palette = resources.pref_colors
    else:
        color_palette = colors
    # draw figure
    if len(species) < 3:
        fig = plt.figure(1, figsize=(8, 4))
    else:
        fig = plt.figure(1, figsize=(8, 8))
    for k in range(total):
        ax = fig.add_subplot(Rows, 2, Position[k])
        # Create a circle for the center of the plot
        my_circle = plt.Circle((0, 0), 0.7, color='white')
        plt.pie(data[0], labels=names, colors=color_palette)
        p = plt.gcf()
        p.gca().add_artist(my_circle)
        plt.title(species[k])
    patches = [mpatches.Patch(color=color_palette[i], label="{:s}".format(
        longnames[i])) for i in range(len(longnames))]
    plt.legend(handles=patches, bbox_to_anchor=(1, 0.5),
               bbox_transform=fig.transFigure, loc="center left", ncol=1)
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)


def drawbarplot(df, output):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib.pyplot as plt
        import seaborn as sns
    # num = len(df.columns) + 1
    sns.set(style="darkgrid")
    fig = plt.figure()
    # colors
    if len(df) > len(resources.pref_colors):
        colorplot = sns.husl_palette(len(df), l=.5).as_hex()
        colorplot = [str(x).upper() for x in colorplot]
    else:
        colorplot = resources.pref_colors[:len(df)]
    ax = sns.barplot(data=df, palette=colorplot)
    plt.xlabel('Genomes')
    plt.ylabel('Secreted Proteins')
    plt.xticks(rotation=90)
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)


def distance2mds(df, distance, type, output):
    import numpy as np
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        from sklearn.metrics.pairwise import pairwise_distances
        from sklearn.manifold import MDS
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    # run distance metric on matrix and then plot using NMDS
    num = len(df.index)
    data = np.array(df).astype(int)
    bc_dm = pairwise_distances(data, metric=distance)
    mds = MDS(n_components=2, metric=False, max_iter=999,
              dissimilarity='precomputed', n_init=10, verbose=0)
    result = mds.fit(bc_dm)
    coords = result.embedding_
    stress = 'stress=' + '{0:.4f}'.format(result.stress_)
    # get axis information and make square plus some padding
    xcoords = abs(maxabs(coords[:, 0])) + 0.1
    ycoords = abs(maxabs(coords[:, 1])) + 0.1
    # setup plot
    fig = plt.figure()
    # colors
    if len(df) > len(resources.pref_colors):
        colorplot = sns.husl_palette(len(df), l=.5).as_hex()
        colorplot = [str(x).upper() for x in colorplot]
    else:
        colorplot = resources.pref_colors[:len(df)]
    for i in range(0, num):
        plt.plot(coords[i, 0], coords[i, 1], 'o', markersize=9,
                 color=colorplot[i], label=df.index.values[i])
    plt.xlabel('NMDS axis 1')
    plt.ylabel('NMDS axis 2')
    plt.ylim(-ycoords, ycoords)
    plt.xlim(-xcoords, xcoords)
    '''
    if num < 13: #if number too large, don't plot
    '''
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title('NMDS analysis of '+type+' domains')
    plt.annotate(stress, xy=(1, 0), xycoords='axes fraction',
                 fontsize=12, ha='right', va='bottom')
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)


def ReciprocalBlast(filelist, protortho, cpus):
    '''
    function to run reciprocal diamond blast for generating proteinortho input
    '''
    # generate dmnd databases for each input
    for x in filelist:
        base = os.path.basename(x)
        cmd = ['diamond', 'makedb', '--in', x, '--db', base+'.dmnd']
        if not checkannotations(os.path.join(protortho, base+'.dmnd')):
            runSubprocess(cmd, protortho, log)
    for p in itertools.permutations(filelist, 2):
        query = p[0]
        target = p[1]
        db = os.path.basename(target)+'.dmnd'
        outname = target+'.vs.'+query+'.bla'
        cmd = ['diamond', 'blastp', '--query', query, '--db', db, '--outfmt', '6',
               '--out', outname, '--evalue', '1e-5', '--more-sensitive', '--threads', str(cpus)]
        if not checkannotations(os.path.join(protortho, outname)):
            runSubprocess4(cmd, protortho, log)
        db = os.path.basename(query)+'.dmnd'
        outname = query+'.vs.'+target+'.bla'
        cmd = ['diamond', 'blastp', '--query', target, '--db', db, '--outfmt', '6',
               '--out', outname, '--evalue', '1e-5', '--more-sensitive', '--threads', str(cpus)]
        if not checkannotations(os.path.join(protortho, outname)):
            runSubprocess4(cmd, protortho, log)
        db = os.path.basename(target)+'.dmnd'
        outname = target+'.vs.'+target+'.bla'
        cmd = ['diamond', 'blastp', '--query', target, '--db', db, '--outfmt', '6',
               '--out', outname, '--evalue', '1e-5', '--more-sensitive', '--threads', str(cpus)]
        if not checkannotations(os.path.join(protortho, outname)):
            runSubprocess4(cmd, protortho, log)
        db = os.path.basename(query)+'.dmnd'
        outname = query+'.vs.'+query+'.bla'
        cmd = ['diamond', 'blastp', '--query', query, '--db', db, '--outfmt', '6',
               '--out', outname, '--evalue', '1e-5', '--more-sensitive', '--threads', str(cpus)]
        if not checkannotations(os.path.join(protortho, outname)):
            runSubprocess4(cmd, protortho, log)


def singletons(poff, name):
    with open(poff, 'r') as input:
        count = 0
        for line in input:
            line = line.replace('\n', '')
            if line.startswith('#'):
                header = line
                species = header.split('\t')[3:]
                i = species.index(name.replace(' ', '_')) + 3
                continue
            col = line.split('\t')
            if col[0] == '1' and col[i] != '*':
                count += 1
        return count


def orthologs(poff, name):
    with open(poff, 'r') as input:
        count = 0
        for line in input:
            line = line.replace('\n', '')
            if line.startswith('#'):
                header = line
                species = header.split('\t')[3:]
                i = species.index(name.replace(' ', '_')) + 3
                continue
            col = line.split('\t')
            if col[0] != '1' and col[i] != '*':
                count += 1
        return count


def iprTSV2dict(file, terms):
    iprDict = {}
    with io.open(file, 'r', encoding="utf-8") as infile:
        for line in infile:
            if line.startswith('ENTRY_AC') or line.startswith('\n'):
                continue
            line = line.rstrip()
            entry, type, name = line.split('\t')
            if not entry in iprDict:
                iprDict[entry] = name
    return iprDict


def iprxml2dict(xmlfile, terms):
    import xml.etree.cElementTree as cElementTree
    iprDict = {}
    for event, elem in cElementTree.iterparse(xmlfile):
        if elem.tag == 'interpro':
            ID = elem.attrib['id']
            if ID in terms:
                for x in elem.getchildren():
                    if x.tag == 'name':
                        description = x.text
                iprDict[ID] = description
                elem.clear()
            else:
                elem.clear()
    return iprDict


def pfam2dict(file):
    pfamDict = {}
    with open(file, 'r') as input:
        for line in input:
            try:
                line = line.decode('utf-8').rstrip()
            except AttributeError:
                line = line.rstrip()
            if line.startswith('PF'):  # just check to be sure
                cols = line.split('\t')
                ID = cols[0]
                desc = cols[4]
                pfamDict[ID] = desc
    return pfamDict


def flipKeyValues(input):
    flipped = {}
    for k, v in list(input.items()):
        for y in v:
            if not y in flipped:
                flipped[y] = k
    return flipped


def dictFlip(input):
    # flip the list of dictionaries
    outDict = {}
    for x in input:
        for k, v in natsorted(iter(x.items())):
            for i in v:
                if i in outDict:
                    outDict[i].append(k)
                else:
                    outDict[i] = [k]
    return outDict


def busco_dictFlip(input):
    # flip the list of dictionaries
    output = []
    for x in input:
        outDict = {}
        for k, v in natsorted(iter(x.items())):
            for i in v:
                if i in outDict:
                    outDict[i].append(k)
                else:
                    outDict[i] = [k]
        output.append(outDict)
    return output


def dictFlipLookup(input, lookup):
    outDict = {}
    for x in input:
        for k, v in natsorted(iter(x.items())):
            # lookup description in another dictionary
            if not lookup.get(k) is None:
                result = k+': '+lookup.get(k)
            else:
                result = k+': No description'
            result = result.encode('utf-8')
            for i in v:
                if i in outDict:
                    outDict[i].append(str(result))
                else:
                    outDict[i] = [str(result)]
    return outDict


def copyDirectory(src, dest, overwrite=False):
    import shutil
    if overwrite:
        if os.path.isdir(dest):
            shutil.rmtree(dest)
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        log.debug('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        log.debug('Directory not copied. Error: %s' % e)


def download_buscos(name, Database):
    if name in resources.busco_links:
        log.info("Downloading %s busco models" % name)
        address = resources.busco_links.get(name)
        filename = address.split('/')[-1]
        if name == 'fungiv1':
            foldername = 'fungi'
        else:
            foldername = filename.split('.')[0]
        cmd = ['wget', '-c', '--tries=0', '--read-timeout=20', address]
        runSubprocess(cmd, '.', log)
        cmd = ['tar', '-zxf', filename]
        runSubprocess(cmd, '.', log)
        copyDirectory(os.path.abspath(foldername),
                      os.path.join(Database, name))
        shutil.rmtree(foldername)
        os.remove(filename)
    else:
        log.error("%s not a valid BUSCO database" % name)
        validBusco = list(resources.busco_links.keys())
        log.error("Valid BUSCO DBs: %s" % (', '.join(validBusco)))
        sys.exit(1)


def fasta2dict(Fasta):
    answer = dict()
    with open(Fasta, 'r') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'fasta')
        for record in SeqRecords:
            if record.id in answer:
                print("WARNING - duplicate key!")
            else:
                answer[record.id] = str(record.seq)
    return answer


def ortho2phylogeny(folder, df, num, dict, cpus, bootstrap, tmpdir, outgroup, sp_file, name, sc_buscos, ml_method):
    import pylab
    from Bio import Phylo
    from Bio.Phylo.Consensus import get_support
    if outgroup:
        # load species fasta ids into dictionary
        OutGroup = {}
        with open(sp_file, 'r') as sp:
            for rec in SeqIO.parse(sp, 'fasta'):
                OutGroup[rec.id] = rec.seq
    # single copy orthologs are in a dataframe, count and then randomly select
    num_species = len(df.columns)
    species = df.columns.values
    if len(df) == 0:
        log.error("0 single copy BUSCO orthologs found, skipping phylogeny")
        return
    if len(df) < int(num):
        number = len(df)
        log.info(
            "Found %i single copy BUSCO orthologs, will use all to infer phylogeny" % (len(df)))
        subsampled = df
    else:
        number = int(num)
        log.info("Found %i single copy BUSCO orthologs, will randomly select %i to infer phylogeny" % (
            len(df), number))
        subsampled = df.sample(n=number)

    if outgroup:  # passed a list to extract from parent script
        busco_list = sc_buscos

    # since you checked for BUSCO id across all previously, loop through first set and print BUSCOs to file
    with open(os.path.join(tmpdir, 'phylogeny.buscos.used.txt'), 'w') as busco_out:
        with open(os.path.join(tmpdir, 'phylogeny.concat.fa'), 'w') as proteinout:
            if outgroup:
                proteinout.write(">%s\n" % name)
                for y in busco_list:
                    proteinout.write("%s" % (OutGroup.get(y)))
                proteinout.write('\n')
            for i in range(0, num_species):
                proteinout.write(">%s\n" % species[i])
                proteins = fasta2dict(os.path.join(folder, species[i]+'.faa'))
                for row in subsampled[species[i]].items():
                    proteinout.write("%s" % proteins.get(row[1]))
                    busco_out.write("%s\t%s\n" % (dict[i].get(row[1]), row[1]))
                proteinout.write('\n')
    cmd = ['mafft', '--anysymbol', '--quiet', os.path.join(tmpdir, 'phylogeny.concat.fa')]
    runSubprocess2(cmd, '.', log, os.path.join(tmpdir, 'phylogeny.mafft.fa'))
    cmd = ['trimal', '-in', os.path.join(tmpdir, 'phylogeny.mafft.fa'), '-out', os.path.join(
        tmpdir, 'phylogeny.trimal.phylip'), '-automated1', '-phylip']
    runSubprocess(cmd, '.', log)
    if ml_method == 'raxml':
        cmd = ['raxmlHPC-PTHREADS', '-T', str(cpus), '-f', 'a', '-m', 'PROTGAMMAAUTO', '-p', '12345',
               '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
        if outgroup:
            cmd = cmd + ['-o', name]
        treefile = os.path.join(tmpdir, 'RAxML_bootstrap.nwk')
        runSubprocess(cmd, tmpdir, log)
        # parse with biopython and draw
        trees = list(Phylo.parse(treefile, 'newick'))
        best = Phylo.read(os.path.join(tmpdir, 'RAxML_bestTree.nwk'), 'newick')
        support_tree = get_support(best, trees)
        Phylo.draw(support_tree, do_show=False)
        pylab.axis('off')
        pylab.savefig(os.path.join(tmpdir, 'ML.phylogeny.pdf'),
                      format='pdf', bbox_inches='tight', dpi=1000)
    else:  # run iqtree as faster and better than raxml in initial testing
        cmd = ['iqtree', '-s', 'phylogeny.trimal.phylip', '-nt', 'AUTO',
               '-ntmax', str(cpus), '-seed', '12345', '-bb', '1000']
        if outgroup:
            cmd = cmd + ['-o', name]
        runSubprocess(cmd, tmpdir, log)
        treefile = os.path.join(tmpdir, 'phylogeny.trimal.phylip.treefile')
        best = Phylo.read(treefile, 'newick')
        Phylo.draw(best, do_show=False)
        pylab.axis('off')
        pylab.savefig(os.path.join(tmpdir, 'ML.phylogeny.pdf'),
                      format='pdf', bbox_inches='tight', dpi=1000)


def getTrainResults(input):
    with open(input, 'r') as train:
        for line in train:
            try:
                line = line.decode('utf-8')
            except AttributeError:
                pass
            line = line.rstrip()
            if line.startswith('nucleotide level'):
                line = line.replace(' ', '')
                values1 = line.split('|')  # get [1] and [2]
            if line.startswith('exon level'):
                line = line.replace(' ', '')  # get [6] and [7]
                values2 = line.split('|')
            if line.startswith('gene level'):
                line = line.replace(' ', '')
                values3 = line.split('|')  # get [6] and [7]
        return (float(values1[1]), float(values1[2]), float(values2[6]), float(values2[7]), float(values3[6]), float(values3[7]))


def count_multi_CDS_genes(input, filterlist):
    # take funannotate annotation dictionary and return number of genes with more than one CDS
    counter = 0
    counter_inList = 0
    for k, v in natsorted(list(input.items())):
        if len(v['CDS'][0]) > 1:
            counter += 1
            if k in filterlist:
                counter_inList += 1
    return len(input), counter, len(filterlist), counter_inList


def selectTrainingModels(input, fasta, genemark_gtf, output):
    from collections import OrderedDict
    '''
    function to take a GFF3 file and filter the gene models so they are non-overalpping
    also sort the models by number of exons, the more the better.
    '''
    def _sortDict(d):
        return (len(d[1]['CDS'][0]))
    # load gene models into funannotate structured dictionary
    gene_inter = defaultdict(InterLap)
    Genes = {}
    Genes = gff2dict(input, fasta, Genes)
    # add to InterLap output proteins
    proteins = 'augustus.training.proteins.fa'
    ignoreList = []
    keeperList = getGenesGTF(genemark_gtf)
    # check number of multi-cds genes
    countGenes, countGenesCDS, countKeeper, countKeeperCDS = count_multi_CDS_genes(
        Genes, keeperList)
    log.debug('{:,} PASA genes; {:,} have multi-CDS; {:,} from filterGeneMark; {:,} have multi-CDS'.format(
        countGenes, countGenesCDS, countKeeper, countKeeperCDS))
    multiCDScheck, keeperCheck = (False,)*2
    if countKeeper >= 200:
        keeperCheck = True
    if keeperCheck:
        if countKeeperCDS >= 200:
            multiCDScheck = True
    else:
        if countGenesCDS >= 200:
            multiCDScheck = True
    log.debug('filterGeneMark GTF filter set to {:}; require genes with multiple CDS set to {:}'.format(
        keeperCheck, multiCDScheck))
    with open(proteins, 'w') as protout:
        for k, v in natsorted(list(Genes.items())):
            if keeperCheck and not k in keeperList:
                ignoreList.append(k)
                continue
            if multiCDScheck and len(v['CDS'][0]) < 2:
                ignoreList.append(k)
                continue
            # add to interlap object and write protein out
            gene_inter[v['contig']].add(
                (v['location'][0], v['location'][1], v['strand'], k, len(v['CDS'][0])))
            protout.write('>%s___%i\n%s\n' %
                          (k, len(v['CDS'][0]), v['protein'][0]))

    # make sure gene models are unique, so do pairwise diamond search @ 80% identity
    cmd = ['diamond', 'makedb', '--in',
           'augustus.training.proteins.fa', '--db', 'aug_training.dmnd']
    runSubprocess4(cmd, '.', log)
    cmd = ['diamond', 'blastp', '--query', 'augustus.training.proteins.fa', '--db', 'aug_training.dmnd', '--more-sensitive', '-o',
           'aug.blast.txt', '-f', '6', 'qseqid', 'sseqid', 'pident', '--query-cover', '80', '--subject-cover', '80', '--id', '80', '--no-self-hits']
    runSubprocess4(cmd, '.', log)
    blast_results = []
    with open('aug.blast.txt', 'r') as blast:
        for line in blast:
            line = line.rstrip()
            line = line.replace('___', '\t')
            blast_results.append(line.split('\t'))
    sortedBlast = natsorted(
        blast_results, key=lambda x: int(x[1]), reverse=True)
    blastignore = []
    for hit in sortedBlast:
        if hit[0] in blastignore or hit[2] in blastignore:
            continue
        if int(hit[1]) >= int(hit[3]):
            if not hit[2] in blastignore:
                blastignore.append(hit[2])
        else:
            if not hit[0] in blastignore:
                blastignore.append(hit[0])
    log.debug('{:,} models fail blast identity threshold'.format(
        len(blastignore)))
    SafeRemove('augustus.training.proteins.fa')
    SafeRemove('aug_training.dmnd')
    SafeRemove('aug.blast.txt')
    # now return cleaned genemark GTF file
    finalIgnoreList = []
    for x in ignoreList:
        if not x in finalIgnoreList:
            finalIgnoreList.append(x)
    for y in blastignore:
        if not y in finalIgnoreList:
            finalIgnoreList.append(y)
    log.debug('{:,} models will be ignored for training Augustus'.format(
        len(finalIgnoreList)))
    GenesPass = {}
    for k, v in natsorted(list(Genes.items())):
        if not k in finalIgnoreList and not k in GenesPass:
            loc = sorted([v['location'][0], v['location'][1]])
            if loc in gene_inter[v['contig']]:
                hits = list(gene_inter[v['contig']].find(loc))
                sortedHits = sorted(
                    hits, key=lambda x: int(x[4]), reverse=True)
                validHits = []
                for y in sortedHits:
                    if not y[3] in finalIgnoreList and y[3] != k:
                        validHits.append(y)
                if len(validHits) > 0:
                    if not validHits[0][3] in GenesPass:
                        GenesPass[validHits[0][3]] = Genes.get(validHits[0][3])
                else:
                    GenesPass[k] = v

    # now sort dictionary number of exons
    sGenes = sorted(iter(GenesPass.items()), key=_sortDict, reverse=True)
    sortedGenes = OrderedDict(sGenes)
    log.info("{:,} of {:,} models pass training parameters".format(
        len(sortedGenes), len(Genes)))
    # x = dict(itertools.islice(sortedGenes.items(), 0, 2500))
    final = {}
    for i, (k, v) in enumerate(natsorted(list(sortedGenes.items()))):
        v['ids'] = ['g_'+str(i+1)+'-T1']
        final['g_'+str(i+1)] = v
    dict2gff3noUTRs(final, output)
    return len(final)


def getGenesGTF(input):
    genes = []
    with open(input, 'r') as infile:
        for line in infile:
            if not line.startswith('\n') or not line.startswith('#'):
                line = line.rstrip()
                info = line.split('\t')[-1]
                attributes = info.split(';')
                ID = None
                for x in attributes:
                    if x.startswith('gene_id'):
                        tmp = x.replace('gene_id ', '')
                        ID = tmp.replace('"', '')
                if ID:
                    if not ID in genes:
                        genes.append(ID)
    return genes


def trainAugustus(AUGUSTUS_BASE, train_species, trainingset,
                  genome, outdir, cpus, num_training, optimize,
                  config_path):
    if which('randomSplit.pl'):
        RANDOMSPLIT = 'randomSplit.pl'
    else:
        RANDOMSPLIT = os.path.join(AUGUSTUS_BASE, 'scripts', 'randomSplit.pl')
    if which('optimize_augustus.pl'):
        OPTIMIZE = 'optimize_augustus.pl'
    else:
        OPTIMIZE = os.path.join(
            AUGUSTUS_BASE, 'scripts', 'optimize_augustus.pl')
    if which('new_species.pl'):
        NEW_SPECIES = 'new_species.pl'
    else:
        NEW_SPECIES = os.path.join(AUGUSTUS_BASE, 'scripts', 'new_species.pl')
    aug_cpus = '--cpus='+str(cpus)
    species = '--species='+train_species
    aug_log = os.path.join(outdir, 'logfiles', 'augustus_training.log')
    TrainSet = os.path.abspath(trainingset)
    onlytrain = '--onlytrain='+TrainSet+'.train'
    testtrain = TrainSet+'.test'
    trainingdir = os.path.join(
        outdir, 'predict_misc', 'tmp_opt_'+train_species)
    myENV = os.environ
    myENV['AUGUSTUS_CONFIG_PATH'] = config_path
    with open(aug_log, 'w') as logfile:
        if not CheckAugustusSpecies(train_species):
            subprocess.call([NEW_SPECIES, '--AUGUSTUS_CONFIG_PATH={:}'.format(
                config_path), species], stdout=logfile, stderr=logfile)
        # run etraining again to only use best models from EVM for training
        p1 = subprocess.Popen(['etraining', species, TrainSet],
                              cwd=os.path.join(outdir, 'predict_misc'),
                              stderr=logfile, stdout=logfile, env=dict(myENV))
        p1.communicate()
        # split off num_training models for testing purposes
        subprocess.call([RANDOMSPLIT, TrainSet, str(num_training)],
                        cwd=os.path.join(outdir, 'predict_misc'))
        if os.path.isfile(os.path.join(outdir, 'predict_misc', TrainSet+'.train')):
            with open(os.path.join(outdir, 'predict_misc', 'augustus.initial.training.txt'), 'w') as initialtraining:
                subprocess.call(['augustus', '--AUGUSTUS_CONFIG_PATH={:}'.format(
                    config_path), species, TrainSet+'.test'], stdout=initialtraining, cwd=os.path.join(outdir, 'predict_misc'))
            train_results = getTrainResults(os.path.join(
                outdir, 'predict_misc', 'augustus.initial.training.txt'))
            trainTable = [['Feature', 'Specificity', 'Sensitivity'],
                          ['nucleotides', '{:.1%}'.format(
                              train_results[0]), '{:.1%}'.format(train_results[1])],
                          ['exons', '{:.1%}'.format(
                              train_results[2]), '{:.1%}'.format(train_results[3])],
                          ['genes', '{:.1%}'.format(
                              train_results[4]), '{:.1%}'.format(train_results[5])]
                          ]
            log.info('Augustus initial training results:')
            train_table = print_table(trainTable, return_str=True)
            sys.stderr.write(train_table)
            if optimize:
                # now run optimization
                subprocess.call([OPTIMIZE, '--AUGUSTUS_CONFIG_PATH={:}'.format(config_path), species, aug_cpus,
                                 onlytrain, testtrain], cwd=os.path.join(outdir, 'predict_misc'), stderr=logfile, stdout=logfile)
                # run etraining again
                p2 = subprocess.Popen(['etraining', species, TrainSet], cwd=os.path.join(
                    outdir, 'predict_misc'), stderr=logfile, stdout=logfile, env=dict(myENV))
                p2.communicate()
                with open(os.path.join(outdir, 'predict_misc', 'augustus.final.training.txt'), 'w') as finaltraining:
                    subprocess.call(['augustus', '--AUGUSTUS_CONFIG_PATH={:}'.format(
                        config_path), species, TrainSet+'.test'], stdout=finaltraining, cwd=os.path.join(outdir, 'predict_misc'))
                train_results = getTrainResults(os.path.join(
                    outdir, 'predict_misc', 'augustus.final.training.txt'))
                trainTable = [['Feature', 'Specificity', 'Sensitivity'],
                              ['nucleotides', '{:.1%}'.format(
                                  train_results[0]), '{:.1%}'.format(train_results[1])],
                              ['exons', '{:.1%}'.format(
                                  train_results[2]), '{:.1%}'.format(train_results[3])],
                              ['genes', '{:.1%}'.format(
                                  train_results[4]), '{:.1%}'.format(train_results[5])]
                              ]
                log.info('Augustus optimized training results:')
                train_table = print_table(trainTable, return_str=True)
                sys.stderr.write(train_table)
                # clean up tmp folder
                shutil.rmtree(trainingdir)
            else:
                if train_results[4] < 0.50:
                    log.info(
                        "Accuracy seems low, you can try to improve by passing the --optimize_augustus option.")
        else:
            log.error("AUGUSTUS training failed, check logfiles")
            sys.exit(1)


def sortList(input, col):
    return natsorted(input, key=operator.itemgetter(col))


def sortHints(input, output):
    data = []
    with open(input, 'r') as infile:
        for line in infile:
            line = line.rstrip()
            data.append(line.split('\t'))
    # replicate this: sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1
    sort1 = sortList(data, 3)
    sort2 = sortList(sort1, 4)
    sort3 = sortList(sort2, 2)
    sort4 = sortList(sort3, 0)
    with open(output, 'w') as sort_out:
        for line in sort4:
            sort_out.write('%s\n' % '\t'.join(line))


def checkgoatools(input):
    with open(input, 'r') as goatools:
        count = -1
        result = False
        headercount = 0
        for line in goatools:
            count += 1
            if line.startswith('GO\tNS') or line.startswith('#'):
                headercount = count
            if line.startswith('GO:'):
                result = True
    return (result, headercount)


def translatemRNA(input, output):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    with open(output, 'w') as outfile:
        with open(input, 'r') as fasta:
            for header, seq in SimpleFastaParser(fasta):
                codon_start = 1
                for x in header.split(' '):
                    if x.startswith('codon_start='):
                        codon_start = int(
                            x.replace('codon_start=', '').rstrip())
                # transcripts should already be in proper orientation
                protSeq = translate(seq, '+', codon_start-1)
                outfile.write('>{:}\n{:}\n'.format(header, protSeq))


def alignMAFFT(input, output):
    FNULL = open(os.devnull, 'w')
    with open(output, 'w') as outfile:
        subprocess.call(['mafft', '--anysymbol', '--quiet', input],
                        stderr=FNULL, stdout=outfile)


def align2Codon(alignment, transcripts, output):
    FNULL = open(os.devnull, 'w')
    with open(output, 'w') as outfile:
        subprocess.call(['perl', os.path.join(parentdir, 'aux_scripts', 'pal2nal.pl'),
                         alignment, transcripts, '-output', 'fasta'], stderr=FNULL, stdout=outfile)
    if getSize(output) < 1:
        os.remove(output)
        log.debug('dNdS Error: pal2nal failed for %s' % alignment)


def counttaxa(input):
    ct = 0
    with open(input, 'r') as tree:
        line = tree.readline()
        ct = line.count(',')+1
    return ct


def getMatchFileName(pattern, directory):
    result = None
    for f in os.listdir(directory):
        if pattern in f:
            result = os.path.join(directory, f)
    return result


def drawPhyMLtree(fasta, tree):
    FNULL = open(os.devnull, 'w')
    fc = countfasta(fasta)
    # need to convert to phylip format
    base = os.path.basename(fasta).split('.')[0]
    dir = os.path.dirname(fasta)
    tmp1 = os.path.join(dir, base+'.draw2tree.phylip')
    subprocess.call(['trimal', '-in', fasta, '-out', tmp1, '-phylip'])
    # draw tree
    subprocess.call(['phyml', '-i', tmp1], stdout=FNULL, stderr=FNULL)
    tmp2 = getMatchFileName(base+'.draw2tree.phylip_phyml_tree', dir)
    # check that num taxa in tree = input
    tc = counttaxa(tmp2)
    if tc != fc:  # something failed...
        log.debug('dNdS Error: phyml tree failed for %s' % fasta)
        # retry
        subprocess.call(['trimal', '-in', fasta, '-out', tmp1, '-phylip'])
        subprocess.call(['phyml', '-i', tmp1], stdout=FNULL, stderr=FNULL)
    # rename and clean
    os.rename(tmp2, tree)
    SafeRemove(tmp1)
    stats = getMatchFileName(base+'.draw2tree.phylip_phyml_stats', dir)
    SafeRemove(stats)


def simplestTreeEver(fasta, tree):
    with open(tree, 'w') as outfile:
        with open(fasta, 'r') as input:
            ids = []
            for rec in SeqIO.parse(input, 'fasta'):
                ids.append(rec.id)
            outfile.write('(%s,%s);' % (ids[0], ids[1]))


def rundNdSexhaustive(folder):
    # setup intermediate files
    tmpdir = os.path.dirname(folder)
    name = os.path.basename(folder)
    transcripts = os.path.join(tmpdir, name+'.transcripts.fa')
    prots = os.path.join(tmpdir, name+'.proteins.fa')
    aln = os.path.join(tmpdir, name+'.aln')
    codon = os.path.join(tmpdir, name+'.codon.aln')
    tree = os.path.join(tmpdir, name+'.tree')
    log = os.path.join(tmpdir, name+'.log')
    finallog = os.path.join(tmpdir, name, name+'.log')
    if not checkannotations(finallog):
        num_seqs = countfasta(transcripts)
        # Translate to protein space
        translatemRNA(transcripts, prots)
        # align protein sequences
        alignMAFFT(prots, aln)
        # convert to codon alignment
        align2Codon(aln, transcripts, codon)
        if checkannotations(codon):
            if num_seqs > 2:
                # now generate a tree using phyml
                drawPhyMLtree(codon, tree)
            else:
                simplestTreeEver(transcripts, tree)
            # now run codeml through ete3
            etecmd = ['ete3', 'evol', '--alg', os.path.abspath(codon), '-t', os.path.abspath(
                tree), '--models', 'M0', 'M1', 'M2', 'M7', 'M8', '-o', name, '--clear_all', '--codeml_param', 'cleandata,1']
            with open(log, 'w') as logfile:
                logfile.write('\n%s\n' % ' '.join(etecmd))
                subprocess.call(etecmd, cwd=tmpdir,
                                stdout=logfile, stderr=logfile)
    # clean up
    for file in os.listdir(tmpdir):
        if file.startswith(name+'.'):
            os.rename(os.path.join(tmpdir, file),
                      os.path.join(tmpdir, name, file))


def rundNdSestimate(folder):
    # setup intermediate files
    tmpdir = os.path.dirname(folder)
    name = os.path.basename(folder)
    transcripts = os.path.join(tmpdir, name+'.transcripts.fa')
    prots = os.path.join(tmpdir, name+'.proteins.fa')
    aln = os.path.join(tmpdir, name+'.aln')
    codon = os.path.join(tmpdir, name+'.codon.aln')
    tree = os.path.join(tmpdir, name+'.tree')
    log = os.path.join(tmpdir, name+'.log')
    finallog = os.path.join(tmpdir, name, name+'.log')
    if not checkannotations(finallog):
        num_seqs = countfasta(transcripts)
        # Translate to protein space
        translatemRNA(transcripts, prots)
        # align protein sequences
        alignMAFFT(prots, aln)
        # convert to codon alignment
        align2Codon(aln, transcripts, codon)
        if checkannotations(codon):
            if num_seqs > 2:
                # now generate a tree using phyml
                drawPhyMLtree(codon, tree)
            else:
                simplestTreeEver(transcripts, tree)
            # now run codeml through ete3
            etecmd = ['ete3', 'evol', '--alg', os.path.abspath(codon), '-t', os.path.abspath(
                tree), '--models', 'M0', '-o', name, '--clear_all', '--codeml_param', 'cleandata,1']
            with open(log, 'w') as logfile:
                logfile.write('\n%s\n' % ' '.join(etecmd))
                subprocess.call(etecmd, cwd=tmpdir,
                                stdout=logfile, stderr=logfile)
    # clean up
    for file in os.listdir(tmpdir):
        if file.startswith(name+'.'):
            os.rename(os.path.join(tmpdir, file),
                      os.path.join(tmpdir, name, file))


def get_subdirs(a_dir):
    return [os.path.join(a_dir, name) for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


def get_subdirs2(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


def parsedNdS(folder):
    results = {}
    hits = get_subdirs2(folder)
    for x in hits:
        finallog = os.path.join(folder, x, x+'.log')
        # parse logfile to get omega
        dnds = 'NA'
        m1m2p = 'NA'
        m7m8p = 'NA'
        if os.path.isfile(finallog):
            with open(finallog, 'r') as input:
                for line in input:
                    line = line.strip()
                    if 'M7' in line and 'M8' in line and '|' in line:
                        m7m8p = line.split('|')[-1].strip()
                        m7m8p = m7m8p.replace('*', '')
                        m7m8p = '{0:.5f}'.format(float(m7m8p))
                    elif 'M1' in line and 'M2' in line and '|' in line:
                        m1m2p = line.split('|')[-1].lstrip()
                        m1m2p = m1m2p.replace('*', '')
                        m1m2p = '{0:.5f}'.format(float(m1m2p))
                    elif line.startswith('- Model M0'):
                        nextline = next(input)
                        dnds = nextline.split('tree: ')[1].rstrip()
        results[x] = (dnds, m1m2p, m7m8p)
    return results


def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out


def getBlastDBinfo(input):
    '''
    function to return a tuple of info using blastdbcmd
    tuple: (name, date, #sequences)
    '''
    cmd = ['blastdbcmd', '-info', '-db', input]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stderr:
        print((stderr.split('\n')[0]))
    results = stdout.split('\n\n')
    results = [x for x in results if x]
    # parse results which are now in list, look for starts with Database and then Date
    Name, Date, NumSeqs = (None,)*3
    for x in results:
        if x.startswith('Database:'):
            hit = x.split('\n\t')
            Name = hit[0].replace('Database: ', '')
            NumSeqs = hit[1].split(' sequences;')[0].replace(',', '')
        if x.startswith('Date:'):
            Date = x.split('\t')[0].replace('Date: ', '')
    return (Name, Date, NumSeqs)


HEADER = '''
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <!-- The above 3 meta tags *must* come first in the head; any other head content must come *after* these tags -->
    <meta name="funannotate comparative genomics output" content="">
    <meta name="Jonathan Palmer" content="">
    <title>Funannotate</title>
    <!-- Bootstrap core CSS -->
    <link href="css/bootstrap.min.css" rel="stylesheet">
    <!-- Custom styles for this template -->
    <link href="css/starter-template.css" rel="stylesheet">
    <script src="js/ie-emulation-modes-warning.js"></script>
  </head>
  <body>
    <nav class="navbar navbar-inverse navbar-fixed-top">
      <div class="container">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
          </button>
          <a class="navbar-brand" href="index.html">Funannotate</a>
        </div>
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li><a href="stats.html">Stats</a></li>
            <li><a href="phylogeny.html">Phylogeny</a></li>
            <li><a href="orthologs.html">Orthologs</a></li>
            <li><a href="interpro.html">InterPro</a></li>
            <li><a href="pfam.html">PFAM</a></li>
            <li><a href="merops.html">Merops</a></li>
            <li><a href="cazy.html">CAZymes</a></li>
            <li><a href="cogs.html">COGs</a></li>
            <li><a href="signalp.html">SignalP</a></li>
            <li><a href="tf.html">TFs</a></li>
            <li><a href="secmet.html">SecMet</a></li>
            <li><a href="go.html">GO</a></li>
            <li><a href="citation.html">Cite</a></li>
          </ul>
        </div><!--/.nav-collapse -->
      </div>
    </nav>
'''
ORTHOLOGS = '''
    <div class="container">
      <div class="table">
        <h2 class="sub-header">Orthologous protein groups</h2>
          <div class="table-responsive">
'''
INDEX = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Funannotate Results</h2>
         <br>
         <p><a href='stats.html'>Genome Summary Stats</a></p>
         <p><a href='phylogeny.html'>Maximum likelihood Phylogeny (RAxML)</a></p>
         <p><a href='merops.html'>MEROPS Protease Stats</a></p>
         <p><a href='cazy.html'>CAZyme carbohydrate activating enzyme Stats</a></p>
         <p><a href='cogs.html'>COGs Stats</a></p>
         <p><a href='signalp.html'>Secreted proteins (SignalP)</a></p>
         <p><a href='interpro.html'>InterProScan Domain Stats</a></p>
         <p><a href='tf.html'>Transcription Factor Summary</a></p>
         <p><a href='secmet.html'>Secondary Metabolism Cluster Summary</a></p>
         <p><a href='pfam.html'>PFAM Domain Stats</a></p>
         <p><a href='go.html'>Gene Ontology Enrichment Analysis</a></p>
         <p><a href='orthologs.html'>Orthologous proteins</a></p>
         <br>
'''
SUMMARY = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Genome Summary Stats</h2>
          <div class="table-responsive">
'''
PHYLOGENY = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">RAxML Maximum Likelihood Phylogeny</h2>
        <a href='phylogeny/ML.phylogeny.pdf'><img src="phylogeny/ML.phylogeny.pdf" height="500" /></a></div>
'''
NOPHYLOGENY = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Number of species too low to generate phylogeny</h2>
'''
MEROPS = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">MEROPS Protease Families per Genome Results</h2>
        <div class='row'>
        <div class="col-sm-7"><a href='merops/MEROPS.graph.pdf'><img src="merops/MEROPS.graph.pdf" height="350" /></a></div>
        <div class="col-sm-5"><a href='merops/MEROPS.heatmap.pdf'><img src="merops/MEROPS.heatmap.pdf" height="500" /></a></div>
        </div>
        <div class="table-responsive">
'''
INTERPRO = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">InterProScan Domains per Genome Results</h2>
        <div class='row'>
        <a href='interpro/InterProScan.nmds.pdf'><img src="interpro/InterProScan.nmds.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''
PFAM = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">PFAM Domains per Genome Results</h2>
        <div class='row'>
        <a href='pfam/PFAM.nmds.pdf'><img src="pfam/PFAM.nmds.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''
SIGNALP = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Secreted Proteins per Genome Results</h2>
        <div class='row'>
        <a href='signalp/signalp.pdf'><img src="signalp/signalp.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''
TF = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Fungal Transcription Factors per Genome Results</h2>
        <div class='row'>
        <a href='tfs/TF.heatmap.pdf'><img src="tfs/TF.heatmap.pdf" height="800" /></a></div>
        <div class="table-responsive">
'''
SECMET = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Secondary Metabolism Clusters per Genome Results</h2>
        <div class='row'>
        <a href='secmet/SM.graph.pdf'><img src="secmet/SM.graph.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''

CAZY = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">CAZyme Families per Genome Results</h2>
        <div class='row'>
        <div class="col-sm-7"><a href='cazy/CAZy.graph.pdf'><img src="cazy/CAZy.graph.pdf" height="350" /></a></div>
        <div class="col-sm-5"><a href='cazy/CAZy.heatmap.pdf'><img src="cazy/CAZy.heatmap.pdf" height="600" /></a></div>
        </div>
        <div class="table-responsive">
'''

COG = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">Clusters of Orthologous Groups (COGs) per Genome Results</h2>
        <div class='row'>
        <a href='cogs/COGS.graph.pdf'><img src="cogs/COGS.graph.pdf" height="500" /></a></div>
        <div class="table-responsive">
'''

GO = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">GO ontology enrichment Results</h2>
        <div class='row'>
'''
MISSING = '''
    <div class="container">
      <div class="starter-template">
        <h2 class="sub-header">These data are missing from annotation.</h2>
'''
CITATION = '''
    <div class="container">
      <div class="starter-template">
        <h3 class="sub-header">If you found Funannotate useful please cite:</h3>
        <p>Palmer JM. 2016. Funannotate: a fungal genome annotation and comparative genomics pipeline. <a href="https://github.com/nextgenusfs/funannotate">https://github.com/nextgenusfs/funannotate</a>.</p>
'''
FOOTER = '''
          </div>
      </div>

    </div><!-- /.container -->


    <!-- Bootstrap core JavaScript
    ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
    <script>window.jQuery || document.write('<script src="js/jquery.min.js"><\/script>')</script>
    <script src="js/bootstrap.min.js"></script>
    <!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->
    <script src="js/ie10-viewport-bug-workaround.js"></script>
  </body>
</html>

'''
HEADER2 = '''
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="funannotate comparative genomics output" content="">
    <meta name="Jonathan Palmer" content="">
    <title>Funannotate</title>
    <link href="css/bootstrap.min.css" rel="stylesheet">
    <link href="css/starter-template.css" rel="stylesheet">
    <script src="js/ie-emulation-modes-warning.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/t/bs/dt-1.10.11/datatables.min.css"/>
    <script type="text/javascript" src="https://cdn.datatables.net/t/bs/dt-1.10.11/datatables.min.js"></script>
  </head>
  <body>
    <nav class="navbar navbar-inverse navbar-fixed-top">
      <div class="container-fluid">
        <div class="navbar-header">
            <span class="sr-only">Toggle navigation</span>
          <a class="navbar-brand" href="index.html">Funannotate</a>
        </div>
        <div class="navbar-header">
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li class="active"><a href="stats.html">Stats</a></li>
            <li><a href="orthologs.html">Orthologs</a></li>
            <li><a href="interpro.html">InterProScan</a></li>
            <li><a href="pfam.html">PFAM</a></li>
            <li><a href="merops.html">Merops</a></li>
            <li><a href="cazy.html">CAZymes</a></li>
            <li><a href="signalp.html">SignalP</a></li>
            <li><a href="go.html">GO ontology</a></li>
            <li><a href="citation.html">Citation</a></li>
            <li class="dropdown">
          <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Genomes <span class="caret"></span></a>
          <ul class="dropdown-menu">
'''
