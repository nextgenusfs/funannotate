'''
Awesome people who have directly contributed to the project:
Jon Palmer - Bug finder & advice on project direction
Mahmut Uludag - Bug finder

Help:       print pybam.wat
Github:     http://github.com/JohnLonginotto/pybam

This code was written by John Longinotto, a PhD student of the Pospisilik Lab at the Max Planck Institute of Immunbiology & Epigenetics, Freiburg.
My PhD is funded by the Deutsches Epigenom Programm (DEEP), and the Max Planck IMPRS Program.
I study Adipose Biology and Circadian Rhythm in mice, although it seems these days I spend most of my time at the computer :-)
'''

import os
import sys
import zlib
import time
import tempfile
import subprocess
from array import array
from struct import unpack

CtoPy = {'A': '<c', 'c': '<b', 'C': '<B', 's': '<h',
         'S': '<H', 'i': '<i', 'I': '<I', 'f': '<f'}
py4py = {'A':  1, 'c':  1, 'C':  1, 's':  2,
         'S':  2, 'i':  4, 'I':  4, 'f':  4}
dna_codes = '=ACMGRSVTWYHKDBN'
cigar_codes = 'MIDNSHP=X'
parse_codes = {
    'sam':                     ' The current alignment in SAM format.',
    'bam':                     ' All the bytes that make up the current alignment ("read"),\n                              still in binary just as it was in the BAM file. Useful\n                              when creating a new BAM file of filtered alignments.',
    'sam_qname':               ' [1st column in SAM] The QNAME (fragment ID) of the alignment.',
    'bam_qname':               ' The original bytes before decoding to sam_qname.',
    'sam_flag':                ' [2nd column in SAM] The FLAG number of the alignment.',
    'bam_flag':                ' The original bytes before decoding to sam_flag.',
    'sam_refID':               ' The chromosome ID (not the same as the name!).\n                              Chromosome names are stored in the BAM header (file_chromosomes),\n                              so to convert refIDs to chromsome names one needs to do:\n                              "my_bam.file_chromosomes[read.sam_refID]" (or use sam_rname)\n                              But for comparisons, using the refID is much faster that using\n                              the actual chromosome name (for example, when reading through a\n                              sorted BAM file and looking for where last_refID != this_refID)\n                              Note that when negative the alignment is not aligned, and thus one\n                              must not perform my_bam.file_chromosomes[read.sam_refID]\n                              without checking that the value is positive first.',
    'sam_rname':               ' [3rd column in SAM] The actual chromosome/contig name for the\n                              alignment. Will return "*" if refID is negative.',
    'bam_refID':               ' The original bytes before decoding to sam_refID.',
    'sam_pos1':                ' [4th column in SAM] The 1-based position of the alignment. Note\n                              that in SAM format values less than 1 are converted to "0" for\n                              "no data" and sam_pos1 will also do this.',
    'sam_pos0':                ' The 0-based position of the alignment. Note that in SAM all\n                              positions are 1-based, but in BAM they are stored as 0-based.\n                              Unlike sam_pos1, negative values are kept as negative values,\n                              essentially giving one the decoded value as it was stored.',
    'bam_pos':                 ' The original bytes before decoding to sam_pos*.',
    'sam_mapq':                ' [5th column in SAM] The Mapping Quality of the current alignment.',
    'bam_mapq':                ' The original bytes before decoding to sam_mapq.',
    'sam_cigar_string':        ' [6th column in SAM] The CIGAR string, as per the SAM format.\n                              Allowed values are "MIDNSHP=X".',
    'sam_cigar_list':          ' A list of tuples with 2 values per tuple:\n                              the number of bases, and the CIGAR operation applied to those\n                              bases. Faster to calculate than sam_cigar_string.',
    'bam_cigar':               ' The original bytes before decoding to sam_cigar_*.',
    'sam_next_refID':          ' The sam_refID of the alignment\'s mate (if any). Note that as per\n                              sam_refID, this value can be negative and is not the actual\n                              chromosome name (see sam_pnext1).',
    'sam_rnext':               ' [7th column in SAM] The chromosome name of the alignment\'s mate.\n                              Value is "*" if unmapped. Note that in a SAM file this value\n                              is "=" if it is the same as the sam_rname, however pybam will\n                              only do this if the user prints the whole SAM entry with "sam".',
    'bam_next_refID':          ' The original bytes before decoding to sam_next_refID.',
    'sam_pnext1':              ' [8th column in SAM] The 1-based position of the alignment\'s mate.\n                              Note that in SAM format values less than 1 are converted to "0"\n                              for "no data", and sam_pnext1 will also do this.',
    'sam_pnext0':              ' The 0-based position of the alignment\'s mate. Note that in SAM all\n                              positions are 1-based, but in BAM they are stored as 0-based.\n                              Unlike sam_pnext1, negative values are kept as negative values\n                              here, essentially giving you the value as it was stored in BAM.',
    'bam_pnext':               ' The original bytes before decoding to sam_pnext0.',
    'sam_tlen':                ' [9th column in SAM] The TLEN value.',
    'bam_tlen':                ' The original bytes before decoding to sam_tlen.',
    'sam_seq':                 ' [10th column in SAM] The SEQ value (DNA sequence of the alignment).\n                              Allowed values are "ACGTMRSVWYHKDBN and =".',
    'bam_seq':                 ' The original bytes before decoding to sam_seq.',
    'sam_qual':                ' [11th column in SAM] The QUAL value (quality scores per DNA base\n                              in SEQ) of the alignment.',
    'bam_qual':                ' The original bytes before decoding to sam_qual.',
    'sam_tags_list':           ' A list of tuples with 3 values per tuple: a two-letter TAG ID, the\n                              type code used to describe the data in the TAG value (see SAM spec.\n                              for details), and the value of the TAG. Note that the BAM format\n                              has type codes like "c" for a number in the range -127 to +127,\n                              and "C" for a number in the range of 0 to 255.\n                              In a SAM file however, all numerical codes appear to just be stored\n                              using "i", which is a number in the range -2147483647 to +2147483647.\n                              sam_tags_list will therefore return the code used in the BAM file,\n                              and not "i" for all numbers.',
    'sam_tags_string':         ' [12th column a SAM] Returns the TAGs in the same format as would be found \n                              in a SAM file (with all numbers having a signed 32bit code of "i").',
    'bam_tags':                ' The original bytes before decoding to sam_tags_*.',
    'sam_bin':                 ' The bin value of the alignment (used for indexing reads).\n                              Please refer to section 5.3 of the SAM spec for how this\n                              value is calculated.',
    'bam_bin':                 ' The original bytes before decoding to sam_bin.',
    'sam_block_size':          ' The number of bytes the current alignment takes up in the BAM\n                              file minus the four bytes used to store the block_size value\n                              itself. Essentially sam_block_size +4 == bytes needed to store\n                              the current alignment.',
    'bam_block_size':          ' The original bytes before decoding to sam_block_size.',
    'sam_l_read_name':         ' The length of the QNAME plus 1 because the QNAME is terminated\n                              with a NUL byte.',
    'bam_l_read_name':         ' The original bytes before decoding to sam_l_read_name.',
    'sam_l_seq':               ' The number of bases in the seq. Useful if you just want to know\n                              how many bases are in the SEQ but do not need to know what those\n                              bases are (which requires more decoding effort).',
    'bam_l_seq':               ' The original bytes before decoding to sam_l_seq.',
    'sam_n_cigar_op':          ' The number of CIGAR operations in the CIGAR field. Useful if one\n                              wants to know how many CIGAR operations there are, but does not\n                              need to know what they are.',
    'bam_n_cigar_op':          ' The original bytes before decoding to sam_n_cigar_op.',
    'file_alignments_read':    ' A running counter of the number of alignments ("reads"),\n                              processed thus far. Note the BAM format does not store\n                              how many reads are in a file, so the usefulness of this\n                              metric is somewhat limited unless one already knows how\n                              many reads are in the file.',
    'file_binary_header':      ' From the first byte in the file, until the first byte of\n                              the first read. The original binary header.',
    'file_bytes_read':         ' A running counter of the bytes read from the file. Note\n                              that as data is read in arbitary chunks, this is literally\n                              the amount of data read from the file/pipe by pybam.',
    'file_chromosome_lengths': ' The binary header of the BAM file includes chromosome names\n                              and chromosome lengths. This is a dictionary of chromosome-name\n                              keys and chromosome-length values.',
    'file_chromosomes':        ' A list of chromosomes from the binary header.',
    'file_decompressor':       ' BAM files are compressed with bgzip. The value here reflects\n                              the decompressor used. "internal" if pybam\'s internal\n                              decompressor is being used, "gzip" or "pigz" if the system\n                              has these binaries installed and pybam can find them.\n                              Any other value reflects a custom decompression command.',
    'file_directory':          ' The directory the input BAM file can be found in. This will be\n                              correct if the input file is specified via a string or python\n                              file object, however if the input is a pipe such as sys.stdin, \n                              then the current working directory will be used.',
    'file_header':             ' The ASCII portion of the BAM header. This is the typical header\n                              users of samtools will be familiar with.',
    'file_name':               ' The file name (base name) of input file if input is a string or\n                              python file object. If input is via stdin this will be "<stdin>"'
}

wat = '''
Main class: pybam.read
Github:     http://github.com/JohnLonginotto/pybam

[ Dynamic Parser Example ]
  for alignment in pybam.read('/my/data.bam'):
      print alignment.sam_seq

[ Static Parser Example ]
  for seq,mapq in pybam.read('/my/data.bam',['sam_seq','sam_mapq']):
       print seq
       print mapq

[ Mixed Parser Example ]
  my_bam = pybam.read('/my/data.bam',['sam_seq','sam_mapq'])
  print my_bam._static_parser_code
  for seq,mapq in my_bam:
       if seq.startswith('ACGT') and mapq > 10:
       print my_bam.sam

[ Custom Decompressor (from file path) Example ]
  my_bam = pybam.read('/my/data.bam.lzma',decompressor='lzma --decompress --stdout /my/data.bam.lzma')

[ Custom Decompressor (from file object) Example ]
  my_bam = pybam.read(sys.stdin,decompressor='lzma --decompress --stdout') # data given to lzma via stdin
    
[ Force Internal bgzip Decompressor ]
  my_bam = pybam.read('/my/data.bam',decompressor='internal')

[ Parse Words (hah) ]'''
wat += '\n'+''.join([('\n===============================================================================================\n\n  ' if code is 'file_alignments_read' or code is 'sam' else '  ') +
                     (code+' ').ljust(25, '-')+description+'\n' for code, description in sorted(parse_codes.items())]) + '\n'


class read:
    '''
    [ Dynamic Parser Example ]
    for alignment in pybam.read('/my/data.bam'):
        print alignment.sam_seq

    [ Static Parser Example ]
    for seq,mapq in pybam.read('/my/data.bam',['sam_seq','sam_mapq']):
        print seq
        print mapq

    [ Mixed Parser Example ]
    my_bam = pybam.read('/my/data.bam',['sam_seq','sam_mapq'])
    print my_bam._static_parser_code
    for seq,mapq in my_bam:
        if seq.startswith('ACGT') and mapq > 10:
            print my_bam.sam

    [ Custom Decompressor (from file path) Example ]
    my_bam = pybam.read('/my/data.bam.lzma',decompressor='lzma --decompress --stdout /my/data.bam.lzma')

    [ Custom Decompressor (from file object) Example ]
    my_bam = pybam.read(sys.stdin,decompressor='lzma --decompress --stdout') # data given to lzma via stdin

    [ Force Internal bgzip Decompressor ]
    my_bam = pybam.read('/my/data.bam',decompressor='internal')

    "print pybam.wat" in the python terminal to see the possible parsable values,
    or visit http://github.com/JohnLonginotto/pybam for the latest info.
    '''

    def __init__(self, f, fields=False, decompressor=False):
        self.file_bytes_read = 0
        self.file_chromosomes = []
        self.file_alignments_read = 0
        self.file_chromosome_lengths = {}

        if fields is not False:
            if type(fields) is not list or len(fields) is 0:
                raise PybamError(
                    '\n\nFields for the static parser must be provided as a non-empty list. You gave a ' + str(type(fields)) + '\n')
            else:
                for field in fields:
                    if field.startswith('sam') or field.startswith('bam'):
                        if field not in parse_codes.keys():
                            raise PybamError('\n\nStatic parser field "' + str(field) + '" from fields ' + str(
                                fields) + ' is not known to this version of pybam!\nPrint "pybam.wat" to see available field names with explinations.\n')
                    else:
                        raise PybamError('\n\nStatic parser field "' + str(field) + '" from fields ' + str(
                            fields) + ' does not start with "sam" or "bam" and thus is not an avaliable field for the static parsing.\nPrint "pybam.wat" in interactive python to see available field names with explinations.\n')

        if decompressor:
            if type(decompressor) is str:
                if decompressor is not 'internal' and '{}' not in decompressor:
                    raise PybamError(
                        '\n\nWhen a custom decompressor is used and the input file is a string, the decompressor string must contain at least one occurence of "{}" to be substituted with a filepath by pybam.\n')
            else:
                raise PybamError('\n\nUser-supplied decompressor must be a string that when run on the command line decompresses a named file (or stdin), to stdout:\ne.g. "lzma --decompress --stdout {}" if pybam is provided a path as input file, where {} is substituted for that path.\nor just "lzma --decompress --stdout" if pybam is provided a file object instead of a file path, as data from that file object will be piped via stdin to the decompression program.\n')

        # First we make a generator that will return chunks of uncompressed data, regardless of how we choose to decompress:
        def generator():
            DEVNULL = open(os.devnull, 'wb')

            # First we need to figure out what sort of file we have - whether it's gzip compressed, uncompressed, or something else entirely!
            if type(f) is str:
                try:
                    self._file = open(f, 'rb')
                except:
                    raise PybamError('\n\nCould not open "' +
                                     str(self._file.name) + '" for reading!\n')
                try:
                    magic = os.read(self._file.fileno(), 4)
                except:
                    raise PybamError('\n\nCould not read from "' +
                                     str(self._file.name) + '"!\n')
            elif type(f) is file:
                self._file = f
                try:
                    magic = os.read(self._file.fileno(), 4)
                except:
                    raise PybamError('\n\nCould not read from "' +
                                     str(self._file.name) + '"!\n')
            else:
                raise PybamError(
                    '\n\nInput file was not a string or a file object. It was: "' + str(f) + '"\n')

            self.file_name = os.path.basename(
                os.path.realpath(self._file.name))
            self.file_directory = os.path.dirname(
                os.path.realpath(self._file.name))

            if magic == 'BAM\1':
                # The user has passed us already unzipped BAM data! Job done :)
                data = 'BAM\1' + self._file.read(35536)
                self.file_bytes_read += len(data)
                self.file_decompressor = 'None'
                while data:
                    yield data
                    data = self._file.read(35536)
                    self.file_bytes_read += len(data)
                self._file.close()
                DEVNULL.close()
                raise StopIteration

            elif magic == "\x1f\x8b\x08\x04":  # The user has passed us compressed gzip/bgzip data, which is typical for a BAM file
                # use custom decompressor if provided:
                if decompressor is not False and decompressor is not 'internal':
                    if type(f) is str:
                        self._subprocess = subprocess.Popen(decompressor.replace(
                            '{}', f),    shell=True, stdout=subprocess.PIPE, stderr=DEVNULL)
                    else:
                        self._subprocess = subprocess.Popen(
                            '{ printf "'+magic+'"; cat; } | ' + decompressor, stdin=self._file, shell=True, stdout=subprocess.PIPE, stderr=DEVNULL)
                    self.file_decompressor = decompressor
                    data = self._subprocess.stdout.read(35536)
                    self.file_bytes_read += len(data)
                    while data:
                        yield data
                        data = self._subprocess.stdout.read(35536)
                        self.file_bytes_read += len(data)
                    self._file.close()
                    DEVNULL.close()
                    raise StopIteration

                # else look for pigz or gzip:
                else:
                    try:
                        self._subprocess = subprocess.Popen(
                            ["pigz"], stdin=DEVNULL, stdout=DEVNULL, stderr=DEVNULL)
                        if self._subprocess.returncode is None:
                            self._subprocess.kill()
                        use = 'pigz'
                    except OSError:
                        try:
                            self._subprocess = subprocess.Popen(
                                ["gzip"], stdin=DEVNULL, stdout=DEVNULL, stderr=DEVNULL)
                            if self._subprocess.returncode is None:
                                self._subprocess.kill()
                            use = 'gzip'
                        except OSError:
                            use = 'internal'

                    if use is not 'internal' and decompressor is not 'internal':
                        if type(f) is str:
                            self._subprocess = subprocess.Popen(
                                [use, '--decompress', '--stdout',       f], stdout=subprocess.PIPE, stderr=DEVNULL)
                        else:
                            self._subprocess = subprocess.Popen(
                                '{ printf "'+magic+'"; cat; } | ' + use + ' --decompress  --stdout', stdin=f, shell=True, stdout=subprocess.PIPE, stderr=DEVNULL)
                        time.sleep(1)
                        if self._subprocess.poll() == None:
                            data = self._subprocess.stdout.read(35536)
                            self.file_decompressor = use
                            self.file_bytes_read += len(data)
                            while data:
                                yield data
                                data = self._subprocess.stdout.read(35536)
                                self.file_bytes_read += len(data)
                            self._file.close()
                            DEVNULL.close()
                            raise StopIteration

                    # Python's gzip module can't read from a stream that doesn't support seek(), and the zlib module cannot read the bgzip format without a lot of help:
                    self.file_decompressor = 'internal'
                    raw_data = magic + self._file.read(65536)
                    self.file_bytes_read = len(raw_data)
                    internal_cache = []
                    blocks_left_to_grab = 50
                    bs = 0
                    checkpoint = 0
                    decompress = zlib.decompress
                    while raw_data:
                        if len(raw_data) - bs < 35536:
                            raw_data = raw_data[bs:] + self._file.read(65536)
                            self.file_bytes_read += len(raw_data) - bs
                            bs = 0
                        magic = raw_data[bs:bs+4]
                        if not magic:
                            break  # a child's heart
                        if magic != "\x1f\x8b\x08\x04":
                            raise PybamError(
                                '\n\nThe input file is not in a format I understand. First four bytes: ' + repr(magic) + '\n')
                        try:
                            more_bs = bs + \
                                unpack("<H", raw_data[bs+16:bs+18])[0] + 1
                            internal_cache.append(decompress(
                                raw_data[bs+18:more_bs-8], -15))
                            bs = more_bs
                        except:  # zlib doesnt have a nice exception for when things go wrong. just "error"
                            header_data = magic + raw_data[bs+4:bs+12]
                            header_size = 12
                            extra_len = unpack("<H", header_data[-2:])[0]
                            while header_size-12 < extra_len:
                                header_data += raw_data[bs+12:bs+16]
                                subfield_id = header_data[-4:-2]
                                subfield_len = unpack(
                                    "<H", header_data[-2:])[0]
                                subfield_data = raw_data[bs +
                                                         16:bs+16+subfield_len]
                                header_data += subfield_data
                                header_size += subfield_len + 4
                                if subfield_id == 'BC':
                                    block_size = unpack("<H", subfield_data)[0]
                            raw_data = raw_data[bs+16+subfield_len:bs +
                                                16+subfield_len+block_size-extra_len-19]
                            # I have left the numbers in verbose, because the above try is the optimised code.
                            crc_data = raw_data[bs+16+subfield_len+block_size -
                                                extra_len-19:bs+16+subfield_len+block_size-extra_len-19+8]
                            bs = bs+16+subfield_len+block_size-extra_len-19+8
                            zipped_data = header_data + raw_data + crc_data
                            # 31 works the same as 47.
                            internal_cache.append(decompress(zipped_data, 47))
                            # Although the following in the bgzip code from biopython, its not needed if you let zlib decompress the whole zipped_data, header and crc, because it checks anyway (in C land)
                            # I've left the manual crc checks in for documentation purposes:
                            '''
                            expected_crc = crc_data[:4]
                            expected_size = unpack("<I", crc_data[4:])[0]
                            if len(unzipped_data) != expected_size: print 'ERROR: Failed to unpack due to a Type 1 CRC error. Could the BAM be corrupted?'; exit()
                            crc = zlib.crc32(unzipped_data)
                            if crc < 0: crc = pack("<i", crc)
                            else:       crc = pack("<I", crc)
                            if expected_crc != crc: print 'ERROR: Failed to unpack due to a Type 2 CRC error. Could the BAM be corrupted?'; exit()
                            '''
                        blocks_left_to_grab -= 1
                        if blocks_left_to_grab == 0:
                            yield ''.join(internal_cache)
                            internal_cache = []
                            blocks_left_to_grab = 50
                    self._file.close()
                    DEVNULL.close()
                    if internal_cache != '':
                        yield ''.join(internal_cache)
                    raise StopIteration

            elif decompressor is not False and decompressor is not 'internal':
                # It wouldn't be safe to just print to the shell four random bytes from the beginning of a file, so instead it's
                # written to a temp file and cat'd. The idea here being that we trust the decompressor string as it was written by
                # someone with access to python, so it has system access anyway. The file/data, however, should not be trusted.
                magic_file = os.path.join(tempfile.mkdtemp(), 'magic')
                with open(magic_file, 'wb') as mf:
                    mf.write(magic)
                if type(f) is str:
                    self._subprocess = subprocess.Popen(decompressor.replace(
                        '{}', f),    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                else:
                    self._subprocess = subprocess.Popen(
                        '{ cat "'+magic_file+'"; cat; } | ' + decompressor, stdin=self._file, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.file_decompressor = decompressor
                data = self._subprocess.stdout.read(35536)
                self.file_bytes_read += len(data)
                while data:
                    yield data
                    data = self._subprocess.stdout.read(35536)
                    self.file_bytes_read += len(data)
                self._file.close()
                DEVNULL.close()
                raise StopIteration
            else:
                raise PybamError(
                    '\n\nThe input file is not in a format I understand. First four bytes: ' + repr(magic) + '\n')

        # At this point, we know that whatever decompression method was used, a call to self._generator will return some uncompressed data.
        self._generator = generator()

        # So lets parse the BAM header:
        header_cache = ''
        while len(header_cache) < 8:
            header_cache += next(self._generator)

        p_from = 0
        p_to = 4
        if header_cache[p_from:p_to] != 'BAM\1':
            raise PybamError('\n\nInput file ' + self.file_name +
                             ' does not appear to be a BAM file.\n')

        # Parse the BAM header:
        p_from = p_to
        p_to += 4
        length_of_header = unpack('<i', header_cache[p_from:p_to])[0]
        p_from = p_to
        p_to += length_of_header
        while len(header_cache) < p_to:
            header_cache += next(self._generator)
        self.file_header = header_cache[p_from:p_to]
        p_from = p_to
        p_to += 4
        while len(header_cache) < p_to:
            header_cache += next(self._generator)
        number_of_reference_sequences = unpack(
            '<i', header_cache[p_from:p_to])[0]

        for _ in range(number_of_reference_sequences):
            p_from = p_to
            p_to += 4
            while len(header_cache) < p_to:
                header_cache += next(self._generator)
            l_name = unpack('<l', header_cache[p_from:p_to])[0]
            p_from = p_to
            p_to += l_name
            while len(header_cache) < p_to:
                header_cache += next(self._generator)
            self.file_chromosomes.append(header_cache[p_from:p_to - 1])
            p_from = p_to
            p_to += 4
            while len(header_cache) < p_to:
                header_cache += next(self._generator)
            self.file_chromosome_lengths[self.file_chromosomes[-1]
                                         ] = unpack('<l', header_cache[p_from:p_to])[0]

        self.file_bytes_read = p_to
        self.file_binary_header = buffer(header_cache[:p_to])
        header_cache = header_cache[p_to:]

        # A quick check to make sure the header of this BAM file makes sense:
        chromosomes_from_header = []
        for line in self.file_header.split('\n'):
            if line.startswith('@SQ\tSN:'):
                chromosomes_from_header.append(line.split('\t')[1][3:])
        if chromosomes_from_header != self.file_chromosomes:
            raise PybamWarn('For some reason the BAM format stores the chromosome names in two locations,\n       the ASCII text header we all know and love, viewable with samtools view -H, and another special binary header\n       which is used to translate the chromosome refID (a number) into a chromosome RNAME when you do bam -> sam.\n\nThese two headers should always be the same, but apparently they are not:\nThe ASCII header looks like: ' +
                            self.file_header + '\nWhile the binary header has the following chromosomes: ' + self.file_chromosomes + '\n')

        # Variable parsing:
        def new_entry(header_cache):
            # we keep a small cache of X bytes of decompressed BAM data, to smoothen out disk access.
            cache = header_cache
            p = 0  # where the next alignment/entry starts in the cache
            while True:
                try:
                    while len(cache) < p + 4:
                        cache = cache[p:] + next(self._generator)
                        p = 0  # Grab enough bytes to parse blocksize
                    self.sam_block_size = unpack('<i', cache[p:p+4])[0]
                    self.file_alignments_read += 1
                    while len(cache) < p + 4 + self.sam_block_size:
                        cache = cache[p:] + next(self._generator)
                        p = 0  # Grab enough bytes to parse entry
                except StopIteration:
                    break
                self.bam = cache[p:p + 4 + self.sam_block_size]
                p = p + 4 + self.sam_block_size
                yield self
        self._new_entry = new_entry(header_cache)

        def compile_parser(self, fields):
            temp_code = ''
            end_of_qname = False
            end_of_cigar = False
            end_of_seq = False
            end_of_qual = False
            dependencies = set(fields)

            if 'bam' in fields:
                fields[fields.index('bam')] = 'self.bam'

            if 'sam_block_size' in fields:
                fields[fields.index('sam_block_size')] = 'self.sam_block_size'

            if 'sam' in dependencies:
                dependencies.update(['sam_qname', 'sam_flag', 'sam_rname', 'sam_pos1', 'sam_mapq', 'sam_cigar_string', 'bam_refID',
                                     'bam_next_refID', 'sam_rnext', 'sam_pnext1', 'sam_tlen', 'sam_seq', 'sam_qual', 'sam_tags_string'])

            if 'sam_tags_string' in dependencies:
                dependencies.update(['sam_tags_list'])

            if 'sam_pos1' in dependencies:
                temp_code += "\n        sam_pos1 = (0 if sam_pos0 < 0 else sam_pos0 + 1)"
                dependencies.update(['sam_pos0'])

            if 'sam_pnext1' in dependencies:
                temp_code += "\n        sam_pnext1 = (0 if sam_pnext0 < 0 else sam_pnext0 + 1)"
                dependencies.update(['sam_pnext0'])

            if 'sam_qname' in dependencies or 'bam_qname' in dependencies:
                temp_code += "\n        _end_of_qname = 36 + sam_l_read_name"
                dependencies.update(['sam_l_read_name'])
                end_of_qname = True

            if 'sam_cigar_string' in dependencies or 'sam_cigar_list' in dependencies or 'bam_cigar' in dependencies:
                if end_of_qname:
                    pass
                else:
                    temp_code += "\n        _end_of_qname = 36 + sam_l_read_name"
                temp_code += "\n        _end_of_cigar = _end_of_qname + (4*sam_n_cigar_op)"
                dependencies.update(['sam_l_read_name', 'sam_n_cigar_op'])
                end_of_cigar = True

            if 'sam_seq' in dependencies or 'bam_seq' in dependencies:
                if end_of_cigar:
                    pass
                elif end_of_qname:
                    temp_code += "\n        _end_of_cigar = _end_of_qname + (4*sam_n_cigar_op)"
                else:
                    temp_code += "\n        _end_of_cigar = 36 + sam_l_read_name + (4*sam_n_cigar_op)"
                temp_code += "\n        _end_of_seq = _end_of_cigar + (-((-sam_l_seq)//2))"
                dependencies.update(
                    ['sam_l_seq', 'sam_n_cigar_op', 'sam_l_read_name'])
                end_of_seq = True

            if 'sam_qual' in dependencies or 'bam_qual' in dependencies:
                if end_of_seq:
                    pass
                elif end_of_cigar:
                    temp_code += "\n        _end_of_seq   = _end_of_cigar + (-((-sam_l_seq)//2))"
                elif end_of_qname:
                    temp_code += "\n        _end_of_seq   = _end_of_qname + (4*sam_n_cigar_op) + (-((-sam_l_seq)//2))"
                else:
                    temp_code += "\n        _end_of_seq   = 36 + sam_l_read_name + (4*sam_n_cigar_op) + (-((-sam_l_seq)//2))"
                temp_code += "\n        _end_of_qual = _end_of_seq + sam_l_seq"
                dependencies.update(
                    ['sam_l_seq', 'sam_n_cigar_op', 'sam_l_read_name'])
                end_of_qual = True

            if 'sam_tags_list' in dependencies or 'bam_tags' in dependencies:
                if end_of_qual:
                    pass
                elif end_of_seq:
                    temp_code += "\n        _end_of_qual = _end_of_seq + sam_l_seq"
                elif end_of_cigar:
                    temp_code += "\n        _end_of_qual = _end_of_cigar + (-((-sam_l_seq)//2)) + sam_l_seq"
                elif end_of_qname:
                    temp_code += "\n        _end_of_qual = _end_of_qname + (4*sam_n_cigar_op) + (-((-sam_l_seq)//2)) + sam_l_seq"
                else:
                    temp_code += "\n        _end_of_qual = 36 + sam_l_read_name + (4*sam_n_cigar_op) + (-((-sam_l_seq)//2)) + sam_l_seq"
                dependencies.update(
                    ['sam_l_seq', 'sam_n_cigar_op', 'sam_l_read_name'])

            if 'sam_rname' in dependencies:
                temp_code += "\n        sam_rname = '*' if sam_refID < 0 else self.file_chromosomes[sam_refID]"
                dependencies.update(['sam_refID'])

            if 'sam_rnext' in dependencies:
                temp_code += "\n        sam_rnext = '*' if sam_next_refID < 0 else self.file_chromosomes[sam_next_refID]"
                dependencies.update(['sam_next_refID'])

            # First we figure out what data from the static portion of the BAM entry we'll need:
            tmp = {}
            tmp['code'] = 'def parser(self):\n    from array import array\n    from struct import unpack\n    for _ in self._new_entry:'
            tmp['last_start'] = None
            tmp['name_list'] = []
            tmp['dtype_list'] = []

            def pack_up(name, dtype, length, end, tmp):
                if name in dependencies:
                    if tmp['last_start'] is None:
                        tmp['last_start'] = end - length
                    tmp['name_list'].append(name)
                    tmp['dtype_list'].append(dtype)
                elif tmp['last_start'] is not None:
                    tmp['code'] += '\n        ' + ', '.join(tmp['name_list']) + ' = unpack("<' + ''.join(
                        tmp['dtype_list']) + '",self.bam[' + str(tmp['last_start']) + ':' + str(end-length) + '])'
                    if len(tmp['dtype_list']) == 1:
                        tmp['code'] += '[0]'
                    tmp['last_start'] = None
                    tmp['name_list'] = []
                    tmp['dtype_list'] = []

            pack_up('sam_refID',       'i', 4, 8, tmp)
            pack_up('sam_pos0',        'i', 4, 12, tmp)
            pack_up('sam_l_read_name', 'B', 1, 13, tmp)
            pack_up('sam_mapq',        'B', 1, 14, tmp)
            pack_up('sam_bin',         'H', 2, 16, tmp)
            pack_up('sam_n_cigar_op',  'H', 2, 18, tmp)
            pack_up('sam_flag',        'H', 2, 20, tmp)
            pack_up('sam_l_seq',       'i', 4, 24, tmp)
            pack_up('sam_next_refID',  'i', 4, 28, tmp)
            pack_up('sam_pnext0',      'i', 4, 32, tmp)
            pack_up('sam_tlen',        'i', 4, 36, tmp)
            # To add anything not yet added.
            pack_up(None,            None, 0, 36, tmp)
            code = tmp['code']
            del tmp

            code += temp_code

            # Fixed-length BAM data (where we just grab the bytes, we dont unpack) can, however, be grabbed individually.
            if 'bam_block_size' in dependencies:
                code += "\n        bam_block_size   = self.bam[0             : 4                ]"
            if 'bam_refID' in dependencies:
                code += "\n        bam_refID        = self.bam[4             : 8                ]"
            if 'bam_pos' in dependencies:
                code += "\n        bam_pos          = self.bam[8             : 12               ]"
            if 'bam_l_read_name' in dependencies:
                code += "\n        bam_l_read_name  = self.bam[12            : 13               ]"
            if 'bam_mapq' in dependencies:
                code += "\n        bam_mapq         = self.bam[13            : 14               ]"
            if 'bam_bin' in dependencies:
                code += "\n        bam_bin          = self.bam[14            : 16               ]"
            if 'bam_n_cigar_op' in dependencies:
                code += "\n        bam_n_cigar_op   = self.bam[16            : 18               ]"
            if 'bam_flag' in dependencies:
                code += "\n        bam_flag         = self.bam[18            : 20               ]"
            if 'bam_l_seq' in dependencies:
                code += "\n        bam_l_seq        = self.bam[20            : 24               ]"
            if 'bam_next_refID' in dependencies:
                code += "\n        bam_next_refID   = self.bam[24            : 28               ]"
            if 'bam_pnext' in dependencies:
                code += "\n        bam_pnext        = self.bam[28            : 32               ]"
            if 'bam_tlen' in dependencies:
                code += "\n        bam_tlen         = self.bam[32            : 36               ]"
            if 'bam_qname' in dependencies:
                code += "\n        bam_qname        = self.bam[36            : _end_of_qname    ]"
            if 'bam_cigar' in dependencies:
                code += "\n        bam_cigar        = self.bam[_end_of_qname : _end_of_cigar    ]"
            if 'bam_seq' in dependencies:
                code += "\n        bam_seq          = self.bam[_end_of_cigar : _end_of_seq      ]"
            if 'bam_qual' in dependencies:
                code += "\n        bam_qual         = self.bam[_end_of_seq   : _end_of_qual     ]"
            if 'bam_tags' in dependencies:
                code += "\n        bam_tags         = self.bam[_end_of_qual  :                  ]"

            if 'sam_qname' in dependencies:
                if 'bam_qname' in dependencies:
                    code += "\n        sam_qname        = bam_qname[:-1]"
                else:
                    code += "\n        sam_qname        = self.bam[36            : _end_of_qname -1 ]"

            if 'sam_cigar_list' in dependencies:
                if 'bam_cigar' in dependencies:
                    code += "\n        sam_cigar_list = [( cig >> 4  , cigar_codes[cig & 0b1111]) for cig in array('I', bam_cigar) ]"
                else:
                    code += "\n        sam_cigar_list = [( cig >> 4  , cigar_codes[cig & 0b1111]) for cig in array('I', self.bam[_end_of_qname : _end_of_cigar]) ]"

            if 'sam_cigar_string'in dependencies:
                if 'bam_cigar' in dependencies:
                    code += "\n        sam_cigar_string = ''.join([        str(cig >> 4) + cigar_codes[cig & 0b1111] for cig     in array('I', bam_cigar)])"
                else:
                    code += "\n        sam_cigar_string = ''.join([        str(cig >> 4) + cigar_codes[cig & 0b1111] for cig     in array('I', self.bam[_end_of_qname : _end_of_cigar]) ])"

            if 'sam_seq' in dependencies:
                if 'bam_seq' in dependencies:
                    code += "\n        sam_seq = ''.join( [ dna_codes[dna >> 4] + dna_codes[dna & 0b1111]   for dna     in array('B', bam_seq)])[:sam_l_seq]"
                else:
                    code += "\n        sam_seq = ''.join( [ dna_codes[dna >> 4] + dna_codes[dna & 0b1111]   for dna     in array('B', self.bam[_end_of_cigar : _end_of_seq])])[:sam_l_seq]"

            if 'sam_qual' in dependencies:
                if 'bam_qual' in dependencies:
                    code += "\n        sam_qual = ''.join( [                   chr(ord(quality) + 33)        for quality in            bam_qual ])"
                else:
                    code += "\n        sam_qual = ''.join( [                   chr(ord(quality) + 33)        for quality in            self.bam[_end_of_seq       : _end_of_qual  ]])"

            if 'sam_tags_list' in dependencies:
                code += '''
        sam_tags_list = []
        offset = _end_of_qual
        while offset != len(self.bam):
            tag_name = self.bam[offset:offset+2]
            tag_type = self.bam[offset+2]
            if tag_type == 'Z':
                offset_end = self.bam.index('\\0',offset+3)+1
                tag_data = self.bam[offset+3:offset_end-1]
            elif tag_type in CtoPy:
                offset_end = offset+3+py4py[tag_type]
                tag_data = unpack(CtoPy[tag_type],self.bam[offset+3:offset_end])[0]
            elif tag_type == 'B':
                offset_end = offset+8+(unpack('<i',self.bam[offset+4:offset+8])[0]*py4py[self.bam[offset+3]])
                tag_data = array(self.bam[offset+3] , self.bam[offset+8:offset_end] )
            else:
                print 'PYBAM ERROR: I dont know how to parse BAM tags in this format: ',repr(tag_type)
                print '             This is simply because I never saw this kind of tag during development.'
                print '             If you could mail the following chunk of text to john at john.uk.com, i will fix this up for everyone :)'
                print repr(tag_type),repr(self.bam[offset+3:end])
                exit()
            sam_tags_list.append((tag_name,tag_type,tag_data))
            offset = offset_end'''

            if 'sam_tags_string' in dependencies:
                code += "\n        sam_tags_string = '\t'.join(A + ':' + ('i' if B in 'cCsSI' else B)  + ':' + ((C.typecode + ',' + ','.join(map(str,C))) if type(C)==array else str(C)) for A,B,C in self.sam_tags_list)"

            if 'sam' in dependencies:
                code += "\n        sam = sam_qname + '\t' + str(sam_flag) + '\t' + sam_rname + '\t' + str(sam_pos1) + '\t' + str(sam_mapq) + '\t' + ('*' if sam_cigar_string == '' else sam_cigar_string) + '\t' + ('=' if bam_refID == bam_next_refID else sam_rnext) + '\t' + str(sam_pnext1) + '\t' + str(sam_tlen) + '\t' + sam_seq + '\t' + sam_qual + '\t' + sam_tags_string"

            code += '\n        yield ' + ','.join([x for x in fields]) + '\n'

            # "code" is the static parser's code as a string (a function called "parser")
            self._static_parser_code = code
            exec_dict = {                   # This dictionary stores things the exec'd code needs to know about, and will store the compiled function after exec()
                'unpack': unpack,
                'array': array,
                'dna_codes': dna_codes,
                'CtoPy': CtoPy,
                'py4py': py4py,
                'cigar_codes': cigar_codes
            }
            # exec() compiles "code" to real code, creating the "parser" function and adding it to exec_dict['parser']
            exec code in exec_dict
            return exec_dict['parser']

        if fields:
            static_parser = compile_parser(self, fields)(self)

            def next_read(): return next(static_parser)
        else:
            def next_read(): return next(self._new_entry)
        self.next = next_read

    def __iter__(self): return self

    def __str__(self): return self.sam

    # Methods to pull out raw bam data from entry (so still in its binary encoding). This can be helpful in some scenarios.
    @property
    def bam_block_size(self): return self.bam[: 4]

    @property
    def bam_refID(self): return self.bam[4: 8]

    @property
    def bam_pos(self): return self.bam[8: 12]

    @property
    def bam_l_read_name(self): return self.bam[12: 13]

    @property
    def bam_mapq(self): return self.bam[13: 14]

    @property
    def bam_bin(self): return self.bam[14: 16]

    @property
    def bam_n_cigar_op(self): return self.bam[16: 18]

    @property
    def bam_flag(self): return self.bam[18: 20]

    @property
    def bam_l_seq(self): return self.bam[20: 24]

    @property
    def bam_next_refID(self): return self.bam[24: 28]

    @property
    def bam_pnext(self): return self.bam[28: 32]

    @property
    def bam_tlen(self): return self.bam[32: 36]

    @property
    def bam_qname(self): return self.bam[36: self._end_of_qname]

    @property
    def bam_cigar(
        self): return self.bam[self._end_of_qname: self._end_of_cigar]

    @property
    def bam_seq(self): return self.bam[self._end_of_cigar: self._end_of_seq]

    @property
    def bam_qual(self): return self.bam[self._end_of_seq: self._end_of_qual]

    @property
    def bam_tags(self): return self.bam[self._end_of_qual:]

    @property
    def sam_refID(self): return unpack('<i', self.bam[4:  8])[0]

    @property
    def sam_pos0(self): return unpack('<i', self.bam[8: 12])[0]

    @property
    def sam_l_read_name(self): return unpack('<B', self.bam[12: 13])[0]

    @property
    def sam_mapq(self): return unpack('<B', self.bam[13: 14])[0]

    @property
    def sam_bin(self): return unpack('<H', self.bam[14: 16])[0]

    @property
    def sam_n_cigar_op(self): return unpack('<H', self.bam[16: 18])[0]

    @property
    def sam_flag(self): return unpack('<H', self.bam[18: 20])[0]

    @property
    def sam_l_seq(self): return unpack('<i', self.bam[20: 24])[0]

    @property
    def sam_next_refID(self): return unpack('<i', self.bam[24: 28])[0]

    @property
    def sam_pnext0(self): return unpack('<i', self.bam[28: 32])[0]

    @property
    def sam_tlen(self): return unpack('<i', self.bam[32: 36])[0]

    @property
    # -1 to remove trailing NUL byte
    def sam_qname(self): return self.bam[36: self._end_of_qname - 1]

    @property
    def sam_cigar_list(self): return [(cig >> 4, cigar_codes[cig & 0b1111]) for cig in array(
        'I', self.bam[self._end_of_qname: self._end_of_cigar])]

    @property
    def sam_cigar_string(self): return ''.join([str(cig >> 4) + cigar_codes[cig & 0b1111]
                                                for cig in array('I', self.bam[self._end_of_qname: self._end_of_cigar])])

    @property
    def sam_seq(self): return ''.join([dna_codes[dna >> 4] + dna_codes[dna & 0b1111] for dna in array('B', self.bam[self._end_of_cigar: self._end_of_seq])])[
        :self.sam_l_seq]  # As DNA is 4 bits packed 2-per-byte, there might be a trailing '0000', so we can either

    @property
    def sam_qual(self): return ''.join(
        [chr(ord(quality) + 33) for quality in self.bam[self._end_of_seq: self._end_of_qual]])

    @property
    def sam_tags_list(self):
        result = []
        offset = self._end_of_qual
        while offset != len(self.bam):
            tag_name = self.bam[offset:offset+2]
            tag_type = self.bam[offset+2]
            if tag_type == 'Z':
                offset_end = self.bam.index('\x00', offset+3)+1
                tag_data = self.bam[offset+3:offset_end-1]
            elif tag_type in CtoPy:
                offset_end = offset+3+py4py[tag_type]
                tag_data = unpack(
                    CtoPy[tag_type], self.bam[offset+3:offset_end])[0]
            elif tag_type == 'B':
                offset_end = offset+8 + \
                    (unpack('<i', self.bam[offset+4:offset+8])
                     [0]*py4py[self.bam[offset+3]])
                tag_data = array(self.bam[offset+3],
                                 self.bam[offset+8:offset_end])
            else:
                print 'PYBAM ERROR: I dont know how to parse BAM tags in this format: ', repr(
                    tag_type)
                print '             This is simply because I never saw this kind of tag during development.'
                print '             If you could mail the following chunk of text to john at john.uk.com, ill fix this up :)'
                print repr(tag_type), repr(self.bam[offset+3:end])
                exit()
            result.append((tag_name, tag_type, tag_data))
            offset = offset_end
        return result

    @property
    def sam_tags_string(self):
        return '\t'.join(A + ':' + ('i' if B in 'cCsSI' else B) + ':' + ((C.typecode + ',' + ','.join(map(str, C))) if type(C) == array else str(C)) for A, B, C in self.sam_tags_list)

    # BONUS methods - methods that mimic how samtools works.
    @property
    def sam_pos1(self): return 0 if self.sam_pos0 < 0 else self.sam_pos0 + 1

    @property
    def sam_pnext1(
        self): return 0 if self.sam_pnext0 < 0 else self.sam_pnext0 + 1

    @property
    def sam_rname(
        self): return '*' if self.sam_refID < 0 else self.file_chromosomes[self.sam_refID]

    @property
    def sam_rnext(
        self): return '*' if self.sam_next_refID < 0 else self.file_chromosomes[self.sam_next_refID]

    @property
    def sam(self): return (
        self.sam_qname + '\t' +
        str(self.sam_flag) + '\t' +
        self.sam_rname + '\t' +
        str(self.sam_pos1) + '\t' +
        str(self.sam_mapq) + '\t' +
        ('*' if self.sam_cigar_string == '' else self.sam_cigar_string) + '\t' +
        ('=' if self.bam_refID == self.bam_next_refID else self.sam_rnext) + '\t' +
        str(self.sam_pnext1) + '\t' +
        str(self.sam_tlen) + '\t' +
        self.sam_seq + '\t' +
        self.sam_qual + '\t' +
        self.sam_tags_string
    )

    # Internal methods - methods used to calculate where variable-length blocks start/end
    @property
    # fixed-length stuff at the beginning takes up 36 bytes.
    def _end_of_qname(self): return self.sam_l_read_name + 36

    @property
    def _end_of_cigar(self): return self._end_of_qname + \
        (4*self.sam_n_cigar_op)   # 4 bytes per n_cigar_op

    @property
    def _end_of_seq(self): return self._end_of_cigar + \
        (-((-self.sam_l_seq)//2))  # {blurgh}

    @property
    def _end_of_qual(self): return self._end_of_seq + \
        self.sam_l_seq            # qual has the same length as seq

    def __del__(self):
        if self._subprocess.returncode is None:
            self._subprocess.kill()
        self._file.close()


class PybamWarn(Exception):
    pass


class PybamError(Exception):
    pass
