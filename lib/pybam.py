import os
import sys
import zlib
import struct
import itertools
import subprocess

##########################################
## Documentation ? What documentation.. ##
#########################################################################################################################################################
##                                                                                                                                                     ##
## Thank you for taking an interest in this program!                                                                                                   ##
## Note that this is really a proof-of-concept rather than a fully-fledged alternative to any of the existing BAM -> Python shared modules like        ##
## pysam and htspython. I just wanted to show that pure python is often as-fast-or-faster than shared C librarys, particularly if you're running pypy. ##
## This is because every time python has to go out of the VM and talk to external things (like htslib), we take a huge performance hit. Its the        ##
## program equivilent of walking to your car to go to the shops, when you might be able to cycle straight there in half the time.                      ##
##                                                                                                                                                     ##
## Particularly if your bicycle has rockets on it because you installed pypy.                                                                          ##
##                                                                                                                                                     ##
## If you would like to contribute to the project in anyway, please do. Better support for tags is greatly needed (ideally we need to find some data   ##
## with some weird tag formats so we can figure out how to parse them), as well as more time spent improving the interface. If you dont like the name  ##
## of something, or the way something looks/feels when you use pybam (or even the name pybam) we can change that :)                                    ##
##                                                                                                                                                     ##
#########################################################################################################################################################

'''
Awesome people who have directly contributed to the project so far:

'''

#############
## bgunzip ##
#########################################################################################################################################################
##                                                                                                                                                     ##
## There are only two parts to pybam. The first is the bgunzip class that will take either pure BAM or bgzip'd BAM data, decompress the blocks to pure ##
## BAM, parse out the BAM header, and from then on act like a generator - where every time you iterate it, it give you a random chunk of BAM data,     ##
## starting from the first read in the file (but blocks end randomly due to the way bgzip was originally designed to be data agnostic).                ##
## So you initialize with a file handle or path, and an optional blocks_at_a_time parameter (how many bgzip blocks to decompress and give you back per ##
## iteration), like so:                                                                                                                                ##
##                                                                                                                                                     ##
## >>>> pure_bam_data = pybam.bgunzip('./ENCFF001LCU.bam.gz')                                                                                          ##
## >>>> pure_bam_data.chromosome_names                                                                                                                 ##
## ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17' ... ##
## >>>> pure_bam_data.chromosome_lengths                                                                                                               ##
## [197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172, 129993255, 121843856, 121257530, 120284312, ... ##
## >>>> pure_bam_data.chromosomes_from_header                                                                                                          ##
## ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17' ... ##
## >>>> pure_bam_data.header_text                                                                                                                      ##
## '@SQ\tSN:chr1\tLN:197195432\tAS:mm9\tSP:mouse\n@SQ\tSN:chr2\tLN:181748087\tAS:mm9\tSP:mouse\n@SQ\tSN:chr3\tLN:159599783\tAS:mm9\tSP:mouse\n@SQ\ ... ##
## >>>> pure_bam_data.original_binary_header                                                                                                           ##
## <read-only buffer for 0x00007fbe9c1344e0, size 1161>                                                                                                ##
##                                                                                                                                                     ##
## You can use the original_binary_header to make a new BAM file, if that's your thing, but I wouldn't recommend it. Use the bam+ format instead.      ##
##                                                                                                                                                     ##
## Now you can iterate it like any other generator:                                                                                                    ##
## >>>> first_50_blocks = next(pure_bam_data)    # 50 because that is the blocks_at_a_time default, or about 3Mb of BAM data                           ##
## >>>> len(first_50_blocks)                                                                                                                           ##
## 3275639                                                                                                                                             ##
## >>>> first_50_blocks[:70]                                                                                                                           ##
## '\x92\x00\x00\x00\x00\x00\x00\x00\xa6\xc9-\x00 \x00\x00\x13\x01\x00\x00\x00$\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\x00\x00\x00\x00SOLEXA1 ... ##
##                                                                                                                                                     ##
## What a time to be alive.                                                                                                                            ##
##                                                                                                                                                     ##
## Truly, none of this could have been done without the already fantastic work by Peter Cock, author of the bgzip format and bgzf.py from Biopython.   ##
## I'm sure he wouldn't mind me posting a link to his fantastic Bioinformatics blog, particularly this  post about bgzip and how it works for any kind ##
## of data, not just BAM data: http://blastedbio.blogspot.de/2011/11/bgzf-blocked-bigger-better-gzip.html                                              ##
## Although I didn't actually use any of Peter's code directly, I didn't need to, because he documents his code SO WELL in bgzf.py that I "got it" on  ##
## the first read-through, and could re-impliment it optimized for PyPy which is what you see here below.                                              ##
#########################################################################################################################################################
class bgunzip:
    def __init__(self,file_handle,blocks_at_a_time=50):
        ## First we init some basic things so that we will fill up in a sec.
        self.blocks_at_a_time = blocks_at_a_time
        self.bytes_read = 0
        self.rows = 0
        self.chromosome_names = []
        self.chromosome_lengths = []
        self.chromosomes_from_header = []
        self.file_handle = file_handle

        # Now we init the generator that, when iterated, will return some chunk of uncompress data.
        self.generator = self._generator()

        # We grab the first chunk, which is more than enough to contain all the header information, and parse the binary header:
        first_chunk = next(self.generator)
        if first_chunk[:4] != 'BAM\1': print 'ERROR: This doesnt look like a bam file to me :/'; exit()
        length_of_header = struct.unpack('<i',first_chunk[4:8])[0]
        self.header_text = struct.unpack(str(length_of_header)+'s',first_chunk[8:8+length_of_header])[0]
        self.number_of_reference_sequences = struct.unpack('<i',first_chunk[8+length_of_header:12+length_of_header])[0]
        rs = length_of_header + 12
        for _ in range(self.number_of_reference_sequences):
            l_name = struct.unpack('<l',first_chunk[rs:rs+4])[0]
            name, l_ref = struct.unpack('<%dsl' % l_name, first_chunk[rs+4:rs+l_name+8]) # disgusting.
            self.chromosome_names.append(name[:-1]) # We dont need the NUL byte.
            self.chromosome_lengths.append(l_ref)
            rs += 8 + l_name
        self.original_binary_header = buffer(first_chunk[:rs])
        self.left_overs = first_chunk[rs:] # bgzip doesnt care about where it drops the division between blocks, so we end up with some
                                           # bytes in the first block that contain read data, not header data.
        for line in self.header_text.split('\n'):
            if line.startswith('@SQ\tSN:'):
                self.chromosomes_from_header.append(line.split('\t')[1][3:]) # God I hope no programs put spaces in the headers instead of tabs..
        if self.chromosomes_from_header != self.chromosome_names:
            print 'ERROR: For some reason, and I have no idea why, the BAM file stores the chromosome name in two locations, the '
            print '       ASCII text header we all know and love, viewable with samtools view -H, and another special binary header'
            print '       which is used to translate the chromosome refID (a number) into a chromosome RNAME when you do bam -> sam.'
            print '       These two headers should always be the same, but apparently someone has messed that up. #YOLO\n'
            print 'Your ASCII header looks like:\n' + self.header_text
            print '\nWhile your binary header has the following chromosomes:'
            print self.chromosome_names
            exit()

        # Now we have the header data stored, but our first call to the generator gave us a block of uncompressed data back that contained way more than
        # just the header data alone. We want this class to be a generator that, on every request, hands back READ data in BAM format uncompressed, starting
        # from the first read, so in order to do that, we need to create a new generator-shim that just regurgitates what was left over from the header parse
        # the first time its called, then uses whatever self.generator would give on all subsequent tries after that.
        # This is how we do that:
        self._iterator = itertools.chain([self.left_overs],self.generator)

    # And this tells python the class is a generator:
    def __iter__(self): return self 
    def next(self): return next(self._iterator)

    def _generator(self):
        DEVNULL = open(os.devnull, 'wb')
        try:
            if type(self.file_handle) == str: p = subprocess.Popen(['pigz','-dc',self.file_handle], stdout=subprocess.PIPE, stderr=DEVNULL)
            elif type(self.file_handle) == file: p = subprocess.Popen(['pigz','-dc'],stdin=self.file_handle, stdout=subprocess.PIPE, stderr=DEVNULL)
            else: print 'ERROR: I do not know how to open and read from "' + str(self.file_handle) + '"'; exit()
            self.file_handle = p.stdout
            #sys.stderr.write('Using pigz!\n')
        except OSError:
            try:
                if type(self.file_handle) == str:    p = subprocess.Popen(['gzip','--stdout','--decompress','--force',self.file_handle]       , stdout=subprocess.PIPE, stderr=DEVNULL)
                elif type(self.file_handle) == file: p = subprocess.Popen(['gzip','--stdout','--decompress','--force'], stdin=self.file_handle, stdout=subprocess.PIPE, stderr=DEVNULL)
                else: print 'ERROR: I do not know how to open and read from "' + str(self.file_handle) + '"'; exit()
                self.file_handle = p.stdout
                #sys.stderr.write('Using gzip!\n')
            except OSError:
                sys.stderr.write('Using internal Python...\n') # We will end up using the python code below. It is faster than the gzip module, but
                                                               # due to how slow Python's zlib module is, it will end up being about 2x slower than pysam.
        data = self.file_handle.read(655360)
        self.bytes_read += 655360
        cache = []
        blocks_left_to_grab = self.blocks_at_a_time
        bs = 0
        checkpoint = 0
        #pool = Pool(processes=3) 
        #def decompress(data): return zlib.decompress(data, 47) # 47 == zlib.MAX_WBITS|32
        decompress = zlib.decompress
        while data:
            if len(data) - bs < 65536:
                data = data[bs:] + self.file_handle.read(35536)
                self.bytes_read += len(data) - bs
                bs = 0

            magic = data[bs:bs+4]
            if not magic: break # a child's heart
            if magic != "\x1f\x8b\x08\x04":
                if magic == 'BAM\1':
                    # The user has passed us already unzipped data, or we're reading from pigz/gzip :)
                    while data:
                        yield data
                        data = self.file_handle.read(35536)
                        self.bytes_read += len(data)
                    raise StopIteration
                elif magic == 'SQLi': print 'OOPS: You have used an SQLite database as your input BAM file!!'; exit()
                else:                 print 'ERROR: The input file is not in a format I understand :('       ; exit()

            try:
                # The gzip format allows compression containers to store metadata about whats inside them. bgzip uses this
                # to encode the virtual file pointers, final decompressed block size, crc checks etc - however we really dont
                # care -- we just want to unzip everything as quickly as possible. So instead of 'following the rules', and parsing this metadata safely,
                # we try to take a short-cut and jump right to the good stuff, and only if that fails we go the long-way-around and parse every bit of metadata properly:

                #cache.append(decompress(data[bs+18:more_bs-8])) ## originally i stored uncompressed data in a list and used map() to decompress in multiple threads. Was not faster.
                #new_zipped_data = data[bs:more_bs]              ## originally i decompressed the data with a standard zlib.decompress(data,47), headers and all. Was slower.
                more_bs = bs + struct.unpack("<H", data[bs+16:bs+18])[0] +1
                cache.append(decompress(data[bs+18:more_bs-8],-15))
                bs = more_bs
            except: ## zlib doesnt have a nice exception for when things go wrong. just "error"
                sys.stderr.write('INFO: Odd bzgip block detected! The author of pybam didnt think this would ever happen... please could you let me know?')
                header_data = magic + data[bs+4:bs+12]
                header_size = 12
                extra_len = struct.unpack("<H", header_data[-2:])[0]
                while header_size-12 < extra_len:
                    header_data += data[bs+12:bs+16]
                    subfield_id = header_data[-4:-2]
                    subfield_len = struct.unpack("<H", header_data[-2:])[0]
                    subfield_data = data[bs+16:bs+16+subfield_len]
                    header_data += subfield_data
                    header_size += subfield_len + 4
                    if subfield_id == 'BC': block_size = struct.unpack("<H", subfield_data)[0]
                raw_data = data[bs+16+subfield_len:bs+16+subfield_len+block_size-extra_len-19]
                crc_data = data[bs+16+subfield_len+block_size-extra_len-19:bs+16+subfield_len+block_size-extra_len-19+8] # I have left the numbers in verbose, because the above try is the optimised code.
                bs = bs+16+subfield_len+block_size-extra_len-19+8
                zipped_data = header_data + raw_data + crc_data
                cache.append(decompress(zipped_data,47)) # 31 works the same as 47.
                # Although the following in the bgzip code from biopython, its not needed if you let zlib decompress the whole zipped_data, header and crc, because it checks anyway (in C land)
                # I've left the manual crc checks in for documentation purposes:
                expected_crc = crc_data[:4]
                expected_size = struct.unpack("<I", crc_data[4:])[0]
                if len(unzipped_data) != expected_size: print 'ERROR: Failed to unpack due to a Type 1 CRC error. Could the BAM be corrupted?'; exit()
                crc = zlib.crc32(unzipped_data)
                if crc < 0: crc = struct.pack("<i", crc)
                else:       crc = struct.pack("<I", crc)
                if expected_crc != crc: print 'ERROR: Failed to unpack due to a Type 2 CRC error. Could the BAM be corrupted?'; exit()

            blocks_left_to_grab -= 1
            if blocks_left_to_grab == 0:
                yield ''.join(cache)
                cache = []
                blocks_left_to_grab = self.blocks_at_a_time
        self.file_handle.close()
        DEVNULL.close()
        if cache != '': yield ''.join(cache)



####################
## compile_parser ##
#########################################################################################################################################################
##                                                                                                                                                     ##
## The second part of pybam is the compile_parser function, which returns a new runtime-generated function, which is what will actually parse you BAM. ##
## We generate this code on-the-fly because doing so is much more efficient that having a whole bunch of "if x: else y" all over the place.            ##
## Actually, the whole point of PyPy is to do exactly this sort of optimisation for you - it figures out that certain paths in your code never happen, ##
## and optimises the code accordingly. However, pybam obviously needs to work on regular python too, so thats why i've done it "manually" here below.  ##
## Some people take offence to the existence of exec(). These people also really really like PEP8, apart from maybe the second paragraph. meh.         ##
##                                                                                                                                                     ##
## You use it like this:                                                                                                                               ##
##                                                                                                                                                     ##
## >>>> pure_bam_data = pybam.bgunzip('./ENCFF001LCU.bam.gz')                                                                                          ##
## >>>>                                                                                                                                                ##
## >>>> parser = pybam.compile_parser(['pos','mapq','qname'])                                                                                          ##
## >>>>                                                                                                                                                ##
## >>>> for read in parser(pure_bam_data):                                                                                                             ##
## ....     print read                                                                                                                                 ##
## ....     break                                                                                                                                      ##
## ....                                                                                                                                                ##
## (3000742, 0, 'SOLEXA1_0001:4:49:11382:21230#0')                                                                                                     ##
##                                                                                                                                                     ##
## So we tell compile_parser how to make our function (that we called parser, but we could have called it anything) by passing it a list of special    ##
## strings. We can see what the function compile_parser came up with by looking at "pybam.code" after compile_paser has run:                           ##
##                                                                                                                                                     ##
## >>>> print pybam.code                                                                                                                               ##
##                                                                                                                                                     ##
## def parser(data_generator):                                                                                                                         ##
##     chunk = next(data_generator)                                                                                                                    ##
##     CtoPy = { 'A':'<c', 'c':'<b', 'C':'<B', 's':'<h', 'S':'<H', 'i':'<i', 'I':'<I' }                                                                ##
##     py4py = { 'A':  1,  'c':  1,  'C':  1,  's':  2,  'S':  2,  'i':  4 , 'I':  4  }                                                                ##
##     dna = '=ACMGRSVTWYHKDBN'                                                                                                                        ##
##     cigar_codes = 'MIDNSHP=X'                                                                                                                       ##
##     from array import array                                                                                                                         ##
##     from struct import unpack                                                                                                                       ##
##     p = 0                                                                                                                                           ##
##     while True:                                                                                                                                     ##
##         try:                                                                                                                                        ##
##             while len(chunk) < p+36: chunk = chunk[p:] + next(data_generator); p = 0                                                                ##
##             block_size,refID,pos,l_read_name,mapq                                                     = unpack('<iiiBB'       ,chunk[p:p+14])       ##
##             while len(chunk) < p + 4 + block_size: chunk = chunk[p:] + next(data_generator); p = 0                                                  ##
##         except StopIteration: break                                                                                                                 ##
##         end = p + block_size + 4                                                                                                                    ##
##         p += 36                                                                                                                                     ##
##         qname = chunk[p:p+l_read_name-1]                                                                                                            ##
##         p = end                                                                                                                                     ##
##         yield pos,mapq,qname                                                                                                                        ##
##                                                                                                                                                     ##
##                                                                                                                                                     ##
## A'int that complicated eh. So now you're probably wondering what special strings exist. Well, for now thats only a few. Many more could be added to ##
## do more abstract things, but right now we've got:                                                                                                   ##
##                                                                                                                                                     ##
#########################################################################################################################################################

def compile_parser(fields):
    deps = set(fields)
    unpack = '''
def parser(data_generator):
    chunk = next(data_generator)
    CtoPy = { 'A':'<c', 'c':'<b', 'C':'<B', 's':'<h', 'S':'<H', 'i':'<i', 'I':'<I' }
    py4py = { 'A':  1,  'c':  1,  'C':  1,  's':  2,  'S':  2,  'i':  4 , 'I':  4  }
    dna = '=ACMGRSVTWYHKDBN'
    cigar_codes = 'MIDNSHP=X'
    from array import array
    from struct import unpack
    p = 0
    while True:
        try:
            while len(chunk) < p+36: chunk = chunk[p:] + next(data_generator); p = 0

            '''

    # Some values require the length of the previous values to be known before they can be grabbed. Add those to our 'dependency' list:
    if  'tags' in deps or 'tags_bam'  in deps: deps.update(['qual_skip','seq_skip','l_seq_bytes','l_seq','cigar_skip','n_cigar_op','qname_skip','l_read_name','fixed_skip'])
    if  'qual' in deps or 'qual_bam'  in deps: deps.update([            'seq_skip','l_seq_bytes','l_seq','cigar_skip','n_cigar_op','qname_skip','l_read_name','fixed_skip'])
    if   'seq' in deps or 'seq_bam'   in deps: deps.update([                       'l_seq_bytes','l_seq','cigar_skip','n_cigar_op','qname_skip','l_read_name','fixed_skip'])
    if 'cigar' in deps or 'cigar_bam' in deps: deps.update([                                                          'n_cigar_op','qname_skip','l_read_name','fixed_skip'])
    if 'qname' in deps or 'qname_bam' in deps: deps.update([                                                                                    'l_read_name','fixed_skip'])

    # Fixed length SAM data. Because its often quicker to unpack more data in a single stuct.unpack call than it is to unpack less data with 2 or more calls to unpack, the logic here
    # is very simple - we unpack everything up until the value desired. It needs to be tested if this is actually holds true in real-world data though.
    if   'tlen'        in deps: unpack += "block_size,refID,pos,l_read_name,mapq,bin_,n_cigar_op,flag,l_seq,next_refID,next_pos,tlen = unpack('<iiiBBHHHiiii',chunk[p:p+36])"    ; fixed_length = 36
    elif 'next_pos'    in deps: unpack += "block_size,refID,pos,l_read_name,mapq,bin_,n_cigar_op,flag,l_seq,next_refID,next_pos      = unpack('<iiiBBHHHiii' ,chunk[p:p+32])"    ; fixed_length = 32
    elif 'next_refID'  in deps: unpack += "block_size,refID,pos,l_read_name,mapq,bin_,n_cigar_op,flag,l_seq,next_refID               = unpack('<iiiBBHHHii'  ,chunk[p:p+28])"    ; fixed_length = 28
    elif 'l_seq'       in deps: unpack += "block_size,refID,pos,l_read_name,mapq,bin_,n_cigar_op,flag,l_seq                          = unpack('<iiiBBHHHi'   ,chunk[p:p+24])"    ; fixed_length = 24
    elif 'flag'        in deps: unpack += "block_size,refID,pos,l_read_name,mapq,bin_,n_cigar_op,flag                                = unpack('<iiiBBHHH'    ,chunk[p:p+20])"    ; fixed_length = 20
    elif 'n_cigar_op'  in deps: unpack += "block_size,refID,pos,l_read_name,mapq,bin_,n_cigar_op                                     = unpack('<iiiBBHH'     ,chunk[p:p+18])"    ; fixed_length = 18
    elif 'bin'         in deps: unpack += "block_size,refID,pos,l_read_name,mapq,bin_                                                = unpack('<iiiBBH'      ,chunk[p:p+16])"    ; fixed_length = 16
    elif 'mapq'        in deps: unpack += "block_size,refID,pos,l_read_name,mapq                                                     = unpack('<iiiBB'       ,chunk[p:p+14])"    ; fixed_length = 14
    elif 'l_read_name' in deps: unpack += "block_size,refID,pos,l_read_name                                                          = unpack('<iiiB'        ,chunk[p:p+13])"    ; fixed_length = 13
    elif 'pos'         in deps: unpack += "block_size,refID,pos                                                                      = unpack('<iii'         ,chunk[p:p+12])"    ; fixed_length = 12
    elif 'refID'       in deps: unpack += "block_size,refID                                                                          = unpack('<ii'          ,chunk[p:p+8 ])"    ; fixed_length = 8
    else:                       unpack += "block_size                                                                                = unpack('<i'           ,chunk[p:p+4 ])[0]" ; fixed_length = 4

    # Fixed-length BAM data (where we just grab the bytes, we dont unpack) can, however, be grabbed individually.
    if 'fixed_bam'       in deps: unpack += "\n            fixed_bam       = chunk[p:p+36]"    # All the fixed data.
    if 'block_size_bam'  in deps: unpack += "\n            block_size_bam  = chunk[p:p+4]"
    if 'refID_bam'       in deps: unpack += "\n            refID_bam       = chunk[p+4:p+8]"
    if 'pos_bam'         in deps: unpack += "\n            pos_bam         = chunk[p+8:p+12]"
    if 'l_read_name_bam' in deps: unpack += "\n            l_read_name_bam = chunk[p+12:p+13]"
    if 'mapq_bam'        in deps: unpack += "\n            mapq_bam        = chunk[p+13:p+14]"
    if 'bin_bam'         in deps: unpack += "\n            bin_bam         = chunk[p+14:p+16]"
    if 'n_cigar_op_bam'  in deps: unpack += "\n            n_cigar_op_bam  = chunk[p+16:p+18]"
    if 'flag_bam'        in deps: unpack += "\n            flag_bam        = chunk[p+18:p+20]"
    if 'l_seq_bam'       in deps: unpack += "\n            l_seq_bam       = chunk[p+20:p+24]"
    if 'next_refID_bam'  in deps: unpack += "\n            next_refID_bam  = chunk[p+24:p+28]"
    if 'next_pos_bam'    in deps: unpack += "\n            next_pos_bam    = chunk[p+28:p+32]"
    if 'tlen_bam'        in deps: unpack += "\n            tlen_bam        = chunk[p+32:p+36]"

    unpack += '''
            while len(chunk) < p + 4 + block_size: chunk = chunk[p:] + next(data_generator); p = 0
        except StopIteration: break
        end = p + block_size + 4
'''

    if 'variable_bam'in deps: unpack += "        variable_bam = chunk[p:end]\n"
    if 'fixed_skip'  in deps: unpack += "        p += 36\n"
    if 'qname_bam'   in deps: unpack += "        qname_bam = chunk[p:p+l_read_name]\n"
    if 'qname'       in deps: unpack += "        qname = chunk[p:p+l_read_name-1]\n"
    if 'qname_skip'  in deps: unpack += "        p += l_read_name\n"
    if 'cigar'       in deps: unpack += "        cigar = [(cigar_codes[cig & 0b1111],cig>>4) for cig in array('I',chunk[p:p+(4*n_cigar_op)])]\n"
    if 'cigar_bam'   in deps: unpack += "        cigar_bam = chunk[p:p+(4*n_cigar_op)]\n"
    if 'cigar_skip'  in deps: unpack += "        p += n_cigar_op*4\n"
    if 'l_seq_bytes' in deps: unpack += "        l_seq_bytes = -((-l_seq)//2)\n"
    if 'seq'         in deps: unpack += "        seq = ''.join([ dna[bit4 >> 4] + dna[bit4 & 0b1111] for bit4 in array('B', chunk[p:p+l_seq_bytes]) ])\n"
    if 'seq_bam'     in deps: unpack += "        seq_bam = chunk[p:p+l_seq_bytes]\n"
    if 'seq_skip'    in deps: unpack += "        p += l_seq_bytes\n"
    if 'qual'        in deps: unpack += "        qual = ''.join([ chr(ord(q)+33) for q in chunk[p:p+l_seq] ])\n"
    if 'qual_bam'    in deps: unpack += "        qual_bam = chunk[p:p+l_seq]\n"
    if 'qual_skip'   in deps: unpack += "        p += l_seq\n"
    if 'tags_bam'    in deps: unpack += "        tags_bam = chunk[p:end]\n"
    if 'tags'        in deps: unpack += '''
        tags = []
        while p != end:
            tag_name = chunk[p:p+2]
            tag_type = chunk[p+2]
            if tag_type == 'Z':
                p_end = chunk.index('\\0',p+3)+1
                tag_data = chunk[p+3:p_end];  # I've opted to keep the NUL byte in the output, although it makes me feel dirty.
            elif tag_type in CtoPy:
                p_end = p+3+py4py[tag_type]
                tag_data = unpack(CtoPy[tag_type],chunk[p+3:p_end])[0]
            else:
                print 'PROGRAMMER ERROR: I dont know how to parse BAM tags in this format: ',repr(tag_type)
                print '                  This is simply because I never saw this kind of tag during development.'
                print '                  If you could mail the following chunk of text to john at john.uk.com, ill fix this up :)'
                print repr(tag_type),repr(chunk[p+3:end])
                exit()
            tags.append((tag_name,tag_type,tag_data))
            p = p_end\n'''
    unpack += '        p = end\n'
    unpack += '        yield ' + ','.join([x for x in fields])
    global code    # To allow user to view
    code = unpack  # code with pybam.code

    exec(unpack)  # "parser" now comes into existance
    return parser

#######################
## THE END / CREDITS ##
#########################################################################################################################################################
##                                                                                                                                                     ##
## Thats it!                                                                                                                                           ##
## This code was written by John Longinotto, a PhD student of the Pospisilik Lab at the Max Planck Institute of Immunbiology & Epigenetics, Freiburg.  ##
## My PhD is funded by the Deutsches Epigenom Programm (DEEP), and the Max Planck IMPRS Program.                                                       ##
## I study Adipose Biology and Circadian Rhythm in mice, although it seems these days I spend most of my time at the computer and not at the bench.    ##
##                                                                                                                                                     ##
#########################################################################################################################################################