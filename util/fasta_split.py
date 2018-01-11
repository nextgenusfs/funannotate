#!/usr/bin/env python

import sys, os, argparse

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='fasta_split.py',
    description='''Script splits fasta file into chunos''',
    epilog="""Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='FASTA file')
parser.add_argument('-n','--num', required=True, type=int, help='Number of "chunks" or files')
parser.add_argument('-o','--out', default='.', help='Output folder')
args=parser.parse_args()

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
      if magic_number == '\x1f\x8b':
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
   def next(self):
      return next(self.f)


def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with gzopen(path) as inf:     
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos

def return_lines(path, linepos, nstart, nstop):
    """return nsamp lines from path where line offsets are in linepos"""
    offsets = linepos[nstart:nstop]
    lines = []
    with gzopen(path) as inf:
        for offset in offsets:
            inf.seek(offset)
            lines.append(inf.readline())
    return lines     

def split_fasta(input, outputdir, chunks):
    #function to return line positions of fasta files for chunking
    fastapos = []
    position = 0
    numseqs = 0
    basename = os.path.basename(input).split('.fa',-1)[0]
    with open(input, 'rU') as infile:
        for line in infile:
            if line.startswith('>'):
                numseqs += 1
                fastapos.append(position)
            position += 1
    splits = []
    n = numseqs / chunks
    num = 0
    for i in range(chunks):
        if i == 0:
            start = 0
            num = n
            lastpos = fastapos[n+1]
        else:
            start = lastpos
            num = num + n
            try:
                lastpos = fastapos[num+1] #find the n+1 seq
            except IndexError:
                lastpos = fastapos[-1]
        splits.append((start, lastpos))
    #check if output folder exists, if not create it
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    #get line positions from file
    linepos = scan_linepos(input)
    #loop through the positions and write output
    for i, x in enumerate(splits):
        num = i+1
        with open(os.path.join(outputdir, basename+'_'+str(num)+'.fasta'), 'w') as output:
            lines = return_lines(input, linepos, x[0], x[1])
            output.write('%s' % ''.join(lines))

split_fasta(args.input, args.out, args.num)
