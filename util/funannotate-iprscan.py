#!/usr/bin/env python

import sys, os, argparse, inspect, multiprocessing, urllib2, subprocess, time
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)

parser=argparse.ArgumentParser(prog='funannotate-iprscan.py',
    description='''Script to run InterProscan locally or with Docker''',
    epilog="""Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', required=True, help='FASTA file')
parser.add_argument('-n','--num', type=int, default=1000, help='Number of files in each chunk')
parser.add_argument('-c','--cpus', default=12, type=int, help='Total number of CPUs')
parser.add_argument('-m','--method', required=True, choices=['local', 'docker'], help='Number of "chunks" or files')
parser.add_argument('--cpus_per_chunk', default=4, type=int, help='Number of CPUs per IPR instance')
parser.add_argument('--iprscan_path', help='Local Path to interproscan.sh')
parser.add_argument('-o','--out', required=True, help='Final output XML file')
args=parser.parse_args()

def which(name):
    try:
        with open(os.devnull) as devnull:
            diff = ['tbl2asn', 'dustmasker', 'mafft']
            if not any(name in x for x in diff):
                subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
            else:
                subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def download(url, name):
    file_name = name
    u = urllib2.urlopen(url)
    f = open(file_name, 'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print("Downloading: {0} Bytes: {1}".format(url, file_size))
    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break
        file_size_dl += len(buffer)
        f.write(buffer)
        p = float(file_size_dl) / file_size
        status = r"{0}  [{1:.2%}]".format(file_size_dl, p)
        status = status + chr(8)*(len(status)+1)
        sys.stdout.write(status)
    f.close()
    
def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with open(path) as inf:     
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos

def return_lines(path, linepos, nstart, nstop):
    """return nsamp lines from path where line offsets are in linepos"""
    offsets = linepos[nstart:nstop]
    lines = []
    with open(path) as inf:
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

def downloadIPRproperties(name, cpus):
    download('https://raw.githubusercontent.com/ebi-pf-team/interproscan/5.22-61.0/core/jms-implementation/support-mini-x86-32/interproscan.properties', 'ipr.tmp')
    with open(name, 'w') as outfile:
        with open('ipr.tmp', 'rU') as infile:
            for line in infile:
                if line.startswith('number.of.embedded.workers='):
                    outfile.write('number.of.embedded.workers=1\n')
                elif line.startswith('maxnumber.of.embedded.workers='):
                    outfile.write('maxnumber.of.embedded.workers=%i\n' % cpus)
                else:
                    outfile.write(line)
    os.remove('ipr.tmp')
    
def runDocker(input):
    cmd = ['docker', 'run', '-u', '$UID:$GROUPS', '--rm', '-v', '$PWD:/dir', '-v', '$PWD:/in', '-v', '$PWD/interproscan.properties:/interproscan-5.22-61.0/interproscan.properties', 'blaxterlab/interproscan:latest', 'interproscan.sh', '-i', '/in/'+input, '-d', '/dir', '-dp', '-f', 'XML', '-goterms', '-pa']
    logfile = os.path.join(tmpdir, input.split('.fasta')[0]+'.log')
    with open(logfile, 'w') as log:
        subprocess.call(cmd, cwd=tmpdir, stdout=log, stderr=log)

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: runDocker(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))

def runMultiProgress(function, inputList, cpus):
    #setup pool
    p = multiprocessing.Pool(cpus)
    #setup results and split over cpus
    tasks = len(inputList)
    results = []
    for i in inputList:
        results.append(p.apply_async(function, [i]))
    #refresh pbar every 5 seconds
    while True:
        incomplete_count = sum(1 for x in results if not x.ready())
        if incomplete_count == 0:
            break
        sys.stdout.write("Progress: %.2f%% \r" % (float(tasks - incomplete_count) / tasks * 100))
        sys.stdout.flush()
        time.sleep(1)
    p.close()
    p.join()

#make temp directory
tmpdir = 'iprscan_' + str(os.getpid())
os.makedirs(tmpdir)

#figure out number of chunks
count = countfasta(args.input)
chunks = count / args.num
threads = args.cpus / args.cpus_per_chunk

#split protein fasta files into chunks
split_fasta(args.input, tmpdir, chunks)

#get list of inputs
file_list = []
for file in os.listdir(tmpdir):
    if file.endswith('.fasta'):
        file_list.append(file)

if args.method == 'docker':
    if not which('docker'):
        print('Docker is not running or not installed, exiting')
        sys.exit(1)
    ipr_properties = os.path.join(tmpdir, 'interproscan.properties')
    downloadIPRproperties(ipr_properties, args.cpus_per_chunk)
    runMultiProgress(safe_run, file_list, threads)
    
final_list = []
for file in os.listdir(tmpdir):
    if file.endswith('.xml'):
        final_list.append(os.path.join(tmpdir, file))
with open(args.out, 'w') as output:
    output.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
    output.write('<protein-matches xmlns="http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5">\n')
    for x in final_list:
        with open(x, 'rU') as infile:
            lines = infile.readlines()
            for line in lines[2:-1]:
                output.write(line)
    output.write('</protein-matches>\n')