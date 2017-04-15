from __future__ import division
import os, subprocess, logging, sys, argparse, inspect, csv, time, re, shutil, datetime, glob, platform, multiprocessing, itertools, hashlib, math
from natsort import natsorted
import warnings
from Bio import SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

#get the working directory, so you can move back into DB folder to find the files you need
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
DB = os.path.join(parentdir, 'DB')
LIB = os.path.join(parentdir, 'lib')
UTIL = os.path.join(parentdir, 'util')
GeneMark2GFF = os.path.join(UTIL, 'genemark_gtf2gff3.pl')

pref_colors=["#CF3C57","#65B23A","#6170DD","#D18738","#D542B5",
"#724A63","#60AABA","#5DB07C","#6C5824","#D74B2B","#6B97D6","#893B2E",
"#B68DB7","#564E91","#ACA13C","#3C6171","#436B33","#D84088",
"#D67A77","#9D55C4","#8B336E","#DA77B9","#D850E5","#B188DF"]

Nogs = {'NOG': 'All organisms (5.0GB)',
'aciNOG': 'Acidobacteria (125.3MB)',
'acidNOG': 'Acidobacteriia (75.4MB)',
'acoNOG': 'Aconoidasida (217.1MB)',
'actNOG': 'Actinobacteria (765.3MB)',
'agaNOG': 'Agaricales (211.1MB)',
'agarNOG': 'Agaricomycetes (236.5MB)',
'apiNOG': 'Apicomplexa (322.7MB)',
'aproNOG': 'Proteobacteria_alpha (638.4MB)',
'aquNOG': 'Aquificae (51.5MB)',
'arNOG': 'Archaea (256.9MB)',
'arcNOG': 'Archaeoglobi (21.8MB)',
'artNOG': 'Arthropoda (725.0MB)',
'arthNOG': 'Arthrodermataceae (111.2MB)',
'ascNOG': 'Ascomycota (1.1GB)',
'aveNOG': 'Aves (186.1MB)',
'bacNOG': 'Bacilli (362.6MB)',
'bactNOG': 'Bacteria (3.3GB)',
'bacteNOG': 'Bacteroidia (199.2MB)',
'basNOG': 'Basidiomycota (356.5MB)',
'bctoNOG': 'Bacteroidetes (508.9MB)',
'biNOG': 'Bilateria (1.7GB)',
'bproNOG': 'Proteobacteria_beta (481.0MB)',
'braNOG': 'Brassicales (275.4MB)',
'carNOG': 'Carnivora (293.5MB)',
'chaNOG': 'Chaetomiaceae (180.9MB)',
'chlNOG': 'Chlorobi (51.3MB)',
'chlaNOG': 'Chlamydiae (39.1MB)',
'chloNOG': 'Chloroflexi (136.8MB)',
'chlorNOG': 'Chloroflexi (75.8MB)',
'chloroNOG': 'Chlorophyta (146.8MB)',
'chorNOG': 'Chordata (1.1GB)',
'chrNOG': 'Chromadorea (392.6MB)',
'cloNOG': 'Clostridia (505.6MB)',
'cocNOG': 'Coccidia (137.4MB)',
'creNOG': 'Crenarchaeota (110.0MB)',
'cryNOG': 'Cryptosporidiidae (105.4MB)',
'cyaNOG': 'Cyanobacteria (254.8MB)',
'cytNOG': 'Cytophagia (164.6MB)',
'debNOG': 'Debaryomycetaceae (145.5MB)',
'defNOG': 'Deferribacteres (41.6MB)',
'dehNOG': 'Dehalococcoidetes (15.0MB)',
'deiNOG': 'Deinococcusthermus (75.4MB)',
'delNOG': 'delta/epsilon (471.4MB)',
'dipNOG': 'Diptera (397.7MB)',
'dotNOG': 'Dothideomycetes (298.2MB)',
'dproNOG': 'Proteobacteria_delta (424.6MB)',
'droNOG': 'Drosophilidae (314.1MB)',
'eproNOG': 'Proteobacteria_epsilon (104.8MB)',
'eryNOG': 'Erysipelotrichi (85.8MB)',
'euNOG': 'Eukaryotes (3.1GB)',
'eurNOG': 'Euryarchaeota (264.7MB)',
'euroNOG': 'Eurotiomycetes (507.2MB)',
'eurotNOG': 'Eurotiales (358.1MB)',
'fiNOG': 'Fishes (641.2MB)',
'firmNOG': 'Firmicutes (728.8MB)',
'flaNOG': 'Flavobacteriia (222.5MB)',
'fuNOG': 'Fungi (1.2GB)',
'fusoNOG': 'Fusobacteria (74.9MB)',
'gproNOG': 'Proteobacteria_gamma (735.0MB)',
'haeNOG': 'Haemosporida (197.1MB)',
'halNOG': 'Halobacteria (106.6MB)',
'homNOG': 'Hominidae (229.9MB)',
'hymNOG': 'Hymenoptera (199.5MB)',
'hypNOG': 'Hypocreales (353.3MB)',
'inNOG': 'Insects (688.9MB)',
'kinNOG': 'Kinetoplastida (259.8MB)',
'lepNOG': 'Lepidoptera (208.0MB)',
'lilNOG': 'Liliopsida (660.0MB)',
'maNOG': 'Mammals (855.5MB)',
'magNOG': 'Magnaporthales (161.3MB)',
'meNOG': 'Animals (1.8GB)',
'metNOG': 'Methanobacteria (38.4MB)',
'methNOG': 'Methanococci (24.5MB)',
'methaNOG': 'Methanomicrobia (99.4MB)',
'necNOG': 'Nectriaceae (200.3MB)',
'negNOG': 'Negativicutes (96.5MB)',
'nemNOG': 'Nematodes (430.0MB)',
'onyNOG': 'Onygenales (282.8MB)',
'opiNOG': 'Opisthokonts (2.8GB)',
'perNOG': 'Peronosporales (154.1MB)',
'plaNOG': 'Planctomycetes (149.3MB)',
'pleNOG': 'Pleosporales (223.4MB)',
'poaNOG': 'Poales (596.3MB)',
'prNOG': 'Primates (448.8MB)',
'proNOG': 'Proteobacteria (1.5GB)',
'rhaNOG': 'Rhabditida (334.5MB)',
'roNOG': 'Rodents (381.4MB)',
'sacNOG': 'Saccharomycetaceae (202.7MB)',
'saccNOG': 'Saccharomycetes (275.9MB)',
'sorNOG': 'Sordariales (296.1MB)',
'sordNOG': 'Sordariomycetes (714.1MB)',
'sphNOG': 'Sphingobacteriia (154.0MB)',
'spiNOG': 'Spirochaetes (121.2MB)',
'spriNOG': 'Supraprimates (635.6MB)',
'strNOG': 'Streptophyta (960.6MB)',
'synNOG': 'Synergistetes (59.5MB)',
'tenNOG': 'Tenericutes (29.9MB)',
'thaNOG': 'Thaumarchaeota (15.3MB)',
'theNOG': 'Thermoplasmata (26.9MB)',
'therNOG': 'Thermotogae (66.5MB)',
'thermNOG': 'Thermococci (31.4MB)',
'treNOG': 'Tremellales (79.9MB)',
'veNOG': 'Vertebrates (1.0GB)',
'verNOG': 'Verrucomicrobia (140.9MB)',
'verrNOG': 'Verrucomicrobiae (73.0MB)',
'virNOG': 'Viridiplantae (1.0GB)'}

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
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'
    
def SafeRemove(input):
    if os.path.isdir(input):
        shutil.rmtree(input)
    elif os.path.isfile(input):
        os.remove(input)
    else:
        return

def runSubprocess(cmd, dir, logfile):
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stdout:
        logfile.debug(stdout)
    if stderr:
        logfile.debug(stderr)

def runSubprocess2(cmd, dir, logfile, output):
    #function where output of cmd is STDOUT, capture STDERR in logfile
    logfile.debug(' '.join(cmd))
    with open(output, 'w') as out:
        proc = subprocess.Popen(cmd, cwd=dir, stdout=out, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        if stderr[0] != None:
            logfile.debug(stderr)

def runSubprocess3(cmd, dir, logfile):
    #function where STDOUT pipes to FNULL, capture STDERR in logfile
    FNULL = open(os.devnull, 'w')
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=FNULL, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if stderr:
        logfile.debug(stderr)
        
def runSubprocess4(cmd, dir, logfile):
    #function where STDOUT and STDERR pipes to FNULL
    FNULL = open(os.devnull, 'w')
    logfile.debug(' '.join(cmd))
    proc = subprocess.Popen(cmd, cwd=dir, stdout=FNULL, stderr=FNULL)


def hashfile(afile, hasher, blocksize=65536):
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.digest()
    
def sha256_check(file1, file2):
    files = [file1, file2]
    output = [(fname, hashfile(open(fname, 'rb'), hashlib.sha256())) for fname in files]
    if output[0][1] == output[1][1]:
        return True
    else:
        return False

def readBlocks(source, pattern):
    buffer = []
    for line in source:
        if line.startswith(pattern):
            if buffer: yield buffer
            buffer = [ line ]
        else:
            buffer.append( line )
    yield buffer

def empty_line_sep(line):
    return line=='\n'

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
    if size < 2: #this is 1 character...
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

def which(name):
    try:
        with open(os.devnull) as devnull:
            diff = ['tbl2asn', 'dustmasker', 'mafft', 'signalp', 'proteinortho5.pl', 'ete3', 'phyml']
            if not any(name in x for x in diff):
                subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
            else:
                if name == 'signalp':
                    subprocess.Popen([name, '-V'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'dustmasker':
                    subprocess.Popen([name, '-version-full'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'tbl2asn':
                    subprocess.Popen([name, '--help'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'raxmlHPC-PTHREADS':
                    subprocess.Popen([name, '-version'], stdout=devnull, stderr=devnull).communicate()
                elif name == 'ete3':
                    subprocess.Popen([name, 'version'], stdout=devnull, stderr=devnull).communicate()
                else:
                    subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def CheckDependencies(input):
    missing = []
    for p in input:
        if which(p) == False:
            missing.append(p)
    if missing != []:
        error = ", ".join(missing)
        log.error("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        sys.exit(1)
        
def checkannotations(input):
    if os.path.isfile(input):
        filesize = getSize(input)
        if int(filesize) < 1:
            return False
        else:
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
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count
    
def get_version():
    version = subprocess.Popen(['funannotate', 'version'], stdout=subprocess.PIPE).communicate()[0].rstrip()
    return version

def checkAugustusFunc(base):
    brakerpass = 0
    buscopass = 0
    version = subprocess.Popen(['augustus', '--version'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0].rstrip()
    version = version.split(' is ')[0]
    bam2hints = which(os.path.join(base, 'bin', 'bam2hints'))
    filterBam = which(os.path.join(base, 'bin', 'filterBam'))
    if bam2hints and filterBam:
        brakerpass = 1
    model = os.path.join(parentdir, 'lib', 'EOG092C0B3U.prfl')
    if not os.path.isfile(model):
        log.error("Testing Augustus Error: installation seems wrong, can't file prfl model")
        sys.exit(1)
    profile = '--proteinprofile='+model
    proteinprofile = subprocess.Popen(['augustus', '--species=anidulans', profile, os.path.join(parentdir, 'lib', 'busco_test.fa')], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0].rstrip()
    if not 'augustus: ERROR' in proteinprofile:
        buscopass = 1
    return (version, brakerpass, buscopass)
    
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
    for i in range(0,cols):
        length = max(map(lambda x: len(x), mylist[i::cols]))
        length += 2
        ljust = map(lambda x: x.ljust(length), mylist[i::cols])
        justify.append(ljust)
    justify = flatten(justify)
    num_lines = len(mylist) / cols
    lines = (' '.join(justify[i::num_lines]) 
             for i in range(0,num_lines))
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
    if cols > len(sobj): cols = len(sobj)
    max_len = max([len(item) for item in sobj])
    if columnwise: cols = int(math.ceil(float(len(sobj)) / float(cols)))
    plist = [sobj[i: i+cols] for i in range(0, len(sobj), cols)]
    if columnwise:
        if not len(plist[-1]) == cols:
            plist[-1].extend(['']*(len(sobj) - len(plist[-1])))
        plist = zip(*plist)
    printer = '\n'.join([
        ''.join([c.ljust(max_len + gap) for c in p])
        for p in plist])
    return printer


def roundup(x):
    return x if x % 100 == 0 else x + 100 - x % 100
    
def maxabs(a, axis=None):
    """Return slice of a, keeping only those values that are furthest away
    from 0 along axis"""
    maxa = a.max(axis=axis)
    mina = a.min(axis=axis)
    p = abs(maxa) > abs(mina) # bool, or indices where +ve values win
    n = abs(mina) > abs(maxa) # bool, or indices where -ve values win
    if axis == None:
        if p: return maxa
        else: return mina
    shape = list(a.shape)
    shape.pop(axis)
    out = np.zeros(shape, dtype=a.dtype)
    out[p] = maxa[p]
    out[n] = mina[n]
    return out
    
def setupLogging(LOGNAME):
    global log
    if 'win32' in sys.platform:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='[%I:%M:%S %p]')
    else:
        stdoutformat = logging.Formatter(colr.GRN+'%(asctime)s'+colr.END+': %(message)s', datefmt='[%I:%M:%S %p]')
    fileformat = logging.Formatter('%(asctime)s: %(message)s')
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

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def countGFFgenes(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if "\tgene\t" in line:
                count += 1
    return count

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
        sys.stdout.write("     Progress: %.2f%% \r" % (float(tasks - incomplete_count) / tasks * 100))
        sys.stdout.flush()
        time.sleep(1)
    p.close()
    p.join()

def cleanProteins(inputList, output):
    #expecting a list of protein fasta files for combining/cleaning headers
    #make sure you aren't duplicated sequences names
    #dropping proteins less than 50 amino acids
    seen = set()
    with open(output, 'w') as out:
        for x in inputList:
            with open(x, 'rU') as input:
                for rec in SeqIO.parse(input, 'fasta'):
                    if len(rec.seq) < 50:
                        continue
                    #explicitly check for swissprot and jgi
                    if rec.id.startswith('sp|') or rec.id.startswith('jgi|'):
                        ID = rec.id.split('|')[-1]
                    else:
                        ID = rec.id
                    #now clean up the shit
                    badshit = [':', ';', '/', '\\', '.', ',', '%']
                    for i in badshit:
                        if i in ID:
                            ID = ID.replace(i, '_')
                    if not ID in seen:
                        seen.add(ID)
                    else:
                        ID = ID+'_1'
                        if not ID in seen:
                            seen.add(ID)
                        else:
                            num = int(ID.split('_')[1])
                            ID = ID.split('_')[0]+str(num+1)
                    out.write('>%s\n%s\n' % (ID, rec.seq))
  
def gb2output(input, output1, output2, output3):
    with open(output1, 'w') as proteins:
        with open(output2, 'w') as transcripts:
            with open(output3, 'w') as scaffolds:
                with open(input, 'rU') as gbk:
                    SeqRecords = SeqIO.parse(gbk, 'genbank')
                    for record in SeqRecords:
                        scaffolds.write(">%s\n%s\n" % (record.id, record.seq))
                        for f in record.features:
                            if f.type == "CDS":
                                proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0].rstrip('*')))
                            if f.type == "mRNA":
                                feature_seq = f.extract(record.seq)
                                transcripts.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))
                                
def checkGenBank(input):
    count = 0
    with open(input, 'rU') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                if f.type == 'CDS':
                    count += 1
    if count == 0:
        return False
    else:
        return True
        
def checkFastaHeaders(input, limit):
    length = 0
    with open(input, 'rU') as fasta:
        for line in fasta:
            if line.startswith('>'):
                line = line.replace('\n', '')
                headlen = len(line) - 1 #subtract one character for fasta carrot 
                if headlen > length:
                    length = headlen
    if length > int(limit):
        return False
    else:
        return True

def BamHeaderTest(genome, mapping):
    import pybam
    #get list of fasta headers from genome
    genome_headers = []
    with open(genome, 'rU') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            if rec.id not in genome_headers:
                genome_headers.append(rec.id)
    #get list of fasta headers from BAM
    bam_headers = []
    with open(mapping, 'rU') as bamin:
        bam_file = pybam.bgunzip(bamin)
        bam_headers = bam_file.chromosomes_from_header
    #now compare lists, basically if BAM headers not in genome headers, then output bad names to logfile and return FALSE
    genome_headers = set(genome_headers)
    diffs = [x for x in bam_headers if x not in genome_headers]
    if len(diffs) > 0:
        log.debug("ERROR: These BAM headers not found in genome FASTA headers\n%s" % ','.join(diffs))
        return False
    else:
        return True
    
def gb2allout(input, GFF, Proteins, Transcripts, DNA):
    #this will not output any UTRs for gene models, don't think this is a problem right now....
    with open(GFF, 'w') as gff:
        gff.write("##gff-version 3\n")
        with open(Proteins, 'w') as proteins:
            with open(Transcripts, 'w') as transcripts:
                with open(DNA, 'w') as scaffolds:
                    with open(input, 'rU') as gbk:
                        for record in SeqIO.parse(gbk, 'genbank'):
                            scaffolds.write(">%s\n%s\n" % (record.id, record.seq))
                            for f in record.features:
                                if f.type == "mRNA":
                                    feature_seq = f.extract(record.seq)
                                    transcripts.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))
                                if f.type == 'CDS':
                                    proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0].rstrip('*')))
                                    chr = record.id
                                    ID = f.qualifiers['locus_tag'][0]
                                    try:
                                        product = f.qualifiers['product'][0]
                                    except KeyError:
                                        product = "hypothetical protein"
                                    start = f.location.nofuzzy_start + 1
                                    end = f.location.nofuzzy_end
                                    strand = f.location.strand
                                    if strand == 1:
                                        strand = '+'
                                    elif strand == -1:
                                        strand = '-'
                                    num_exons = len(f.location.parts)
                                    current_phase = 0
                                    gff.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n" % (chr, start, end, strand, ID))
                                    gff.write("%s\tGenBank\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, ID, product))
                                    if num_exons < 2: #only a single exon
                                        ex_start = str(f.location.nofuzzy_start + 1)
                                        ex_end = str(f.location.nofuzzy_end)
                                        gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon1;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ID))
                                        gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t0\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ID))
                                    else: #more than 1 exon, so parts are actually in correct orientation, so loop through
                                        for i in range(0,num_exons):
                                            ex_start = str(f.location.parts[i].nofuzzy_start + 1)
                                            ex_end = str(f.location.parts[i].nofuzzy_end)
                                            ex_num = i + 1
                                            gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon%i;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ex_num, ID))
                                            gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, current_phase, ID, ID))
                                            current_phase = (current_phase - (int(ex_end) - int(ex_start) + 1)) % 3
                                            if current_phase == 3:
                                                current_phase = 0

                                if f.type == 'tRNA':
                                    ID = f.qualifiers['locus_tag'][0]
                                    start = f.location.nofuzzy_start
                                    end = f.location.nofuzzy_end
                                    strand = f.location.strand
                                    if strand == 1:
                                        strand = '+'
                                    elif strand == -1:
                                        strand = '-'
                                    try:
                                        product = f.qualifiers['product'][0]
                                    except KeyError:
                                        product = "tRNA-XXX"
                                    chr = record.id
                                    gff.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n" % (chr, start, end, strand, ID))
                                    gff.write("%s\tGenBank\ttRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, ID, product))
                                    gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon1;Parent=%s-T1\n" % (chr, start, end, strand, ID, ID))

def runGMAP(transcripts, genome, cpus, intron, tmpdir, output):
    #first build genome database
    build_log = os.path.join(tmpdir, 'gmap-build.log')
    with open(build_log, 'w') as logfile:
        subprocess.call(['gmap_build', '-D', tmpdir, '-d', 'genome', '-k', '13', genome], stdout = logfile, stderr = logfile)
    #now map transcripts
    map_log = os.path.join(tmpdir, 'gmap-map.log')
    with open(map_log, 'w') as logfile:
        with open(output, 'w') as out:
            subprocess.call(['gmap', '--cross-species', '-f', '3', '-K', str(intron), '-n', '1', '-t', str(cpus), '-B', '5', '-D', tmpdir, '-d', 'genome', transcripts], stdout = out, stderr = logfile)
    
def runBUSCO(input, DB, cpus, tmpdir, output):
    #run busco in protein mapping mode
    BUSCO = os.path.join(UTIL, 'funannotate-BUSCO2.py')
    cmd = [BUSCO, '-i', input, '-m', 'proteins', '-l', DB, '-o', 'busco', '-c', str(cpus), '-f']
    runSubprocess(cmd, tmpdir, log)
    #now parse output and write to annotation file
    with open(output, 'w') as out:
        with open(os.path.join(tmpdir, 'run_busco', 'full_table_busco.tsv'), 'rU') as busco:
            for line in busco:
                if line.startswith('#'):
                    continue
                col = line.split('\t')
                if col[1] == 'Complete' or col[1] == 'Duplicated': #if diploid these should show up, but problematic for drawing trees....
                    if col[2].endswith('-T1'):
                        ID = col[2]
                    else:
                        ID = col[2]+'-T1'
                    out.write("%s\tnote\tBUSCO:%s\n" % (ID, col[0]))

def dupBUSCO2gff(ID, base_folder, locationID):
    hmmerfolder = os.path.join(base_folder, 'hmmer_output')
    geneID = ''
    AugFile = ''
    GFFfile = os.path.join(base_folder, 'augustus_output', 'gffs', ID+'.gff')
    if geneID == '':
        for file in os.listdir(hmmerfolder):
            if file.startswith(ID):
                with open(os.path.join(hmmerfolder, file), 'rU') as hmmer:
                    for line in hmmer:
                        if not line.startswith('#'):
                            longID = line.split()[0]
                            longID = longID.replace(']', '')
                            partsID = longID.split('[')
                            if locationID == partsID[1]:
                                geneID = partsID[0]
                                AugFile = os.path.join(base_folder, 'augustus_output', 'predicted_genes', file)
                                break
    #so now should have gene name, get the GFF from augustus
    with open(GFFfile, 'w') as gffout:
        with open(AugFile, 'rU') as augustus:
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
     
                  
def parseBUSCO2genome(input, ploidy, ContigSizes, output):
    #input is BUSCO output, ploidy is integer, ContigSizes is dictionary, output is a bedfile, function returns dictionary
    busco_complete = {}
    hits = {}
    with open(output, 'w') as bedfile:
        with open(input, 'rU') as buscoinput:
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
                            hits[cols[0]] = (ID,score,contig,start,end,length)
                    if ploidy > 1:
                        if cols[1] == 'Duplicated':
                            if not cols[0] in hits:
                                hits[cols[0]] = (ID,score,contig,start,end,length)
                                dupBUSCO2gff(cols[0], os.path.dirname(input), ID)
                            else:
                                oldscore = float(hits.get(cols[0])[1])
                                if float(score) > oldscore:
                                    hits[cols[0]] = (ID,score,contig,start,end,length)
                                    dupBUSCO2gff(cols[0], os.path.dirname(input), ID)
            for k,v in natsorted(hits.items()):
                #validate locations for bedfile, move 100 bp in each direction for bedfile
                start = int(v[3]) - 100
                if start < 1: #negative no good
                    start = 1
                end = int(v[4]) + 100
                if end > ContigSizes.get(contig): #check it doesn't go past contig length
                    end = ContigSizes.get(contig)
                bedfile.write('%s\t%i\t%i\t%s\n' % (contig,start,end,k))         
                busco_complete[k] = v[0]
    return busco_complete

def SwissProtBlast(input, cpus, evalue, tmpdir, output):
    #run blastp against uniprot
    blast_tmp = os.path.join(tmpdir, 'uniprot.xml')
    blastdb = os.path.join(DB, 'uniprot')
    cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus), '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input]
    runSubprocess(cmd, '.', log)
    #parse results
    with open(output, 'w') as out:
        with open(blast_tmp, 'rU') as results:
            for qresult in SearchIO.parse(results, "blast-xml"):
                hits = qresult.hits
                qlen = qresult.seq_len
                ID = qresult.id
                num_hits = len(hits)
                if num_hits > 0:
                    length = hits[0].hsps[0].aln_span
                    pident = hits[0].hsps[0].ident_num / float(length)
                    if pident < 0.6:
                        continue
                    diff = length / float(qlen)
                    if diff < 0.6:
                        continue
                    description = hits[0].description.split("=")
                    hdescript = description[0].replace(' OS','')
                    name = description[2].replace(' PE','').upper()
                    #need to do some filtering here of certain words
                    bad_words = ['(Fragment)', 'homolog', 'homolog,']
                    descript = hdescript.split(' ') #turn string into array, splitting on spaces
                    final_desc = [x for x in descript if x not in bad_words]
                    final_desc = ' '.join(final_desc)
                    #okay, print out annotations for GAG
                    if ID.endswith('-T1'):
                        out.write("%s\tprot_desc\t%s\n" % (ID,final_desc))
                        geneID = ID.replace('-T1','')
                    else:
                        mrnaID = ID + '-T1'
                        out.write("%s\tprot_desc\t%s\n" % (mrnaID,final_desc))

def RepeatBlast(input, cpus, evalue, tmpdir, output):
    #run blastp against repeats
    blast_tmp = os.path.join(tmpdir, 'repeats.xml')
    blastdb = os.path.join(DB,'REPEATS')
    cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus), '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input]
    runSubprocess(cmd, '.', log)
    #parse results   
    with open(output, 'w') as out:
        with open(blast_tmp, 'rU') as results:
            for qresult in SearchIO.parse(results, "blast-xml"):
                hits = qresult.hits
                qlen = qresult.seq_len
                ID = qresult.id
                num_hits = len(hits)
                if num_hits > 0:
                    length = 0
                    for i in range(0,len(hits[0].hsps)):
                        length += hits[0].hsps[i].aln_span
                    pident = hits[0].hsps[0].ident_num / float(length)
                    out.write("%s\t%s\t%f\t%s\n" % (ID, hits[0].id, pident, hits[0].hsps[0].evalue))

def MEROPSBlast(input, cpus, evalue, tmpdir, output):
    #run blastp against merops
    blast_tmp = os.path.join(tmpdir, 'merops.xml')
    blastdb = os.path.join(DB,'MEROPS')
    cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus), '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input]
    runSubprocess(cmd, '.', log)
    #parse results
    with open(output, 'w') as out:
        with open(blast_tmp, 'rU') as results:
            for qresult in SearchIO.parse(results, "blast-xml"):
                hits = qresult.hits
                qlen = qresult.seq_len
                ID = qresult.id
                num_hits = len(hits)
                if num_hits > 0:
                    if hits[0].hsps[0].evalue > evalue:
                        continue
                    sseqid = hits[0].id
                    family = hits[0].description
                    #okay, print out annotations for GAG
                    if not ID.endswith('-T1'):
                        ID = ID + '-T1'
                    out.write("%s\tnote\tMEROPS:%s\n" % (ID,sseqid))

def eggnog2dict(annotations):
    #load in annotation dictionary
    EggNog = {}
    with open(annotations, 'rU') as input:
        reader = csv.reader(input, delimiter='\t')
        for line in reader:
            EggNog[line[1]] = line[5]
    return EggNog

def runEggNog(file, HMM, annotations, cpus, evalue, tmpdir, output):
    #kind of hacky, but hmmersearch doesn't allow me to get sequence length from hmmer3-text, only domtbl, but then I can't get other values, so read seqlength into dictionary for lookup later.
    SeqLength = {}
    with open(file, 'rU') as proteins:
        SeqRecords = SeqIO.parse(proteins, 'fasta')
        for rec in SeqRecords:
            length = len(rec.seq)
            SeqLength[rec.id] = length
    #run hmmerscan
    eggnog_out = os.path.join(tmpdir, 'eggnog.txt')
    cmd = ['hmmsearch', '-o', eggnog_out, '--cpu', str(cpus), '-E', str(evalue), HMM, file]
    runSubprocess3(cmd, '.', log)
    #now parse results
    Results = {}
    with open(output, 'w') as out:
        with open(eggnog_out, 'rU') as results:
            for qresult in SearchIO.parse(results, "hmmer3-text"):
                query_length = qresult.seq_len #length of HMM model
                hits = qresult.hits
                num_hits = len(hits)
                if num_hits > 0:
                    for i in range(0,num_hits):
                        if round(hits[i].domain_exp_num) != hits[i].domain_obs_num: #make sure # of domains is nearly correct
                            continue
                        query = hits[i].id
                        #get total length from dictionary
                        seq_length = SeqLength.get(query)
                        lower = seq_length * 0.75
                        upper = seq_length * 1.25
                        aln_length = 0
                        num_hsps = len(hits[i].hsps)
                        for x in range(0,num_hsps):
                            aln_length += hits[i].hsps[x].aln_span
                        if aln_length < lower or aln_length > upper: #make sure most of the protein aligns to the model, 50% flex
                            continue
                        if not query.endswith('-T1'):
                            query = query + '-T1'
                        hit = hits[i].query_id.split(".")[1]
                        score = hits[i].bitscore
                        evalue = hits[i].evalue
                        if not query in Results:
                            Results[query] = (hit, score, evalue, seq_length, aln_length)
                        else:
                            OldScore = Results.get(query)[1]
                            if score > OldScore:
                                Results[query] = (hit, score, evalue, seq_length, aln_length)
            #look up descriptions in annotation dictionary
            for k, v in Results.items():
                out.write("%s\tnote\tEggNog:%s\n" % (k, v[0]))

def PFAMsearch(input, cpus, evalue, tmpdir, output):
    #run hmmerscan
    HMM = os.path.join(DB, 'Pfam-A.hmm')
    pfam_out = os.path.join(tmpdir, 'pfam.txt')
    pfam_filtered = os.path.join(tmpdir, 'pfam.filtered.txt')
    cmd = ['hmmsearch', '--domtblout', pfam_out, '--cpu', str(cpus), '-E', str(evalue), HMM, input]
    runSubprocess3(cmd, '.', log)
    #now parse results
    with open(output, 'w') as out:
        with open(pfam_filtered, 'w') as filtered:
            with open(pfam_out, 'rU') as results:
                for qresult in SearchIO.parse(results, "hmmsearch3-domtab"):
                    hits = qresult.hits
                    num_hits = len(hits)
                    if num_hits > 0:
                        for i in range(0,num_hits):
                            hit_evalue = hits[i].evalue
                            if hit_evalue > evalue:
                                continue
                            query = hits[i].id
                            pfam = qresult.accession.split('.')[0]
                            hmmLen = qresult.seq_len
                            hmm_aln = int(hits[i].hsps[0].hit_end) - int(hits[i].hsps[0].hit_start)
                            coverage = hmm_aln / float(hmmLen)
                            if coverage < 0.50: #coverage needs to be at least 50%
                                continue
                            hit = hits[i].query_id
                            #description = hits[i].description
                            if not query.endswith('-T1'):
                                query = query + '-T1'
                            filtered.write("%s\t%s\t%s\t%f\n" % (query, pfam, hit_evalue, coverage))
                            out.write("%s\tdb_xref\tPFAM:%s\n" % (query, pfam))



def dbCANsearch(input, cpus, evalue, tmpdir, output):
    CAZY = {'CBM': 'Carbohydrate-binding module', 'CE': 'Carbohydrate esterase','GH': 'Glycoside hydrolase', 'GT': 'Glycosyltransferase', 'PL': 'Polysaccharide lyase', 'AA': 'Auxillary activities'}
    #run hmmerscan
    HMM = os.path.join(DB, 'dbCAN.hmm')
    dbCAN_out = os.path.join(tmpdir, 'dbCAN.txt')
    dbCAN_filtered = os.path.join(tmpdir, 'dbCAN.filtered.txt')
    cmd = ['hmmscan', '--domtblout', dbCAN_out, '--cpu', str(cpus), '-E', str(evalue), HMM, input]
    runSubprocess3(cmd, '.', log)
    #now parse results
    with open(output, 'w') as out:
        with open(dbCAN_filtered, 'w') as filtered:
            filtered.write("#HMM_family\tHMM_len\tQuery_ID\tQuery_len\tE-value\tHMM_start\tHMM_end\tQuery_start\tQuery_end\tCoverage\n")
            with open(dbCAN_out, 'rU') as results:
                for qresult in SearchIO.parse(results, "hmmscan3-domtab"):
                    query_length = qresult.seq_len
                    hits = qresult.hits
                    num_hits = len(hits)
                    if num_hits > 0:
                        for i in range(0,num_hits):
                            hit_evalue = hits[i].evalue
                            if hit_evalue > evalue:
                                continue
                            hit = hits[i].id
                            hmmLen = hits[i].seq_len
                            hmm_aln = int(hits[i].hsps[0].hit_end) - int(hits[i].hsps[0].hit_start)
                            coverage = hmm_aln / float(hmmLen)
                            if coverage < 0.45:
                                continue
                            query = hits[i].query_id
                            filtered.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\n" % (hit, hmmLen, query, query_length, hit_evalue, hits[i].hsps[0].hit_start, hits[i].hsps[0].hit_end, hits[i].hsps[0].query_start, hits[i].hsps[0].query_end, coverage))
                            #get type of hit for writing the annotation note
                            type = ''.join(i for i in hit if not i.isdigit())
                            descript = CAZY.get(type)
                            if not query.endswith('-T1'):
                                query = query + '-T1'
                            out.write("%s\tnote\tCAZy:%s\n" % (query, hit))

def batch_iterator(iterator, batch_size):
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def fasta2chunks(input, chunks, tmpdir, output):
    #split the input fasta file into 20 chunks to process
    with open(input, 'rU') as seqs:
        SeqCount = countfasta(input)
        SeqRecords = SeqIO.parse(seqs, 'fasta')
        chunks = SeqCount / int(chunks)
        #divide into chunks, store in tmp file
        folder = os.path.join(tmpdir, output)
        if not os.path.exists(folder):
            os.makedirs(folder)
        else:
            shutil.rmtree(folder)
            os.makedirs(folder)
        for i, batch in enumerate(batch_iterator(SeqRecords, chunks)) :
            filename = "chunk_%i.fa" % (i+1)
            tmpout = os.path.join(folder, filename)
            handle = open(tmpout, "w")
            count = SeqIO.write(batch, handle, "fasta")
            handle.close()

def signalP(input, tmpdir, output):
    #split input file into chunks, 20 should mean < 200 proteins per chunk
    fasta2chunks(input, 40, tmpdir, 'signalp_tmp')
    for file in os.listdir(os.path.join(tmpdir, 'signalp_tmp')):
        if file.startswith('chunk'):
            file = os.path.join(tmpdir, 'signalp_tmp', file)
            tmp_out = file.replace('.fa', '.signalp.out')
            cmd = ['signalp', '-t', 'euk', '-f', 'short', file]
            runSubprocess2(cmd, '.', log, tmp_out)
    #now concatenate all outputs
    if os.path.isfile(output):
        os.remove(output)            
    with open(output, 'a') as finalout:
        for file in os.listdir(os.path.join(tmpdir, 'signalp_tmp')):
            if file.endswith('.signalp.out'):
                file = os.path.join(tmpdir, 'signalp_tmp', file)
                with open(file) as infile:
                    finalout.write(infile.read())
    #cleanup tmp directory
    shutil.rmtree(os.path.join(tmpdir, 'signalp_tmp'))
    
def parsePhobiusSignalP(phobius, sigP, membrane_annot, secretome_annot):
    #give directory of annotate_misc, first get phobius results
    '''
    This is what phobius results look like
    ID	TM	SP	Prediction
    VE00_00001	0	0	o
    VE00_00002	2	0	i198-219o283-301i
    VE00_00003	0	0	o
    VE00_00004	0	Y	n8-18c23/24o
    VE00_00005	12	0	i49-69o89-107i119-138o144-167i179-200o212-234i280-299o319-341i348-366o378-398i410-430o442-465i
    '''
    pSecDict = {}
    pTMDict = {}
    sigpDict = {}
    with open(phobius, 'rU') as input1:
        for line in input1:
            line = line.replace('\n', '')
            if line.startswith('ID\t'):
                continue
            cols = line.split('\t')
            geneID = cols[0]
            if not geneID.endswith('-T1'):
                geneID = geneID + '-T1'
            if int(cols[1]) > 0: #then found TM domain
                annot = cols[3]
                if '/' in annot:
                    annotation = annot.split('/')[-1]
                if not geneID in pTMDict:
                    pTMDict[geneID] = 'TransMembrane:'+cols[1]+' ('+annot+')'
            if cols[2] == 'Y': #then sig pep discovered
                location = cols[3].split('/')[0]
                clevage = location.split('c')[-1]
                if not geneID in pSecDict:
                    pSecDict[geneID] = clevage
    if sigP: #will be passed FALSE if signalP data missing
        #parse signalp output and turn into annotation file
        with open(sigP, 'rU') as results:
            for line in results:
                line = line.replace('\n', '')
                if line.startswith('#'):
                    continue
                col = line.split(' ') #not tab delimited
                col = filter(None, col) #clean up empty spaces
                if col[9] == 'Y': #then there is signal peptide
                    ID = col[0]
                    if not ID.endswith('-T1'):
                        ID = ID + '-T1'
                    end = int(col[2]) - 1
                    #save as secreted only if also found in phobius
                    if ID in pSecDict:
                        sigpDict[ID] = end
    else:
        sigpDict = pSecDict
    #write annotation files
    with open(membrane_annot, 'w') as memout:
        for k,v in natsorted(pTMDict.items()):
            memout.write("%s\tnote\t%s\n" % (k, v))
    with open(secretome_annot, 'w') as secout:
         for k,v in natsorted(sigpDict.items()):   
            secout.write("%s\tnote\tSECRETED:SignalP(1-%s)\n" % (k, v))
                
def RepeatModelMask(input, cpus, tmpdir, output, debug):
    log.info("Loading sequences and soft-masking genome")
    outdir = os.path.join(tmpdir, 'RepeatModeler')
    input = os.path.abspath(input)
    output = os.path.abspath(output)
    #lets run RepeatModeler here to get repeat library
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    log.info("Soft-masking: building RepeatModeler database")
    with open(debug, 'a') as debug_log:
        subprocess.call(['BuildDatabase', '-name', 'Repeats', input], cwd=outdir, stdout = debug_log, stderr=debug_log)
    log.info("Soft-masking: generating repeat library using RepeatModeler")
    with open(debug, 'a') as debug_log:
        subprocess.call(['RepeatModeler', '-e', 'ncbi', '-database', 'Repeats', '-pa', str(cpus)], cwd=outdir, stdout = debug_log, stderr=debug_log)
    #find name of folder
    for i in os.listdir(outdir):
        if i.startswith('RM_'):
            RP_folder = i
    library = os.path.join(tmpdir, 'repeatmodeler.lib.fa')
    library = os.path.abspath(library)
    try:
        os.rename(os.path.join(outdir, RP_folder, 'consensi.fa.classified'), library)
    except OSError:
        pass
    #now soft-mask the genome for gene predictors
    outdir2 = os.path.join(tmpdir, 'RepeatMasker')
    if os.path.isdir(outdir2):
        shutil.rmtree(outdir2)
    os.makedirs(outdir2)
    if not os.path.isfile(library):
        log.info("Soft-masking: running RepeatMasker with default library (RepeatModeler found 0 models)")
        with open(debug, 'a') as debug_log:
            subprocess.call(['RepeatMasker', '-e', 'ncbi', '-pa', str(cpus), '-xsmall', '-dir','.', input], cwd=outdir2, stdout=debug_log, stderr = debug_log)
    else:
        log.info("Soft-masking: running RepeatMasker with custom library")
        with open(debug, 'a') as debug_log:
            subprocess.call(['RepeatMasker', '-e', 'ncbi', '-lib', library, '-pa', str(cpus), '-xsmall', '-dir', '.', input], cwd=outdir2, stdout=debug_log, stderr = debug_log)
    for file in os.listdir(outdir2):
        if file.endswith('.masked'):
            shutil.copyfile(os.path.join(outdir2, file), output)
        if file.endswith('.out'):
            rm_gff3 = os.path.join(tmpdir, 'repeatmasker.gff3')
            cmd = ['rmOutToGFF3.pl', file]
            runSubprocess2(cmd, outdir2, log, rm_gff3)

def RepeatMask(input, library, cpus, tmpdir, output, debug):
    FNULL = open(os.devnull, 'w')
    outdir = os.path.join(tmpdir, 'RepeatMasker')
    #now soft-mask the genome for gene predictors
    log.info("Soft-masking: running RepeatMasker with custom library")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    with open(debug, 'a') as debug_log:
        subprocess.call(['RepeatMasker', '-e', 'ncbi', '-lib', os.path.abspath(library), '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input], stderr = debug_log, stdout=debug_log, cwd = tmpdir)
    for file in os.listdir(outdir):
        if file.endswith('.masked'):
            os.rename(os.path.join(outdir, file), output)
        if file.endswith('.out'):
            rm_gff3 = os.path.join(tmpdir, 'repeatmasker.gff3')
            cmd = ['rmOutToGFF3.pl', file]
            runSubprocess2(cmd, outdir, log, rm_gff3)
    
def CheckAugustusSpecies(input):
    #get the possible species from augustus
    augustus_list = []
    for i in os.listdir(os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], 'species')):
        if not i.startswith('.'):
            augustus_list.append(i)
    augustus_list = set(augustus_list)
    if input in augustus_list:
        return True
    else:
        return False

def SortRenameHeaders(input, output):
    #sort records and write temp file
    with open(output, 'w') as out:
        with open(input, 'rU') as input:
            records = list(SeqIO.parse(input, 'fasta'))
            records.sort(cmp=lambda x,y: cmp(len(y),len(x)))
            counter = 1
            for rec in records:
                rec.name = ''
                rec.description = ''
                rec.id = 'scaffold_' + str(counter)
                counter +=1
            SeqIO.write(records, out, 'fasta')

def RunGeneMarkES(input, cpus, tmpdir, output, fungus):
    #make directory to run script from
    outdir = os.path.join(tmpdir, 'genemark')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    contigs = os.path.abspath(input)
    log.info("Running GeneMark-ES on assembly")
    if fungus:
        cmd = ['gmes_petap.pl', '--ES', '--fungus', '--soft_mask', '5000', '--cores', str(cpus), '--sequence', contigs]
    else:
        cmd = ['gmes_petap.pl', '--ES', '--soft_mask', '5000', '--cores', str(cpus), '--sequence', contigs]
    runSubprocess3(cmd, outdir, log)
    #rename results and grab mod file
    try:
        os.rename(os.path.join(outdir,'output','gmhmm.mod'), os.path.join(tmpdir, 'gmhmm.mod'))
    except OSError:
        log.error("GeneMark-ES failed, likely input was not sufficient for training, provide a gmhmm.mod file and re-run")
        sys.exit(1)
    #convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join(outdir, 'genemark.gtf')
    log.info("Converting GeneMark GTF file to GFF3")
    with open(output, 'w') as out:
        subprocess.call([GeneMark2GFF, gm_gtf], stdout = out)
    log.info('Found {0:,}'.format(countGFFgenes(output)) +' gene models')
        
def RunGeneMark(input, mod, cpus, tmpdir, output, fungus):
    #make directory to run script from
    outdir = os.path.join(tmpdir, 'genemark')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    contigs = os.path.abspath(input)
    mod = os.path.abspath(mod)
    log.info("Running GeneMark-ES on assembly")
    if fungus:
        cmd = ['gmes_petap.pl', '--ES', '--soft_mask', '5000', '--ini_mod', mod, '--fungus', '--cores', str(cpus), '--sequence', contigs]
    else:
        cmd = ['gmes_petap.pl', '--ES', '--soft_mask', '5000', '--ini_mod', mod, '--cores', str(cpus), '--sequence', contigs]
    runSubprocess3(cmd, outdir, log)
    #convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join(outdir, 'genemark.gtf')
    log.info("Converting GeneMark GTF file to GFF3")
    with open(output, 'w') as out:
        subprocess.call([GeneMark2GFF, gm_gtf], stdout = out)
    log.info('Found {0:,}'.format(countGFFgenes(output)) +' gene models')

def MemoryCheck():
    import psutil
    mem = psutil.virtual_memory()
    RAM = int(mem.total)
    return round(RAM / 1024000000)

def systemOS():
    if sys.platform == 'darwin':
        system_os = 'MacOSX '+ platform.mac_ver()[0]
    elif sys.platform == 'linux':
        linux_version = platform.linux_distribution()
        system_os = linux_version[0]+ ' '+linux_version[1]
    else:
        system_os = sys.platform
    return system_os

def SystemInfo():
    system_os = systemOS()
    python_vers = str(sys.version_info[0])+'.'+str(sys.version_info[1])+'.'+str(sys.version_info[2])   
    log.info("OS: %s, %i cores, ~ %i GB RAM. Python: %s" % (system_os, multiprocessing.cpu_count(), MemoryCheck(), python_vers))    

def runtRNAscan(input, tmpdir, output):
    tRNAout = os.path.join(tmpdir, 'tRNAscan.out')
    if os.path.isfile(tRNAout): # tRNAscan can't overwrite file, so check first
        os.remove(tRNAout)
    cmd = ['tRNAscan-SE', '-o', tRNAout, input]
    runSubprocess(cmd, '.', log)
    trna2gff = os.path.join(UTIL, 'trnascan2gff3.pl')
    with open(output, 'w') as out:
        subprocess.call(['perl', trna2gff, '--input', tRNAout], stdout = out)
    log.info('Found {0:,}'.format(countGFFgenes(output)) +' tRNA gene models')


def gb2smurf(input, prot_out, smurf_out):
    with open(smurf_out, 'w') as smurf:
        with open(prot_out, 'w') as proteins:
            with open(input, 'rU') as gbk:
                SeqRecords = SeqIO.parse(gbk, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        name = re.sub('[^0-9]','', record.name)
                        if f.type == "CDS":
                            proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0].rstrip('*')))
                            locus_tag = f.qualifiers.get("locus_tag", ["No ID"])[0]
                            product_name = f.qualifiers.get("product", ["No Description"])[0]
                            mystart = f.location.start
                            myend = f.location.end
                            strand = f.location.strand
                            if strand == 1:
                                smurf.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(mystart), int(myend), product_name))
                            else:
                                smurf.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(myend), int(mystart), product_name))
                            
def RemoveBadModels(proteins, gff, length, repeats, BlastResults, tmpdir, output):
    #first run bedtools to intersect models where 90% of gene overlaps with repeatmasker region
    repeat_temp = os.path.join(tmpdir, 'genome.repeats.to.remove.gff')
    cmd = ['bedtools', 'intersect', '-f', '0.9', '-a', gff, '-b', repeats]
    runSubprocess2(cmd, '.', log, repeat_temp)
    #now remove those proteins that do not have valid starts, less then certain length, and have internal stops
    remove = []
    reason = {}
    #parse the results from bedtools and add to remove list
    with open(repeat_temp, 'rU') as input:
        for line in input:
            if "\tgene\t" in line:
                ninth = line.split('ID=')[-1]
                ID = ninth.split(";")[0]
                remove.append(ID)
                if not ID in reason:
                    reason[ID] = 'remove_reason=repeat_overlap;'
    #parse the results from BlastP search of transposons
    with open(BlastResults, 'rU') as input:
        for line in input:
            col = line.split('\t')
            remove.append(col[0])
            if not col[0] in reason:
                ID = col[0].replace('evm.model.', 'evm.TU.')
                reason[ID] = 'remove_reason=repeat_match;'
            else:
                ID = col[0].replace('evm.model.', 'evm.TU.')
                reason[ID] = 'remove_reason=repeat_overalap|repeat_match;'
 
    #I'm only seeing these models with GAG protein translations, so maybe that is a problem? skip for now
    with open(proteins, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)[:-1]
            '''
            if not Seq.startswith('M'):
                remove.append(rec.id)
                if not rec.id in reason:
                    ID = rec.id.replace('evm.model.', 'evm.TU.')
                    reason[ID] = 'remove_reason=bad_start;'
            '''        
            if len(Seq) < int(length):
                remove.append(rec.id)
                if not rec.id in reason:
                    ID = rec.id.replace('evm.model.', 'evm.TU.')
                    reason[ID] = 'remove_reason=seq_too_short;'
            if 'XX' in Seq:
                remove.append(rec.id)
                if not rec.id in reason:
                    ID = rec.id.replace('evm.model.', 'evm.TU.')
                    reason[ID] = 'remove_reason=model_span_gap;'
    remove = [w.replace('evm.TU.','') for w in remove]
    remove = [w.replace('evm.model.','') for w in remove]
    remove = set(remove)
    if len(remove) > 0:
        remove_match = re.compile(r'\b\.(?:%s)[\.;]\b' % '|'.join(remove))
        with open(output, 'w') as out:
            with open(os.path.join(tmpdir, 'bad_models.gff'), 'w') as out2:
                with open(gff, 'rU') as GFF:
                    for line in GFF:
                        if '\tstart_codon\t' in line:
                            continue
                        if '\tstop_codon\t' in line:
                            continue
                        if not remove_match.search(line):
                            line = re.sub(';Name=.*$', ';', line) #remove the Name attribute as it sticks around in GBK file
                            out.write(line)           
                        else:
                            if "\tgene\t" in line:
                                bad_ninth = line.split('ID=')[-1]
                                bad_ID = bad_ninth.split(";")[0]                
                                bad_reason = reason.get(bad_ID)
                                if bad_reason:
                                    line = line.replace('\n', ';'+bad_reason+'\n')
                                else:
                                    log.debug("%s was removed in removeBadModels function for unknown reason, please check manually" % bad_ID)
                                    line = line.replace('\n', ';remove_reason=unknown;\n')
                            out2.write(line)
    else: #if nothing to remove, just print out GFF
        with open(output, 'w') as out:
            with open(gff, 'rU') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue
                    line = re.sub(';Name=.*$', ';', line) #remove the Name attribute as it sticks around in GBK file
                    out.write(line)
                               
def CleantRNAtbl(GFF, TBL, output):
    #clean up genbank tbl file from gag output
    #try to read through GFF file, make dictionary of tRNA genes and products
    TRNA = {}
    with open(GFF, 'rU') as gff:
        for line in gff:
            if '\ttRNA\t' in line:
                cols = line.split('\t')
                ID = cols[8].split(';')[0].replace('ID=', '')
                ID = ID.replace('-T1', '')
                product = cols[8].split('product=')[-1].replace('\n', '')
                TRNA[ID] = product
    #print TRNA           
    with open(output, 'w') as out:
        with open(TBL, 'rU') as input:
            for line in input:
                if line.startswith('\t\t\tlocus_tag\t'):
                    out.write(line)
                    geneID = line.split('locus_tag\t')[-1].replace('\n', '')
                    if geneID in TRNA:
                        if 'tRNA-Xxx' == TRNA.get(geneID):
                            out.write("\t\t\tpseudo\n")       
                elif line.startswith("\t\t\tproduct\ttRNA-Xxx"):
                    out.write(line)
                    out.write("\t\t\tpseudo\n")
                    input.next()
                    input.next()
                elif line.startswith("\t\t\tproduct\ttRNA-"):
                    out.write(line)
                    input.next()
                    input.next()
                else:
                    out.write(line)

def ParseErrorReport(input, Errsummary, val, Discrep, output, keep_stops):
    errors = []
    gapErrors = []
    remove = []
    with open(Errsummary) as summary:
        for line in summary:
            if 'ERROR' in line:
                if 'SEQ_DESCR.OrganismIsUndefinedSpecies' in line: #there are probably other errors you are unaware of....
                    pass
                elif 'SEQ_FEAT.MissingTrnaAA' in line:
                    pass
                elif 'SEQ_INST.TerminalNs' in line:
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
    #parse the discrepency report and look for overlapping genes, so far, all have been tRNA's in introns, so just get those for now.
    with open(Discrep, 'rU') as discrep:
        for line in discrep:
            if line.startswith('DiscRep_ALL:FIND_OVERLAPPED_GENES::'): #skip one line and then move through next lines until line starts with nothing
                num = line.split(' ')[0]
                num = num.split('::')[-1]
                num = int(num)
                for i in range(num):
                    gene = discrep.next().split('\t')[1]
                    tRNA = gene + '_tRNA'
                    exon = gene + '_exon'
                    remove.append(gene)
                    remove.append(tRNA)
                    remove.append(exon)              
    if len(errors) < 1 and len(remove) < 1: #there are no errors, then just remove stop/start codons and move on
        with open(output, 'w') as out:
            with open(input, 'rU') as GFF:
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
                #this is only picking up tRNAs right now, which "probably" is all that it needs to.....but u never know
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
        remove = set(remove)
        remove_match = re.compile(r'\b(?:%s)+\b' % '|'.join(remove))
        with open(output, 'w') as out:
            with open(input, 'rU') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue                 
                    if not remove_match.search(line):
                        out.write(line)
                        
def ParseAntiSmash(input, tmpdir, output, annotations):
    log.info("Now parsing antiSMASH results, finding SM clusters")
    global bbDomains, bbSubType, BackBone
    BackBone = {}; SMCOGs = {}; bbSubType = {}; bbDomains = {}; smProducts = {}
    backboneCount = 0; clusterCount = 0; cogCount = 0
    #parse antismash genbank to get clusters in bed format and slice the record for each cluster prediction
    with open(output, 'w') as antibed:
        with open(input, 'rU') as input:
            SeqRecords = SeqIO.parse(input, 'genbank')
            for record in SeqRecords:
                for f in record.features:
                    if f.type == "source":
                        record_start = f.location.start
                        record_end = f.location.end
                    if f.type == "cluster":
                        clusterCount += 1
                        chr = record.id
                        start = f.location.start
                        end = f.location.end
                        clusternum = f.qualifiers.get("note")[0].replace("Cluster number: ", "")
                        antibed.write("%s\t%s\t%s\tCluster_%s\t0\t+\n" % (chr, start, end, clusternum))
                    Domains = []
                    if f.type == "CDS":
                        ID = f.qualifiers.get('locus_tag')[0]                    
                        if f.qualifiers.get('sec_met'):            
                            for k, v in f.qualifiers.items():
                                if k == 'sec_met':
                                    for i in v:
                                        if i.startswith('Type:'):
                                            type = i.replace('Type: ', '')
                                            backboneCount += 1
                                            BackBone[ID] = type
                                        if i.startswith('NRPS/PKS subtype:'):
                                            subtype = i.replace('NRPS/PKS subtype: ', '')
                                            bbSubType[ID] = subtype
                                        if i.startswith('NRPS/PKS Domain:'):
                                            doms = i.replace('NRPS/PKS Domain: ', '')
                                            doms = doms.split('. ')[0]
                                            Domains.append(doms)
                                bbDomains[ID] = Domains
                        for k,v in f.qualifiers.items():
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
                            
    log.info("Found %i clusters, %i biosynthetic enyzmes, and %i smCOGs predicted by antiSMASH" % (clusterCount, backboneCount, cogCount))
    #now generate the annotations to add to genome
    with open(annotations, 'w') as out:
        #add product annotations - use bbSubType --> BackBone
        for k, v in BackBone.items():
            if k in bbSubType:
                if not k.endswith('-T1'):
                    ID = k + '-T1'
                else:
                    ID = k
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
            out.write("%s\tproduct\t%s\n" % (ID, hit))          
        #add annots from smProducts
        for k, v in smProducts.items():
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            out.write("%s\tproduct\t%s\n" % (ID, v))               
        #add smCOGs into note section
        for k, v in SMCOGs.items():
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            out.write("%s\tnote\t%s\n" % (ID, v))
              
def GetClusterGenes(input, GFF, output, annotations):
    global dictClusters
    #pull out genes in clusters from GFF3, load into dictionary
    cmd = ['bedtools', 'intersect','-wo', '-a', input, '-b', GFF]
    runSubprocess2(cmd, '.', log, output)
    dictClusters = {}
    with open(output, 'rU') as input:
        for line in input:
            cols = line.split('\t')
            if cols[8] != 'gene':
                continue
            gene = cols[14].replace('ID=', '')
            ID = cols[3]
            if ID not in dictClusters:
                dictClusters[ID] = [gene]
            else:
                dictClusters[ID].append(gene)
    with open(annotations, 'w') as annotout: 
        for k, v in dictClusters.items():
            for i in v:
                if not i.endswith('-T1'):
                    ID = i + ('-T1')
                else:
                    ID = i
                annotout.write("%s\tnote\tantiSMASH:%s\n" % (ID, k))

def splitFASTA(input, outputdir):
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    with open(input, 'rU') as InputFasta:
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
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            lengths.append(len(record.seq))
            GeeCee.append(str(record.seq))
            for f in record.features:
                if f.type == "source":
                    organism = f.qualifiers.get("organism", ["???"])[0]
                    isolate = f.qualifiers.get("isolate", ["???"])[0]
                    if isolate == "???":
                        isolate = f.qualifiers.get("strain", ["???"])[0]
                if f.type == "CDS":
                    Prots += 1
                if f.type == "gene":
                    Genes += 1
                    if Genes == 1:
                        locus_tag = f.qualifiers.get("locus_tag")[0].split('_')[0]
                if f.type == "tRNA":
                    tRNA += 1
    log.info("working on %s genome" % organism)
    GenomeSize = sum(lengths)
    LargestContig = max(lengths)
    ContigNum = len(lengths)
    AvgContig = int(round(GenomeSize / ContigNum))
    pctGC = round(GC("".join(GeeCee)), 2)
    
    #now get N50
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
    #return values in a list
    return [organism, isolate, locus_tag, "{0:,}".format(GenomeSize)+' bp', "{0:,}".format(LargestContig)+' bp', "{0:,}".format(AvgContig)+' bp', "{0:,}".format(ContigNum), "{0:,}".format(N50)+' bp', "{:.2f}".format(pctGC)+'%', "{0:,}".format(Genes), "{0:,}".format(Prots), "{0:,}".format(tRNA)]

def MEROPS2dict(input):
    dict = {}
    with open(input, 'rU') as fasta:
        for line in fasta:
            if line.startswith('>'):
                cols = line.split(' ')
                ID = cols[0].replace('>', '')
                family = cols[1].replace('\n', '')
                dict[ID] = family
    return dict

def getEggNogfromNote(input):
    dict = {}
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    ID = f.qualifiers['locus_tag'][0]
                    for k,v in f.qualifiers.items():
                        if k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith('EggNog:'):
                                    hit = i.replace('EggNog:', '')
                                    if not ID in dict:
                                        dict[ID] = hit
    return dict
                                      
def getStatsfromNote(input, word):
    dict = {}
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    ID = f.qualifiers['locus_tag'][0]
                    for k,v in f.qualifiers.items():
                        if k == 'note':
                            notes = v[0].split('; ')
                            for i in notes:
                                if i.startswith(word+':'):
                                    hit = i.replace(word+':', '')
                                    if hit.startswith('MER'): #change to family name
                                        meropsDict = MEROPS2dict(os.path.join(parentdir, 'DB', 'merops_formatted.fa'))
                                        hit = meropsDict.get(hit)
                                    if not hit in dict:
                                        dict[hit] = [ID]
                                    else:
                                        dict[hit].append(ID)
    return dict
    
def getSMBackbones(input):
    dict = {'NRPS': 0, 'PKS': 0, 'Hybrid': 0}
    with open(input, 'rU') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                if f.type == 'CDS':
                    product = f.qualifiers['product'][0]
                    if not product == 'hypothetical protein':
                        ID = f.qualifiers['locus_tag'][0]
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
            with open(input, 'rU') as gbk:
                SeqRecords = SeqIO.parse(gbk, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        if f.type == 'CDS':
                            ID = f.qualifiers['locus_tag'][0]
                            GOS = []
                            for k,v in f.qualifiers.items():
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
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == 'CDS':
                    ID = f.qualifiers['locus_tag'][0]
                    for k,v in f.qualifiers.items():
                        if k == 'db_xref':
                            for i in v:
                                if i.startswith(word+':'):
                                    hit = i.replace(word+':', '')
                                    if not hit in dict:
                                        dict[hit] = [ID]
                                    else:
                                        dict[hit].append(ID)
    return dict


def convert2counts(input):
    import pandas as pd
    Counts = []
    for i in range(0,len(input)):
        dict = {}
        for k,v in input[i].items():
            dict[k] = len(v)
        Counts.append(dict)
    df = pd.DataFrame(Counts)
    df.fillna(0, inplace=True) #fill in zeros for missing data
    return df

def gb2proteinortho(input, folder, name):
    history = []
    gffOut = os.path.join(folder, name+'.gff')
    FastaOut = os.path.join(folder, name+'.faa')
    Transcripts = os.path.join(folder, name+'.transcripts.fa')
    with open(gffOut, 'w') as gff:
        with open(FastaOut, 'w') as fasta:
            with open(Transcripts, 'w') as transcripts:
                with open(input, 'rU') as input:
                    SeqRecords = SeqIO.parse(input, 'genbank')
                    for record in SeqRecords:
                        for f in record.features:
                            if f.type == "mRNA":
                                feature_seq = f.extract(record.seq)
                                transcripts.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))
                            if f.type == 'CDS':
                                locusID = f.qualifiers['locus_tag'][0]
                                try:  #saw in a genome downloaded from Genbank that a few models don't have protID?  
                                    protID = f.qualifiers['protein_id'][0]
                                except KeyError:
                                    protID = 'ncbi:'+locusID+'-T1'                       
                                start = f.location.nofuzzy_start
                                end = f.location.nofuzzy_end
                                strand = f.location.strand
                                if strand == 1:
                                    strand = '+'
                                elif strand == -1:
                                    strand = '-'
                                translation = f.qualifiers['translation'][0].rstrip('*')
                                product = f.qualifiers['product'][0]
                                chr = record.id
                                if '.' in chr:
                                    chr = chr.split('.')[0]
                                if not protID in history:
                                    history.append(protID)
                                    gff.write("%s\tNCBI\tCDS\t%s\t%s\t.\t%s\t.\tID=%s;Alias=%s;Product=%s;\n" % (chr, start, end, strand, locusID, protID, product))
                                    fasta.write(">%s\n%s\n" % (locusID, translation))

def drawStackedBar(panda, type, labels, ymax, output):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from stackedBarGraph import StackedBarGrapher as StackedBarGrapher
    #stackedbargraph from summary data
    SBG = StackedBarGrapher()
    #labels
    d_labels = panda.index.values
    #y-ticks
    ticks = np.linspace(0,ymax,6)
    ticks = list(ticks)
    nums = [ int(x) for x in ticks ]
    vals = [ str(x) for x in nums ]
    yticks = [nums,vals]
    #colors
    color_palette = sns.hls_palette(len(panda.columns), l=.4, s=.8).as_hex()
    color_palette = [ str(x).upper() for x in color_palette ]
    #set up plot
    sns.set_style('darkgrid')
    sns.set_context('paper')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    YLabel = "Number of "+type
    SBG.stackedBarPlot(ax,panda,color_palette,xLabels=panda.index.values,endGaps=True,gap=0.25,xlabel="Genomes",ylabel=YLabel,yTicks=yticks) 
    plt.title(type+" summary")
    #get the legend
    legends = [] 
    i = 0 
    for column in panda.columns: 
        legends.append(mpatches.Patch(color=color_palette[i], label=panda.columns.values[i]+ ": " + labels.get(panda.columns.values[i]))) 
        i+=1 
    lgd = ax.legend(handles=legends, fontsize=6, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    plt.ylim([0,ymax]) 
    #set the font size - i wish I knew how to do this proportionately.....but setting to something reasonable.
    for item in ax.get_xticklabels():
        item.set_fontsize(8)
    #setup the plot
    fig.subplots_adjust(bottom=0.4)
    fig.savefig(output, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig) 

def drawHeatmap(df, color, output, labelsize, annotate):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    #get size of table
    width = len(df.columns) / 2
    height = len(df.index) / 4
    fig, ax = plt.subplots(figsize=(width,height))
    cbar_ax = fig.add_axes(shrink=0.4)
    if annotate:
        sns.heatmap(df,linewidths=0.5, cmap=color, ax=ax, fmt="d", annot_kws={"size": 4}, annot=True)
    else:
        sns.heatmap(df,linewidths=0.5, cmap=color, ax=ax, annot=False)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    for item in ax.get_xticklabels():
        item.set_fontsize(8)
    for item in ax.get_yticklabels():
        item.set_fontsize(int(labelsize))
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)

def drawbarplot(df, output):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib.pyplot as plt
        import seaborn as sns
    #num = len(df.columns) + 1
    sns.set(style="darkgrid")
    fig = plt.figure()
    #colors
    colorplot = sns.husl_palette(len(df), l=.5).as_hex()
    #colorplot = sns.hls_palette(len(df), l=.4, s=.8).as_hex()
    colorplot = [ str(x).upper() for x in colorplot ]
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
    #run distance metric on matrix and then plot using NMDS
    num = len(df.index)
    data = np.array(df).astype(int)
    bc_dm = pairwise_distances(data, metric=distance)
    mds = MDS(n_components=2, metric=False, max_iter=999, dissimilarity='precomputed', n_init=10, verbose=0)
    result = mds.fit(bc_dm)
    coords = result.embedding_
    stress = 'stress=' + '{0:.4f}'.format(result.stress_)
    #get axis information and make square plus some padding
    xcoords = abs(maxabs(coords[:,0])) + 0.1
    ycoords = abs(maxabs(coords[:,1])) + 0.1
    #setup plot
    fig = plt.figure()
    #colors
    colorplot = sns.husl_palette(len(df), l=.5).as_hex()
    colorplot = [ str(x).upper() for x in colorplot ]
    for i in range(0,num):
        plt.plot(coords[i,0], coords[i,1], 'o', markersize=14, color=colorplot[i], label=df.index.values[i])
    plt.xlabel('NMDS axis 1')
    plt.ylabel('NMDS axis 2')
    plt.ylim(-ycoords,ycoords)
    plt.xlim(-xcoords,xcoords)
    '''
    if num < 13: #if number too large, don't plot
    '''
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.title('NMDS analysis of '+type+' domains')
    plt.annotate(stress, xy=(1,0), xycoords='axes fraction', fontsize=12, ha='right', va='bottom')
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)

def singletons(poff, name):
    with open(poff, 'rU') as input:
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
    with open(poff, 'rU') as input:
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
    with open(file, 'rU') as input:
        for line in input:
            if line.startswith('PF'): #just check to be sure
                line = line.replace('\n', '')
                cols = line.split('\t')
                ID = cols[0]
                desc = cols[4]
                pfamDict[ID] = desc
    return pfamDict

def dictFlip(input):
    #flip the list of dictionaries
    outDict = {}
    for x in input:  
        for k,v in natsorted(x.iteritems()):
            for i in v:
                if i in outDict:
                    outDict[i].append(k)
                else:
                    outDict[i] = [k]
    return outDict

def busco_dictFlip(input):
    #flip the list of dictionaries
    output = []
    for x in input:
        outDict = {}
        for k,v in natsorted(x.iteritems()):
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
        for k,v in natsorted(x.iteritems()):
            #lookup description in another dictionary
            if not lookup.get(k) is None:
                result = k+': '+lookup.get(k)
            else:
                result = k+': No description'
            for i in v:
                if i in outDict:
                    outDict[i].append(str(result))
                else:
                    outDict[i] = [str(result)]
    return outDict

def copyDirectory(src, dest):
    import shutil
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)

buscoTree='eukaryota (303)\n\tmetazoa (978)\n\t\tnematoda (982)\n\t\tarthropoda (1066)\n\t\t\tinsecta (1658)\n\t\t\tendopterygota (2442)\n\t\t\thymenoptera (4415)\n\t\t\tdiptera (2799)\n\t\tvertebrata (2586)\n\t\t\tactinopterygii (4584)\n\t\t\ttetrapoda (3950)\n\t\t\taves (4915)\n\t\t\tmammalia (4104)\n\t\teuarchontoglires (6192)\n\t\t\tlaurasiatheria (6253)\n\tfungi (290)\n\t\tdikarya (1312)\n\t\t\tascomycota (1315)\n\t\t\t\tpezizomycotina (3156)\n\t\t\t\t\teurotiomycetes (4046)\n\t\t\t\t\tsordariomycetes (3725)\n\t\t\t\t\tsaccharomycetes (1759)\n\t\t\t\t\t\tsaccharomycetales (1711)\n\t\t\tbasidiomycota (1335)\n\t\tmicrosporidia (518)\n\tembryophyta (1440)\n\tprotists (215)\n\t\talveolata_stramenophiles (234)\n'
        
busco_links = {
'fungiv1': 'http://busco.ezlab.org/v1/files/fungi_buscos.tar.gz',
'fungi': 'http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz',
'microsporidia': 'http://busco.ezlab.org/v2/datasets/microsporidia_odb9.tar.gz',
'dikarya': 'http://busco.ezlab.org/v2/datasets/dikarya_odb9.tar.gz',
'ascomycota': 'http://busco.ezlab.org/v2/datasets/ascomycota_odb9.tar.gz',
'pezizomycotina' :'http://busco.ezlab.org/v2/datasets/pezizomycotina_odb9.tar.gz',
'eurotiomycetes' : 'http://busco.ezlab.org/v2/datasets/eurotiomycetes_odb9.tar.gz',
'sordariomycetes' : 'http://busco.ezlab.org/v2/datasets/sordariomyceta_odb9.tar.gz',
'saccharomycetes' : 'http://busco.ezlab.org/v2/datasets/saccharomyceta_odb9.tar.gz',
'saccharomycetales' : 'http://busco.ezlab.org/v2/datasets/saccharomycetales_odb9.tar.gz',
'basidiomycota' : 'http://busco.ezlab.org/v2/datasets/basidiomycota_odb9.tar.gz',
'eukaryota' : 'http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz',
'protists' : 'http://busco.ezlab.org/v2/datasets/protists_ensembl.tar.gz',
'alveolata_stramenophiles' : 'http://busco.ezlab.org/v2/datasets/alveolata_stramenophiles_ensembl.tar.gz',
'metazoa' : 'http://busco.ezlab.org/v2/datasets/metazoa_odb9.tar.gz',
'nematoda' : 'http://busco.ezlab.org/v2/datasets/nematoda_odb9.tar.gz',
'arthropoda' : 'http://busco.ezlab.org/v2/datasets/arthropoda_odb9.tar.gz',
'insecta' : 'http://busco.ezlab.org/v2/datasets/insecta_odb9.tar.gz',
'endopterygota' : 'http://busco.ezlab.org/v2/datasets/endopterygota_odb9.tar.gz',
'hymenoptera' : 'http://busco.ezlab.org/v2/datasets/hymenoptera_odb9.tar.gz',
'diptera' : 'http://busco.ezlab.org/v2/datasets/diptera_odb9.tar.gz',
'vertebrata' : 'http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz',
'actinopterygii' : 'http://busco.ezlab.org/v2/datasets/actinopterygii_odb9.tar.gz',
'tetrapoda' : 'http://busco.ezlab.org/v2/datasets/tetrapoda_odb9.tar.gz',
'aves' : 'http://busco.ezlab.org/v2/datasets/aves_odb9.tar.gz',
'mammalia' : 'http://busco.ezlab.org/v2/datasets/mammalia_odb9.tar.gz',
'euarchontoglires' : 'http://busco.ezlab.org/v2/datasets/euarchontoglires_odb9.tar.gz',
'lauraiatheria' : 'http://busco.ezlab.org/v2/datasets/laurasiatheria_odb9.tar.gz',
'embryophyta' : 'http://busco.ezlab.org/v2/datasets/embryophyta_odb9.tar.gz'}
   
def download_buscos(name):
    if name in busco_links:
        log.info("Downloading %s busco models" % name)
        address = busco_links.get(name)
        filename = address.split('/')[-1]
        if name == 'fungiv1':
            foldername = 'fungi'
        else:
            foldername = filename.split('.')[0]
        cmd = ['wget', '-c', '--tries=0', '--read-timeout=20', address]
        runSubprocess(cmd, '.', log)
        cmd = ['tar', '-zxf', filename]
        runSubprocess(cmd, '.', log)
        copyDirectory(os.path.abspath(foldername), os.path.join(parentdir, 'DB', name))
        shutil.rmtree(foldername)
        os.remove(filename)
    else:
        log.error("%s not a valid BUSCO database" % name)
        validBusco = list(busco_links.keys())
        log.error("Valid BUSCO DBs: %s" % (', '.join(validBusco)))
        sys.exit(1)
    
def fasta2dict(Fasta):
    answer = dict()
    with open(Fasta, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'fasta')
        for record in SeqRecords:
            if record.id in answer:
                print "WARNING - duplicate key!"
            else:
                answer[record.id] = str(record.seq)
    return answer 

def ortho2phylogeny(folder, df, num, dict, cpus, bootstrap, tmpdir, outgroup, sp_file, name, sc_buscos):
    import random, pylab
    from Bio import Phylo
    from Bio.Phylo.Consensus import get_support
    if outgroup:
        #load species fasta ids into dictionary
        OutGroup = {}
        with open(sp_file, 'rU') as sp:
            for rec in SeqIO.parse(sp, 'fasta'):
                OutGroup[rec.id] = rec.seq
    #single copy orthologs are in a dataframe, count and then randomly select
    num_species = len(df.columns)
    species = df.columns.values
    if len(df) == 0:
        log.error("0 single copy BUSCO orthologs found, skipping phylogeny")
        return
    if len(df) < int(num):
        number = len(df)
        log.info("Found %i single copy BUSCO orthologs, will use all to infer phylogeny" % (len(df)))
        subsampled = df
    else:
        number = int(num) 
        log.info("Found %i single copy BUSCO orthologs, will randomly select %i to infer phylogeny" % (len(df), number))
        subsampled = df.sample(n=number)

    if outgroup: #passed a list to extract from parent script
        busco_list = sc_buscos

    #since you checked for BUSCO id across all previously, loop through first set and print BUSCOs to file
    with open(os.path.join(tmpdir, 'phylogeny.buscos.used.txt'), 'w') as busco_out:                
        with open(os.path.join(tmpdir, 'phylogeny.concat.fa'), 'w') as proteinout:
            if outgroup:
                proteinout.write(">%s\n" % name)
                for y in busco_list:
                    proteinout.write("%s" % (OutGroup.get(y)))
                proteinout.write('\n')
            for i in range(0,num_species):
                proteinout.write(">%s\n" % species[i])
                proteins = fasta2dict(os.path.join(folder, species[i]+'.faa'))
                for row in subsampled[species[i]].iteritems():
                    proteinout.write("%s" % proteins.get(row[1]))
                    busco_out.write("%s\t%s\n" % (dict[i].get(row[1]), row[1]))
                proteinout.write('\n')
    cmd = ['mafft', '--quiet', os.path.join(tmpdir,'phylogeny.concat.fa')]
    runSubprocess2(cmd, '.', log, os.path.join(tmpdir,'phylogeny.mafft.fa'))
    cmd = ['trimal', '-in', os.path.join(tmpdir,'phylogeny.mafft.fa'), '-out', os.path.join(tmpdir, 'phylogeny.trimal.phylip'), '-automated1', '-phylip']
    runSubprocess(cmd, '.', log)
    if int(cpus) == 1:
        if not outgroup:
            cmd = ['raxmlHPC-PTHREADS', '-f', 'a', '-m', 'PROTGAMMAAUTO', '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
        else:
            cmd = ['raxmlHPC-PTHREADS', '-f', 'a', '-m', 'PROTGAMMAAUTO', '-o', name, '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
    else:
        if not outgroup:
            cmd = ['raxmlHPC-PTHREADS', '-T', str(cpus), '-f', 'a', '-m', 'PROTGAMMAAUTO', '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
        else:
            cmd = ['raxmlHPC-PTHREADS', '-T', str(cpus), '-f', 'a', '-m', 'PROTGAMMAAUTO', '-o', name, '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk']
    #run RAxML
    runSubprocess(cmd, tmpdir, log)
    #parse with biopython and draw
    trees = list(Phylo.parse(os.path.join(tmpdir, 'RAxML_bootstrap.nwk'), 'newick'))
    best = Phylo.read(os.path.join(tmpdir,'RAxML_bestTree.nwk'), 'newick')
    support_tree = get_support(best, trees)
    Phylo.draw(support_tree, do_show=False)
    pylab.axis('off')
    pylab.savefig(os.path.join(tmpdir, 'RAxML.phylogeny.pdf'), format='pdf', bbox_inches='tight', dpi=1000) 

def getTrainResults(input):  
    with open(input, 'rU') as train:
        for line in train:
            if line.startswith('nucleotide level'):
                line = line.replace(' ', '')
                values1 = line.split('|') #get [1] and [2]
            if line.startswith('exon level'):
                line = line.replace(' ', '') #get [6] and [7]
                values2 = line.split('|')
            if line.startswith('gene level'):
                line = line.replace(' ', '')
                values3 = line.split('|') #get [6] and [7]
        return (values1[1], values1[2], values2[6], values2[7], values3[6], values3[7])

def trainAugustus(AUGUSTUS_BASE, train_species, trainingset, genome, outdir, cpus, optimize):
    RANDOMSPLIT = os.path.join(AUGUSTUS_BASE, 'scripts', 'randomSplit.pl')
    OPTIMIZE = os.path.join(AUGUSTUS_BASE, 'scripts', 'optimize_augustus.pl')
    NEW_SPECIES = os.path.join(AUGUSTUS_BASE, 'scripts', 'new_species.pl')
    aug_cpus = '--cpus='+str(cpus)
    species = '--species='+train_species
    aug_log = os.path.join(outdir, 'logfiles', 'augustus_training.log')
    trainingdir = 'tmp_opt_'+train_species
    with open(aug_log, 'w') as logfile:
        if not CheckAugustusSpecies(train_species):
            subprocess.call(['perl', NEW_SPECIES, species], stdout = logfile, stderr = logfile)
        #run etraining again to only use best models from EVM for training
        subprocess.call(['etraining', species, trainingset], stderr = logfile, stdout = logfile)
        subprocess.call(['perl', RANDOMSPLIT, trainingset, '200']) #split off 200 models for testing purposes
        if os.path.isfile(os.path.join(outdir, 'predict_misc', 'busco.training.gb.train')):
            with open(os.path.join(outdir, 'predict_misc', 'augustus.initial.training.txt'), 'w') as initialtraining:
                subprocess.call(['augustus', species, trainingset+'.test'], stdout=initialtraining)
            train_results = getTrainResults(os.path.join(outdir, 'predict_misc', 'augustus.initial.training.txt'))
            log.info('Initial training: '+'{0:.2%}'.format(float(train_results[4]))+' genes predicted exactly and '+'{0:.2%}'.format(float(train_results[2]))+' of exons predicted exactly')
            if optimize:
                #now run optimization
                subprocess.call(['perl', OPTIMIZE, species, aug_cpus, '--onlytrain='+trainingset+'.train', trainingset+'.test'], stderr = logfile, stdout = logfile)
                #run etraining again
                subprocess.call(['etraining', species, trainingset], stderr = logfile, stdout = logfile)
                with open(os.path.join(outdir, 'predict_misc', 'augustus.final.training.txt'), 'w') as finaltraining:
                    subprocess.call(['augustus', species, trainingset+'.test'], stdout=finaltraining)
                train_results = getTrainResults(os.path.join(outdir, 'predict_misc', 'augustus.final.training.txt'))
                log.info('Optimized training: '+'{0:.2%}'.format(float(train_results[4]))+' genes predicted exactly and '+'{0:.2%}'.format(float(train_results[2]))+' of exons predicted exactly')
                #clean up tmp folder
                shutil.rmtree(trainingdir)
            else:
                if float(train_results[4]) < 0.50:
                    log.info("Accuracy seems low, you can try to improve by passing the --optimize_augustus option.")
        else:
            log.error("AUGUSTUS training failed, check logfiles")
            sys.exit(1)

def checkgoatools(input):
    with open(input, 'rU') as goatools:
        count = -1
        result = False
        headercount = 0
        for line in goatools:
            count += 1
            if line.startswith('GO\tNS'):
                header = line.replace('\n', '')
                headercount = count
            if line.startswith('GO:'):
                result = True
    return (result, headercount)

def translatemRNA(input, output):
    with open(output, 'w') as outfile:
        with open(input, 'rU') as fasta:
            for rec in SeqIO.parse(fasta, 'fasta'):
                rec.seq = rec.seq.translate()
                SeqIO.write(rec, outfile, 'fasta')

def alignMAFFT(input, output):
    FNULL = open(os.devnull, 'w')
    with open(output, 'w') as outfile:
        subprocess.call(['mafft', '--quiet', input], stderr = FNULL, stdout = outfile)

def align2Codon(alignment, transcripts, output):
    FNULL = open(os.devnull, 'w')
    with open(output, 'w') as outfile:
        subprocess.call(['perl', os.path.join(UTIL,'pal2nal.pl'), alignment, transcripts, '-output', 'fasta'], stderr=FNULL, stdout = outfile)
    if getSize(output) < 1:
    	os.remove(output)
    	log.debug('dNdS Error: pal2nal failed for %s' % alignment)

def counttaxa(input):
    ct = 0
    with open(input, 'rU') as tree:
        line = tree.readline()
        ct = line.count(',')+1
    return ct

def drawPhyMLtree(fasta, tree):
    FNULL = open(os.devnull, 'w')
    fc = countfasta(fasta)
    #need to convert to phylip format
    base = fasta.split('.')[0]
    tmp1 = base+'.draw2tree.phylip'
    subprocess.call(['trimal', '-in', fasta, '-out', tmp1, '-phylip'])
    #draw tree
    subprocess.call(['phyml', '-i', tmp1], stdout = FNULL, stderr = FNULL)
    tmp2 = base+'.draw2tree.phylip_phyml_tree.txt'
    #check that num taxa in tree = input
    tc = counttaxa(tmp2)
    if tc != fc: #something failed...
        log.debug('dNdS Error: phyml tree failed for %s' % fasta)
        #retry
        subprocess.call(['trimal', '-in', fasta, '-out', tmp1, '-phylip'])
        subprocess.call(['phyml', '-i', tmp1], stdout = FNULL, stderr = FNULL)
    #rename and clean
    os.rename(tmp2, tree)
    os.remove(tmp1)
    os.remove(base+'.draw2tree.phylip_phyml_stats.txt')

def simplestTreeEver(fasta, tree):
    with open(tree, 'w') as outfile:
        with open(fasta, 'rU') as input:
            ids = []
            for rec in SeqIO.parse(input, 'fasta'):
                ids.append(rec.id)
            outfile.write('(%s,%s);' % (ids[0], ids[1]))

def rundNdSexhaustive(folder):
    FNULL = open(os.devnull, 'w')
    #setup intermediate files
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
        #Translate to protein space
        translatemRNA(transcripts, prots)
        #align protein sequences
        alignMAFFT(prots, aln)
        #convert to codon alignment
        align2Codon(aln, transcripts, codon)
        if checkannotations(codon):
            if num_seqs > 2:
                #now generate a tree using phyml
                drawPhyMLtree(codon, tree)
            else:
                simplestTreeEver(transcripts, tree)
            #now run codeml through ete3
            etecmd = ['ete3', 'evol', '--alg', os.path.abspath(codon), '-t', os.path.abspath(tree), '--models', 'M0', 'M1', 'M2', 'M7', 'M8', '-o', name, '--clear_all', '--codeml_param', 'cleandata,1']
            with open(log, 'w') as logfile:
                logfile.write('\n%s\n' % ' '.join(etecmd))
                subprocess.call(etecmd, cwd = tmpdir, stdout = logfile, stderr = logfile)
    #clean up
    for file in os.listdir(tmpdir):
        if file.startswith(name+'.'):
            os.rename(os.path.join(tmpdir, file), os.path.join(tmpdir, name, file))            
              

def rundNdSestimate(folder):
    FNULL = open(os.devnull, 'w')
    #setup intermediate files
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
        #Translate to protein space
        translatemRNA(transcripts, prots)
        #align protein sequences
        alignMAFFT(prots, aln)
        #convert to codon alignment
        align2Codon(aln, transcripts, codon)
        if checkannotations(codon):
            if num_seqs > 2:
                #now generate a tree using phyml
                drawPhyMLtree(codon, tree)
            else:
                simplestTreeEver(transcripts, tree)
            #now run codeml through ete3
            etecmd = ['ete3', 'evol', '--alg', os.path.abspath(codon), '-t', os.path.abspath(tree), '--models', 'M0', '-o', name, '--clear_all', '--codeml_param', 'cleandata,1']
            with open(log, 'w') as logfile:
                logfile.write('\n%s\n' % ' '.join(etecmd))
                subprocess.call(etecmd, cwd = tmpdir, stdout = logfile, stderr = logfile)
    #clean up
    for file in os.listdir(tmpdir):
        if file.startswith(name+'.'):
            os.rename(os.path.join(tmpdir, file), os.path.join(tmpdir, name, file))            

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
        #parse logfile to get omega
        dnds = 'NA'
        m1m2p = 'NA'
        m7m8p = 'NA'
        if os.path.isfile(finallog):
            with open(finallog, 'rU') as input:
                for line in input:
                    line = line.replace('\n', '')
                    if line.startswith('                       M7 |                       M8 | '):
                        m7m8p = line.split('|')[-1].lstrip()
                        m7m8p = m7m8p.replace('*','')
                    if line.startswith('                       M7 |                       M2 | '):
                        m1m2p = line.split('|')[-1].lstrip()
                        m1m2p = m1m2p.replace('*','')
                    if line.startswith(' - Model M0'):
                        nextline = next(input)             
                        dnds = nextline.split('tree: ')[1].rstrip() 
        results[x] =  (dnds, m1m2p, m7m8p)
    return results     


def chunkIt(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0
  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg
  return out


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
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="index.html">Funannotate</a>
        </div>
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li><a href="stats.html">Stats</a></li>
            <li><a href="phylogeny.html">Phylogeny</a></li>
            <li><a href="orthologs.html">Orthologs</a></li>
            <li><a href="interpro.html">InterProScan</a></li>
            <li><a href="pfam.html">PFAM</a></li>
            <li><a href="merops.html">Merops</a></li>
            <li><a href="cazy.html">CAZymes</a></li>
            <li><a href="signalp.html">SignalP</a></li>
            <li><a href="tf.html">TFs</a></li>
            <li><a href="secmet.html">SecMet</a></li>
            <li><a href="go.html">GO ontology</a></li>
            <li><a href="citation.html">Citation</a></li>
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
        <a href='phylogeny/RAxML.phylogeny.pdf'><img src="phylogeny/RAxML.phylogeny.pdf" height="500" /></a></div>
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
