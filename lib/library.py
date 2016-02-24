from __future__ import division
import os, subprocess, logging, sys, argparse, inspect, csv, time, re, shutil, datetime, glob
from natsort import natsorted
import warnings
from Bio import SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO

#get the working directory, so you can move back into DB folder to find the files you need
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
DB = os.path.join(parentdir, 'DB')
UTIL = os.path.join(parentdir, 'util')
GeneMark2GFF = os.path.join(UTIL, 'genemark_gtf2gff3.pl')

pref_colors=["#CF3C57","#65B23A","#6170DD","#D18738","#D542B5",
"#724A63","#60AABA","#5DB07C","#6C5824","#D74B2B","#6B97D6","#893B2E",
"#B68DB7","#564E91","#ACA13C","#3C6171","#436B33","#D84088",
"#D67A77","#9D55C4","#8B336E","#DA77B9","#D850E5","#B188DF"]

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
    
def checkInternet():
    import urllib2
    try:
        response=urllib2.urlopen('http://74.125.224.72/', timeout=1)
        return True
    except urllib2.URLError as err: pass
    return False

def get_parent_dir(directory):
    return os.path.dirname(directory)

def getSize(filename):
    st = os.stat(filename)
    return st.st_size

def multipleReplace(text, wordDict):
    for key in wordDict:
        text = text.replace(key, wordDict[key])
    return text

def which(name):
    try:
        with open(os.devnull) as devnull:
            diff = ['tbl2asn', 'dustmasker', 'proteinortho5.pl', 'mafft']
            if not any(name in x for x in diff):
                subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
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

def update_progress(progress):
    barLength = 30 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\r IPR progress: [{0}] {1:.2f}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


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
    
def runBUSCO(input, cpus, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #run busco in protein mapping mode
    BUSCO = os.path.join(UTIL, 'funannotate-BUSCO.py')
    proteins = input.split('/')[-1]
    subprocess.call([BUSCO, '-in', proteins, '-m', 'ogs', '-l', os.path.join(DB, 'fungi'), '-o', 'busco', '-c', str(cpus), '-f'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)
    #now parse output and write to annotation file
    with open(output, 'w') as out:
        with open(os.path.join(tmpdir, 'run_busco', 'full_table_busco'), 'rU') as busco:
            for line in busco:
                col = line.split('\t')
                if col[0].startswith('#'):
                    continue
                if col[1] == 'Complete':
                    if col[2].endswith('-T1'):
                        ID = col[2]
                    else:
                        ID = col[2]+'-T1'
                    out.write("%s\tnote\tBUSCO:%s\n" % (ID, col[0]))   

def SwissProtBlast(input, cpus, evalue, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #run blastp against uniprot
    blast_tmp = os.path.join(tmpdir, 'uniprot.xml')
    blastdb = os.path.join(DB,'uniprot')
    subprocess.call(['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus), '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input], stdout = FNULL, stderr = FNULL)
    #parse results
    with open(output, 'w') as output:
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
                        output.write("%s\tprot_desc\t%s\n" % (ID,final_desc))
                        geneID = ID.replace('-T1','')
                    else:
                        mrnaID = ID + '-T1'
                        output.write("%s\tprot_desc\t%s\n" % (mrnaID,final_desc))


def MEROPSBlast(input, cpus, evalue, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #run blastp against merops
    blast_tmp = os.path.join(tmpdir, 'merops.xml')
    blastdb = os.path.join(DB,'MEROPS')
    subprocess.call(['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus), '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input], stdout = FNULL, stderr = FNULL)
    #parse results
    with open(output, 'w') as output:
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
                    output.write("%s\tnote\tMEROPS:%s\n" % (ID,sseqid))

def eggnog2dict():
    #load in annotation dictionary
    EggNog = {}
    with open(os.path.join(DB,'FuNOG.annotations.tsv'), 'rU') as input:
        reader = csv.reader(input, delimiter='\t')
        for line in reader:
            EggNog[line[1]] = line[5]
    return EggNog

def runEggNog(file, cpus, evalue, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #kind of hacky, but hmmersearch doesn't allow me to get sequence length from hmmer3-text, only domtbl, but then I can't get other values, so read seqlength into dictionary for lookup later.
    SeqLength = {}
    with open(file, 'rU') as proteins:
        SeqRecords = SeqIO.parse(proteins, 'fasta')
        for rec in SeqRecords:
            length = len(rec.seq)
            SeqLength[rec.id] = length
    #run hmmerscan
    HMM = os.path.join(DB, 'fuNOG_4.5.hmm')
    eggnog_out = os.path.join(tmpdir, 'eggnog.txt')
    subprocess.call(['hmmsearch', '-o', eggnog_out, '--cpu', str(cpus), '-E', str(evalue), HMM, file], stdout = FNULL, stderr = FNULL)
    #now parse results
    Results = {}
    with open(output, 'w') as output:
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
                output.write("%s\tnote\tEggNog:%s\n" % (k, v[0]))


def PFAMsearch(input, cpus, evalue, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #run hmmerscan
    HMM = os.path.join(DB, 'Pfam-A.hmm')
    pfam_out = os.path.join(tmpdir, 'pfam.txt')
    pfam_filtered = os.path.join(tmpdir, 'pfam.filtered.txt')
    subprocess.call(['hmmsearch', '--domtblout', pfam_out, '--cpu', str(cpus), '-E', str(evalue), HMM, input], stdout = FNULL, stderr = FNULL)
    #now parse results
    with open(output, 'w') as output:
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
                            output.write("%s\tdb_xref\tPFAM:%s\n" % (query, pfam))



def dbCANsearch(input, cpus, evalue, tmpdir, output):
    CAZY = {'CBM': 'Carbohydrate-binding module', 'CE': 'Carbohydrate esterase','GH': 'Glycoside hydrolase', 'GT': 'Glycosyltransferase', 'PL': 'Polysaccharide lyase', 'AA': 'Auxillary activities'}
    FNULL = open(os.devnull, 'w')
    #run hmmerscan
    HMM = os.path.join(DB, 'dbCAN.hmm')
    dbCAN_out = os.path.join(tmpdir, 'dbCAN.txt')
    dbCAN_filtered = os.path.join(tmpdir, 'dbCAN.filtered.txt')
    subprocess.call(['hmmscan', '--domtblout', dbCAN_out, '--cpu', str(cpus), '-E', str(evalue), HMM, input], stdout = FNULL, stderr = FNULL)
    #now parse results
    with open(output, 'w') as output:
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
                            output.write("%s\tnote\tCAZy:%s\n" % (query, hit))

def RepeatModelMask(input, cpus, tmpdir, output):
    log.info("Loading sequences and soft-masking genome")
    FNULL = open(os.devnull, 'w')
    input = os.path.abspath(input)
    output = os.path.abspath(output)
    #lets run RepeatModeler here to get repeat library
    if not os.path.exists('RepeatModeler'):
        os.makedirs('RepeatModeler')
    log.info("Soft-masking: building RepeatModeler database")
    subprocess.call(['BuildDatabase', '-name', tmpdir, input], cwd='RepeatModeler', stdout = FNULL, stderr = FNULL)
    log.info("Soft-masking: generating repeat library using RepeatModeler")
    subprocess.call(['RepeatModeler', '-database', tmpdir, '-pa', str(cpus)], cwd='RepeatModeler', stdout = FNULL, stderr = FNULL)
    #find name of folder
    for i in os.listdir('RepeatModeler'):
        if i.startswith('RM_'):
            RP_folder = i
    library = os.path.join(tmpdir, 'repeatmodeler.lib.fa')
    library = os.path.abspath(library)
    try:
        os.rename(os.path.join('RepeatModeler', RP_folder, 'consensi.fa.classified'), library)
    except OSError:
        pass
    #now soft-mask the genome for gene predictors
    if not os.path.isdir('RepeatMasker'):
            os.makedirs('RepeatMasker')
    if not os.path.isfile(library):
        log.info("Soft-masking: running RepeatMasker with default library (Repeat Modeler found 0 models)")
        subprocess.call(['RepeatMasker', '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input], stdout=FNULL, stderr=FNULL)
    else:
        log.info("Soft-masking: running RepeatMasker with custom library")
        subprocess.call(['RepeatMasker', '-lib', library, '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input], stdout=FNULL, stderr=FNULL)
    for file in os.listdir('RepeatMasker'):
        if file.endswith('.masked'):
            os.rename(os.path.join('RepeatMasker', file), output)
        if file.endswith('.out'):
            rm_gff3 = os.path.join(tmpdir, 'repeatmasker.gff3')
            with open(rm_gff3, 'w') as output:
                subprocess.call(['rmOutToGFF3.pl', file], cwd='RepeatMasker', stdout = output, stderr = FNULL)

def RepeatMask(input, library, cpus, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #now soft-mask the genome for gene predictors
    log.info("Soft-masking: running RepeatMasker with custom library")
    if not os.path.isdir('RepeatMasker'):
        os.makedirs('RepeatMasker')
    subprocess.call(['RepeatMasker', '-lib', library, '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input], stdout=FNULL, stderr=FNULL)
    for file in os.listdir('RepeatMasker'):
        if file.endswith('.masked'):
            os.rename(os.path.join('RepeatMasker', file), output)
        if file.endswith('.out'):
            rm_gff3 = os.path.join(tmpdir, 'repeatmasker.gff3')
            with open(rm_gff3, 'w') as output:
                subprocess.call(['rmOutToGFF3.pl', file], cwd='RepeatMasker', stdout = output, stderr = FNULL)
    
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
    with open(output, 'w') as output:
        with open(input, 'rU') as input:
            records = list(SeqIO.parse(input, 'fasta'))
            records.sort(cmp=lambda x,y: cmp(len(y),len(x)))
            counter = 1
            for rec in records:
                rec.name = ''
                rec.description = ''
                rec.id = 'scaffold_' + str(counter)
                counter +=1
            SeqIO.write(records, output, 'fasta')

def RunGeneMarkES(input, cpus, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #make directory to run script from
    if not os.path.isdir('genemark'):
        os.makedirs('genemark')
    contigs = os.path.abspath(input)
    log.info("Running GeneMark-ES on assembly")
    log.debug("gmes_petap.pl --ES --fungus --cores %i --sequence %s" % (cpus, contigs))
    subprocess.call(['gmes_petap.pl', '--ES', '--fungus', '--soft_mask', '5000', '--cores', str(cpus), '--sequence', contigs], cwd='genemark', stdout = FNULL, stderr = FNULL)
    os.rename(os.path.join('genemark','output','gmhmm.mod'), os.path.join(tmpdir, 'gmhmm.mod'))
    #convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join('genemark', 'genemark.gtf')
    log.info("Converting GeneMark GTF file to GFF3")
    with open(output, 'w') as gff:
        subprocess.call([GeneMark2GFF, gm_gtf], stdout = gff)
        
def RunGeneMark(input, mod, cpus, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #make directory to run script from
    if not os.path.isdir('genemark'):
        os.makedirs('genemark')
    contigs = os.path.abspath(input)
    mod = os.path.abspath(mod)
    log.info("Running GeneMark-ES on assembly")
    log.debug("gmes_petap.pl --ES --ini_mod %s  --fungus --cores %i --sequence %s" % (mod, cpus, contigs))
    subprocess.call(['gmes_petap.pl', '--ES', '--soft_mask', '5000', '--ini_mod', mod, '--fungus', '--cores', str(cpus), '--sequence', contigs], cwd='genemark', stdout = FNULL, stderr = FNULL)
    #convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join('genemark', 'genemark.gtf')
    log.info("Converting GeneMark GTF file to GFF3")
    with open(output, 'w') as gff:
        subprocess.call([GeneMark2GFF, gm_gtf], stdout = gff)

def MemoryCheck():
    from psutil import virtual_memory
    mem = virtual_memory()
    RAM = int(mem.total)
    return round(RAM / 1024000000)

def runtRNAscan(input, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    tRNAout = os.path.join(tmpdir, 'tRNAscan.out')
    if os.path.isfile(tRNAout): # tRNAscan can't overwrite file, so check first
        os.remove(tRNAout)
    subprocess.call(['tRNAscan-SE', '-o', tRNAout, input], stdout = FNULL, stderr = FNULL)
    trna2gff = os.path.join(UTIL, 'trnascan2gff3.pl')
    with open(output, 'w') as output:
        subprocess.call(['perl', trna2gff, '--input', tRNAout], stdout = output, stderr = FNULL)

def gb2smurf(input, prot_out, smurf_out):
    with open(smurf_out, 'w') as smurf:
        with open(prot_out, 'w') as proteins:
            with open(input, 'rU') as gbk:
                SeqRecords = SeqIO.parse(gbk, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        name = re.sub('[^0-9]','', record.name)
                        if f.type == "CDS":
                            proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0]))
                            locus_tag = f.qualifiers.get("locus_tag", ["No ID"])[0]
                            product_name = f.qualifiers.get("product", ["No Description"])[0]
                            mystart = f.location.start
                            myend = f.location.end
                            strand = f.location.strand
                            if strand == 1:
                                smurf.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(mystart), int(myend), product_name))
                            else:
                                smurf.write("%s\t%s\t%s\t%s\t%s\n" % (locus_tag, name.lstrip("0"), int(myend), int(mystart), product_name))
                            

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
                                proteins.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], f.qualifiers['translation'][0]))
                            if f.type == "mRNA":
                                feature_seq = f.extract(record.seq)
                                transcripts.write(">%s\n%s\n" % (f.qualifiers['locus_tag'][0], feature_seq))

def RemoveBadModels(proteins, gff, length, repeats, tmpdir, Output):
    #first run bedtools to intersect models where 90% of gene overlaps with repeatmasker region
    FNULL = open(os.devnull, 'w')
    repeat_temp = os.path.join(tmpdir, 'genome.repeats.to.remove.gff')
    with open(repeat_temp, 'w') as repeat_out:
        subprocess.call(['bedtools', 'intersect', '-f', '0.9', '-a', gff, '-b', repeats], stdout = repeat_out, stderr = FNULL)
    #now remove those proteins that do not have valid starts, less then certain length, and have internal stops
    remove = []
    #parse the results from bedtools and add to remove list
    with open(repeat_temp, 'rU') as input:
        for line in input:
            if "\tgene\t" in line:
                ninth = line.split('ID=')[-1]
                ID = ninth.split(";")[0]
                remove.append(ID)
        
    #I'm only seeing these models with GAG protein translations, so maybe that is a problem? skip for now
    with open(proteins, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)[:-1]
            if not Seq.startswith('M'):
                remove.append(rec.id)
            if len(Seq) < int(length):
                remove.append(rec.id)
            if 'XX' in Seq:
                remove.append(rec.id)
    
    remove = [w.replace('evm.TU.','') for w in remove]
    remove = [w.replace('evm.model.','') for w in remove]
    remove = set(remove)
    remove_match = re.compile(r'\b(?:%s)[\.;]+\b' % '|'.join(remove))
    with open(Output, 'w') as output:
        with open(os.path.join(tmpdir, 'bad_models.gff'), 'w') as output2:
            with open(gff, 'rU') as GFF:
                for line in GFF:
                    if '\tstart_codon\t' in line:
                        continue
                    if '\tstop_codon\t' in line:
                        continue
                    if not remove_match.search(line):
                        line = re.sub(';Name=.*$', ';', line) #remove the Name attribute as it sticks around in GBK file
                        output.write(line)           
                    else:
                        output2.write(line)

def CleantRNAtbl(GFF, TBL, Output):
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
    with open(Output, 'w') as output:
        with open(TBL, 'rU') as input:
            for line in input:
                if line.startswith('\t\t\tlocus_tag\t'):
                    output.write(line)
                    geneID = line.split('locus_tag\t')[-1].replace('\n', '')
                    if geneID in TRNA:
                        if 'tRNA-Xxx' == TRNA.get(geneID):
                            output.write("\t\t\tpseudo\n")       
                elif line.startswith("\t\t\tproduct\ttRNA-Xxx"):
                    output.write(line)
                    output.write("\t\t\tpseudo\n")
                    input.next()
                    input.next()
                elif line.startswith("\t\t\tproduct\ttRNA"):
                    output.write(line)
                    input.next()
                    input.next()
                else:
                    output.write(line)

def ParseErrorReport(input, Errsummary, val, Discrep, output):
    errors = []
    gapErrors = []
    remove = []
    with open(Errsummary) as summary:
        for line in summary:
            if 'ERROR' in line:
                if 'SEQ_DESCR.OrganismIsUndefinedSpecies' in line: #there are probably other errors you are unaware of....
                    pass
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
    if len(errors) < 1: #there are no errors, then just remove stop/start codons and move on
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
    with open(annotations, 'w') as output:
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
            output.write("%s\tproduct\t%s\n" % (ID, hit))          
        #add annots from smProducts
        for k, v in smProducts.items():
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            output.write("%s\tproduct\t%s\n" % (ID, v))               
        #add smCOGs into note section
        for k, v in SMCOGs.items():
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            output.write("%s\tnote\t%s\n" % (ID, v))
              
def GetClusterGenes(input, GFF, Output, annotations):
    global dictClusters
    #pull out genes in clusters from GFF3, load into dictionary
    with open(Output, 'w') as output:
        subprocess.call(['bedtools', 'intersect','-wo', '-a', input, '-b', GFF], stdout = output)
    dictClusters = {}
    with open(Output, 'rU') as input:
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
    with open(annotations, 'w') as output: 
        for k, v in dictClusters.items():
            for i in v:
                if not i.endswith('-T1'):
                    ID = i + ('-T1')
                else:
                    ID = i
                output.write("%s\tnote\tantiSMASH:%s\n" % (ID, k))

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
    with open(input, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            lengths.append(len(record.seq))
            GeeCee.append(str(record.seq))
            for f in record.features:
                if f.type == "source":
                    organism = f.qualifiers.get("organism", ["???"])[0]
                    isolate = f.qualifiers.get("isolate", ["???"])[0]
                if f.type == "CDS":
                    Prots += 1
                if f.type == "gene":
                    Genes += 1
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
    return [organism, isolate, "{0:,}".format(GenomeSize)+' bp', "{0:,}".format(LargestContig)+' bp', "{0:,}".format(AvgContig)+' bp', "{0:,}".format(ContigNum), "{0:,}".format(N50)+' bp', "{:.2f}".format(pctGC)+'%', "{0:,}".format(Genes), "{0:,}".format(Prots), "{0:,}".format(tRNA)]

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
    with open(gffOut, 'w') as gff:
        with open(FastaOut, 'w') as fasta:
            with open(input, 'rU') as input:
                SeqRecords = SeqIO.parse(input, 'genbank')
                for record in SeqRecords:
                    for f in record.features:
                        if f.type == 'CDS':
                            protID = f.qualifiers['protein_id'][0]
                            locusID = f.qualifiers['locus_tag'][0]
                            start = f.location.nofuzzy_start
                            end = f.location.nofuzzy_end
                            strand = f.location.strand
                            if strand == 1:
                                strand = '+'
                            elif strand == -1:
                                strand = '-'
                            translation = f.qualifiers['translation'][0]
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
    YLabel = "Number of "+type+" families"
    SBG.stackedBarPlot(ax,panda,color_palette,xLabels=panda.index.values,endGaps=True,gap=0.25,xlabel="Genomes",ylabel=YLabel,yTicks=yticks) 
    plt.title(type+" family summary")
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

def drawHeatmap(df, color, output, annotate):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import matplotlib.pyplot as plt
        import seaborn as sns
    #get size of table
    width = len(df.columns) / 2
    height = len(df.index) / 4
    fig, ax = plt.subplots(figsize=(width,height))
    cbar_ax = fig.add_axes(shrink=0.4)
    if annotate:
        sns.heatmap(df,linewidths=0.5, cmap=color, ax=ax, annot=True)
    else:
        sns.heatmap(df,linewidths=0.5, cmap=color, ax=ax, annot=False)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    for item in ax.get_xticklabels():
        item.set_fontsize(8)
    for item in ax.get_yticklabels():
        item.set_fontsize(4)
    fig.savefig(output, format='pdf', dpi=1000, bbox_inches='tight')
    plt.close(fig)
 
def distance2mds(df, distance, type, output):
    import numpy as np
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        from sklearn.metrics.pairwise import pairwise_distances
        from sklearn.manifold import MDS
        import seaborn as sns
        import matplotlib.pyplot as plt
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
    if num < 13:
        for i in range(0,num):
            plt.plot(coords[i,0], coords[i,1], 'o', markersize=14, color=pref_colors[i], label=df.index.values[i])
    else:
        for i in range(0,num):
            plt.plot(coords[i,0], coords[i,1], 'o', markersize=14, color='Blue', label=df.index.values[i])    
    plt.xlabel('NMDS axis 1')
    plt.ylabel('NMDS axis 2')
    plt.ylim(-ycoords,ycoords)
    plt.xlim(-xcoords,xcoords)
    if num < 13: #if number too large, don't plot
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
                header = line.replace('.faa', '')
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
                header = line.replace('.faa', '')
                species = header.split('\t')[3:]
                i = species.index(name.replace(' ', '_')) + 3
                continue
            col = line.split('\t')
            if col[0] != '1' and col[i] != '*':
                count += 1
        return count

def iprxml2dict(xmlfile):
    from xml.dom import minidom
    iprDict = {}
    xmldoc = minidom.parse(xmlfile)
    iprlist = xmldoc.getElementsByTagName('interpro')
    for i in iprlist:
        ID = i.attributes['id'].value
        desc = i.getElementsByTagName('name')[0]
        description = desc.firstChild.data
        iprDict[ID] = description
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

def dictFlipLookup(input, lookup):
    outDict = {}
    for x in input:
        for k,v in natsorted(x.iteritems()):
            #lookup description in another dictionary
            result = k+': '+lookup.get(k)
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

def ortho2phylogeny(poff, num, dict, cpus, bootstrap, tmpdir):
    import random, pylab
    from Bio import Phylo
    from Bio.Phylo.Consensus import get_support
    FNULL = open(os.devnull, 'w')
    #get folder name for poff
    folder = os.path.dirname(poff) 
    with open(poff, 'rU') as input:
        count = 0
        sco = {}
        rando = {}
        for line in input:
            if line.startswith('#'):
                header = line.replace('\n', '')
                species = header.split('\t')[3:]
                num_species = line.count('\t') - 2
            line = line.replace('\n', '')
            col = line.split('\t')
            if col[0] == str(num_species) and col[1] == str(num_species) and col[2] == '1':
                count += 1
                ID = 'ortho'+str(count)
                prots = col[3:]
                busco_check = []
                for i in prots:
                    if i in dict:
                        busco_check.append(dict.get(i))
                if len(prots) == len(busco_check): #check that all hits are buscos
                    if len(set(busco_check)) == 1: #check that all busco hits are identical
                        sco[ID] = prots #finally write to dictionary if all checks match
        if len(sco) < int(num):
            num = len(sco)
        else:
            num = int(num)  
        for key in random.sample(sco.keys(), num):
            rando[key] = sco.get(key)
        final = []
        for k,v in sorted(rando.items()):
            final.append(v)
        test = [list(x) for x in zip(*final)] #transpose list
        #since you checked for BUSCO id across all previously, loop through first set and print BUSCOs to file
        with open(os.path.join(tmpdir, 'phylogeny.buscos.used.txt'), 'w') as busco_out:                
            with open(os.path.join(tmpdir, 'phylogeny.concat.fa'), 'w') as proteinout:
                for i in range(0,num_species):
                    proteinout.write(">%s\n" % species[i].split('.faa')[0])
                    proteins = fasta2dict(os.path.join(folder, species[i]))
                    for x in test[i]:
                        proteinout.write("%s" % proteins.get(x))
                        busco_out.write("%s\t%s\n" % (dict.get(x), x))
                    proteinout.write('\n')
        
        with open(os.path.join(tmpdir,'phylogeny.mafft.fa'), 'w') as output:
            subprocess.call(['mafft', os.path.join(tmpdir,'phylogeny.concat.fa')], stdout = output, stderr = FNULL)
        subprocess.call(['trimal', '-in', os.path.join(tmpdir,'phylogeny.mafft.fa'), '-out', os.path.join(tmpdir, 'phylogeny.trimal.phylip'), '-automated1', '-phylip'], stderr = FNULL, stdout = FNULL)
        if int(cpus) == 1:
            subprocess.call(['raxmlHPC-PTHREADS', '-f', 'a', '-m', 'PROTGAMMAAUTO', '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk'], cwd = tmpdir)
        else:
            subprocess.call(['raxmlHPC-PTHREADS', '-T', str(cpus), '-f', 'a', '-m', 'PROTGAMMAAUTO', '-p', '12345', '-x', '12345', '-#', str(bootstrap), '-s', 'phylogeny.trimal.phylip', '-n', 'nwk'], cwd = tmpdir)
    
        #parse with biopython and draw
        trees = list(Phylo.parse(os.path.join(tmpdir, 'RAxML_bootstrap.nwk'), 'newick'))
        best = Phylo.read(os.path.join(tmpdir,'RAxML_bestTree.nwk'), 'newick')
        support_tree = get_support(best, trees)
        Phylo.draw(support_tree, do_show=False)
        pylab.axis('off')
        pylab.savefig(os.path.join(tmpdir, 'RAxML.phylogeny.pdf'), format='pdf', bbox_inches='tight', dpi=1000) 

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
          </button>
          <a class="navbar-brand" href="index.html">Funannotate</a>
        </div>
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li class="active"><a href="stats.html">Stats</a></li>
            <li><a href="phylogeny.html">Phylogeny</a></li>
            <li><a href="orthologs.html">Orthologs</a></li>
            <li><a href="interpro.html">InterProScan</a></li>
            <li><a href="pfam.html">PFAM</a></li>
            <li><a href="merops.html">Merops</a></li>
            <li><a href="cazy.html">CAZymes</a></li>
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
         <p><a href='interpro.html'>InterProScan Domain Stats</a></p>
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
            <li><a href="go.html">GO ontology</a></li>
            <li><a href="citation.html">Citation</a></li>
            <li class="dropdown">
          <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Genomes <span class="caret"></span></a>
          <ul class="dropdown-menu">
'''