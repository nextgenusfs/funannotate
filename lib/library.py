import os, subprocess, logging, sys, argparse, inspect, csv, time, re, shutil
import warnings
from BCBio import GFF
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

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'
    
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
            diff = ['tbl2asn', 'dustmasker']
            if not any(name in x for x in diff):
                subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
            else:
                subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def line_count(fname):
    with open(fname) as f:
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
                    #species = description[1].replace(' GN','')
                    name = description[2].replace(' PE','').upper()
                    #hdescript = 'hypothetical protein similar to ' + name
                    #okay, print out annotations for GAG
                    if ID.endswith('-T1'):
                        output.write("%s\tprot_desc\t%s\n" % (ID,hdescript))
                        geneID = ID.replace('-T1','')
                        #output.write("%s\tname\t%s\n" % (geneID,name))
                    else:
                        #output.write("%s\tname\t%s\n" % (ID,name))
                        mrnaID = ID + '-T1'
                        output.write("%s\tprot_desc\t%s\n" % (mrnaID,hdescript))

def MEROPSBlast(input, cpus, evalue, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #run blastp against uniprot
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
    #load in annotation dictionary
    EggNog = {}
    with open(os.path.join(DB,'FuNOG.annotations.tsv'), 'rU') as input:
        reader = csv.reader(input, delimiter='\t')
        for line in reader:
            EggNog[line[1]] = line[5]
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
                description = EggNog.get(v[0])
                final_result = v[0] + ': ' + description
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

def fCEGMA(input, cpus, evalue, tmpdir, gff, output):
    FNULL = open(os.devnull, 'w')
    #now run hmmsearch against fCEGMA models
    fCEGMA_HMM = os.path.join(parentdir, 'DB', 'fCEGMA.hmm')
    temp_out = os.path.join(tmpdir, 'fCEGMA.hmmsearch.txt')
    subprocess.call(['hmmsearch', '-o', temp_out, '-E', str(evalue), '--cpu', str(cpus), fCEGMA_HMM, input], stdout = FNULL, stderr = FNULL)
    #now parse results, getting only high quality matches
    keep = {}
    with open(output, 'w') as output:
        with open(temp_out, 'rU') as results:
            for qresult in SearchIO.parse(results, "hmmer3-text"):
                hits = qresult.hits
                model = qresult.id
                #here we just want the best hit for each model
                if len(hits) > 0:
                    beste = hits[0].evalue
                    if beste >= evalue:
                        continue
                    model_length = qresult.seq_len
                    hit_start = hits[0].hsps[0].hit_start
                    hit_end = hits[0].hsps[0].hit_end
                    hit_aln = hit_end - hit_start
                    coverage = hit_aln / float(model_length)
                    if coverage < 0.9:
                        continue
                    hit = hits[0].id
                    if hit not in keep:
                        keep[hit] = model
                    output.write("%s\t%s\t%s\t%s\t%s\t%s\t%f\n" % (hit, model, beste, model_length, hit_start, hit_end, coverage))
    #loop through genemark GFF3 and pull out genes, rename to 'MODEL', and then slice to just get those that pass.
    for key, value in keep.items():
        keep[key] = value + "-T1"
        new_key = key.replace("_t", "_g")
        keep[new_key] = value
    import re
    pattern = re.compile(r'\b(' + '|'.join(keep.keys()) + r')\b')
    gff_out = os.path.join(tmpdir, 'training.gff3')
    with open(gff_out, 'w') as output:
        with open(gff, 'rU') as input:
            for line in input:
                line = pattern.sub(lambda x: keep[x.group()], line)
                if 'MODEL' in line:
                    output.write(line)

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
    os.rename(os.path.join('RepeatModeler', RP_folder, 'consensi.fa.classified'), library)

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

def CheckDependencies(input):
    missing = []
    for p in input:
        if which(p) == False:
            missing.append(p)
    if missing != []:
        error = ", ".join(missing)
        log.error("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        sys.exit(1)

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
    subprocess.call(['gmes_petap.pl', '--ES', '--fungus', '--cores', str(cpus), '--sequence', contigs], cwd='genemark', stdout = FNULL, stderr = FNULL)
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
    subprocess.call(['gmes_petap.pl', '--ES', '--ini_mod', mod, '--fungus', '--cores', str(cpus), '--sequence', contigs], cwd='genemark', stdout = FNULL, stderr = FNULL)
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
                    geneID = line.split('locus_tag\t')[-1].replace('\n', '')
                    if geneID in TRNA:
                        if 'tRNA-Xxx' == TRNA.get(geneID):
                            output.write(line)
                            output.write("\t\t\tpseudo\n")
                        else:
                            output.write(line)
                    
                if line.startswith("\t\t\tproduct\ttRNA-Xxx"):
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
    global BackBone, SMCOGs, bbSubType, bbDomains, Offset, smProducts
    BackBone = {}; SMCOGs = {}; bbSubType = {}; bbDomains = {}; Offset = {}; smProducts = {}
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
                        #get cluster + 15 kb on each side just to be safe for writing output files later on
                        sub_start = start - 15000
                        sub_end = end + 15000
                        if sub_start < 1:
                            sub_start = 1
                        if sub_end > record_end:
                            sub_end = record_end
                        sub_record = record[sub_start:sub_end]
                        Offset['Cluster_'+clusternum] = sub_start
                        sub_record_name = os.path.join(tmpdir, 'Cluster_'+clusternum+'.gbk')
                        with open(sub_record_name, 'w') as clusterout:
                            SeqIO.write(sub_record, clusterout, 'genbank')
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
            else:
                hit = v
            if hit == 'terpene':
                hit = 'terpene cyclase'
            elif hit == 'other':
                hit = 'putative secondary metabolism biosynthetic enzyme'
            output.write("%s\tproduct\t%s" % (ID, hit))          
        #add annots from smProducts
        for k, v in smProducts.items()
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            output.write("%s\tproduct\t%s" % (ID, v))               
        #add smCOGs into note section
        for k, v in SMCOGs.items():
            if not k.endswith('-T1'):
                ID = k + '-T1'
            else:
                ID = k
            output.write("%s\tnote\t%s" % (ID, v))
              
def GetClusterGenes(input, GFF, Output, annotations):
    global dictClusters
    #pull out genes in clusters from GFF3, load into dictionary
    with open(Output, 'w') as output:
        subprocess.call(['bedtools', 'intersect','-wo', '-a', input, '-b', GFF], stdout = output)
    dictClusters = {}
    with open(GenesInClusters, 'rU') as input:
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
                    ID = i + ('-T1)
                else:
                    ID = i
                output.write("%s\tnote\tantiSMASH:%s" % (ID, k))
            
        

