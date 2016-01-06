import os, subprocess, logging, sys, argparse, inspect, csv, time, re
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
UTIL = os.path.join(parentdir, 'Util')

class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

def multipleReplace(text, wordDict):
    for key in wordDict:
        text = text.replace(key, wordDict[key])
    return text

def which(name):
    try:
        with open(os.devnull) as devnull:
            if not name == 'tbl2asn':
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
                    #okay, print out annotations for GAG
                    if ID.endswith('-T1'):
                        output.write("%s\tproduct\t%s\n" % (ID,hdescript))
                        geneID = ID.replace('-T1','')
                        output.write("%s\tname\t%s\n" % (geneID,name))
                    else:
                        output.write("%s\tname\t%s\n" % (ID,name))
                        mrnaID = ID + '-T1'
                        output.write("%s\tproduct\t%s\n" % (mrnaID,hdescript))

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
                    output.write("%s\tnote\tMEROPS:%s %s\n" % (ID,sseqid,family))


def runEggNog(file, cpus, evalue, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #run hmmerscan
    HMM = os.path.join(DB, 'fuNOG_4.5.hmm')
    eggnog_out = os.path.join(tmpdir, 'eggnog.txt')
    subprocess.call(['hmmscan', '-o', eggnog_out, '--cpu', str(cpus), '-E', str(evalue), HMM, file], stdout = FNULL, stderr = FNULL)
    #load in annotation dictionary
    EggNog = {}
    with open(os.path.join(DB,'fuNOG.annotations.tsv'), 'rU') as input:
        reader = csv.reader(input, delimiter='\t')
        for line in reader:
            EggNog[line[1]] = line[5]
    #now parse results
    with open(output, 'w') as output:
        with open(eggnog_out, 'rU') as results:
            for qresult in SearchIO.parse(results, "hmmer3-text"):
                query_length = qresult.seq_len
                lower = query_length * 0.50
                upper = query_length * 1.50
                hits = qresult.hits
                num_hits = len(hits)
                if num_hits > 0:
                    for i in range(0,num_hits):
                        if hits[i].domain_exp_num != hits[i].domain_obs_num: #make sure # of domains is correct
                            continue
                        aln_length = 0
                        num_hsps = len(hits[i].hsps)
                        for x in range(0,num_hsps):
                            aln_length += hits[i].hsps[x].aln_span
                        if aln_length < lower or aln_length > upper: #make sure most of the protein aligns to the model
                            continue
                        hit = hits[i].id.split(".")[1]
                        query = hits[i].query_id
                        #look up descriptions in annotation dictionary
                        description = EggNog.get(hit)
                        final_result = hit + ': ' + description
                        output.write("%s\tnote\t%s\n" % (query, final_result))
                        break

def PFAMsearch(input, cpus, evalue, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    #run hmmerscan
    HMM = os.path.join(DB, 'Pfam-A.hmm')
    pfam_out = os.path.join(tmpdir, 'pfam.txt')
    pfam_filtered = os.path.join(tmpdir, 'pfam.filtered.txt')
    subprocess.call(['hmmscan', '--domtblout', pfam_out, '--cpu', str(cpus), '-E', str(evalue), HMM, input], stdout = FNULL, stderr = FNULL)
    #now parse results
    with open(output, 'w') as output:
        with open(pfam_filtered, 'w') as filtered:
            with open(pfam_out, 'rU') as results:
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
                            pfam = hits[i].accession.split('.')[0]
                            hmmLen = hits[i].seq_len
                            hmm_aln = int(hits[i].hsps[0].hit_end) - int(hits[i].hsps[0].hit_start)
                            coverage = hmm_aln / float(hmmLen)
                            if coverage < 0.50: #coverage needs to be at least 50%
                                continue
                            query = hits[i].query_id
                            description = hits[i].description
                            filtered.write("%s\t%s\t%s\t%s\t%f\n" % (query, pfam, description, hit_evalue, coverage))
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
                            output.write("%s\tnote\t%s enzyme from CAZy family %s\n" % (query, descript, hit))

def fCEGMA(input, cpus, evalue, tmpdir, gff, output):
    FNULL = open(os.devnull, 'w')
    #now run hmmsearch against fCEGMA models
    fCEGMA_HMM = os.path.join(DB, 'fCEGMA.hmm')
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
    #os.rename(os.path.join('RepeatModeler', RP_folder, 'consensi.fa.classified'), library)

    #now soft-mask the genome for gene predictors
    log.info("Soft-masking: running RepeatMasker with custom library")
    if not os.path.exists('RepeatMasker'):
        os.makedirs('RepeatMasker')
    subprocess.call(['RepeatMasker', '-lib', library, '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input], stdout=FNULL, stderr=FNULL)
    for file in os.listdir('RepeatMasker'):
        if file.endswith('.masked'):
            os.rename(os.path.join('RepeatMasker', file), os.path.join(tmpdir, output))
        if file.endswith('.out'):
            rm_gff3 = output.split('.softmasked.fa')[0]
            rm_gff3 = rm_gff3 + '.repeatmasked.gff3'
            rm_gff3 = os.path.join(tmpdir, rm_gff3)
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
    if not os.path.exists('genemark'):
        os.makedirs('genemark')
    contigCount = countfasta(input)
    log.info('Loading genome assembly: ' + '{0:,}'.format(contigCount) + ' contigs')
    contigs = os.path.abspath(input)
    log.info("Running GeneMark-ES on assembly")
    log.debug("gmes_petap.pl --ES --fungus --cores %i --sequence %s" % (cpus, contigs))
    subprocess.call(['gmes_petap.pl', '--ES', '--fungus', '--cores', str(cpus), '--sequence', contigs], cwd='genemark', stdout = FNULL, stderr = FNULL)
    os.rename(os.path.join('genemark','output','gmhmm.mod'), os.path.join(tmpdir, 'gmhmm.mod'))
    #convert genemark gtf to gff3 so GAG can interpret it
    gm_gtf = os.path.join('genemark', 'genemark.gtf')
    log.info("Converting GeneMark GTF file to GFF3")
    with open(output, 'w') as gff:
        subprocess.call(['genemark_gtf2gff3', gm_gtf], stdout = gff)

def MemoryCheck():
    from psutil import virtual_memory
    mem = virtual_memory()
    RAM = int(mem.total)
    return round(RAM / 1024000000)

def runtRNAscan(input, tmpdir, output):
    FNULL = open(os.devnull, 'w')
    tRNAout = os.path.join(tmpdir, 'tRNAscan.out')
    subprocess.call(['tRNAscan-SE', '-o', tRNAout, input], stdout = FNULL, stderr = FNULL)
    trna2gff = os.path.join(UTIL, 'trnascan2gff3.pl')
    with open(output, 'w') as output:
        subprocess.call(['perl', trna2gff, tRNAout], stdout = output, stderr = FNULL)

def runEVM(genome, predictions, proteins, weights, tmpdir, cpus, output):
    import multiprocessing, shutil
    from itertools import izip_longest
    #make sure absolute path is specified for all inputs to avoid problems
    genome = os.path.abspath(genome)
    predictions = os.path.abspath(predictions)
    proteins = os.path.abspath(proteins)
    WEIGHTS = os.path.abspath(weights)

    #define some commands used below
    perl = 'perl'
    EVM = os.environ['EVM_HOME']
    Partition = os.path.join(EVM, 'EvmUtils', 'partition_EVM_inputs.pl')
    Commands = os.path.join(EVM, 'EvmUtils', 'write_EVM_commands.pl')
    Execute = os.path.join(EVM, 'EvmUtils', 'execute_EVM_commands.pl')
    Combine = os.path.join(EVM, 'EvmUtils', 'recombine_EVM_partial_outputs.pl')
    Convert = os.path.join(EVM, 'EvmUtils', 'convert_EVM_outputs_to_GFF3.pl')

    def grouper(n, iterable, fillvalue=None):
        "Collect data into fixed-length chunks or blocks"
        # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
        args = [iter(iterable)] * n
        return izip_longest(fillvalue=fillvalue, *args)

    def worker(input):
        logfile = input + '.log'
        with open(logfile, 'w') as output:
            subprocess.call([perl, Execute, input], cwd = tmpdir, stdout = output, stderr = output)

    def safe_run(*args, **kwargs):
        """Call run(), catch exceptions."""
        try: worker(*args, **kwargs)
        except Exception as e:
            print("error: %s run(*%r, **%r)" % (e, args, kwargs))

    FNULL = open(os.devnull, 'w')
    #split partitions
    log.info("Running EVidence Modeler: setting up EVM partitions")
    subprocess.call([perl, Partition, '--genome', genome, '--gene_predictions', predictions, '--protein_alignments', proteins, '--min_intron_length', '10', '--segmentSize', '100000', '--overlapSize', '10000', '--partition_listing', 'partitions_list.out'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)

    #generate commands
    log.info("generating EVM command list")
    with open('commands.list', 'w') as output:
        subprocess.call([perl, Commands, '--genome', genome, '--gene_predictions', predictions, '--protein_alignments', proteins, '--weights', WEIGHTS, '--min_intron_length', '10', '--output_file_name', 'evm.out', '--partitions', 'partitions_list.out'], cwd = tmpdir, stdout = output, stderr = FNULL)

    #count total lines
    log.info("Running EVM commands in parallel, using %i CPUs" % (int(cpus)))
    commands = os.path.join(tmpdir, 'commands.list')
    num_lines = sum(1 for line in open(commands))
    n = int(round(num_lines / cpus))

    with open(commands, 'rU') as f:
        for i, g in enumerate(grouper(n, f, fillvalue=''), 1):
            with open(os.path.join(tmpdir, 'split_{0}'.format(i * n)+'.cmds'), 'w') as fout:
                fout.writelines(g)

    #now launch a process for each split file
    files = ((os.path.join(tmpdir, f))
                for f in os.listdir(tmpdir) if f.endswith('cmds'))
    p = multiprocessing.Pool(cpus)
    rs = p.map_async(safe_run, files)
    p.close()
    while (True):
        if (rs.ready()): break
        #remaining = rs._number_left
        #print "Waiting for", remaining, "tasks to complete..."
        #time.sleep(30)

    #now combine the paritions
    log.info("Combining EVM partitioned outputs")
    subprocess.call([perl, Combine, '--partitions', 'partitions_list.out', '--output_file_name', 'evm.out'], cwd = tmpdir, stdout = FNULL, stderr = FNULL)

    #now convert to GFF3
    log.info("Converting EVM output to GFF3")
    subprocess.call([perl, Convert, '--partitions', 'partitions_list.out', '--output', 'evm.out', '--genome', genome], cwd = tmpdir, stdout = FNULL, stderr = FNULL)

    #now concatenate all GFF3 files together for a genome then
    log.info("Now collecting all EVM results")
    with open(output, 'w') as output:
        for root, dirs, files in os.walk(tmpdir):
            for file in files:
                if file == 'evm.out.gff3':
                    filename = os.path.join(root,file)
                    with open(filename, 'rU') as readfile:
                        shutil.copyfileobj(readfile, output)

    #remove all the folders in this directory, all
    dirs = [os.path.join(tmpdir,o) for o in os.listdir(tmpdir) if os.path.isdir(os.path.join(tmpdir,o))]
    for i in dirs:
        shutil.rmtree(i)
    for file in os.listdir(tmpdir):
        if file.endswith('cmds') or file.endswith('log'):
            os.remove(file)


def RemoveBadModels(proteins, gff, length, repeats, tmpdir, output):
    #first run bedtools to intersect models where 90% of gene overlaps with repeatmasker region
    FNULL = open(os.devnull, 'w')
    repeat_temp = os.path.join(tmpdir, 'genome.repeats.filtered.gff')
    with open(repeat_temp, 'w') as repeat_out:
        subprocess.call(['bedtools', 'intersect', '-f', '0.9', '-v', '-a', gff, '-b', repeats], stdout = repeat_out, stderr = FNULL)
    #now remove those proteins that do not have valid starts, less then certain length, and have internal stops
    remove = []
    with open(proteins, 'rU') as input:
        SeqRecords = SeqIO.parse(input, 'fasta')
        for rec in SeqRecords:
            Seq = str(rec.seq)[:-1]
            if '*' in Seq:
                remove.append(rec.id)
            if not Seq.startswith('M'):
                remove.append(rec.id)
            if len(Seq) < int(length):
                remove.append(rec.id)
            if 'XX' in Seq:
                remove.append(rec.id)
    remove = [w.replace('evm.model.','') for w in remove]
    remove = set(remove)
    remove_match = re.compile(r'\b(?:%s)[\.;]+\b' % '|'.join(remove))
    with open(output, 'w') as output:
        with open(repeat_out, 'rU') as gff:
            for line in gff:
                if not remove_match.search(line):
                    output.write(line)


def runMaker(input, tmpdir, repeats, mod, species, proteins, transcripts, alt, shortname):
    FNULL = open(os.devnull, 'w')
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    subprocess.call(['maker', '-CTL'], cwd = tmpdir, stdout = FNULL, stderr = FNULL) #create
    #edit maker control file
    os.rename(os.path.join(tmpdir,'maker_opts.ctl'), os.path.join(tmpdir, 'maker_opts.ctl.bak'))
    with open(os.path.join(tmpdir,'maker_opts.ctl'), 'w') as output:
        with open(os.path.join(tmpdir,'maker_opts.ctl.bak'), 'rU') as input:
            for line in input:
                if line.startswith('genome='):
                    line.split(' ', 1)
                    newline = line[0] + input + ' ' + line[1]
                    output.write(newline)+'\n'
                elif line.startswith('protein'):
                    line.split(' ', 1)
                    newline = line[0] + proteins + ' ' + line[1]
                    output.write(newline)+'\n'
                    continue
                if alt == True:
                    if line.startswith('altest'):
                        line.split(' ', 1)
                        newline = line[0] + proteins + ' ' + line[1]
                        output.write(newline)+'\n'
                        continue
                elif alt == False:
                    if line.startswith('est'):
                        line.split(' ', 1)
                        newline = line[0] + proteins + ' ' + line[1]
                        output.write(newline)+'\n'
                        continue

