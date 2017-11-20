#!/usr/bin/env python
from __future__ import division

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time, csv, glob, re
from natsort import natsorted
import warnings
from Bio import SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO
#import funannotate library
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

IPR2ANNOTATE = os.path.join(parentdir, 'util', 'iprscan2annotations.py')

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-functional.py', usage="%(prog)s [options] -i folder --eggnog emapper.annotations --iprscan proteins.xml --cpus 12",
    description='''Script that adds functional annotation to a genome.''',
    epilog="""Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', help='Folder from funannotate predict.')
parser.add_argument('--genbank', help='Annotated genome in GenBank format')
parser.add_argument('--fasta', help='Genome in FASTA format')
parser.add_argument('--proteins', help='Proteins in FASTA format')
parser.add_argument('--transcripts', help='Transcripts in FASTA format')
parser.add_argument('--gff', help='GFF3 annotation file')
parser.add_argument('-o','--out', help='Basename of output files')
parser.add_argument('--sbt', default='SBT', help='Basename of output files')
parser.add_argument('-s','--species', help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('-t','--tbl2asn', help='Custom parameters for tbl2asn, example: linkage and gap info')
parser.add_argument('-a','--annotations', help='Custom annotations, tsv 3 column file')
parser.add_argument('--isolate', help='Isolate name (e.g. Af293)')
parser.add_argument('--strain', help='Strain name (e.g. CEA10)')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to use')
parser.add_argument('--iprscan', help='IPR5 XML file or folder of pre-computed InterProScan results')
parser.add_argument('--antismash', help='antiSMASH results in genbank format')
parser.add_argument('--force', action='store_true', help='Over-write output folder')
parser.add_argument('--AUGUSTUS_CONFIG_PATH', help='Path to Augustus config directory, $AUGUSTUS_CONFIG_PATH')
parser.add_argument('--phobius', help='Phobius results')
parser.add_argument('--eggnog', help='EggNog Mapper annotations')
parser.add_argument('--busco_db', default='dikarya', help='BUSCO model database')
parser.add_argument('--p2g', help='NCBI p2g file from previous annotation')
parser.add_argument('-d','--database', help='Path to funannotate database, $FUNANNOTATE_DB')
args=parser.parse_args()

#functions
def PfamHmmer(input):
    HMM = os.path.join(FUNDB, 'Pfam-A.hmm')
    base = os.path.basename(input).split('.fa')[0]
    pfam_out = os.path.join(os.path.dirname(input), base+'.pfam.txt')
    cmd = ['hmmsearch', '--domtblout', pfam_out, '--cpu', '1', '-E', '1e-50', HMM, input]
    lib.runSubprocess3(cmd, '.', lib.log)

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: PfamHmmer(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))
        
def combineHmmerOutputs(inputList, output):
    #function to combine multiple HMMER runs with proper header/footer so biopython can read
    allHeadFoot = []
    with open(inputList[0], 'rU') as infile:
        for line in infile:
            if line.startswith('#'):
                allHeadFoot.append(line)
    with open(output, 'w') as out:
        for x in allHeadFoot[:3]:
            out.write(x)
        for file in inputList:
            with open(file, 'rU') as resultin:
                for line in resultin:
                    if line.startswith('#') or line.startswith('\n'): 
                        continue
                    out.write(line)
        for y in allHeadFoot[3:]:
            out.write(y)
 
def multiPFAMsearch(inputList, cpus, evalue, tmpdir, output):
    #run hmmerscan multithreaded by running at same time
    #input is a list of files, run multiprocessing on them
    pfam_results = os.path.join(os.path.dirname(tmpdir), 'pfam.txt')
    pfam_filtered = os.path.join(os.path.dirname(tmpdir), 'pfam.filtered.txt')
    lib.runMultiNoProgress(safe_run, inputList, cpus)
    
    #now grab results and combine, kind of tricky as there are header and footers for each
    resultList = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if os.path.isfile(os.path.join(tmpdir, f)) and f.endswith('.pfam.txt')]
    combineHmmerOutputs(resultList, pfam_results)

    #now parse results
    with open(output, 'w') as out:
        with open(pfam_filtered, 'w') as filtered:
            with open(pfam_results, 'rU') as results:
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

def dbCANHmmer(input):
    HMM = os.path.join(FUNDB, 'dbCAN.hmm')
    base = os.path.basename(input).split('.fa')[0]
    outfiles = os.path.join(os.path.dirname(input), base+'.dbcan.txt')
    cmd = ['hmmscan', '--domtblout', outfiles, '--cpu', '1', '-E', '1e-17', HMM, input]
    lib.runSubprocess3(cmd, '.', lib.log)

def safe_run2(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: dbCANHmmer(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))

def dbCANsearch(inputList, cpus, evalue, tmpdir, output):
    CAZY = {'CBM': 'Carbohydrate-binding module', 'CE': 'Carbohydrate esterase','GH': 'Glycoside hydrolase', 'GT': 'Glycosyltransferase', 'PL': 'Polysaccharide lyase', 'AA': 'Auxillary activities'}
    #run hmmerscan
    dbCAN_out = os.path.join(tmpdir, 'dbCAN.txt')
    dbCAN_filtered = os.path.join(tmpdir, 'dbCAN.filtered.txt')
    lib.runMultiNoProgress(safe_run2, inputList, cpus)
    #cmd = ['hmmscan', '--domtblout', dbCAN_out, '--cpu', str(cpus), '-E', str(evalue), HMM, input]
    #lib.runSubprocess3(cmd, '.', lib.log)
    
    #now grab results
    resultList = [os.path.join(tmpdir, f) for f in os.listdir(tmpdir) if os.path.isfile(os.path.join(tmpdir, f)) and f.endswith('.dbcan.txt')]
    combineHmmerOutputs(resultList, dbCAN_out)
    
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

def MEROPSBlast(input, cpus, evalue, tmpdir, output, diamond=True):
    #run blastp against merops
    blast_tmp = os.path.join(tmpdir, 'merops.xml')
    if diamond:
        blastdb = os.path.join(FUNDB,'merops.dmnd')
        cmd = ['diamond', 'blastp', '--sensitive', '--query', input, '--threads', str(cpus), '--out', blast_tmp, '--db', blastdb, '--evalue', str(evalue), '--max-target-seqs', '1', '--outfmt', '5']
    else:
        blastdb = os.path.join(FUNDB,'MEROPS')
        cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus), '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input]
    lib.runSubprocess(cmd, '.', lib.log)
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

def SwissProtBlast(input, cpus, evalue, tmpdir, GeneDict, diamond=True):
    #run blastp against uniprot
    blast_tmp = os.path.join(tmpdir, 'uniprot.xml')
    if diamond:
        blastdb = os.path.join(FUNDB,'uniprot.dmnd')
        cmd = ['diamond', 'blastp', '--sensitive', '--query', input, '--threads', str(cpus), '--out', blast_tmp, '--db', blastdb, '--evalue', str(evalue), '--max-target-seqs', '1', '--outfmt', '5']
    else:
        blastdb = os.path.join(FUNDB, 'uniprot')
        cmd = ['blastp', '-db', blastdb, '-outfmt', '5', '-out', blast_tmp, '-num_threads', str(cpus), '-max_target_seqs', '1', '-evalue', str(evalue), '-query', input]
    if not lib.checkannotations(blast_tmp):
        lib.runSubprocess(cmd, '.', lib.log)
    #parse results
    counter = 0
    total = 0
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
                name = name.replace('-', '')
                passname = None
                if not '_' in name and not ' ' in name and not '.' in name and number_present(name) and len(name) > 2 and not morethanXnumbers(name, 3):
                    passname = name
                #need to do some filtering here of certain words
                bad_words = ['(Fragment)', 'homolog', 'homolog,', 'AltName:']
                descript = hdescript.split(' ') #turn string into array, splitting on spaces
                final_desc = [x for x in descript if x not in bad_words]
                final_desc = ' '.join(final_desc)
                total += 1
                #add to GeneDict
                if passname:
                    counter += 1
                    if not ID in GeneDict:
                        GeneDict[ID] = [{'name': passname, 'product': final_desc}]
                    else:
                        GeneDict[ID].append({'name': passname, 'product': final_desc})
    lib.log.info('{:,} valid gene/product annotations from {:,} total'.format(counter, total))
    
    
def number_present(s):
    return any(i.isdigit() for i in s)

def morethanXnumbers(s, num):
    count = 0
    for i in s:
        if number_present(i):
            count += 1
    if count >= num:
        return True
    else:
        return False
    
def capfirst(x):
    return x[0].upper() + x[1:]
    
def item2index(inputList, item):
    #return the index of an item in the input list
    item_index = None
    for x in inputList:
        if item in x:
            item_index = inputList.index(x)
    return item_index

def getEggNogHeaders(input):
    IDi, DBi, OGi, Genei, COGi, Desci = (None,)*6
    with open(input, 'rU') as infile:
        for line in infile:
            line = line.replace('\n', '')
            if line.startswith('#query_name'): #this is HEADER
                headerCols = line.split('\t')
                IDi = item2index(headerCols, 'query_name')
                Genei = item2index(headerCols, 'predicted_gene_name')
                DBi = item2index(headerCols, 'Annotation_tax_scope')
                OGi = item2index(headerCols, 'OGs')
                COGi = item2index(headerCols, 'COG cat')
                Desci = item2index(headerCols, 'eggNOG annot')
                break
    return IDi, DBi, OGi, Genei, COGi, Desci
    
def parseEggNoggMapper(input, output, GeneDict):
    Definitions = {}
    #indexes from header file
    IDi, DBi, OGi, Genei, COGi, Desci = getEggNogHeaders(input)
    #take annotations file from eggnog-mapper and create annotations
    with open(output, 'w') as out:
        with open(input, 'rU') as infile:
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
                        NOG = 'ENOG41'+ x.split('@')[0]
                Gene = ''
                if cols[Genei] != '':
                    if not '_' in cols[Genei] and not '.' in cols[Genei] and number_present(cols[Genei]) and len(cols[Genei]) > 2 and not morethanXnumbers(cols[Genei], 3):
                        Gene = cols[Genei]
                Description = cols[Desci].split('. ')[0]
                if not ID.endswith('-T1'):
                    ID = ID+'-T1'
                if NOG == '':
                    continue
                if not NOG in Definitions:
                    Definitions[NOG] = Description
                out.write("%s\tnote\tEggNog:%s\n" % (ID, NOG))
                if cols[COGi] != '':
                    out.write("%s\tnote\tCOG:%s\n" % (ID, cols[COGi].replace(' ','')))
                if Gene != '':
                    product = Gene.lower()+'p'
                    product = capfirst(product)                  
                    #out.write("%s\tname\t%s\n" % (ID.split('-T1')[0], Gene))
                    #out.write("%s\tproduct\t%s\n" % (ID, product))
                    #if Description != '':
                    #    out.write("%s\tnote\t%s\n" % (ID, Description))
                    GeneID = ID.split('-T1')[0]
                    if not GeneID in GeneDict:
                        GeneDict[GeneID] = [{'name': Gene, 'product': Description}]
                    else:
                        GeneDict[GeneID].append({'name': Gene, 'product': Description})                     
    return Definitions


#start here rest of script
#create log file
log_name = 'funannotate-functional.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.SystemInfo()

#get version of funannotate
version = lib.get_version()
lib.log.info("Running %s" % version)

#check dependencies
if args.antismash:
    programs = ['hmmscan', 'hmmsearch', 'gag.py', 'diamond', 'bedtools']
else:
    programs = ['hmmscan', 'hmmsearch', 'gag.py', 'diamond']
lib.CheckDependencies(programs)

#setup funannotate DB path
if args.database:
    FUNDB = args.database
else:
    try:
        FUNDB = os.environ["FUNANNOTATE_DB"]
    except KeyError:
        lib.log.error('Funannotate database not properly configured, run funannotate setup.')
        sys.exit(1)

#check database sources, so no problems later
sources = [os.path.join(FUNDB, 'Pfam-A.hmm.h3p'), os.path.join(FUNDB, 'dbCAN.hmm.h3p'), os.path.join(FUNDB,'merops.dmnd'), os.path.join(FUNDB,'uniprot.dmnd')]
if not all([os.path.isfile(f) for f in sources]):
    lib.log.error('Database files not found in %s, run funannotate database and/or funannotate setup' % FUNDB)
    sys.exit(1)

#write versions of Databases used to logfile
versDB = {}
if not lib.checkannotations(os.path.join(FUNDB, 'funannotate-db-info.txt')):
	lib.log.error('Database not properly configured, run funannotate database and/or funannotate setup')
	sys.exit(1)
with open(os.path.join(FUNDB, 'funannotate-db-info.txt'), 'rU') as dbfile:
	for line in dbfile:
		line = line.strip()
		name, type, file, version, date, num_records = line.split('\t')
		versDB[name] = version

#check Augustus config path as BUSCO needs it to validate species to use
if args.AUGUSTUS_CONFIG_PATH:
    AUGUSTUS = args.AUGUSTUS_CONFIG_PATH
else:
    try:
        AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
    except KeyError:
        lib.log.error("$AUGUSTUS_CONFIG_PATH variable not found. You can use the --AUGUSTUS_CONFIG_PATH argument to specify a path at runtime.")
        sys.exit(1)
        
if not os.path.isdir(os.path.join(AUGUSTUS, 'species')):
    lib.log.error("Augustus species folder not found at %s, exiting" % (os.path.join(AUGUSTUS, 'species')))
    sys.exit(1)

#take care of some preliminary checks
if args.sbt == 'SBT':
    SBT = os.path.join(parentdir, 'lib', 'test.sbt')
    lib.log.info("No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'")
else:
    SBT = args.sbt
    
#check other input files
if not os.path.isfile(SBT):
    lib.log.error("SBT file not found, exiting")
    sys.exit(1)
if args.antismash:
    if not os.path.isfile(args.antismash):
        lib.log.error("Antismash GBK file not found, exiting")
        sys.exit(1)

#check buscos, download if necessary
if not os.path.isdir(os.path.join(FUNDB, args.busco_db)):
    lib.download_buscos(args.busco_db, FUNDB)

#need to do some checks here of the input
genbank, Scaffolds, Protein, Transcripts, GFF = (None,)*5

if not args.input:
    #did not parse folder of funannotate results, so need either gb + gff or fasta + proteins, + gff and also need to have args.out for output folder
    if not args.out:
        lib.log.error("If you are not providing funannotate predict input folder, then you need to provide an output folder (--out)")
        sys.exit(1)
    else:
        outputdir = args.out
        #create outputdir and subdirs
        if not os.path.isdir(outputdir):
            os.makedirs(outputdir)
            os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            os.makedirs(os.path.join(outputdir, 'annotate_results'))
            os.makedirs(os.path.join(outputdir, 'logfiles')) 
    if not args.genbank:
        if not args.fasta or not args.proteins or not args.gff:
            lib.log.error("You did not specifiy the apropriate input files, either: \n1) GenBank \n2) Genome FASTA + Protein FASTA + GFF3")
            sys.exit(1)
        else:
            Scaffolds = args.fasta
            Proteins = args.proteins
            Transcripts = args.transcripts
            GFF = args.gff
    else:
        #create output directories
        if not os.path.isdir(outputdir):
            os.makedirs(outputdir)
            os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            os.makedirs(os.path.join(outputdir, 'annotate_results'))
            os.makedirs(os.path.join(outputdir, 'logfiles'))
        else:
            lib.log.error("Output directory %s already exists, will use any existing data.  If this is not what you want, exit, and provide a unique name for output folder" % (outputdir))
            if not os.path.isdir(os.path.join(outputdir, 'annotate_misc')):
                os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            if not os.path.isdir(os.path.join(outputdir, 'annotate_results')):
                os.makedirs(os.path.join(outputdir, 'annotate_results'))
            if not os.path.isdir(os.path.join(outputdir, 'logfiles')):
                os.makedirs(os.path.join(outputdir, 'logfiles'))                
        genbank = args.genbank
        Scaffolds = os.path.join(outputdir, 'annotate_misc', 'genome.scaffolds.fasta')
        Proteins = os.path.join(outputdir, 'annotate_misc', 'genome.proteins.fasta')
        Transcripts = os.path.join(outputdir, 'annotate_misc', 'genome.transcripts.fasta')
        GFF = os.path.join(outputdir, 'annotate_misc', 'genome.gff3')
        lib.log.info("Checking GenBank file for annotation")
        if not lib.checkGenBank(genbank):
            lib.log.error("Found no annotation in GenBank file, exiting")
            sys.exit(1)
        lib.gb2allout(genbank, GFF, Proteins, Transcripts, Scaffolds)
    
else:
    #should be a folder, with funannotate files, thus store results there, no need to create output folder
    if not os.path.isdir(args.input):
        lib.log.error("%s directory does not exist" % args.input)
        sys.exit(1)
    if os.path.isdir(os.path.join(args.input, 'update_results')): #funannotate results 1) in update folder or 2) in predict folder
        inputdir = os.path.join(args.input, 'update_results')
        outputdir = args.input
    elif os.path.isdir(os.path.join(args.input, 'predict_results')):
        inputdir = os.path.join(args.input, 'predict_results')
        outputdir = args.input
    else:
        inputdir = os.path.join(args.input) #here user specified the predict_results folder, or it is a custom folder

    #get files that you need
    for file in os.listdir(inputdir):
        if file.endswith('.gbk'):
            genbank = os.path.join(inputdir, file)
        if file.endswith('.gff3'):
            GFF = os.path.join(inputdir, file)
    
    #now create the files from genbank input file for consistency in gene naming, etc
    if not genbank or not GFF:
        lib.log.error("Properly formatted 'funannotate predict' files do no exist in this directory")
        sys.exit(1)
    else:
        if 'predict_results' in inputdir or 'update_results' in inputdir: #if user gave predict_results folder, then set output to up one directory
            outputdir = lib.get_parent_dir(inputdir)
        else:
            if not args.out:
                outputdir = inputdir #output the results in the input directory
            else:
                outputdir = args.out
                if not os.path.isdir(outputdir):
                    os.makedirs(outputdir)
        #create output directories
        if not os.path.isdir(os.path.join(outputdir, 'annotate_misc')):
            os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            os.makedirs(os.path.join(outputdir, 'annotate_results'))
        else:
            lib.log.error("Output directory %s already exists, will use any existing data.  If this is not what you want, exit, and provide a unique name for output folder" % (outputdir))
        lib.log.info("Parsing input files")
        Scaffolds = os.path.join(outputdir, 'annotate_misc', 'genome.scaffolds.fasta')
        Proteins = os.path.join(outputdir, 'annotate_misc','genome.proteins.fasta')
        Transcripts = os.path.join(outputdir, 'annotate_misc', 'genome.transcripts.fasta')
        lib.gb2output(genbank, Proteins, Transcripts, Scaffolds)

#make sure logfiles directory is present, will need later
if not os.path.isdir(os.path.join(outputdir, 'logfiles')):
    os.makedirs(os.path.join(outputdir, 'logfiles'))

#get absolute path for all input so there are no problems later, not using Transcripts yet could be error? so take out here
Scaffolds, Proteins, GFF = [os.path.abspath(i) for i in [Scaffolds, Proteins, GFF]] #suggestion via GitHub

#get organism and isolate from GBK file
organism, strain, isolate, accession, WGS_accession, gb_gi, version = (None,)*7
if genbank:
    organism, strain, isolate, accession, WGS_accession, gb_gi, version = lib.getGBKinfo(genbank)
    #since can't find a way to propage the WGS_accession, writing to a file and then parse here
    if os.path.isfile(os.path.join(outputdir, 'update_results', 'WGS_accession.txt')):
        with open(os.path.join(outputdir, 'update_results', 'WGS_accession.txt'), 'rU') as infile:
            for line in infile:
                line = line.replace('\n', '')
                WGS_accession = line

#if command line species/strain/isolate passed, over-write detected 
#check if organism/species/isolate passed at command line, if so, overwrite what you detected.
if args.species:
    organism = args.species
if args.strain:
    strain = args.strain
if args.isolate:
    isolate = args.isolate
if not organism:
    lib.log.error("No GenBank species and no species name given will cause problems downstream, please pass a name to -s,--species")
    sys.exit(1)
if strain:
    organism_name = organism+'_'+strain
elif isolate:
    organism_name = organism+'_'+isolate
else:
    organism_name = organism
organism_name = organism_name.replace(' ', '_')
    
lib.log.info("Adding Functional Annotation to %s, NCBI accession: %s" % (organism, WGS_accession))
lib.log.info("Annotation consists of: {:,} gene models".format(lib.countGFFgenes(GFF)))

############################################################################
#start workflow here
ProtCount = lib.countfasta(Proteins)
lib.log.info('{0:,}'.format(ProtCount) + ' protein records loaded')
if ProtCount < 1:
    lib.log.error("There are no gene models in this genbank file")
    sys.exit(1)

#create tmpdir folder and split proteins into X CPUs to run with HMMER3 searches
protDir = os.path.join(outputdir, 'annotate_misc', 'split_prots')
if not os.path.isdir(protDir):
    os.makedirs(protDir)
lib.fasta2chunks(Proteins, args.cpus, os.path.join(outputdir, 'annotate_misc'), 'split_prots')
splitProts = [os.path.join(protDir, f) for f in os.listdir(protDir) if os.path.isfile(os.path.join(protDir, f))]

#run PFAM-A search
lib.log.info("Running HMMer search of PFAM version %s" % versDB.get('pfam'))
pfam_results = os.path.join(outputdir, 'annotate_misc', 'annotations.pfam.txt')
if not lib.checkannotations(pfam_results):
    multiPFAMsearch(splitProts, args.cpus, 1e-50, protDir, pfam_results)
num_annotations = lib.line_count(pfam_results)
lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')

#initiate Gene Name/Product dictionary
GeneProducts = {}

#run SwissProt Blast search
lib.log.info("Running Diamond blastp search of UniProt DB version %s" % versDB.get('uniprot'))
blast_out = os.path.join(outputdir, 'annotate_misc', 'annotations.swissprot.txt')
SwissProtBlast(Proteins, args.cpus, 1e-5, os.path.join(outputdir, 'annotate_misc'), GeneProducts)
#num_annotations = lib.line_count(blast_out)
#lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')

#Check for EggNog annotations, parse if present
eggnog_out = os.path.join(outputdir, 'annotate_misc', 'annotations.eggnog.txt')
eggnog_result = os.path.join(outputdir, 'annotate_misc', 'eggnog.emapper.annotations')
if args.eggnog:
    if os.path.isfile(eggnog_result):
        os.remove(eggnog_result)
    shutil.copyfile(args.eggnog, eggnog_result)
if not lib.checkannotations(eggnog_result):
    if lib.which('emapper.py'): #eggnog installed, so run it
        lib.log.info("Running Eggnog-mapper")
        cmd = ['emapper.py', '-m', 'diamond', '-i', Proteins, '-o', 'eggnog', '--cpu', str(args.cpus)]
        lib.runSubprocess(cmd, os.path.join(outputdir, 'annotate_misc'), lib.log)
    else:
    	lib.log.info("Install eggnog-mapper or use webserver to improve functional annotation: https://github.com/jhcepas/eggnog-mapper")  
if lib.checkannotations(eggnog_result):
    lib.log.info("Parsing EggNog Annotations")
    EggNog = parseEggNoggMapper(eggnog_result, eggnog_out, GeneProducts)
    num_annotations = lib.line_count(eggnog_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' COG and EggNog annotations added')            

else:
    lib.log.error("No Eggnog-mapper results found.")
    EggNog = {}

#combine the results from UniProt and Eggnog to parse Gene names and product descriptions
#load curated list
lib.log.info("Combining UniProt/EggNog gene and product names using Gene2Product version %s" % versDB.get('gene2product'))
CuratedNames = {}
with open(os.path.join(FUNDB, 'ncbi_cleaned_gene_products.txt'), 'rU') as input:
    for line in input:
        line = line.strip()
        if line.startswith('#'):
            continue
        ID, product = line.split('\t')
        if not ID in CuratedNames:
            CuratedNames[ID] = product

GeneSeen = {}
NeedCurating = {}
NotInCurated = {}
for k,v in natsorted(GeneProducts.items()):
    GeneName = None
    GeneProduct = None
    for x in v:
        if x['name'] in CuratedNames:
            GeneProduct = CuratedNames.get(x['name'])
            GeneName = x['name']
        elif x['name'].lower() in CuratedNames:
            GeneProduct = CuratedNames.get(x['name'].lower())
            GeneName = x['name']    
    if not GeneName: #taking first one will default to swissprot if products for both
        GeneName = v[0]['name']
        GeneProduct = v[0]['product']
        if not GeneName in NotInCurated:
            NotInCurated[GeneName] = GeneProduct
    #now attempt to clean the product name
    rep = {'potential': 'putative', 'possible': 'putative', 'probable': 'putative', 'predicted': 'putative', 
           'uncharacterized': 'putative', 'uncharacterised': 'putative', 'homolog': '', 'EC': '', 'COG': '', 
           'inactivated': '', 'related': '', 'family': '', 'gene': 'protein', 'homologue': ''}
    # replace words in dictionary, from https://stackoverflow.com/questions/6116978/python-replace-multiple-strings
    rep = dict((re.escape(k), v) for k, v in rep.iteritems())
    pattern = re.compile("|".join(rep.keys()))
    GeneProduct = pattern.sub(lambda m: rep[re.escape(m.group(0))], GeneProduct)
    if 'By similarity' in GeneProduct or 'Required for' in GeneProduct or 'nvolved in' in GeneProduct or 'protein '+GeneName == GeneProduct: #some eggnog descriptions are paragraphs....
        if not GeneName in NeedCurating:
            NeedCurating[GeneName] = GeneProduct
        GeneProduct = GeneName.lower()+'p'
        GeneProduct = capfirst(GeneProduct)
    #make sure not multiple spaces
    GeneProduct = ' '.join(GeneProduct.split())
    #if gene name in product, convert to lowercase
    if GeneName in GeneProduct:
        GeneProduct = GeneProduct.replace(GeneName, GeneName.lower())
    if not GeneName in GeneSeen:
        GeneSeen[GeneName] = [(k,GeneProduct)]
    else:
        GeneSeen[GeneName].append((k,GeneProduct))

#finally output the annotations
#which genes are duplicates, need to append numbers to those gene names and then finally output annotations
Gene2ProdFinal = {}
with open(os.path.join(outputdir, 'annotate_misc', 'annotations.genes-products.txt'), 'w') as gene_annotations:
    for key,value in natsorted(GeneSeen.items()):
        if len(value) > 1:
            for i in range(0, len(value)):
                gene_annotations.write("%s\tname\t%s_%i\n" % (value[i][0], key, i+1))
                gene_annotations.write("%s-T1\tproduct\t%s\n" % (value[i][0], value[i][1]))
                Gene2ProdFinal[value[i][0]] = (key+'_'+str(i+1), value[i][1])
        else:
            gene_annotations.write("%s\tname\t%s\n" % (value[0][0], key))
            gene_annotations.write("%s-T1\tproduct\t%s\n" % (value[0][0], value[0][1]))
            Gene2ProdFinal[value[0][0]] = (key, value[0][1])      
num_annotations = int(lib.line_count(os.path.join(outputdir, 'annotate_misc', 'annotations.genes-products.txt')) / 2)
lib.log.info('{:,} gene name and product description annotations added'.format(num_annotations))

#run MEROPS Blast search
lib.log.info("Running Diamond blastp search of MEROPS version %s" % versDB.get('merops'))
blast_out = os.path.join(outputdir, 'annotate_misc', 'annotations.merops.txt')
if not lib.checkannotations(blast_out):
    MEROPSBlast(Proteins, args.cpus, 1e-5, os.path.join(outputdir, 'annotate_misc'), blast_out)
num_annotations = lib.line_count(blast_out)
lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')

#run dbCAN search
dbCAN_out = os.path.join(outputdir, 'annotate_misc', 'annotations.dbCAN.txt')
lib.log.info("Annotating CAZYmes using HMMer search of dbCAN version %s" % versDB.get('dbCAN'))
if not lib.checkannotations(dbCAN_out):
    dbCANsearch(splitProts, args.cpus, 1e-17, protDir, dbCAN_out)
num_annotations = lib.line_count(dbCAN_out)
lib.log.info('{:,} annotations added'.format(num_annotations))

#run BUSCO OGS search
busco_out = os.path.join(outputdir, 'annotate_misc', 'annotations.busco.txt')
lib.log.info("Annotating proteins with BUSCO %s models" % args.busco_db)
buscoDB = os.path.join(FUNDB, args.busco_db)
if not lib.checkannotations(busco_out):
    lib.runBUSCO(Proteins, buscoDB, args.cpus, os.path.join(outputdir, 'annotate_misc'), busco_out)
num_annotations = lib.line_count(busco_out)
lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')

#run Phobius if local is installed, otherwise you will have to use funannotate remote
phobius_out = os.path.join(outputdir, 'annotate_misc', 'phobius.results.txt')
phobiusLog = os.path.join(outputdir, 'logfiles', 'phobius.log')
if args.phobius:
    if os.path.isfile(phobius_out):
        os.remove(phobius_out)
    shutil.copyfile(args.phobius, phobius_out)
if not lib.checkannotations(phobius_out):
    if lib.which('phobius.pl'):
        if not lib.checkannotations(phobius_out):
            lib.log.info("Predicting secreted and transmembrane proteins using Phobius")
            subprocess.call([os.path.join(parentdir, 'util', 'phobius-multiproc.py'), '-i', Proteins, '-o', phobius_out, '-l', phobiusLog])        
    else:
        lib.log.info("Skipping phobius predictions, try funannotate remote -m phobius")
else:
    if lib.checkannotations(phobius_out):
        lib.log.info("Found phobius pre-computed results")
        
#run signalP if installed, have to manually install, so test if exists first, then run it if it does, parse results
signalp_out = os.path.join(outputdir, 'annotate_misc', 'signalp.results.txt')
secreted_out = os.path.join(outputdir, 'annotate_misc', 'annotations.secretome.txt')
membrane_out = os.path.join(outputdir, 'annotate_misc', 'annotations.transmembrane.txt')
if lib.which('signalp'):
    lib.log.info("Predicting secreted proteins with SignalP")
    if not lib.checkannotations(signalp_out):
        lib.signalP(Proteins, os.path.join(outputdir, 'annotate_misc'), signalp_out)
    if lib.checkannotations(phobius_out):
        lib.parsePhobiusSignalP(phobius_out, signalp_out, membrane_out, secreted_out)
    else:
        lib.parseSignalP(signalp_out, secreted_out)
else:
    if not args.phobius:
        lib.log.info("Skipping secretome: neither SignalP nor Phobius installed")
    else:
        lib.log.info("SignalP not installed, secretome prediction less accurate using only Phobius")
        lib.parsePhobiusSignalP(phobius_out, False, membrane_out, secreted_out)
if lib.checkannotations(secreted_out):
    num_secreted = lib.line_count(secreted_out)
else:
    num_secreted = 0
if lib.checkannotations(membrane_out):
    num_mem = lib.line_count(membrane_out)
else:
    num_mem = 0
lib.log.info('{0:,}'.format(num_secreted) + ' secretome and '+ '{0:,}'.format(num_mem) + ' transmembane annotations added')

#interproscan
IPRCombined = os.path.join(outputdir, 'annotate_misc', 'iprscan.xml')
IPR_terms = os.path.join(outputdir, 'annotate_misc', 'annotations.iprscan.txt')
if args.iprscan:
    if os.path.isfile(IPRCombined):
        os.remove(IPRCombined)
    shutil.copyfile(args.iprscan, IPRCombined)
if not lib.checkannotations(IPRCombined):
    lib.log.error("InterProScan error, %s is empty, or no XML file passed via --iprscan. Functional annotation will be lacking." % IPRCombined)
else:
    lib.log.info("Parsing InterProScan5 XML file")
    if os.path.isfile(IPR_terms):
        os.remove(IPR_terms)
    cmd = [sys.executable, IPR2ANNOTATE, IPRCombined, IPR_terms]
    lib.runSubprocess(cmd, '.', lib.log)

#check if antiSMASH data is given, if so parse and reformat for annotations and cluster textual output
antismash_input = os.path.join(outputdir, 'annotate_misc', 'antiSMASH.results.gbk')
if args.antismash:
    if os.path.isfile(antismash_input):
        os.remove(antismash_input)
    shutil.copyfile(args.antismash, antismash_input)
if lib.checkannotations(antismash_input): #result found
    AntiSmashFolder = os.path.join(outputdir, 'annotate_misc', 'antismash')
    AntiSmashBed = os.path.join(AntiSmashFolder,'clusters.bed')
    GFF2clusters = os.path.join(AntiSmashFolder,'secmet.clusters.txt')
    AntiSmash_annotations = os.path.join(outputdir, 'annotate_misc', 'annotations.antismash.txt')
    Cluster_annotations = os.path.join(outputdir, 'annotate_misc', 'annotations.antismash.clusters.txt')
    if not os.path.isdir(AntiSmashFolder):
        os.makedirs(AntiSmashFolder)
    lib.ParseAntiSmash(antismash_input, AntiSmashFolder, AntiSmashBed, AntiSmash_annotations) #results in several global dictionaries
    lib.GetClusterGenes(AntiSmashBed, GFF, GFF2clusters, Cluster_annotations) #results in dictClusters dictionary

#if custom annotations passed, parse here
if args.annotations:
    lib.log.info("Parsing custom annotations from %s" % args.annotations)
    shutil.copyfile(args.annotations, os.path.join(outputdir, 'annotate_misc', 'annotations.custom.txt'))
    num_annotations = lib.line_count(os.path.join(outputdir, 'annotate_misc', 'annotations.custom.txt'))
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    
#now bring all annotations together and annotated genome using gag, remove any duplicate annotations
ANNOTS = os.path.join(outputdir, 'annotate_misc', 'all.annotations.txt')
GeneNames = lib.getGeneBasename(Proteins)
total_annotations = 0
filtered_annotations = 0
lines_seen = set()
with open(ANNOTS, 'w') as output:
    for file in os.listdir(os.path.join(outputdir, 'annotate_misc')):
        if file.startswith('annotations'):
            file = os.path.join(outputdir, 'annotate_misc', file)
            with open(file) as input:
                for line in input:
                    if "\tproduct\ttRNA-" in line: #make sure no collisions when downstream tRNA filtering
                        line.replace('\ttRNA-', '\ttRNA ')
                    total_annotations += 1
                    if not line.startswith(tuple(GeneNames)):
                        continue
                    if line not in lines_seen:
                        output.write(line)
                        lines_seen.add(line)
                        filtered_annotations += 1
ANNOTS = os.path.abspath(ANNOTS)
diff_annotations = total_annotations - filtered_annotations
lib.log.info("Found " + '{0:,}'.format(diff_annotations) + " duplicated annotations, adding " + '{0:,}'.format(filtered_annotations) + ' valid annotations')

#launch gag
GAG = os.path.join(outputdir, 'annotate_misc', 'gag')
lib.log.info("Adding annotations to GFF using GAG")
cmd = ['gag.py', '-f', Scaffolds, '-g', GFF, '-a', ANNOTS, '-o', GAG]
lib.runSubprocess(cmd, '.', lib.log)

#fix the tbl file for tRNA genes
lib.log.info("Fixing tRNA annotations in GenBank tbl file")
original = os.path.join(outputdir, 'annotate_misc','gag', 'genome.tbl')
tmp_tbl = os.path.join(outputdir, 'annotate_misc','gag', 'genome.tbl.original')
os.rename(original, tmp_tbl)
lib.CleantRNAtbl(GFF, tmp_tbl, original)

#if this is reannotation, then need to fix tbl file to track gene changes
if WGS_accession:
    os.rename(original, os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.tbl.bak'))
    p2g = {}
    #see if p2g file is present
    p2gfile = None
    if os.path.isfile(os.path.join(outputdir, 'update_results', 'ncbi.p2g')):
        p2gfile = os.path.join(outputdir, 'update_results', 'ncbi.p2g')
    else:
        if args.p2g:
            p2gfile = args.p2g
    if p2gfile:
        with open(p2gfile, 'rU') as input:
            for line in input:
                cols = line.split('\t')
                if not cols[0] in p2g:
                    p2g[cols[0]] = cols[1]
        with open(original, 'w') as outfile:
            with open(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.tbl.bak'), 'rU') as infile:
                for line in infile:
                    line = line.replace('\n', '')
                    if line.startswith('\t\t\tprotein_id') or line.startswith('\t\t\ttranscript_id'):
                        ID = line.rsplit('|',1)[-1].replace('_mrna', '')
                        type = 'prot'
                        if 'transcript_id' in line:
                            type = 'transcript'
                        if not ID in p2g:
                            if type == 'prot':
                                outfile.write('\t\t\tprotein_id\tgnl|%s|%s\n' % (WGS_accession, ID))
                            elif type == 'transcript':
                                outfile.write('\t\t\ttranscript_id\tgnl|%s|%s_mrna\n' % (WGS_accession, ID))
                        else:
                            p2gID = p2g.get(ID)
                            if type == 'prot':
                                outfile.write('\t\t\tprotein_id\tgnl|%s|%s|gb|%s\n' % (WGS_accession, ID, p2gID))
                            elif type == 'transcript':
                                outfile.write('\t\t\ttranscript_id\tgnl|%s|%s_mrna\n' % (WGS_accession, ID))       
                    else:
                        outfile.write('%s\n' % line)  
    else:
        lib.log.error("Detected NCBI reannotation, but couldn't locate p2g file, please pass via --p2g")
        os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.tbl.bak'), original)
        

#launch tbl2asn to create genbank submission files
shutil.copyfile(os.path.join(GAG, 'genome.fasta'), os.path.join(GAG, 'genome.fsa'))
discrep = 'discrepency.report.txt'
lib.log.info("Converting to final Genbank format, good luck!")
if not version:
    annot_version = 1
else:
    annot_version = version
lib.runtbl2asn(GAG, SBT, discrep, organism, args.isolate, args.strain, args.tbl2asn, annot_version)

#parse discrepancy report to see which names/product descriptions failed/passed
BadProducts = lib.getFailedProductNames(discrep, Gene2ProdFinal) #return list of tuples of (GeneName, GeneProduct)
Gene2ProductPassed = os.path.join(outputdir, 'annotate_results', 'Gene2Products.new-names-passed.txt')
with open(Gene2ProductPassed, 'w') as prodpassed:
    prodpassed.write('#Name\tDescription\n')
    for key, value in natsorted(NotInCurated.items()):
        if not key in BadProducts:
            prodpassed.write('%s\t%s\n' % (key, value))
Gene2ProductHelp = os.path.join(outputdir, 'annotate_results', 'Gene2Products.need-curating.txt')
with open(Gene2ProductHelp, 'w') as needhelp:
    needhelp.write('#Name\tDescription\tError-message\n')
    for key, value in natsorted(NeedCurating.items()):
        needhelp.write('%s\t%s\tProduct defline failed funannotate checks\n' % (key, value))
    for key, value in natsorted(BadProducts.items()):
        needhelp.write('%s\t%s\tProduct defline failed tbl2asn checks\n' % (key, value))

#collected output files and rename accordingly
ResultsFolder = os.path.join(outputdir, 'annotate_results')
os.rename(discrep, os.path.join(ResultsFolder, organism_name+'.discrepency.report.txt'))
final_gbk = os.path.join(ResultsFolder, organism_name+'.gbk')
final_proteins = os.path.join(ResultsFolder, organism_name+'.proteins.fa')
final_transcripts = os.path.join(ResultsFolder, organism_name+'.transcripts.fa')
final_fasta = os.path.join(ResultsFolder, organism_name+'.scaffolds.fa')
final_annotation = os.path.join(ResultsFolder, organism_name+'.annotations.txt')
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.gbf'), final_gbk)
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.gff'), os.path.join(ResultsFolder, organism_name+'.gff3'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.tbl'), os.path.join(ResultsFolder, organism_name+'.tbl'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.sqn'), os.path.join(ResultsFolder, organism_name+'.sqn'))
lib.gb2output(final_gbk, final_proteins, final_transcripts, final_fasta)

#write AGP output so all files in correct directory
lib.log.info("Creating AGP file and corresponding contigs file")
agp2fasta = os.path.join(parentdir, 'util', 'fasta2agp.pl')
AGP = os.path.join(ResultsFolder, organism_name+'.agp')
cmd = ['perl', agp2fasta, organism_name+'.scaffolds.fa']
lib.runSubprocess2(cmd, ResultsFolder, lib.log, AGP)

#write secondary metabolite clusters output using the final genome in gbk format
if lib.checkannotations(antismash_input): 
    lib.log.info("Cross referencing SM cluster hits with MIBiG database version %s" % versDB.get('mibig'))
    #do a blast best hit search against MIBiG database for cluster annotation, but looping through gene cluster hits
    AllProts = []
    for k, v in lib.dictClusters.items():
        for i in v:
            if not i in AllProts:
                AllProts.append(i)
    AllProts = set(AllProts)
    mibig_fasta = os.path.join(AntiSmashFolder, 'smcluster.proteins.fasta')
    mibig_blast = os.path.join(AntiSmashFolder, 'smcluster.MIBiG.blast.txt')
    mibig_db = os.path.join(FUNDB, 'mibig.dmnd')
    with open(mibig_fasta, 'w') as output:
        with open(Proteins, 'rU') as input:
            SeqRecords = SeqIO.parse(Proteins, 'fasta')
            for record in SeqRecords:
                if record.id in AllProts:
                    SeqIO.write(record, output, 'fasta')
    cmd = ['diamond', 'blastp', '--sensitive', '--query', mibig_fasta, '--threads', str(args.cpus), '--out', mibig_blast, '--db', mibig_db, '--max-hsps', '1', '--evalue', '0.001', '--max-target-seqs', '1', '--outfmt', '6']
    #cmd = ['blastp', '-query', mibig_fasta, '-db', mibig_db, '-num_threads', str(args.cpus), '-max_target_seqs', '1', '-max_hsps', '1', '-evalue', '0.001', '-outfmt', '6', '-out', mibig_blast]
    lib.runSubprocess(cmd, '.', lib.log)
    #now parse blast results to get {qseqid: hit}
    MIBiGBlast = {}
    with open(mibig_blast, 'rU') as input:
        for line in input:
            cols = line.split('\t')
            ID = cols[0]
            hit = cols[1].split('|')
            desc = hit[5]
            cluster = hit[0]
            db_ref = hit[6]
            evalue = cols[10]
            pident = cols[2]
            result = (desc, cluster, db_ref, pident, evalue)
            MIBiGBlast[ID] = result
            
    lib.log.info("Creating tab-delimited SM cluster output")

    #load in antismash cluster bed file to slice record
    slicing = []
    with open(AntiSmashBed, 'rU') as antibed:
        for line in antibed:
            cols = line.split('\t')
            cluster = (cols[0],cols[3],cols[1],cols[2]) #chr, cluster, start, stop in a tuple
            slicing.append(cluster)
    Offset = {}
    #Get each cluster + 15 Kb in each direction to make sure you can see the context of the cluster
    with open(os.path.join(ResultsFolder, organism_name+'.gbk'), 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == "source":
                    record_start = f.location.start
                    record_end = f.location.end
            for slice in slicing:
                if record.id == slice[0]:
                    sub_start = int(slice[2]) - 15000
                    sub_stop = int(slice[3]) + 15000
                    if sub_start < 1:
                        sub_start = 1
                    if sub_stop > record_end:
                        sub_stop = record_end
                    sub_record = record[sub_start:sub_stop]
                    cluster_name = slice[1]
                    sub_record_name = os.path.join(AntiSmashFolder, cluster_name+'.gbk')
                    Offset[cluster_name] = sub_start
                    with open(sub_record_name, 'w') as clusterout:
                        SeqIO.write(sub_record, clusterout, 'genbank')

    #okay, now loop through each cluster 
    for file in os.listdir(AntiSmashFolder):
        if file.endswith('.gbk'):
            base = file.replace('.gbk', '')
            outputName = os.path.join(AntiSmashFolder, base+'.secmet.cluster.txt')
            file = os.path.join(AntiSmashFolder, file)
            with open(outputName, 'w') as output:
                output.write("#%s\n" % base)
                output.write("#GeneID\tChromosome:start-stop\tStrand\tClusterPred\tBackbone Enzyme\tBackbone Domains\tProduct\tsmCOGs\tEggNog\tInterPro\tPFAM\tGO terms\tNotes\tMIBiG Blast\tProtein Seq\tDNA Seq\n")
                with open(file, 'rU') as input:
                    SeqRecords = SeqIO.parse(input, 'genbank')
                    for record in SeqRecords:
                        for f in record.features:
                            if f.type == "CDS":
                                name = f.qualifiers["locus_tag"][0]
                                prot_seq = f.qualifiers['translation'][0]
                                start = f.location.nofuzzy_start
                                actualStart = int(start) + int(Offset.get(base)) + 1 #account for python numbering shift?
                                end = f.location.nofuzzy_end
                                actualEnd = int(end) + int(Offset.get(base))
                                strand = f.location.strand
                                if strand == 1:
                                    strand = '+'
                                    DNA_seq = record.seq[start:end]
                                elif strand == -1:
                                    strand = '-'
                                    DNA_seq = record.seq[start:end].reverse_complement()
                                chr = record.id
                                product = f.qualifiers["product"][0]
                                #now get the info out of the note and db_xref fields, need to clear each field for each record
                                note = []
                                goTerms = []
                                pFAM = []
                                IPR = []  
                                eggnogDesc = 'NA'
                                location = 'flanking'
                                cog = '.'                  
                                for k,v in f.qualifiers.items():
                                    if k == 'note':
                                        #multiple notes are split with a semi colon
                                        items = v[0].split('; ')
                                        for i in items:
                                            if i.startswith('EggNog:'):
                                                eggnogID = i.replace('EggNog:', '')
                                                eggnogDesc = EggNog.get(eggnogID)
                                            elif i.startswith('GO_'):
                                                goterm = i.split(': ', 1)[-1]
                                                goTerms.append(goterm)
                                            elif i.startswith('SMCOG'):
                                                cog = i
                                            elif i.startswith('antiSMASH:'):
                                                location = 'cluster'
                                            else:
                                                note.append(i)
                                    if k == 'db_xref':
                                        for i in v:
                                            if i.startswith('InterPro:'):
                                                r = i.replace('InterPro:', '')
                                                IPR.append(r)
                                            if i.startswith('PFAM:'):
                                                p = i.replace('PFAM:', '')
                                                pFAM.append(p)
                                if name in lib.bbDomains:
                                    domains = ";".join(lib.bbDomains.get(name))
                                else:
                                    domains = '.'
                                if name in lib.bbSubType:
                                    enzyme = lib.bbSubType.get(name)
                                else:
                                    if name in lib.BackBone:
                                        enzyme = lib.BackBone.get(name)
                                    else:
                                        enzyme = '.'
                                if name in MIBiGBlast:
                                    mibigTup = MIBiGBlast.get(name)
                                    mibig = mibigTup[0]+' from '+mibigTup[1]+' ('+mibigTup[2]+':pident='+mibigTup[3]+', evalue='+mibigTup[4]+')'
                                    mibig = str(mibig)
                                else:
                                    mibig = '.'
                                if IPR:
                                    IP = ";".join(IPR)
                                else:
                                    IP = '.'
                                if pFAM:
                                    PF = ";".join(pFAM)
                                else:
                                    PF = '.'
                                if goTerms:
                                    GO = ";".join(goTerms)
                                else:
                                    GO = '.'            
                                if note:
                                    No = ";".join(note)
                                else:
                                    No = '.'              
                                output.write("%s\t%s:%i-%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, chr, actualStart, actualEnd, strand, location, enzyme, domains, product, cog, eggnogDesc, IP, PF, GO, No, mibig, prot_seq, DNA_seq))
                                                             
    #now put together into a single file
    finallist = []
    ClustersOut = os.path.join(ResultsFolder, organism_name+'.clusters.txt')
    for file in os.listdir(AntiSmashFolder):
        if file.endswith('secmet.cluster.txt'):
            file = os.path.join(AntiSmashFolder, file)
            finallist.append(file)
    with open(ClustersOut, 'w') as output:
        for file in natsorted(finallist):
            with open(file, 'rU') as input:
                output.write(input.read())
                output.write('\n\n')

#write tsv annotation table
lib.log.info("Writing genome annotation table.")
lib.annotationtable(final_gbk, FUNDB, final_annotation)

#final wrap up message
lib.log.info("Funannotate annotate has completed successfully!\n\n\
We need YOUR help to improve gene names/product descriptions:\n\
   {:,} gene/product names did not pass, see {:}\n\
   {:,} gene/product names passed but are not in Database, see {:}\n\n\
Please consider contributing a PR at https://github.com/nextgenusfs/gene2product\n".format(len(NeedCurating)+len(BadProducts),Gene2ProductHelp,len(NotInCurated),Gene2ProductPassed))
print "-------------------------------------------------------"
#move logfile to logfiles directory
if os.path.isfile(log_name):
    if not os.path.isdir(os.path.join(outputdir, 'logfiles')):
        os.makedirs(os.path.join(outputdir, 'logfiles'))
    os.rename(log_name, os.path.join(outputdir, 'logfiles', log_name))
