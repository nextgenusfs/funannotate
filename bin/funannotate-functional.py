#!/usr/bin/env python
from __future__ import division

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time, csv, glob
from Bio import SeqIO
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

RUNIPRSCAN = os.path.join(parentdir, 'util', 'runIPRscan.py')

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-functional.py', usage="%(prog)s [options] -i genome.fasta -g genome.gff -o test -e youremail@mail.edu",
    description='''Script that adds functional annotation to a genome.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', help='Folder from funannotate predict.')
parser.add_argument('--genbank', help='Annotated genome in GenBank format')
parser.add_argument('--fasta', help='Genome in FASTA format')
parser.add_argument('--proteins', help='Proteins in FASTA format')
parser.add_argument('--transcripts', help='Transcripts in FASTA format')
parser.add_argument('--gff', help='GFF3 annotation file')
parser.add_argument('-o','--out', help='Basename of output files')
parser.add_argument('-e','--email', help='Email address for IPRSCAN server')
parser.add_argument('--sbt', default='SBT', help='Basename of output files')
parser.add_argument('-s','--species', help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('--isolate', help='Isolate/strain name (e.g. Af293)')
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use')
parser.add_argument('--iprscan', help='Folder of pre-computed InterProScan results (1 xml per protein)')
parser.add_argument('--antismash', help='antiSMASH results in genbank format')
parser.add_argument('--skip_iprscan', action='store_true', help='skip InterProScan remote query')
parser.add_argument('--force', action='store_true', help='Over-write output folder')
args=parser.parse_args()

def runIPRpython(Input):
    base = Input.split('/')[-1]
    base = base.split('.fa')[0]
    OUTPATH = os.path.join(IPROUT, base)
    subprocess.call([sys.executable, RUNIPRSCAN, '--goterms','--email', args.email, '--outfile', OUTPATH, '--input', Input], stderr=FNULL, stdout=FNULL)
    #now rename output files, just keep xml and tsv files?
    time.sleep(3) #make sure there is time for all files to show up
    os.rename(OUTPATH+'.xml.xml', OUTPATH+'.xml')
    os.rename(OUTPATH+'.tsv.txt', OUTPATH+'.tsv')
    os.rename(OUTPATH+'.svg.svg', OUTPATH+'.svg')
    os.rename(OUTPATH+'.sequence.txt', OUTPATH+'.fa')
    os.remove(OUTPATH+'.gff.txt')
    os.remove(OUTPATH+'.htmltarball.html.tar.gz')
    os.remove(OUTPATH+'.log.txt')
    os.remove(OUTPATH+'.out.txt')

#create log file
log_name = 'funnannotate-functional.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.log.info("Operating system: %s, %i cores, ~ %i GB RAM" % (sys.platform, multiprocessing.cpu_count(), lib.MemoryCheck()))

#check dependencies
if args.antismash:
    programs = ['hmmscan', 'hmmsearch', 'blastp', 'gag.py','bedtools']
else:
    programs = ['hmmscan', 'hmmsearch', 'blastp', 'gag.py']
lib.CheckDependencies(programs)

#take care of some preliminary checks
if args.sbt == 'SBT':
    SBT = os.path.join(parentdir, 'lib', 'test.sbt')
    lib.log.info("No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'")
else:
    SBT = args.sbt

if not args.iprscan and not args.email:
    lib.log.error("To run InterProScan you need to specify an email address to identify yourself to the online service")
    os._exit(1)

#need to do some checks here of the input
if not args.input:
    #did not parse folder of funannotate results, so need either gb + gff or fasta + proteins, + gff and also need to have args.out for output folder
    if not args.out:
        lib.log.error("If you are not providing funannotate predict input folder, then you need to provide an output folder (--out)")
        os._exit(1)
    else:
        outputdir = args.out
    if not args.genbank or not args.gff:
        if not args.fasta or not args.proteins or not args.gff or not args.transcripts:
            lib.log.error("You did not specifiy the apropriate input files, either: \n1) GenBank + GFF3\n2) Genome FASTA + Protein FASTA + Transcript FASTA + GFF3")
            os._exit(1)
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
        else:
            lib.log.error("Output directory %s already exists, will use any existing data.  If this is not what you want, exit, and provide a unique name for output folder" % (outputdir))
            if not os.path.isdir(os.path.join(outputdir, 'annotate_misc')):
                os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            if not os.path.isdir(os.path.join(outputdir, 'annotate_results')):
                os.makedirs(os.path.join(outputdir, 'annotate_results'))
        genbank = args.genbank
        Scaffolds = os.path.join(outputdir, 'annotate_misc', 'genome.scaffolds.fasta')
        Proteins = os.path.join(outputdir, 'annotate_misc', 'genome.proteins.fasta')
        Transcripts = os.path.join(outputdir, 'annotate_misc', 'genome.transcripts.fasta')
        GFF = args.gff
        lib.gb2output(genbank, Proteins, Transcripts, Scaffolds)
    
else:
    #should be a folder, with funannotate files, thus store results there, no need to create output folder
    if not os.path.isdir(args.input):
        lib.log.error("%i directory does not exist" % args.input)
        os._exit(1)
    if os.path.isdir(os.path.join(args.input, 'predict_results')): #funannotate results should be here
        inputdir = os.path.join(args.input, 'predict_results')
        outputdir = args.input
    else:
        inputdir = os.path.join(args.input) #here user specified the predict_results folder, or it is all together wrong, find out in next few lines
    #get files that you need
    for file in os.listdir(inputdir):
        if file.endswith('.gbk'):
            genbank = os.path.join(inputdir, file)
        if file.endswith('.gff3'):
            GFF = os.path.join(inputdir, file)
    
    #now create the files from genbank input file for consistency in gene naming, etc
    if not genbank or not GFF:
        lib.log.error("Properly formatted 'funannotate predict' files do no exist in this directory")
        os._exit(1)
    else:
        if inputdir == outputdir: #I'm assuming means that user supplied the predict folder (probably bad to assume...)
            outputdir = lib.get_parent_dir(inputdir)
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

#get absolute path for all input so there are no problems later
for i in Scaffolds, Proteins, Transcripts, GFF:
    i = os.path.abspath(i)
        

#get organism and isolate from GBK file
if not args.species:
    with open(genbank, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == "source":
                    organism = f.qualifiers.get("organism", ["???"])[0]
                    if not args.isolate:
                        isolate = f.qualifiers.get("isolate", ["???"])[0]
                    else:
                        isolate = args.isolate
                    break
else:
    organism = args.species
    if not args.isolate:
        isolate = '???'
    else:
        isolate = args.isolate

############################################################################
#start workflow here
ProtCount = lib.countfasta(Proteins)
lib.log.info('{0:,}'.format(ProtCount) + ' protein records loaded')  

#run interpro scan
IPROUT = os.path.join(outputdir, 'annotate_misc', 'iprscan')
PROTS = os.path.join(outputdir, 'annotate_misc', 'protein_tmp')
for i in IPROUT,PROTS:
    if not os.path.exists(i):
        os.makedirs(i)
 
#run PFAM-A search
lib.log.info("Running HMMer search of PFAM domains")
pfam_results = os.path.join(outputdir, 'annotate_misc', 'annotations.pfam.txt')
if not os.path.isfile(pfam_results):
    lib.PFAMsearch(Proteins, args.cpus, 1e-50, os.path.join(outputdir, 'annotate_misc'), pfam_results)
num_annotations = lib.line_count(pfam_results)
lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
#run SwissProt Blast search
lib.log.info("Running Blastp search of UniProt DB")
blast_out = os.path.join(outputdir, 'annotate_misc', 'annotations.swissprot.txt')
if not os.path.isfile(blast_out):
    lib.SwissProtBlast(Proteins, args.cpus, 1e-5, os.path.join(outputdir, 'annotate_misc'), blast_out)
num_annotations = lib.line_count(blast_out)
lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
#run MEROPS Blast search
lib.log.info("Running Blastp search of MEROPS protease DB")
blast_out = os.path.join(outputdir, 'annotate_misc', 'annotations.merops.txt')
if not os.path.isfile(blast_out):
    lib.MEROPSBlast(Proteins, args.cpus, 1e-5, os.path.join(outputdir, 'annotate_misc'), blast_out)
num_annotations = lib.line_count(blast_out)
lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
#run dbCAN search
dbCAN_out = os.path.join(outputdir, 'annotate_misc', 'annotations.dbCAN.txt')
lib.log.info("Annotating CAZYmes using dbCAN")
if not os.path.isfile(dbCAN_out):
    lib.dbCANsearch(Proteins, args.cpus, 1e-17, os.path.join(outputdir, 'annotate_misc'), dbCAN_out)
num_annotations = lib.line_count(dbCAN_out)
lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
#run EggNog search
eggnog_out = os.path.join(outputdir, 'annotate_misc', 'annotations.eggnog.txt')
lib.log.info("Annotating proteins with EggNog 4.5 database")
if not os.path.isfile(eggnog_out):
    lib.runEggNog(Proteins, args.cpus, 1e-10, os.path.join(outputdir, 'annotate_misc'), eggnog_out)
num_annotations = lib.line_count(eggnog_out)
lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')

internet = lib.checkInternet()
if not internet:
    lib.log.info("No internet connection detected, skipping InterProScan search")
    
if not args.iprscan or args.skip_iprscan or not internet:
    #now run interproscan
    #split input into individual files
    lib.splitFASTA(Proteins, PROTS)

    #now iterate over list using pool and up to 25 submissions at a time
    proteins = []
    for file in os.listdir(PROTS):
        if file.endswith('.fa'):
            file = os.path.join(PROTS, file)
            proteins.append(file)
        
    num_files = len(glob.glob1(IPROUT,"*.xml"))
    num_prots = len(proteins)
    lib.log.info("Now running InterProScan search remotely using EBI servers on " + '{0:,}'.format(num_prots) + ' proteins')
    while (num_files < num_prots):
        #build in a check before running (in case script gets stopped and needs to restart
        finished = []
        for file in os.listdir(IPROUT):
            if file.endswith('.xml'):
                base = file.split('.xml')[0]
                fasta_file = os.path.join(PROTS, base+'.fa')
                finished.append(fasta_file)

        finished = set(finished)
        runlist = [x for x in proteins if x not in finished]
        #start up the list
        p = multiprocessing.Pool(25) #max searches at a time for IPR server
        rs = p.map_async(runIPRpython, runlist)
        p.close()
        while (True):
            if (rs.ready()): break
            num_files = len(glob.glob1(IPROUT,"*.xml"))
            pct = num_files / num_prots
            lib.update_progress(pct)
            time.sleep(10)
        num_files = len(glob.glob1(IPROUT,"*.xml"))
        pct = num_files / num_prots
        lib.update_progress(pct)
else:
    IPROUT = args.iprscan

if not args.skip_iprscan:
    #now collect the results from InterProscan, then start to reformat results
    lib.log.info("InterProScan has finished, now pulling out annotations from results")
    IPR_terms = os.path.join(outputdir, 'annotate_misc', 'annotations.iprscan.txt')
    if not os.path.isfile(IPR_terms):
        IPR2TSV = os.path.join(parentdir, 'util', 'ipr2tsv.py')
        with open(IPR_terms, 'w') as output:
            subprocess.call([sys.executable, IPR2TSV, IPROUT], stdout = output, stderr = FNULL)
    GO_terms = os.path.join(outputdir, 'annotate_misc', 'annotations.GO.txt')
    if not os.path.isfile(GO_terms):
        IPR2GO = os.path.join(parentdir, 'util', 'ipr2go.py')
        OBO = os.path.join(parentdir, 'DB', 'go.obo')
        with open(GO_terms, 'w') as output:
            subprocess.call([sys.executable, IPR2GO, OBO, IPROUT], stdout = output, stderr = FNULL)
        
#check if antiSMASH data is given, if so parse and reformat for annotations and cluster textual output
if args.antismash:
    AntiSmashFolder = os.path.join(outputdir, 'annotate_misc', 'antismash')
    AntiSmashBed = os.path.join(AntiSmashFolder, 'annotate_misc', 'clusters.bed')
    GFF2clusters = os.path.join(AntiSmashFolder, 'annotate_misc', 'secmet.clusters.txt')
    AntiSmash_annotations = os.path.join(outputdir, 'annotate_misc', 'annotations.antismash.txt')
    Cluster_annotations = os.path.join(outputdir, 'annotate_misc', 'annotations.antismash.clusters.txt')
    if not os.path.isdir(AntiSmashFolder):
        os.makedirs(AntiSmashFolder)
    lib.ParseAntiSmash(args.antismash, AntiSmashFolder, AntiSmashBed, AntiSmash_annotations) #results in several global dictionaries
    lib.GetClusterGenes(AntiSmashBed, GFF, GFF2clusters, Cluster_annotations) #results in dictClusters dictionary
     
#now bring all annotations together and annotated genome using gag
ANNOTS = os.path.join(outputdir, 'annotate_misc', 'all.annotations.txt')
with open(ANNOTS, 'w') as output:
    for file in os.listdir(os.path.join(outputdir, 'annotate_misc')):
        if file.startswith('annotations'):
            file = os.path.join(outputdir, 'annotate_misc', file)
            with open(file) as input:
                output.write(input.read())
ANNOTS = os.path.abspath(ANNOTS)

#launch gag
GAG = os.path.join(outputdir, 'annotate_misc', 'gag')
lib.log.info("Adding annotations to GFF using GAG")
subprocess.call(['gag.py', '-f', Scaffolds, '-g', GFF, '-a', ANNOTS, '-o', GAG], stdout = FNULL, stderr = FNULL)

#fix the tbl file for tRNA genes
lib.log.info("Fixing tRNA annotations in GenBank tbl file")
original = os.path.join(outputdir, 'annotate_misc','gag', 'genome.tbl')
tmp_tbl = os.path.join(outputdir, 'annotate_misc','gag', 'genome.tbl.original')
os.rename(original, tmp_tbl)
lib.CleantRNAtbl(GFF, tmp_tbl, original)

#write to GBK file
if not isolate == '???':
    ORGANISM = "[organism=" + organism + "] " + "[isolate=" + isolate + "]"
    baseOUTPUT = organism + '_' + isolate
else:
    ORGANISM = "[organism=" + organism + "]"
    baseOUTPUT = organism
#remove any spaces from baseoutput 
baseOUTPUT = baseOUTPUT.replace(' ', '_')

#launch tbl2asn to create genbank submission files
shutil.copyfile(os.path.join(GAG, 'genome.fasta'), os.path.join(GAG, 'genome.fsa'))
discrep = 'discrepency.report.txt'
lib.log.info("Converting to final Genbank format, good luck!.....")
subprocess.call(['tbl2asn', '-p', GAG, '-t', SBT, '-M', 'n', '-Z', discrep, '-a', 'r10u', '-l', 'paired-ends', '-j', ORGANISM, '-V', 'b', '-c', 'fx'], stdout = FNULL, stderr = FNULL)

#collected output files and rename accordingly
ResultsFolder = os.path.join(outputdir, 'annotate_results')
os.rename(discrep, os.path.join(ResultsFolder, baseOUTPUT+'.discrepency.report.txt'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.gbf'), os.path.join(ResultsFolder, baseOUTPUT+'.gbk'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.gff'), os.path.join(ResultsFolder, baseOUTPUT+'.gff3'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.tbl'), os.path.join(ResultsFolder, baseOUTPUT+'.tbl'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.sqn'), os.path.join(ResultsFolder, baseOUTPUT+'.sqn'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag', 'genome.fasta'), os.path.join(ResultsFolder, baseOUTPUT+'.scaffolds.fa'))
shutil.rmtree(PROTS)

#write AGP output so all files in correct directory
lib.log.info("Creating AGP file and corresponding contigs file")
agp2fasta = os.path.join(parentdir, 'util', 'fasta2agp.pl')
AGP = os.path.join(ResultsFolder, baseOUTPUT+'.agp')
with open(AGP, 'w') as output:
    subprocess.call(['perl', agp2fasta, baseOUTPUT+'.scaffolds.fa'], cwd = ResultsFolder, stdout = output, stderr = FNULL)


#write secondary metabolite clusters output using the final genome in gbk format
if args.antismash:
    lib.log.info("Cross referencing SM cluster hits with MIBiG database")
    #do a blast best hit search against MIBiG database for cluster annotation, but looping through gene cluster hits
    AllProts = []
    for k, v in lib.dictClusters.items():
        for i in v:
            if not i in AllProts:
                AllProts.append(i)
    AllProts = set(AllProts)
    mibig_fasta = os.path.join(AntiSmashFolder, 'smcluster.proteins.fasta')
    mibig_blast = os.path.join(AntiSmashFolder, 'smcluster.MIBiG.blast.txt')
    mibig_db = os.path.join(parentdir, 'DB', 'MIBiG')
    with open(mibig_fasta, 'w') as output:
        with open(Proteins, 'rU') as input:
            SeqRecords = SeqIO.parse(Proteins, 'fasta')
            for record in SeqRecords:
                if record.id in AllProts:
                    SeqIO.write(record, output, 'fasta')
    subprocess.call(['blastp', '-query', mibig_fasta, '-db', mibig_db, '-num_threads', str(args.cpus), '-max_target_seqs', '1', '-max_hsps', '1', '-evalue', '0.001', '-outfmt', '6', '-out', mibig_blast])
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
    #load in EggNog annotations to get descriptions for table
    EggNog = {}
    with open(os.path.join(parentdir, 'DB','FuNOG.annotations.tsv'), 'rU') as input:
        reader = csv.reader(input, delimiter='\t')
        for line in reader:
            EggNog[line[1]] = line[5]
    #load in antismash cluster bed file to slice record
    slicing = []
    with open(AntiSmashBed, 'rU') as antibed:
        for line in antibed:
            cols = line.split('\t')
            cluster = (cols[0],cols[3],cols[1],cols[2]) #chr, cluster, start, stop in a tuple
            slicing.append(cluster)
    Offset = {}
    #Get each cluster + 15 Kb in each direction to make sure you can see the context of the cluster
    with open(os.path.join(ResultsFolder, baseOUTPUT+'.gbk'), 'rU') as gbk:
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
    ClustersOut = os.path.join(ResultsFolder, baseOUTPUT+'.clusters.txt')
    for file in os.listdir(AntiSmashFolder):
        if file.endswith('secmet.cluster.txt'):
            file = os.path.join(AntiSmashFolder, file)
            finallist.append(file)
    with open(ClustersOut, 'w') as output:
        for file in natsorted(finallist):
            with open(file, 'rU') as input:
                output.write(input.read())
                output.write('\n\n')
#move logfile to logfiles directory
os.rename(log_name, os.path.join(outputdir, 'logfiles', log_name))
os._exit(1)
    

