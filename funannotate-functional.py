#!/usr/bin/env python

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time, csv
from Bio import SeqIO
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib


#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-functional.py', usage="%(prog)s [options] -i genome.fasta -g genome.gff -o test -e youremail@mail.edu",
    description='''Script that adds functional annotation to a genome.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', help='Results folder from funannotate predict.')
parser.add_argument('--genbank', help='Annotated genome in GenBank format')
parser.add_argument('--fasta', help='Genome in FASTA format')
parser.add_argument('--proteins', help='Proteins in FASTA format')
parser.add_argument('--transcripts', help='Transcripts in FASTA format')
parser.add_argument('--gff', help='GFF3 annotation file')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('-e','--email', help='Email address for IPRSCAN server')
parser.add_argument('--sbt', default='SBT', help='Basename of output files')
parser.add_argument('-s','--species', help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('--isolate', help='Isolate/strain name (e.g. Af293)')
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use')
parser.add_argument('--iprscan', help='Folder of pre-computed InterProScan results (1 xml per protein)')
parser.add_argument('--antismash', help='antiSMASH results in genbank format')
parser.add_argument('--force', action='store_true', help='Over-write output folder')
args=parser.parse_args()

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

#get executable for runiprscan
PATH2JAR = os.path.join(currentdir, 'util', 'RunIprScan-1.1.0', 'RunIprScan.jar')
RUNIPRSCAN_PATH = os.path.join(currentdir, 'util', 'RunIprScan-1.1.0')

#check dependencies
if args.antismash:
    programs = ['hmmscan', 'hmmsearch', 'blastp', 'gag.py', 'java', 'bedtools']
else:
    programs = ['hmmscan', 'hmmsearch', 'blastp', 'gag.py', 'java']
lib.CheckDependencies(programs)

#create temp folder to house intermediate files
if not os.path.exists(args.out):
    os.makedirs(args.out)
else:
    lib.log.error("Output directory %s already exists, will use any existing data.  If this is not what you want, exit, and provide a unique name for output folder" % (args.out))

#need to do some checks here of the input
if not args.input:
    #did not parse folder of funannotate results, so need either gb + gff or fasta + proteins, + gff.
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
        Scaffolds = os.path.join(args.out, 'genome.scaffolds.fasta')
        Proteins = os.path.join(args.out, 'genome.proteins.fasta')
        Transcripts = os.path.join(args.out, 'genome.transcripts.fasta')
        GFF = args.gff
        lib.gb2output(genbank, Proteins, Transcripts, Scaffolds)
    
else:
    #should be a folder, with funannotate files
    if not os.path.isdir(args.input):
        lib.log.error("%i directory does not exist" % args.input)
        os._exit(1)
    for file in os.listdir(args.input):
        if file.endswith('.gbk'):
            genbank = os.path.join(args.input, file)
        if file.endswith('.gff3'):
            GFF = os.path.join(args.input, file)
    
    #now create the files from genbank input file for consistency in gene naming, etc
    if not genbank or not GFF:
        lib.log.error("Properly formatted 'funannotate predict' files do no exist in this directory")
        os._exit(1)
    else:
        lib.log.info("Parsing input files")
        Scaffolds = os.path.join(args.out, 'genome.scaffolds.fasta')
        Proteins = os.path.join(args.out, 'genome.proteins.fasta')
        Transcripts = os.path.join(args.out, 'genome.transcripts.fasta')
        lib.gb2output(genbank, Proteins, Transcripts, Scaffolds)

#get absolute path for all input so there are no problems later
for i in Scaffolds, Proteins, Transcripts, GFF:
    i = os.path.abspath(i)

if not args.iprscan and not args.email:
    lib.log.error("In order to run InterProScan you need to supply a valid email address to identify yourself to the online service")
    os._exit(1)
            
#take care of some preliminary checks
IPROUT = os.path.join(args.out, 'iprscan')
if args.sbt == 'SBT':
    SBT = os.path.join(currentdir, 'lib', 'test.sbt')
    lib.log.info("You did not specify an NCBI SBT file, will use default, however you might want to create one and pass it here under the '--sbt' argument")
else:
    SBT = args.sbt

#get organism and isolate from GBK file
if not args.species:
    with open(genbank, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == "source":
                    organism = f.qualifiers.get("organism", ["???"])[0]
                    if not args.isolate:
                        isolate = f.qualifiers.get("strain", ["???"])[0]
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

#run interpro scan, in background hopefully....
if not os.path.exists(os.path.join(args.out, 'iprscan')):
    os.makedirs(os.path.join(args.out, 'iprscan'))

if not args.iprscan: #here run the routine of IPRscan in the background
    lib.log.info("Starting RunIprScan and running in background")
    p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #while RunIprScan is running in background, run more functional annotation methods
    while p.poll() is None:
        #run PFAM-A search
        lib.log.info("Running HMMer search of PFAM domains")
        pfam_results = os.path.join(args.out, 'annotations.pfam.txt')
        if not os.path.isfile(pfam_results):
            lib.PFAMsearch(Proteins, args.cpus, 1e-50, args.out, pfam_results)
        num_annotations = lib.line_count(pfam_results)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:   #run it again to recover any that did not work
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run SwissProt Blast search
        lib.log.info("Running Blastp search of UniProt DB")
        blast_out = os.path.join(args.out, 'annotations.swissprot.txt')
        if not os.path.isfile(blast_out):
            lib.SwissProtBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
        num_annotations = lib.line_count(blast_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run MEROPS Blast search
        lib.log.info("Running Blastp search of MEROPS protease DB")
        blast_out = os.path.join(args.out, 'annotations.merops.txt')
        if not os.path.isfile(blast_out):
            lib.MEROPSBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
        num_annotations = lib.line_count(blast_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run EggNog search
        eggnog_out = os.path.join(args.out, 'annotations.eggnog.txt')
        lib.log.info("Annotating proteins with EggNog 4.5 database")
        if not os.path.isfile(eggnog_out):
            lib.runEggNog(Proteins, args.cpus, 1e-10, args.out, eggnog_out)
        num_annotations = lib.line_count(eggnog_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run dbCAN search
        dbCAN_out = os.path.join(args.out, 'annotations.dbCAN.txt')
        lib.log.info("Annotating CAZYmes using dbCAN")
        if not os.path.isfile(dbCAN_out):
            lib.dbCANsearch(Proteins, args.cpus, 1e-17, args.out, dbCAN_out)
        num_annotations = lib.line_count(dbCAN_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("Waiting for RunIprScan to complete")
            p.wait()
    
    #lets do a loop over the IPRresults and run until it is complete, sometimes this can get stuck and stop downloading results. final check of IPRresults, i.e. the number of proteins should equal number of files in iprscan folder, let run for 1 hour, check again, relaunch, etc.
    num_ipr = len([name for name in os.listdir(IPROUT) if os.path.isfile(os.path.join(IPROUT, name))])
    if num_ipr < ProtCount:
        lib.log.info("Number of IPR xml files (%i) does not equal number of proteins (%i), will run RunIprScan until complete" % (num_ipr, ProtCount))
        p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=FNULL, stderr=FNULL)
        while (num_ipr < ProtCount):
            time.sleep(60)
            num_ipr = len([name for name in os.listdir(IPROUT) if os.path.isfile(os.path.join(IPROUT, name))])      
        lib.log.info("Number of proteins (%i) is less than or equal to number of XML files (%i)" % (ProtCount, num_ipr))
        p.terminate()
        
    else:
        lib.log.info("Number of proteins (%i) is less than or equal to number of XML files (%i)" % (ProtCount, num_ipr))
else:   
    #check that remaining searches have been done, if not, do them.
    #run PFAM-A search
    lib.log.info("Running HMMer search of PFAM domains")
    pfam_results = os.path.join(args.out, 'annotations.pfam.txt')
    if not os.path.isfile(pfam_results):
        lib.PFAMsearch(Proteins, args.cpus, 1e-50, args.out, pfam_results)
    num_annotations = lib.line_count(pfam_results)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    #run SwissProt Blast search
    lib.log.info("Running Blastp search of UniProt DB")
    blast_out = os.path.join(args.out, 'annotations.swissprot.txt')
    if not os.path.isfile(blast_out):
        lib.SwissProtBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
    num_annotations = lib.line_count(blast_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    #run MEROPS Blast search
    lib.log.info("Running Blastp search of MEROPS protease DB")
    blast_out = os.path.join(args.out, 'annotations.merops.txt')
    if not os.path.isfile(blast_out):
        lib.MEROPSBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
    num_annotations = lib.line_count(blast_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    #run EggNog search
    eggnog_out = os.path.join(args.out, 'annotations.eggnog.txt')
    lib.log.info("Annotating proteins with EggNog 4.5 database")
    if not os.path.isfile(eggnog_out):
        lib.runEggNog(Proteins, args.cpus, 1e-10, args.out, eggnog_out)
    num_annotations = lib.line_count(eggnog_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    #run dbCAN search
    dbCAN_out = os.path.join(args.out, 'annotations.dbCAN.txt')
    lib.log.info("Annotating CAZYmes using dbCAN")
    if not os.path.isfile(dbCAN_out):
        lib.dbCANsearch(Proteins, args.cpus, 1e-17, args.out, dbCAN_out)
    num_annotations = lib.line_count(dbCAN_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    
#now collect the results from InterProscan, then start to reformat results
lib.log.info("RunIprScan has finished, now pulling out annotations from results")
IPR_terms = os.path.join(args.out, 'annotations.iprscan.txt')
if not os.path.isfile(IPR_terms):
    IPR2TSV = os.path.join(RUNIPRSCAN_PATH, 'ipr2tsv.py')
    with open(IPR_terms, 'w') as output:
        subprocess.call([sys.executable, IPR2TSV, IPROUT], stdout = output, stderr = FNULL)
GO_terms = os.path.join(args.out, 'annotations.GO.txt')
if not os.path.isfile(GO_terms):
    IPR2GO = os.path.join(RUNIPRSCAN_PATH, 'ipr2go.py')
    OBO = os.path.join(currentdir, 'DB', 'go.obo')
    with open(GO_terms, 'w') as output:
        subprocess.call([sys.executable, IPR2GO, OBO, IPROUT], stdout = output, stderr = FNULL)
        
#check if antiSMASH data is given, if so parse and reformat for annotations and cluster textual output
if args.antismash:
    AntiSmashFolder = os.path.join(args.out, 'antismash')
    AntiSmashBed = os.path.join(AntiSmashFolder, 'clusters.bed')
    GFF2clusters = os.path.join(AntiSmashFolder, 'secmet.clusters.txt')
    AntiSmash_annotations = os.path.join(args.out, 'annotations.antismash.txt')
    Cluster_annotations = os.path.join(args.out, 'annotations.antismash.clusters.txt')
    if not os.path.isdir(AntiSmashFolder):
        os.makedirs(AntiSmashFolder)
    lib.ParseAntiSmash(args.antismash, AntiSmashFolder, AntiSmashBed, AntiSmash_annotations) #results in several global dictionaries
    lib.GetClusterGenes(AntiSmashBed, GFF, GFF2clusters, Cluster_annotations) #results in dictClusters dictionary
     
#now bring all annotations together and annotated genome using gag
ANNOTS = os.path.join(args.out, 'all.annotations.txt')
with open(ANNOTS, 'w') as output:
    for file in os.listdir(args.out):
        if file.startswith('annotations'):
            file = os.path.join(args.out, file)
            with open(file) as input:
                output.write(input.read())
ANNOTS = os.path.abspath(ANNOTS)

#launch gag
GAG = os.path.join(args.out, 'gag')
lib.log.info("Adding annotations to GFF using GAG")
subprocess.call(['gag.py', '-f', Scaffolds, '-g', GFF, '-a', ANNOTS, '-o', GAG], stdout = FNULL, stderr = FNULL)

#fix the tbl file for tRNA genes
lib.log.info("Fixing tRNA annotations in GenBank tbl file")
original = os.path.join(args.out, 'gag', 'genome.tbl')
tmp_tbl = os.path.join(args.out, 'gag', 'genome.tbl.original')
os.rename(original, tmp_tbl)
lib.CleantRNAtbl(GFF, tmp_tbl, original)

#write to GBK file
if not isolate == '???':
    ORGANISM = "[organism=" + organism + "] " + "[isolate=" + isolate + "]"
else:
    ORGANISM = "[organism=" + organism + "]"

shutil.copyfile(os.path.join(GAG, 'genome.fasta'), os.path.join(GAG, 'genome.fsa'))
discrep = 'discrepency.report.txt'
lib.log.info("Converting to final Genbank format, good luck!.....")
subprocess.call(['tbl2asn', '-p', GAG, '-t', SBT, '-M', 'n', '-Z', discrep, '-a', 'r10u', '-l', 'paired-ends', '-j', ORGANISM, '-V', 'b', '-c', 'fx'], stdout = FNULL, stderr = FNULL)

#rename and organize output files
os.rename(discrep, os.path.join(args.out, discrep))
OutputGBK = organism + '_' + isolate + '.gbk'
os.rename(os.path.join(args.out, 'gag', 'genome.gbf'), OutputGBK)

#write secondary metabolite clusters output using the final genome in gbk format
if args.antismash:
    #do a blast best hit search against MIBiG database for cluster annotation, but looping through gene cluster hits
    AllProts = []
    for k, v in lib.dictClusters.items():
        for i in v:
            if not i in AllProts:
                AllProts.append(i)
    AllProts = set(AllProts)
    mibig_fasta = os.path.join(AntiSmashFolder, 'smcluster.proteins.fasta')
    mibig_blast = os.path.join(AntiSmashFolder, 'smcluster.MIBiG.blast.txt')
    mibig_db = os.path.join(currentdir, 'DB', 'MIBiG')
    with open(mibig_fasta, 'w') as output:
        with open(Proteins, 'rU') as input:
            SeqRecords = SeqIO.parse(Proteins, 'fasta')
            for record in SeqRecords:
                if record.id in AllProts:
                    SeqIO.write(record, output, 'fasta')
    subprocess.call(['blastp', '-query', mibig_fasta, '-db', mibig_db, '-num_threads', args.cpus, '-max_target_seqs', '1', '-outformat', '6', '-out', mibig_blast])
    #now parse blast results to get {qseqid: hit}
    MIBiGBlast = {}
    with open(mibig_blast, 'rU') as input:
        for line in input:
            cols = line.split('\t')
            ID = cols[0]
            hit = cols[2].split('|')
            desc = hit[5]
            cluster = hit[0]
            db_ref = hit[-1]
            evalue = cols[10]
            pident = cols[2]
            result = (desc, cluster, db_ref, pident, evalue)
            MIBiGBlast[ID] = result
            
    lib.log.info("Creating tab-delimited SM cluster output")
    #load in EggNog annotations to get descriptions for table
    EggNog = {}
    with open(os.path.join(currentdir, 'DB','FuNOG.annotations.tsv'), 'rU') as input:
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
    with open(OutputGBK, 'rU') as gbk:
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
                                output.write("%s\t%s:%i-%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, chr, actualStart, actualEnd, strand, location, enzyme, domains, product, cog, eggnogDesc, IP, PF, GO, No, mibig, prot_seq, DNA_seq))
                                                             
    #now put together into a single file
    finallist = []
    ClustersOut = organism + '_' + isolate + '.clusters.txt'
    for file in os.listdir(AntiSmashFolder):
        if file.endswith('secmet.cluster.txt'):
            file = os.path.join(AntiSmashFolder, file)
            finallist.append(file)
    with open(ClustersOut, 'w') as output:
        for file in natsorted(finallist):
            with open(file, 'rU') as input:
                output.write(input.read())
                output.write('\n\n')

os._exit(1)
    

