#!/usr/bin/env python
from __future__ import division
import sys, os, subprocess, inspect, shutil, argparse, fnmatch
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib
from lib.interlap import InterLap
from collections import defaultdict
from natsort import natsorted
import numpy as np

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='funannotate-update.py', usage="%(prog)s [options] -i genome.gbk -l left.fq.gz -r right.fg.gz",
    description = '''Script is a wrapper for automated Trinity/PASA reannotation.''',
    epilog = """Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--input', required=True, help='Genome in GBK format or funannotate folder')
parser.add_argument('-l', '--left', nargs='+', help='Left (R1) FASTQ Reads')
parser.add_argument('--left_norm', help='Left (R1) FASTQ Reads')
parser.add_argument('--right_norm', help='Right (R2) normalized FASTQ Reads')
parser.add_argument('--single_norm', help='single normalized FASTQ Reads')
parser.add_argument('-r', '--right', nargs='+', help='Right (R2) FASTQ Reads')
parser.add_argument('-s', '--single', nargs='+', help='Single ended FASTQ Reads')
parser.add_argument('-o', '--out', help='Basename of output files')
parser.add_argument('--species', help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('-c', '--coverage', default=50, type=int, help='Depth to normalize reads to')
parser.add_argument('--isolate', help='Isolate name (e.g. Af293)')
parser.add_argument('--strain', help='Strain name (e.g. CEA10)')
parser.add_argument('--trinity', help='Trinity genome guided FASTA results')
parser.add_argument('--pasa_gff', help='PASA GFF')
parser.add_argument('--pasa_alignment_overlap', default='30.0', help='PASA --stringent_alingment_overlap')
parser.add_argument('--pasa_config', help='PASA assembly configuration file')
parser.add_argument('--memory', default='50G', help='RAM to use for Jellyfish/Trinity')
parser.add_argument('--no_normalize_reads', action='store_true', help='skip normalization')
parser.add_argument('--no_trimmomatic', action='store_true', help='skip quality trimming via trimmomatic')
parser.add_argument('--no_antisense_filter', action='store_true', help='skip antisense filtering')
parser.add_argument('--jaccard_clip', action='store_true', help='Turn on jaccard_clip for dense genomes')
parser.add_argument('--kallisto', help='Kallisto abundances table')
parser.add_argument('--name', help='Shortname for genes, perhaps assigned by NCBI, eg. VC83_')
parser.add_argument('--max_intronlen', default=3000, help='Maximum intron length for gene models')
parser.add_argument('--min_protlen', default=50, type=int, help='Minimum amino acid length for valid gene model')
parser.add_argument('--stranded', default = 'no', choices=['RF','FR','F','R','no'], help='RNA seq strandedness')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to use')
parser.add_argument('-t','--tbl2asn', default='-l paired-ends', help='Parameters for tbl2asn, linkage and gap info')
parser.add_argument('--sbt', default='SBT', help='Basename of output files')
parser.add_argument('--p2g', help='NCBI p2g file from previous annotation')
parser.add_argument('--PASAHOME', help='Path to PASA home directory, $PASAHOME')
parser.add_argument('--TRINITYHOME', help='Path to Trinity config directory, $TRINITYHOME')
parser.add_argument('--SeqCenter', default='CFMR', help='Sequencing center for GenBank tbl file')
parser.add_argument('--SeqAccession', default='12345', help='Sequencing accession number')
parser.add_argument('--alt_transcripts', default='0.10', help='Threshold to keep alt-transcripts, percent highest expression')
args=parser.parse_args()

#functions
def worker(input):
    logfile = os.path.join(tmpdir, 'Trinity-gg.log')
    cmd = input.split(' ')
    #make sure no empty items
    cmd = [x for x in cmd if x]
    with open(logfile, 'a') as output:
        subprocess.call(cmd, cwd=os.path.join(tmpdir, 'trinity_gg'), stdout = output, stderr = output)

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: worker(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

def gbk2pasa(input, gffout, trnaout, fastaout, spliceout, exonout, proteinsout):
    LocusTags = []
    multiExon = {}
    geneMapping = {}
    with open(gffout, 'w') as gff:
        gff.write("##gff-version 3\n")
        with open(trnaout, 'w') as trna:
            with open(fastaout, 'w') as fasta:
                with open(proteinsout, 'w') as proteins:
                    with open(input, 'rU') as gbk:
                        for record in SeqIO.parse(gbk, 'genbank'):
                            fasta.write(">%s\n%s\n" % (record.id, record.seq))
                            for f in record.features:
                                if f.type == 'CDS':
                                    chr = record.id
                                    ID = f.qualifiers['locus_tag'][0]
                                    if not ID in LocusTags:
                                        LocusTags.append(ID)
                                    if not ID in geneMapping:
                                        geneMapping[ID] = 1
                                    else:
                                        geneMapping[ID] += 1
                                    protID = f.qualifiers['protein_id'][0]
                                    protSeq = f.qualifiers['translation'][0]
                                    product = f.qualifiers['product'][0]
                                    proteins.write('>%s|%s\n%s\n' % (ID,protID,protSeq))
                                    start = f.location.nofuzzy_start + 1
                                    end = f.location.nofuzzy_end
                                    strand = f.location.strand
                                    if strand == 1:
                                        strand = '+'
                                    elif strand == -1:
                                        strand = '-'
                                    num_exons = len(f.location.parts)
                                    TranNum = geneMapping.get(ID)
                                    current_phase = int(f.qualifiers['codon_start'][0]) - 1 #need to adjust NCBI to GFF3 notation
                                    gff.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n" % (chr, start, end, strand, ID))
                                    gff.write("%s\tGenBank\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s-T%i;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, TranNum, ID, product))
                                    if num_exons < 2: #only a single exon
                                        ex_start = str(f.location.nofuzzy_start + 1)
                                        ex_end = str(f.location.nofuzzy_end)
                                        gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T%i.exon1;Parent=%s-T%i\n" % (chr, ex_start, ex_end, strand, ID, TranNum, ID, TranNum))
                                        gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T%i.cds;Parent=%s-T%i\n" % (chr, ex_start, ex_end, strand, current_phase, ID, TranNum, ID, TranNum))
                                    else: #more than 1 exon, so parts sub_features
                                        splices = []
                                        for i in range(0,num_exons):
                                            ex_start = str(f.location.parts[i].nofuzzy_start + 1)
                                            ex_end = str(f.location.parts[i].nofuzzy_end)
                                            ex_num = i + 1
                                            gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T%i.exon%i;Parent=%s-T%i\n" % (chr, ex_start, ex_end, strand, ID, TranNum, ex_num, ID, TranNum))
                                            gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T%i.cds;Parent=%s-T%i\n" % (chr, ex_start, ex_end, strand, current_phase, ID, TranNum, ID, TranNum))
                                            current_phase = (current_phase - (int(ex_end) - int(ex_start) + 1)) % 3
                                            if current_phase == 3:
                                                current_phase = 0
                                            #ss and exons are 0-based position, so 1 less than GFF
                                            exons_start = int(ex_start) - 1
                                            exons_end = int(ex_end) -1                                            
                                            #add to exon dictionary
                                            if not ID in multiExon:
                                                multiExon[ID] = [chr, strand, [(exons_start,exons_end)]]
                                            else:
                                                multiExon[ID][2].append((exons_start,exons_end))
                                #here lets enforce new NCBI tRNA rules about length to prevent tbl2asn errors later on     
                                if f.type == 'tRNA':
                                    ID = f.qualifiers['locus_tag'][0]
                                    if not ID in LocusTags:
                                        LocusTags.append(ID)
                                    if not ID in geneMapping:
                                        geneMapping[ID] = 1
                                    else:
                                        geneMapping[ID] += 1
                                    start = f.location.nofuzzy_start
                                    end = f.location.nofuzzy_end
                                    strand = f.location.strand
                                    if strand == 1:
                                        strand = '+'
                                        length = abs(int(end) - int(start))
                                    elif strand == -1:
                                        strand = '-'
                                        length = abs(int(start) - int(end))
                                    try:
                                        product = f.qualifiers['product'][0]
                                    except KeyError:
                                        product = "tRNA-XXX"
                                    chr = record.id
                                    if length < 50 or length > 150:
                                        continue
                                    else:
                                        trna.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n" % (chr, start, end, strand, ID))
                                        trna.write("%s\tGenBank\ttRNA\t%s\t%s\t.\t%s\t.\tID=%s-T%i;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, TranNum, ID, product))
                                        trna.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T%i.exon1;Parent=%s-T%i\n" % (chr, start, end, strand, ID, TranNum, ID, TranNum))
    #parse splice sites and write to file
    with open(exonout, 'w') as exon:
        with open(spliceout, 'w') as splicer:
            for k,v in natsorted(multiExon.items()):
                sortedList = sorted(v[2], key=lambda tup: tup[0])
                for y in sortedList:
                    exon.write("%s\t%i\t%i\t%s\n" % (v[0], y[0], y[1], v[1])) 
                splices = []
                for i in range(1, len(sortedList)):
                    splices.append((sortedList[i-1][1], sortedList[i][0]))
                for x in splices:
                    splicer.write('%s\t%i\t%i\t%s\n' % (v[0], x[0], x[1], v[1]))
    
    #finally lets return the base locus tag name and the last number
    lastTag = natsorted(LocusTags)[-1]
    if '_' in lastTag:
        tag, count = lastTag.split('_')
        tag = tag+'_'
    else:
        for i,c in enumerate(lastTag):
            if c.isdigit():
                tag = lastTag[:i]
                count = lastTag[i:]
                break
    justify = len(count)
    return tag, count, justify

def Funzip(input, output, cpus):
    '''
    function to unzip as fast as it can, pigz -> bgzip -> gzip
    '''
    if lib.which('pigz'):
        cmd = ['pigz', '--decompress', '-c', '-p', str(cpus), input]
    elif lib.which('bgzip'):
        cmd = ['bgzip', '--decompress', '-c', '-@', str(cpus), input]
    else:
        cmd = ['gzip', '--decompress', '-c', input]
    try:
        lib.runSubprocess2(cmd, '.', lib.log, output)
    except NameError:
        with open(output, 'w') as outfile:
            subprocess.call(cmd, stdout=outfile)

def Fzip(input, output, cpus):
    '''
    function to zip as fast as it can, pigz -> bgzip -> gzip
    '''
    if lib.which('pigz'):
        cmd = ['pigz', '-c', '-p', str(cpus), input]
    elif lib.which('bgzip'):
        cmd = ['bgzip', '-c', '-@', str(cpus), input]
    else:
        cmd = ['gzip', '-c', input]
    try:
        lib.runSubprocess2(cmd, '.', lib.log, output)
    except NameError:
        with open(output, 'w') as outfile:
            subprocess.call(cmd, stdout=outfile)

def Fzip_inplace(input, cpus):
    '''
    function to zip as fast as it can, pigz -> bgzip -> gzip
    '''
    if lib.which('pigz'):
        cmd = ['pigz', '-f', '-p', str(cpus), input]
    elif lib.which('bgzip'):
        cmd = ['bgzip', '-f', '-@', str(cpus), input]
    else:
        cmd = ['gzip', '-f', input]
    try:
        lib.runSubprocess(cmd, '.', lib.log)
    except NameError:
        subprocess.call(cmd)

def choose_aligner():
    '''
    function to choose alignment method for mapping reads to transcripts to determine
    orientation of the trinity transcripts. rapmap -> bowtie2 -> hisat2
    note hisat2 is probably not ideal for this, but should work okay.
    '''
    aligner = ''
    if lib.which('rapmap'):
        aligner = 'rapmap'
    elif lib.which('bowtie2'):
        aligner = 'bowtie2'
    else:
        aligner = 'hisat2'
    return aligner

def runTrimmomaticPE(left, right):
    '''
    function is wrapper for Trinity trimmomatic
    '''
    #create tmpdir
    folder = os.path.join(tmpdir, 'trimmomatic')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    lib.log.info("Adapter and Quality trimming PE reads with Trimmomatic")
    left_paired = os.path.join(folder, 'trimmed_left.fastq')
    left_single = os.path.join(folder, 'trimmed_left.unpaired.fastq')
    right_paired = os.path.join(folder, 'trimmed_right.fastq')
    right_single = os.path.join(folder, 'trimmed_right.unpaired.fastq')
    TRIMMOMATIC_DIR = os.path.join(TRINITY, 'trinity-plugins', 'Trimmomatic-0.36')
    cmd = ['java', '-jar', os.path.join(TRIMMOMATIC_DIR, 'trimmomatic.jar'), 'PE', '-threads', str(args.cpus), '-phred33', 
            left, right, left_paired, left_single, right_paired, right_single,
            'ILLUMINACLIP:'+os.path.join(TRIMMOMATIC_DIR,'adapters','TruSeq3-PE.fa')+':2:30:10', 'SLIDINGWINDOW:4:5', 'LEADING:5', 'TRAILING:5', 'MINLEN:25']
    lib.runSubprocess(cmd, '.', lib.log)
    for x in [left_paired, left_single, right_paired, right_single]:
        Fzip_inplace(x, args.cpus)
    trim_left = os.path.join(folder, 'trimmed_left.fastq.gz')
    trim_right = os.path.join(folder, 'trimmed_right.fastq.gz')
    return trim_left, trim_right

def runTrimmomaticSE(reads):
    '''
    function is wrapper for Trinity trimmomatic
    '''
    #create tmpdir
    folder = os.path.join(tmpdir, 'trimmomatic')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    lib.log.info("Adapter and Quality trimming SE reads with Trimmomatic")
    output = os.path.join(folder, 'trimmed_single.fastq')
    TRIMMOMATIC_DIR = os.path.join(TRINITY, 'trinity-plugins', 'Trimmomatic-0.36')
    cmd = ['java', '-jar', os.path.join(TRIMMOMATIC_DIR, 'trimmomatic.jar'), 'SE', '-threads', str(args.cpus), '-phred33', 
            reads, output, 'ILLUMINACLIP:'+os.path.join(TRIMMOMATIC_DIR,'adapters','TruSeq3-SE.fa')+':2:30:10', 'SLIDINGWINDOW:4:5', 'LEADING:5', 'TRAILING:5', 'MINLEN:25']
    lib.runSubprocess(cmd, '.', lib.log)
    Fzip_inplace(output, args.cpus)
    trim_single = os.path.join(folder, 'trimmed_single.fastq.gz')
    return trim_single

def runNormalization(readTuple, memory):
    '''
    function is wrapper for Trinity read normalization
    have to run normalization separately for PE versus single
    '''
    left_norm, right_norm, single_norm = (None,)*3
    SENormalLog = os.path.join(tmpdir, 'trinity_normalization.SE.log')
    PENormalLog = os.path.join(tmpdir, 'trinity_normalization.PE.log')
    if args.stranded != 'no':
        cmd = [os.path.join(TRINITY, 'util', 'insilico_read_normalization.pl'), '--PARALLEL_STATS', '--JM', memory, '--max_cov', str(args.coverage), '--seqType', 'fq', '--output', os.path.join(tmpdir, 'normalize'), '--CPU', str(args.cpus), '--SS_lib_type', args.stranded]
    else:
        cmd = [os.path.join(TRINITY, 'util', 'insilico_read_normalization.pl'), '--PARALLEL_STATS', '--JM', memory, '--max_cov', str(args.coverage), '--seqType', 'fq', '--output', os.path.join(tmpdir, 'normalize'), '--CPU', str(args.cpus)]
    if readTuple[2]: #single reads present, so run normalization just on those reads
        cmd = cmd + ['--single', readTuple[2]]
        lib.runSubprocess2(cmd, '.', lib.log, SENormalLog)
        single_norm = os.path.join(tmpdir, 'normalize', 'single.norm.fq')
    if readTuple[0] and readTuple[1]:
        cmd = cmd + ['--pairs_together', '--left', readTuple[0], '--right', readTuple[1]]
        left_norm = os.path.join(tmpdir, 'normalize', 'left.norm.fq')
        right_norm = os.path.join(tmpdir, 'normalize', 'right.norm.fq') 
        lib.runSubprocess2(cmd, '.', lib.log, PENormalLog)

    return left_norm, right_norm, single_norm

def concatenateReads(input, output):
    '''
    Since I can't seem to get the comma separated lists to work with subprocess modules, just 
    concatenate FASTQ files in order and use a single file, input should be a list of FASTQ files
    using system cat here so that gzipped files are concatenated correctly
    '''
    cmd = ['cat']
    cmd = cmd + input
    lib.runSubprocess2(cmd, '.', lib.log, output)

def runTrinityGG(genome, readTuple, hisatexons, hisatss, output):
    '''
    function will run genome guided Trinity. First step will be to run hisat2 to align reads
    to the genome, then pass that BAM file to Trinity to generate assemblies
    '''
    #build hisat2 index, using exons and splice sites
    lib.log.info("Starting Trinity genome guided")
    lib.log.info("Building Hisat2 genome index, incorporating exons and splice-sites")
    cmd = ['hisat2-build', '--exon', hisatexons, '--ss', hisatss, genome, os.path.join(tmpdir, 'hisat2.genome')]
    lib.runSubprocess4(cmd, '.', lib.log)
    #align reads using hisat2
    lib.log.info("Aligning reads to genome using Hisat2")
    hisat2bam = os.path.join(tmpdir, 'hisat2.coordSorted.bam')
    #use bash wrapper for samtools piping for SAM -> BAM -> sortedBAM
    bamthreads = (args.cpus + 2 // 2) // 2 #use half number of threads for bam compression threads
    if args.stranded != 'no' and not readTuple[2]:
        hisat2cmd = ['hisat2', '-p', str(args.cpus), '--max-intronlen', str(args.max_intronlen), '--dta', '-x', os.path.join(tmpdir, 'hisat2.genome'), '--rna-strandness', args.stranded]
    else:
        hisat2cmd = ['hisat2', '-p', str(args.cpus), '--max-intronlen', str(args.max_intronlen), '--dta', '-x', os.path.join(tmpdir, 'hisat2.genome')]
    if readTuple[0] and readTuple[1]:
        hisat2cmd = hisat2cmd + ['-1', readTuple[0], '-2', readTuple[1]]
    if readTuple[2]:
        hisat2cmd = hisat2cmd + ['-U', readTuple[2]]
        
    cmd = [os.path.join(parentdir, 'util', 'sam2bam.sh'), " ".join(hisat2cmd), str(bamthreads), hisat2bam]
    lib.runSubprocess(cmd, '.', lib.log)

    #now launch Trinity genome guided
    TrinityLog = os.path.join(tmpdir, 'Trinity-gg.log')
    lib.log.info("Running genome-guided Trinity, logfile: %s" % TrinityLog)
    lib.log.info("Clustering of reads from BAM and preparing assembly commands")
    jaccard_clip = []
    if args.jaccard_clip:
        jaccard_clip = ['--jaccard_clip']
    if args.stranded != 'no' and not readTuple[2]:
        cmd = ['Trinity', '--SS_lib_type', args.stranded, '--no_distributed_trinity_exec', '--genome_guided_bam', hisat2bam, '--genome_guided_max_intron', str(args.max_intronlen), '--CPU', str(args.cpus), '--max_memory', args.memory, '--output', os.path.join(tmpdir, 'trinity_gg')]
    else:
        cmd = ['Trinity', '--no_distributed_trinity_exec', '--genome_guided_bam', hisat2bam, '--genome_guided_max_intron', str(args.max_intronlen), '--CPU', str(args.cpus), '--max_memory', args.memory, '--output', os.path.join(tmpdir, 'trinity_gg')]
    cmd = cmd + jaccard_clip
    lib.runSubprocess2(cmd, '.', lib.log, TrinityLog)
    commands = os.path.join(tmpdir, 'trinity_gg', 'trinity_GG.cmds')

    #this will create all the Trinity commands, will now run these in parallel using multiprocessing in Python (seems to be much faster than Parafly on my system)
    file_list = []
    with open(commands, 'rU') as cmdFile:
        for line in cmdFile:
            line = line.replace('\n', '')
            line = line.replace('--no_distributed_trinity_exec', '') #don't think this should be appended to every command....
            line = line.replace('"', '') #don't need these double quotes
            file_list.append(line)
    lib.log.info("Assembling "+"{0:,}".format(len(file_list))+" Trinity clusters using %i CPUs" % (args.cpus-1))
    lib.runMultiProgress(safe_run, file_list, args.cpus-1)

    #collected output files and clean
    outputfiles = os.path.join(tmpdir, 'trinity_gg', 'trinity_output_files.txt')
    with open(outputfiles, 'w') as fileout:
        for filename in find_files(os.path.join(tmpdir, 'trinity_gg'), '*inity.fasta'):
            fileout.write('%s\n' % filename)
    #now grab them all using Trinity script
    cmd = [os.path.join(TRINITY, 'util', 'support_scripts', 'GG_partitioned_trinity_aggregator.pl'), 'Trinity_GG']
    lib.runSubprocess5(cmd, '.', lib.log, outputfiles, output)

def polyAclip(input, output):
    #now parse the input fasta file removing records in list
    totalclipped = 0
    with open(output, 'w') as outfile:
        with open(input, 'rU') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                #now check for poly-A tail
                Seq = str(record.seq)
                if Seq.endswith('AAAAAAA'): #then count number of As
                    count = 0
                    for base in reversed(Seq):
                        if base == 'A':
                            count += 1
                        else:
                            break
                    Seq = Seq[:-count]
                    totalclipped += 1
                outfile.write(">%s\n%s\n" % (record.description, Seq))
    lib.log.info("Clipped {:,} poly-A tails from transcripts".format(totalclipped))

def removeAntiSense(input, readTuple, output):
    '''
    function will map reads to the input transcripts, determine strandedness, and then filter
    out transcripts that were assembled in antisense orientation. idea here is that the antisense
    transcripts, while potentially valid, aren't going to help update the gene models and perhaps
    could hurt the annotation effort?
    '''
    lib.log.info("Running anti-sense filtering of Trinity transcripts")
    bamthreads = (args.cpus + 2 // 2) // 2 #use half number of threads for bam compression threads
    aligner = choose_aligner()
    if aligner == 'hisat2':
        lib.log.info("Building Hisat2 index of "+"{0:,}".format(lib.countfasta(input))+" trinity transcripts")
        cmd = ['hisat2-build', input, os.path.join(tmpdir, 'hisat2.transcripts')]
        lib.runSubprocess4(cmd, '.', lib.log)

        #now launch the aligner
        lib.log.info("Aligning reads to trinity transcripts with Hisat2")
        hisat2cmd = ['hisat2', '-p', str(args.cpus), '-k', '50', '--max-intronlen', str(args.max_intronlen), '-x', os.path.join(tmpdir, 'hisat2.transcripts')]
        if readTuple[2]:
            hisat2cmd = hisat2cmd + ['-U', readTuple[2]]
        if readTuple[0] and readTuple[1]:
            hisat2cmd = hisat2cmd + ['-1', readTuple[0], '-2', readTuple[1]]
        bowtie2bam = os.path.join(tmpdir, 'hisat2.transcripts.coordSorted.bam')
        cmd = [os.path.join(parentdir, 'util', 'sam2bam.sh'), " ".join(hisat2cmd), str(bamthreads), bowtie2bam]
        lib.runSubprocess4(cmd, '.', lib.log)

    elif aligner == 'bowtie2':
        #using bowtie2
        lib.log.info("Building Bowtie2 index of "+"{0:,}".format(lib.countfasta(input))+" trinity transcripts")
        cmd = ['bowtie2-build', input, os.path.join(tmpdir, 'bowtie2.transcripts')]
        lib.runSubprocess4(cmd, '.', lib.log)

        #now launch the subprocess commands in order
        lib.log.info("Aligning reads to trinity transcripts with Bowtie2")
        bowtie2cmd = ['bowtie2', '-p', str(args.cpus), '-k', '50', '--local', '--no-unal', '-x', os.path.join(tmpdir, 'bowtie2.transcripts')]
        if readTuple[2]:
            bowtie2cmd = bowtie2cmd + ['-U', readTuple[2]]
        if readTuple[0] and readTuple[1]:
            bowtie2cmd = bowtie2cmd + ['-1', readTuple[0], '-2', readTuple[1]]
        bowtie2bam = os.path.join(tmpdir, 'bowtie2.transcripts.coordSorted.bam')
        cmd = [os.path.join(parentdir, 'util', 'sam2bam.sh'), " ".join(bowtie2cmd), str(bamthreads), bowtie2bam]
        lib.runSubprocess4(cmd, '.', lib.log)
        
    elif aligner == 'rapmap':
        #using bowtie2
        lib.log.info("Building RapMap index of "+"{0:,}".format(lib.countfasta(input))+" trinity transcripts")
        cmd = ['rapmap', 'quasiindex', '-t', input, '-i', os.path.join(tmpdir, 'rapmap_index')]
        lib.runSubprocess4(cmd, '.', lib.log)
        #now launch the subprocess commands in order
        lib.log.info("Aligning reads to trinity transcripts with RapMap")
        rapmapcmd = ['rapmap', 'quasimap', '-t', str(args.cpus), '-i', os.path.join(tmpdir, 'rapmap_index'), '-1', readTuple[0], '-2', readTuple[1]]
        bowtie2bam = os.path.join(tmpdir, 'rapmap.transcripts.coordSorted.bam')
        cmd = [os.path.join(parentdir, 'util', 'sam2bam.sh'), " ".join(rapmapcmd), str(bamthreads), bowtie2bam]
        lib.runSubprocess(cmd, '.', lib.log)        

    #now run Trinity examine strandeness tool
    lib.log.info("Examining strand specificity")
    cmd = [os.path.join(TRINITY, 'util', 'misc', 'examine_strand_specificity.pl'), bowtie2bam, os.path.join(tmpdir, 'strand_specific')]
    lib.runSubprocess(cmd, '.', lib.log)
    #parse output dat file and get list of transcripts to remove
    removeList = []
    with open(os.path.join(tmpdir, 'strand_specific.dat'), 'rU') as infile:
        for line in infile:
            line = line.replace('\n', '')
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            if args.stranded == 'RF': #then we want to keep negative ratios in cols[4]
                if not cols[4].startswith('-'):
                    removeList.append(cols[0])
            elif args.stranded == 'FR': #keep + values
                if cols[4].startswith('-'):
                    removeList.append(cols[0])
    
    #now parse the input fasta file removing records in list
    with open(output, 'w') as outfile:
        with open(input, 'rU') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                if not record.id in removeList:
                    outfile.write(">%s\n%s\n" % (record.description, str(record.seq)))
    lib.log.info("Removing %i antisense transcripts" % (len(removeList)))

def getPASAinformation(DBname, folder, genome):
    '''
    function to dump GFF from existing PASA database, compare genome headers to what is in PASA
    DB to make sure genome is same, return True if good to go, else error out
    '''
    #run some checks of the data to make sure it is same assembly
    mysqlDB, mysqlUser, mysqlPass = (None,)*3
    with open(os.path.join(PASA, 'pasa_conf', 'conf.txt'), 'rU') as pasaconf:
        for line in pasaconf:
            line = line.replace('\n', '')
            if line.startswith('MYSQLSERVER='):
                mysqlDB = line.split('=')[-1]
            if line.startswith('MYSQL_RW_USER='):
                mysqlUser = line.split('=')[-1]
            if line.startswith('MYSQL_RW_PASSWORD='):
                mysqlPass = line.split('=')[-1]
    pasaExistingGFF = os.path.join(folder, 'existing_pasa.gff3')
    cmd = [os.path.join(PASA, 'scripts', 'pasa_asmbl_genes_to_GFF3.dbi'), '-M', DBname+':'+mysqlDB, '-p', mysqlUser+':'+mysqlPass]
    lib.runSubprocess2(cmd, folder, lib.log, pasaExistingGFF)
    #now get number of genes and list of contigs
    pasaContigs = []
    geneCount = 0
    with open(pasaExistingGFF, 'rU') as infile:
        for line in infile:
            if line.startswith('\n'):
                continue
            cols = line.split('\t')
            if not cols[0] in pasaContigs:
                pasaContigs.append(cols[0])
            if cols[2] == 'gene':
                geneCount += 1
    #now get fasta headers from genome
    genomeContigs = []
    with open(genome, 'rU') as fasta:
        for line in fasta:
            if line.startswith('>'):
                line = line.replace('\n', '')
                line = line.replace('>', '')
                if not line in genomeContigs:
                    genomeContigs.append(line)
    #now make sure PASA headers in genome
    genomeContigs = set(genomeContigs)
    for contig in pasaContigs:
        if not contig in genomeContigs:
            return False
    lib.log.info("Existing PASA database contains {:,} gene models, validated FASTA headers match".format(geneCount))
    return True
    
 
def runPASA(genome, transcripts, stranded, intronlen, cpus, previousGFF, dbname, output, configFile):
    '''
    function will run PASA align assembly, followed by 2 rounds of comparison to update
    annotations for preexisting gene models
    '''
    #pasa cpus are IO bound? and for gmap and blat are split, so this sould be 1/2 of funannotate cpus
    if cpus > 2:
        pasa_cpus = cpus / 2
    else:
        pasa_cpus = 2
    pasa_cpus = int(pasa_cpus)
    #create tmpdir
    folder = os.path.join(tmpdir, 'pasa')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    #get config files and edit
    alignConfig = os.path.join(folder, 'alignAssembly.txt')
    annotConfig = os.path.join(folder, 'annotCompare.txt')
    #check if config file is passed, if so, get databasename and copy to assembly config file
    DataBaseName = dbname.replace('-', '_') #dashes will get stripped in MySQL
    if configFile:
        with open(configFile, 'rU') as infile:
            for line in infile:
                line = line.replace('\n', '')
                if line.startswith('MYSQLDB='):
                    DataBaseName = line.split('=')[-1]
        shutil.copyfile(configFile, alignConfig)
        #check existing database
        if not getPASAinformation(DataBaseName, folder, genome):
            lib.log.error("Headers in PASA database, do not match those in FASTA, exiting.")
            sys.exit(1)
        #finally need to index the genome using cdbfasta so lookups can be done
        cmd = [os.path.join(PASA, 'bin', 'cdbfasta'), genome]
        lib.runSubprocess(cmd, '.', lib.log)
    else:
        #create new config file from template
        with open(alignConfig, 'w') as config1:
            with open(os.path.join(PASA, 'pasa_conf', 'pasa.alignAssembly.Template.txt'), 'rU') as template1:
                for line in template1:
                    line = line.replace('<__MYSQLDB__>', DataBaseName)
                    config1.write(line)
        #now run PASA alignment step
        lib.log.info("Running PASA alignment step using "+"{0:,}".format(lib.countfasta(transcripts))+" transcripts")
        if stranded == 'no':
            cmd = [os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl'), '-c', os.path.abspath(alignConfig), '-r', '-C', '-R', '-g', os.path.abspath(genome), '--ALIGNERS', 'blat,gmap', '-t', os.path.abspath(transcripts), '--stringent_alignment_overlap', args.pasa_alignment_overlap, '--TRANSDECODER', '--MAX_INTRON_LENGTH', str(intronlen), '--CPU', str(pasa_cpus )]
        else:
            cmd = [os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl'), '-c', os.path.abspath(alignConfig), '-r', '-C', '-R', '-g', os.path.abspath(genome), '--ALIGNERS', 'blat,gmap', '-t', os.path.abspath(transcripts), '--transcribed_is_aligned_orient', '--stringent_alignment_overlap', args.pasa_alignment_overlap, '--TRANSDECODER', '--MAX_INTRON_LENGTH', str(intronlen), '--CPU', str(pasa_cpus )]
        lib.runSubprocess(cmd, folder, lib.log)
        
    #generate comparison template file
    with open(annotConfig, 'w') as config2:
        with open(os.path.join(PASA, 'pasa_conf', 'pasa.annotationCompare.Template.txt'), 'rU') as template2:
            for line in template2:
                line = line.replace('<__MYSQLDB__>', DataBaseName)
                config2.write(line)
    
    #now run Annotation comparisons
    lib.log.info("Running PASA annotation comparison step 1")
    cmd = [os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl'), '-c', os.path.abspath(annotConfig), '-g', os.path.abspath(genome), '-t', os.path.abspath(transcripts), '-A', '-L', '--annots_gff3', os.path.abspath(previousGFF), '--CPU', str(pasa_cpus)]
    lib.runSubprocess(cmd, folder, lib.log)
    round1GFF = None
    for file in os.listdir(folder):
        if not file.endswith('.gff3'):
            continue
        if 'gene_structures_post_PASA_updates' in file:     
            round1GFF = os.path.join(folder, file)
    if not round1GFF:
        lib.log.error("PASA failed, check log, exiting")
        sys.exit(1)
    #run round 2 comparison
    lib.log.info("Running PASA annotation comparison step 2")
    cmd = [os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl'), '-c', os.path.abspath(annotConfig), '-g', os.path.abspath(genome), '-t', os.path.abspath(transcripts), '-A', '-L', '--annots_gff3', os.path.abspath(round1GFF), '--CPU', str(pasa_cpus)]
    lib.runSubprocess(cmd, folder, lib.log)
    round2GFF = None
    for file in os.listdir(folder):
        if not file.endswith('.gff3'):
            continue
        if file == os.path.basename(round1GFF):
            continue
        if 'gene_structures_post_PASA_updates' in file: 
            round2GFF = os.path.join(folder, file)
    if not round2GFF:
        lib.log.error("PASA failed, check log, exiting")
        sys.exit(1)
    lib.log.debug("copying final PASA GFF3 to output: %s" % round2GFF)
    #grab final result
    shutil.copyfile(round2GFF, output)


def runKallisto(input, fasta, readTuple, stranded, cpus, output):
    '''
    function takes GFF3 output from PASA compare, extracts transcripts, and then calculates TPM
    using Kallisto to idenitfy the best scoring gene model for each locus, the left and right
    these should be the adapter cleaned non-normalized Illumina reads
    '''
    lib.log.info("Using Kallisto TPM data to determine which PASA gene models to select at each locus")
    #convert GFF to transcripts
    folder = os.path.join(tmpdir, 'getBestModel')
    os.makedirs(folder)
    PASAtranscripts = os.path.join(folder, 'transcripts.fa')
    cmd = [os.path.join(PASA, 'misc_utilities', 'gff3_file_to_proteins.pl'), input, fasta, 'cDNA']
    lib.log.info("Building Kallisto index")
    lib.runSubprocess2(cmd, '.', lib.log, PASAtranscripts)
    #generate kallisto index
    cmd = ['kallisto', 'index', '-i', os.path.join(folder, 'bestModel'), PASAtranscripts]
    lib.runSubprocess(cmd, '.', lib.log)
    #use kallisto to map reads to index
    #base command
    cmd = ['kallisto', 'quant', '-i', os.path.join(folder, 'bestModel'), '-o', os.path.join(folder, 'kallisto'), '--plaintext', '-t', str(cpus)]
    #parse the strand information
    if stranded == 'RF':
        strandcmd = ['--rf-stranded']
    elif stranded == 'FR':
        strandcmd = ['--fr-stranded']
    else:
        strandcmd = []
    #adapt command for input, i.e. single or PE ends -> what do you do if you have both?
    if readTuple[2] and not readTuple[0] and not readTuple[1]: #single, not just using estimated lengths and SD, I think this is okay? can make this an option otherwise
        cmd = cmd + ['--single', '-l', '200', '-s', '20', readTuple[2]]
    elif readTuple[0] and readTuple[1]:
        cmd = cmd + strandcmd + [readTuple[0], readTuple[1]]
    lib.log.info("Mapping reads using pseudoalignment in Kallisto")
    lib.runSubprocess(cmd, '.', lib.log)

    #modify kallisto ouput to map gene names to each mRNA ID so you know what locus they have come from
    mRNADict = {}
    #since mRNA is unique, parse the transcript file which has mRNAID geneID in header
    with open(PASAtranscripts, 'rU') as transin:
        for line in transin:
            if line.startswith('>'):
                line = line.rstrip()
                line = line.replace('>', '')
                cols = line.split(' ')
                mRNAID = cols[0]
                geneID = cols[1]
                location = cols[-1]
                if not mRNAID in mRNADict:
                    mRNADict[mRNAID] = (geneID, location)

    #now make new tsv file with #mRNAID geneID location TPM
    with open(output, 'w') as outfile:
        outfile.write("#mRNA-ID\tgene-ID\tLocation\tTPM\n")
        with open(os.path.join(folder, 'kallisto', 'abundance.tsv'), 'rU') as infile:
            for line in infile:
                if line.startswith('targed_id'):
                    continue
                line = line.rstrip()
                cols = line.split('\t')
                if cols[0] in mRNADict:
                    geneHit = mRNADict.get(cols[0])
                    geneID = geneHit[0]
                    location = geneHit[1]
                    outfile.write('%s\t%s\t%s\t%s\n' % (cols[0], geneID, location, cols[4]))

def getBestModels(input, fasta, abundances, alt_transcripts, outfile):
    #function to parse PASA results and generate GFF3; supports multiple transcripts
    if float(alt_transcripts) == 0:
        lib.log.info("Parsing Kallisto results. Keeping all alt-splicing transcripts at each locus.")
    elif float(alt_transcripts) < 1:
        lib.log.info("Parsing Kallisto results. Keeping alt-splicing transcipts if expressed at least {0:.1f}% of highest transcript per locus.".format(float(alt_transcripts)*100))
    else:
        lib.log.info("Parsing Kallisto results. Keeping best transcript at each locus.")
    bestModels = {}
    locations = {}
    with open(abundances, 'rU') as tpms:
        for line in tpms:
            line = line.rstrip()
            if line.startswith('#') or line.startswith('target_id'):
                continue
            transcriptID, geneID, Loc, TPM = line.split('\t')
            if not Loc in locations:
                locations[Loc] = geneID
                geneLocus = geneID
            else:
                geneLocus = locations.get(geneID)
            if not geneLocus in bestModels:
                bestModels[geneLocus] = [(transcriptID, float(TPM))]
            else:
                bestModels[geneLocus].append((transcriptID,float(TPM)))
    #now we have geneID dictionary containing list of tuples of of transcripts
    #loop through each locus grabbing transcript IDs of those models to keep
    #use best expression value * alt_transcript threshold to filter models
    #alt_transcript == 1 would be then only keeping best hit
    extractList = []
    ExpValues = {}
    for k,v in natsorted(bestModels.items()):
        if len(v) < 2:
            extractList.append(v[0][0])
            if not v[0][0] in ExpValues:
                ExpValues[v[0][0]] = v[0][1]
        else:
            sortedTranscripts = sorted(v, key=lambda tup: tup[1], reverse=True)
            ExpThreshold = sortedTranscripts[0][1] * float(alt_transcripts)
            for hit in sortedTranscripts:
                if hit[1] >= ExpThreshold:
                    extractList.append(hit[0])
                    if not hit[0] in ExpValues:
                        ExpValues[hit[0]] = hit[1]
    #now go through the PASA GFF file and generate filtered GFF3 file composed of extractList
    extractList = set(extractList)
    with open(outfile, 'w') as output:
        output.write("##gff-version 3\n")
        with open(input, 'rU') as gff:
            for line in gff:
                if line.startswith("#") or line.startswith('\n'):
                    continue
                line = line.rstrip()
                cols = line.split('\t')
                gffID = cols[8].split(';Parent')[0].replace('ID=', '')
                if 'gene' in cols[2]:
                    continue
                elif 'mRNA' in cols[2]:
                    if gffID in extractList:
                        geneID = cols[8].split(';Name=')[0]
                        geneID = 'ID='+ geneID.split(';Parent=')[-1]
                        mRNAID = cols[8].split(';Name=')[0]
                        expression = ExpValues.get(gffID)
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};\n'.format(cols[0], 'PASA', 'gene', cols[3], cols[4], cols[5], cols[6], cols[7], geneID))
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};note=TPM:{:0.2f};\n'.format(cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], mRNAID, expression))
                elif '_prime_UTR' in cols[2]:
                    utrID = gffID.split('.utr')[0]
                    if utrID in extractList:
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};\n'.format(cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
                elif 'exon' in cols[2]:
                    exonID = gffID.split('.exon')[0]
                    if exonID in extractList:
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};\n'.format(cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
                elif 'CDS' in cols[2]:
                    cdsID = gffID.split('cds.')[-1]
                    if cdsID in extractList:
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};\n'.format(cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
    lib.log.info('Wrote {:,} transcripts derived from {:,} protein coding loci.'.format(len(extractList),len(bestModels)))     
                                
def GFF2tblCombined(evm, genome, trnascan, prefix, genenumber, justify, SeqCenter, SeqRefNum, tblout):
    from collections import OrderedDict
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    '''
    function to take GFF3 annotation to produce a GBK tbl file, support multiple transcripts per locus.
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['start'])
    
    #nested function to do some sorting and return values
    def _locusSort(input):
        #input is a list of dictionaries
        sortedList = []
        startPos = None
        stopPos = None
        pStart = True
        pStop = True
        sortedInput = sorted(input, key=lambda x: x['TPM'], reverse=True)
        gStrand = sortedInput[0]['strand']     
        for i in sortedInput:
            if not i['proper_start']:
                pStart = False
            if not i['proper_stop']:
                pStop = False
            if not startPos:
                startPos = i['start']
            else:
                if i['start'] < startPos:
                    startPos = i['start']               
            if not stopPos:
                stopPos = i['end']
            else:
                if i['end'] > stopPos:
                    stopPos = i['end']
        return sortedInput, startPos, stopPos, pStart, pStop, gStrand
    
    #make sure genenumber is integer 
    genenumber = int(genenumber)
    
    #generate genome length dictionary used for feature tbl generation
    scaffLen = {}
    with open(genome, 'rU') as fastain:
        for record in SeqIO.parse(fastain, 'fasta'):
            if not record.id in scaffLen:
                scaffLen[record.id] = len(record.seq)
    #setup interlap database for genes on each chromosome
    gene_inter = defaultdict(InterLap)
    Genes = {}
    Transcripts = {}
    idParent = {}
    with open(evm, 'rU') as infile:
        for line in infile:
            if line.startswith('\n') or line.startswith('#'):
                continue
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            start = int(start)
            end = int(end)
            if feature == 'gene':
                ID = attributes.split(';')[0].replace('ID=', '')
                if not ID in Genes:
                    Genes[ID] = {'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'transcript_ids': []}
                    gene_inter[contig].add((start, end, strand, ID))   
                else:
                    if not args.alt_transcripts == '1':
                        lib.log.debug("Duplicate Gene IDs found: %s" % ID)      
            else:
                Product = 'hypothetical protein'
                info = attributes.split(';')
                ID,Parent,Note = (None,)*3
                for x in info:
                    if x.startswith('ID='):
                        ID = x.replace('ID=', '')
                    elif x.startswith('Parent='):
                        Parent = x.replace('Parent=', '')
                    elif x.startswith('note='):
                        Note = x.replace('note=', '')
                    elif x.startswith('product='):
                        Product = x.replace('product=', '')
                if not ID or not Parent:
                    lib.log.error("Error, can't find ID or Parent. Malformed GFF file.")
                    sys.exit(1)
                if feature == 'mRNA': #add to Transcript dictionary
                    if not ID in Transcripts:
                        Transcripts[ID] = {'type': 'mRNA', 'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'mRNA': [], 'CDS': [], 'proper_start': False, 'proper_stop': False, 'phase': [], 'product': Product, 'note': None, 'TPM': None, 'codon_start': ''}
                        Genes[Parent]['transcript_ids'].append(ID)
                    if Note:
                        Transcripts[ID]['note'] = Note
                        if 'TPM' in Note:
                            expval = Note.split(':')[-1]
                            Transcripts[ID]['TPM'] = float(expval)
                elif feature == 'exon':
                    if Parent in Transcripts:
                        Transcripts[Parent]['mRNA'].append((start,end))
                elif feature == 'CDS':
                    if Parent in Transcripts:
                        Transcripts[Parent]['CDS'].append((start, end))
                        Transcripts[Parent]['phase'].append(phase)
    #now load tRNA predictions
    with open(trnascan, 'rU') as trnain:
        for line in trnain:
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            start = int(start)
            end = int(end)
            if feature == 'gene':
                ID = attributes.split(';')[0].replace('ID=', '')
                if not ID in Genes:
                    gene_inter[contig].add((start, end, strand, ID))
                    Genes[ID] = {'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'transcript_ids': []}
            elif feature == 'tRNA':
                info = attributes.split(';')
                ID = info[0].replace('ID=', '')
                Parent = info[1].replace('Parent=', '')
                Product = info[2].replace('product=', '')
                if Product == 'tRNA-OTHER':
                    Product = 'tRNA-Xxx'
                Genes[Parent]['product'] = Product
                Genes[Parent]['transcript_ids'].append(ID)
                if not ID in Transcripts:
                    Transcripts[ID] = {'type': 'tRNA', 'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'mRNA': [], 'CDS': [], 'proper_start': False, 'proper_stop': False, 'phase': [], 'product': Product, 'note': None, 'TPM': None, 'codon_start': ''}
            elif feature == 'exon':
                info = attributes.split(';')
                ID = info[0].replace('ID=', '')
                Parent = info[1].replace('Parent=', '')
                if Parent in Transcripts:
                    Transcripts[Parent]['mRNA'].append((start, end))
                    Transcripts[Parent]['CDS'].append((start, end))
                    Transcripts[Parent]['phase'].append('0')
    #now sort dictionary by contig and location, rename using prefix
    sGenes = sorted(Genes.iteritems(), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    scaff2genes = {}
    SeqRecords = SeqIO.index(genome, 'fasta')
    inter = defaultdict(InterLap)
    skipList = []
    dropped = 0
    keeper = 0
    tooShort = 0
    internalStop = 0
    lib.log.info("Renaming gene models and filtering for any that are completely contained (overlapping).")
    for k,v in sortedGenes.items():
        GoodModel = True
        #check if gene model completely contained inside another one on same strand
        if args.alt_transcripts == '1':
            loc = sorted([v['start'],v['end']])
            if loc in gene_inter[v['contig']]:
                for hit in list(gene_inter[v['contig']].find(loc)):
                    if hit[3] != k and hit[2] == v['strand']:  #same strand but diff gene
                        sortedhit = sorted([hit[0],hit[1]])
                        if loc[0] >= sortedhit[0] and loc[1] <= sortedhit[1]: #then this gene is fully contained, skip it
                            #if two gene models have exact same start stop they will both be removed, not really what I want, so run check
                            if loc[0] == sortedhit[0] and loc[1] == sortedhit[1]: #exact same, then choose which has higher TPM
                                if k in ExpressionValues and hit[3] in ExpressionValues:
                                    currExp = ExpressionValues.get(k)
                                    oldExp = ExpressionValues.get(hit[3])
                                    if currExp < oldExp:
                                        GoodModel = False
                                    else:
                                        skipList.append(hit[3])
                            else:
                                GoodModel = False
        if not GoodModel:
            dropped += 1
            continue
        #rename gene locus here
        keeper += 1  
        #renaming scheme, leave if startswith locustag
        if k.startswith('novel_gene') or k.startswith('temp_gene') or k.startswith('split_gene'):
            genenumber += 1
            locusTag = prefix+str(genenumber).zfill(justify)            
        elif k.startswith(prefix) and '_'+prefix in k: #means models were merged
            locusTag = k.split('_'+prefix)[0]
        else:
            locusTag = k   
        #translate to protein space, drop if less than minimum
        #translate to protein sequence, construct cDNA Seq object, translate
        #if passes then add to final output dictionary
        #get orientation and order exons/CDS
        for y in v['transcript_ids']:
            if not y in Transcripts:
                continue
            #set some variables
            protSeq = None
            model = Transcripts.get(y)
            #sort the exons and CDSs
            if model['strand'] == '+':
                sortedExons = sorted(model['mRNA'], key=lambda tup: tup[0])
                sortedCDS = sorted(model['CDS'], key=lambda tup: tup[0])
            else:
                sortedExons = sorted(model['mRNA'], key=lambda tup: tup[0], reverse=True)
                sortedCDS = sorted(model['CDS'], key=lambda tup: tup[0], reverse=True)    
            #now translate
            if model['strand'] == '+' and model['type'] == 'mRNA':
                cdsSeq = ''
                for s in sortedCDS:
                    singleCDS = SeqRecords[model['contig']][s[0]-1:s[1]]
                    cdsSeq += (str(singleCDS.seq))
                mySeq = Seq(cdsSeq, IUPAC.ambiguous_dna)
                protSeq = mySeq.translate(cds=False, table=1)
            elif model['strand'] == '-' and model['type'] == 'mRNA':
                cdsSeq = ''
                for s in sorted(sortedCDS, key=lambda tup: tup[0]):
                    singleCDS = SeqRecords[model['contig']][s[0]-1:s[1]]
                    cdsSeq += (str(singleCDS.seq))
                mySeq = Seq(cdsSeq, IUPAC.ambiguous_dna)
                mySeq = mySeq.reverse_complement()
                protSeq = mySeq.translate(cds=False, table=1)
            if protSeq and len(protSeq) - 1 < 50:
                tooShort += 1
                continue
            if protSeq and '*' in protSeq[:-1]:
                internalStop += 1
                continue
            #add the properstart/stop info
            if model['type'] == 'mRNA':
                if protSeq.endswith('*'):
                    protStop = True
                else:
                    protStop = False
                if protSeq.startswith('M'):
                    protStart = True
                else:
                    protStart = False
                model['proper_start'] = protStart
                model['proper_stop'] = protStop             
        
                #get the codon_start by getting first CDS phase + 1
                indexStart = [x for x, y in enumerate(model['CDS']) if y[0] == sortedCDS[0][0]]
                codon_start = int(model['phase'][indexStart[0]]) + 1
                model['codon_start'] = codon_start
                
            model['mRNA'] = sortedExons
            model['CDS'] = sortedCDS
            
            #add transcript to output dictionary
            if not locusTag in renamedGenes:
                renamedGenes[locusTag] = [model]
            else:
                renamedGenes[locusTag].append(model)
        
        #add the locus to scaff2genes so everything stays in order
        if locusTag in renamedGenes:
            if not v['contig'] in scaff2genes:
                scaff2genes[v['contig']] = [locusTag]
            else:
                scaff2genes[v['contig']].append(locusTag)  
            
    lib.log.info('Writing {:,} loci to TBL format: dropped {:,} overlapping, {:,} too short, and {:,} frameshift gene models'.format(len(renamedGenes),dropped,tooShort,internalStop))
    #now have scaffolds dict and gene dict, loop through scaff dict printing tbl
    with open(tblout, 'w') as tbl:
        for k,v in natsorted(scaff2genes.items()):
            tbl.write('>Feature %s\n' % k)
            tbl.write('1\t%s\tREFERENCE\n' % scaffLen.get(k))
            tbl.write('\t\t\t%s\t%s\n' % (SeqCenter, SeqRefNum))
            for genes in v: #now loop through each gene on the scaffold
                if genes in skipList:
                    continue
                transcriptInfo = renamedGenes.get(genes) #this is now a list of dictionaries
                #need a method to sort this list based on TPM value in each dictionary
                #also need to see if any of these transcripts are partial, as that has to be
                #annotated in the gene portion of TBL as well (at least I think it does)
                sorted_transcriptInfo, gene_start, gene_stop, propStart, propStop, geneStrand = _locusSort(transcriptInfo)
                #check for partial models
                if geneStrand == '+':
                    if propStart:
                        ps = ''
                    else:
                        ps = '<'
                    if propStop:
                        pss = ''
                    else:
                        pss = '>'
                    tbl.write('%s%i\t%s%i\tgene\n' % (ps, gene_start, pss, gene_stop))
                    tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                else:
                    if propStart:
                        ps = ''
                    else:
                        ps = '>'
                    if propStop:
                        pss = ''
                    else:
                        pss = '<'
                    tbl.write('%s%i\t%s%i\tgene\n' % (pss, gene_stop, ps, gene_start))
                    tbl.write('\t\t\tlocus_tag\t%s\n' % genes)                                         

                #now write the mRNA and CDSs
                partialStart,partialStop = ('',)*2
                for tran_num,geneInfo in enumerate(sorted_transcriptInfo):
                    if geneInfo['type'] == 'mRNA':
                        if geneInfo['strand'] == '+':
                            for num, exon in enumerate(geneInfo['mRNA']):
                                if num == 0 and num == len(geneInfo['mRNA']) - 1: #single exon, so slightly differnt method
                                    tbl.write('%s%s\t%s%s\tmRNA\n' % (ps, exon[0], pss, exon[1]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tmRNA\n' % (ps, exon[0], exon[1]))
                                elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                    tbl.write('%s\t%s%s\n' % (exon[0], pss, exon[1]))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,tran_num+1))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,tran_num+1))
                            if geneInfo['TPM']:
                                tbl.write('\t\t\tnote\tKallisto TPM expression: {:.2f}\n'.format(geneInfo['TPM']))
                            for num, cds in enumerate(geneInfo['CDS']):
                                if num == 0 and num == len(geneInfo['CDS']) - 1: #single exon, so slightly differnt method
                                    tbl.write('%s%s\t%s%s\tCDS\n' % (ps, cds[0], pss, cds[1]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tCDS\n' % (ps, cds[0], cds[1]))
                                elif num == len(geneInfo['CDS']) - 1: #this is last one
                                    tbl.write('%s\t%s%s\n' % (cds[0], pss, cds[1]))
                                else:
                                    tbl.write('%s\t%s\n' % (cds[0], cds[1]))
                            tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'])
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,tran_num+1))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,tran_num+1))
                            if geneInfo['TPM']:
                                tbl.write('\t\t\tnote\tKallisto TPM expression: {:.2f}\n'.format(geneInfo['TPM']))                                   
                        else: #means this is on crick strand                         
                            for num, exon in enumerate(geneInfo['mRNA']):
                                if num == 0 and num == len(geneInfo['mRNA']) - 1: #single exon, so slightly differnt method
                                    tbl.write('%s%s\t%s%s\tmRNA\n' % (pss, exon[1], ps, exon[0]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tmRNA\n' % (pss, exon[1], exon[0]))
                                elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                    tbl.write('%s\t%s%s\n' % (exon[1], ps, exon[0]))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[1], exon[0]))                 
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,tran_num+1))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,tran_num+1))
                            if geneInfo['TPM']:
                                tbl.write('\t\t\tnote\tKallisto TPM expression: {:.2f}\n'.format(geneInfo['TPM']))
                            for num, cds in enumerate(geneInfo['CDS']):
                                if num == 0 and num == len(geneInfo['CDS']) - 1: #single exon, so slightly differnt method
                                    tbl.write('%s%s\t%s%s\tCDS\n' % (pss, cds[1], ps, cds[0]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tCDS\n' % (pss, cds[1], cds[0]))
                                elif num == len(geneInfo['CDS']) - 1: #this is last one
                                    tbl.write('%s\t%s%s\n' % (cds[1], ps, cds[0]))
                                else:
                                    tbl.write('%s\t%s\n' % (cds[1], cds[0]))
                            tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'])
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,tran_num+1))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,tran_num+1))
                            if geneInfo['TPM']:
                                tbl.write('\t\t\tnote\tKallisto TPM expression: {:.2f}\n'.format(geneInfo['TPM']))
                    elif geneInfo['type'] == 'tRNA':
                        if geneInfo['strand'] == '+':
                            for num, exon in enumerate(geneInfo['mRNA']):
                                if num == 0:
                                    tbl.write('<%s\t>%s\ttRNA\n' % (exon[0], exon[1]))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tbl.write('\t\t\tpseudo\n')
                            if geneInfo['note'] and geneInfo['note'] != '':
                                tbl.write('\t\t\tnote\t%s\n' % geneInfo['note'])                                    
                        else:
                            for num, exon in enumerate(geneInfo['mRNA']):
                                if num == 0:
                                    tbl.write('<%s\t>%s\ttRNA\n' % (exon[1], exon[0]))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[1], exon[0]))
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tbl.write('\t\t\tpseudo\n')
                            if geneInfo['note'] and geneInfo['note'] != '':
                                tbl.write('\t\t\tnote\t%s\n' % geneInfo['note']) 

def gbk2interlap(input):
    '''
    function to parse GBK file, construct scaffold/gene interlap dictionary, transcript
    dictionary, and protein dictionary
    '''
    inter = defaultdict(InterLap)
    Transcripts = {}
    Proteins = {}
    Genes = {}
    with open(input, 'rU') as filein:
        for record in SeqIO.parse(filein, 'genbank'):
            Contig = record.id
            for f in record.features:
                #reset for every feature
                ID,start,end,strand,num_parts,exons,cds,protID,transcriptID,protSeq,locusTag,Parent,genenum = (None,)*13
                if not f.type in ['gene', 'mRNA', 'tRNA', 'CDS']:
                    continue
                #get info from features
                locusTag, ID, Parent = lib.getID(f, f.type)
                start = int(f.location.nofuzzy_start)
                end = int(f.location.nofuzzy_end)
                strand = f.location.strand
                if strand == 1:
                    direction = '+'
                else:
                    direction = '-'
                num_parts = len(f.location.parts)           
                if f.type == 'gene':
                    inter[Contig].add((start,end,locusTag))
                elif f.type == 'mRNA' or f.type == 'tRNA':
                    exonTuples = []
                    if num_parts < 2: #only single exon
                        exonTuples.append((int(start),int(end)))
                    else: #more than 1 exon, so loop through
                        for i in range(0, num_parts):
                            ex_start = f.location.parts[i].nofuzzy_start
                            ex_end = f.location.parts[i].nofuzzy_end
                            exonTuples.append((int(ex_start),int(ex_end)))
                    #now we want to sort the positions I think...
                    sortedExons = sorted(exonTuples, key=lambda tup: tup[0])
                    #gbk file has no transcript_id field, so I guess make one up
                    if not locusTag in Genes:
                        if ID:
                            transcriptID = ID
                        else:
                            transcriptID = locusTag+'-T1'
                        Genes[locusTag] = {'transcript_id': [transcriptID], 'protein_id': [], 'type': f.type}
                    else:
                        if ID:
                            transcriptID = ID
                        else:
                            genenum = len(Genes[locusTag]['transcript_id'])+1
                            transcriptID = locusTag+'-T'+str(genenum)
                        Genes[locusTag]['transcript_id'].append(transcriptID)
                    #now add relevent info to Transcripts dictionary
                    if not transcriptID in Transcripts:
                        Transcripts[transcriptID] = {'mRNA': sortedExons, 'strand': direction, 'contig': Contig, 'location': (start, end)}
                elif f.type == 'CDS':
                    protSeq = f.qualifiers['translation'][0]
                    if not locusTag in Genes:
                        if ID:
                            protID = ID
                        else:
                            protID = locusTag+'-1'
                        Genes[locusTag] = {'protein_id': [protID], 'transcript_id': [], 'type': None}
                    else:
                        if ID:
                            protID = ID
                        else:
                            genenum = len(Genes[locusTag]['protein_id'])+1
                            protID = locusTag+'-T'+str(genenum)
                        Genes[locusTag]['protein_id'].append(protID)
                    cdsTuples = []
                    if num_parts < 2: #only single CDS
                        cdsTuples.append((int(start),int(end)))
                    else:
                        for i in range(0, num_parts):
                            ex_start = f.location.parts[i].nofuzzy_start
                            ex_end = f.location.parts[i].nofuzzy_end
                            cdsTuples.append((int(ex_start),int(ex_end)))
                    sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0])
                    if not protID in Proteins:
                        Proteins[protID] = {'CDS': sortedCDS, 'strand': direction, 'contig': Contig, 'location': (start, end), 'seq': protSeq}
    return inter, Genes, Transcripts, Proteins
    
def merge_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z
    
def matchCDS2mRNA_AED(locus, genes, cds, mrna):
    '''
    since there is currently no way to get transcript_id in genebank flat files, need to be able
    to determine which mRNA feature corresponds to which CDS feature if there are multiple
    transcripts per locus. Return a list of combined transcript/cds dictionary.
    '''
    #locus is a string, transcripts/proteins be found in genes dictionary
    #exon and cds strings are in the cds and mrna dictionaries respectivly
    #want to return a list of combined dictionaries with correct CDS and mRNA mapped
    #correct IDs are in the protein_id field of the genes dictionary
    combinedList = []
    #get the data out of dictionaries for the locus
    mrnaLoci = genes[locus]['transcript_id']
    cdsLoci = genes[locus]['protein_id']
    #if no mrnaLoci, then it is a feature that script is ignoring for now
    if len(mrnaLoci) < 1:
        return False
    #if no cds hits, then it is a tRNA model, no reason to pair up, just return None type for all cds
    if len(cdsLoci) < 1:
        combinedList = []
        for i in mrnaLoci:
            combined = merge_dicts({'CDS': '', 'strand': mrna[i]['strand'], 'contig': mrna[i]['contig'], 'location': mrna[i]['location'], 'seq': '', 'protein_id': '', 'transcript_id': i}, mrna[i])
            combinedList.append(combined)
        return combinedList
    else:
        #check for multiple transcripts, if not, combine transcripts/proteins and move on
        if len(mrnaLoci) > 1:
            #want to loop through the CDS data, looking for compatible exons, calculate AED for each, take best hit?
            seen = []
            bestMatch = {}
            for prot in cdsLoci:
                sortCDS = sorted(cds[prot]['CDS'], key=lambda tup: tup[0])
                aeds = []
                for transcript in mrnaLoci:
                    if not transcript in seen:
                        sortMRNA = sorted(mrna[transcript]['mRNA'], key=lambda tup: tup[0])
                        if len(sortMRNA) >= len(sortCDS):
                            calcAED = getAED(sortCDS, sortMRNA)
                            aeds.append((prot, transcript,float(calcAED)))
                if len(aeds) == 1: #only a single hit, combine and add to seen list
                    if not prot in bestMatch:
                        bestMatch[prot] = transcript
                    seen.append(transcript)
                elif len(aeds) > 1:
                    sortedAED = sorted(aeds, key=lambda tup: tup[2])
                    if not prot in bestMatch:
                        bestMatch[prot] = sortedAED[0][1]
                        seen.append(sortedAED[0][1])
                else:
                    lib.log.error("Multiple transcript CDS not matched to mRNA for %s" % prot)
            for k,v in bestMatch.items():
                combined = merge_dicts(cds[k], mrna[v])
                combined['protein_id'] = prot
                combined['transcript_id'] = transcript
                combinedList.append(combined)
        else:
            combined = merge_dicts(cds[cdsLoci[0]], mrna[mrnaLoci[0]])
            combined['protein_id'] = cdsLoci[0]
            combined['transcript_id'] = mrnaLoci[0]
            combinedList.append(combined)
    return combinedList

def matchCDS2mRNA(locus, genes, cds, mrna):
    '''
    simple function that assumes the GBK mRNA/CDS features are in order, i.e. when we loop through
    each feature, we assume the first mRNA we see corresponds to the first CDS we see.
    '''

    combinedList = []
    #get the data out of dictionaries for the locus
    mrnaLoci = genes[locus]['transcript_id']
    cdsLoci = genes[locus]['protein_id']
    #if no mrnaLoci, then it is a feature that script is ignoring for now
    if len(mrnaLoci) < 1:
        return False
    #if no cds hits, then it is a tRNA model, no reason to pair up, just return None type for all cds
    if len(cdsLoci) < 1:
        combinedList = []
        for i in mrnaLoci:
            combined = merge_dicts({'CDS': '', 'strand': mrna[i]['strand'], 'contig': mrna[i]['contig'], 'location': mrna[i]['location'], 'seq': '', 'protein_id': '', 'transcript_id': i}, mrna[i])
            combinedList.append(combined)
        return combinedList
    else:
        #check for multiple transcripts, if not, combine transcripts/proteins and move on
        if len(mrnaLoci) > 1:
            if not len(mrnaLoci) == len(cdsLoci):
                print("Error: number of transcripts does not equal coding sequences")
                print locus, mrnaLoci, cdsLoci
                sys.exit(1)
            for i,loc in enumerate(mrnaLoci):
                mRNAdict = mrna[loc]
                CDSdict = cds[cdsLoci[i]]           
                combined = merge_dicts(CDSdict, mRNAdict)
                combined['protein_id'] = cdsLoci[i]
                combined['transcript_id'] = loc
                combinedList.append(combined)
        else:
            combined = merge_dicts(cds[cdsLoci[0]], mrna[mrnaLoci[0]])
            combined['protein_id'] = cdsLoci[0]
            combined['transcript_id'] = mrnaLoci[0]
            combinedList.append(combined)
    return combinedList


def message(loc1, loc2, cdsAED, mrnaAED):
    msg = ''
    if not cdsAED or cdsAED == '':
        cds = 0
    else:
        cds = float(cdsAED)
    mrna = float(mrnaAED)
    pos = loc1 == loc2
    if pos and cds == 0 and mrna == 0:
        msg = 'no change'
    elif pos and cds > 0 and mrna == 0:
        msg = 'CDS changed'
    elif pos and cds == 0 and mrna > 0:
        msg = 'mRNA changed'
    elif pos and cds > 0 and mrna > 0:
        msg = 'mRNA and CDS changed'
    elif not pos and cds == 0 and mrna == 0:
        msg = 'locus coordinates changed'
    elif not pos and cds > 0 and mrna == 0:
        msg = 'locus coordinates changed; CDS changed'
    elif not pos and cds == 0 and mrna > 0:
        msg = 'locus coordinates changed; UTR added'
    elif not pos and cds > 0 and mrna > 0:
        msg = 'locus coordinates changed; mRNA and CDS changed'
    return msg

def outputCounter(message, no_change, UTR_added, yardSale, exonChange, modelChangeNotProt):
    if message == 'locus coordinates changed; mRNA and CDS changed' or message == 'locus coordinates changed; CDS changed' or message == 'mRNA and CDS changed' or message == 'CDS changed':
        yardSale += 1  
    elif message == 'no change':
        no_change += 1
    elif message == 'mRNA changed':
        exonChange += 1
    elif message == 'locus coordinates changed; UTR added':
        UTR_added += 1
    elif message == 'locus coordinates changed':
        modelChangeNotProt += 1
    return no_change, UTR_added, yardSale, exonChange, modelChangeNotProt

def compareAnnotations2(old, new, output):
    '''
    function takes two GenBank annotated genomes and compares gene models
    output is a tsv file for each locus and a description of what is different
    can handle multiple transcripts per locus
    '''
    result = {}
    global no_change, UTR_added, yardSale, exonChange, modelChangeNotProt, dropped, added, total_transcripts, total_genes
    no_change, UTR_added, yardSale, exonChange, modelChangeNotProt, dropped, added, total_transcripts, total_genes = (0,)*9
    lib.log.info("Parsing GenBank files...comparing annotation")
    oldInter, oldGenes, oldTranscripts, oldProteins = gbk2interlap(old)
    newInter, newGenes, newTranscripts, newProteins = gbk2interlap(new)
    #do the simple stuff first, find models that were deleted
    for contig in oldInter:
        for gene in oldInter[contig]:
            if not gene in newInter[contig]: #these models are removed
                dropped += 1
                if not gene[2] in oldGenes:
                    continue
                #parse results, get list of dictionaries
                pairedList = matchCDS2mRNA(gene[2], oldGenes, oldProteins, oldTranscripts)
                if not pairedList:
                    continue
                #populate output dictionary with results
                for x in pairedList:
                    if not gene[2] in result:
                        result[gene[2]] = {'contig': x['contig'], 'old_num_transcripts': len(pairedList), 'old_location': x['location'], 'num_transcripts': len(pairedList), 'strand': x['strand'], 'mRNA': [x['mRNA']], 'location': x['location'], 'CDS': [x['CDS']], 'message': 'gene model removed', 'cdsAED': '1.000', 'exonAED': '1.000', 'transcript_id': [x['transcript_id']], 'protein_id':[x['protein_id']], 'seq': [x['seq']]}
                    else:
                        result[gene[2]]['mRNA'].append(x['mRNA'])
                        result[gene[2]]['CDS'].append(x['CDS'])
                        result[gene[2]]['seq'].append(x['seq'])
                        result[gene[2]]['transcript_id'].append(x['transcript_id'])
                        result[gene[2]]['protein_id'].append(x['protein_id'])

    #now go through the updated annotation, comparing to old annot
    for contig in newInter:
        for gene in newInter[contig]:
            if not gene in oldInter[contig]: #means this is a new model, so add it
                added += 1
                total_genes += 1
                if not gene[2] in newGenes:
                    continue
                #parse results, get list of dictionaries
                pairedList = matchCDS2mRNA(gene[2], newGenes, newProteins, newTranscripts)
                if not pairedList:
                    continue
                #populate output dictionary with results
                for x in pairedList:
                    total_transcripts += len(pairedList)
                    if not gene[2] in result:
                        result[gene[2]] = {'contig': x['contig'], 'old_num_transcripts': 0, 'old_location': x['location'], 'num_transcripts': len(pairedList), 'strand': x['strand'], 'mRNA': [x['mRNA']], 'location': x['location'], 'CDS': [x['CDS']], 'message': 'new gene model', 'cdsAED': '0.000', 'exonAED': '0.000', 'transcript_id': [x['transcript_id']], 'protein_id':[x['protein_id']], 'seq': [x['seq']]}
                    else:
                        result[gene[2]]['mRNA'].append(x['mRNA'])
                        result[gene[2]]['CDS'].append(x['CDS'])
                        result[gene[2]]['seq'].append(x['seq'])
                        result[gene[2]]['transcript_id'].append(x['transcript_id'])
                        result[gene[2]]['protein_id'].append(x['protein_id'])

            else: #means this is existing model, and need to do some comparisons
                hitList = list(oldInter[contig].find(gene))
                #there might be some overlapping transcripts, so enforce locus name
                hit = None
                for z in hitList:
                    if gene[2] == z[2]:
                        hit = z
                if not hit:
                    #there is no real hit, so this a new gene
                    pairedList = matchCDS2mRNA(gene[2], newGenes, newProteins, newTranscripts)
                    if not pairedList:
                        continue
                    total_transcripts += len(pairedList)
                    #populate output dictionary with results
                    for x in pairedList:
                        added += 1
                        total_genes += 1
                        if not gene[2] in result:
                            result[gene[2]] = {'contig': x['contig'], 'old_num_transcripts': 0, 'old_location': x['location'], 'num_transcripts': len(pairedList), 'strand': x['strand'], 'mRNA': [x['mRNA']], 'location': x['location'], 'CDS': [x['CDS']], 'message': 'new gene model', 'cdsAED': '0.000', 'exonAED': '0.000', 'transcript_id': [x['transcript_id']], 'protein_id':[x['protein_id']], 'seq': [x['seq']]}
                        else:
                            result[gene[2]]['mRNA'].append(x['mRNA'])
                            result[gene[2]]['CDS'].append(x['CDS'])
                            result[gene[2]]['seq'].append(x['seq'])
                            result[gene[2]]['transcript_id'].append(x['transcript_id'])
                            result[gene[2]]['protein_id'].append(x['protein_id'])                   
                else:
                    #since we may have multiple transcripts from hit as well as new annotation we need to be aware of that
                    #also, tRNA annotations do not exist in Proteins dictionary, so process them differently
                    #get the reference hits, pull out CDS and mRNA for pairwiseAED calculation
                    total_genes += 1
                    #initiation variables that could be empty
                    hitCDSs, CDSs = (None,)*2
                    
                    #get old model information
                    hit_pairedList = matchCDS2mRNA(gene[2], oldGenes, oldProteins, oldTranscripts)
                    if not hit_pairedList: #overlap with a non tRNA or mRNA feature, by default these are dropped....
                        pairedList = matchCDS2mRNA(gene[2], newGenes, newProteins, newTranscripts)
                        if not pairedList:
                            continue
                        #populate output dictionary with results
                        total_transcripts += len(pairedList)
                        for x in pairedList:
                            if not gene[2] in result:
                                result[gene[2]] = {'contig': x['contig'], 'old_num_transcripts': 0, 'num_transcripts': len(pairedList), 'old_location': x['location'], 'strand': x['strand'], 'mRNA': [x['mRNA']], 'location': x['location'], 'CDS': [x['CDS']], 'message': 'new gene model', 'cdsAED': '0.000', 'exonAED': '0.000', 'transcript_id': [x['transcript_id']], 'protein_id':[x['protein_id']], 'seq': [x['seq']]}
                            else:
                                result[gene[2]]['mRNA'].append(x['mRNA'])
                                result[gene[2]]['CDS'].append(x['CDS'])
                                result[gene[2]]['seq'].append(x['seq'])
                                result[gene[2]]['transcript_id'].append(x['transcript_id'])
                                result[gene[2]]['protein_id'].append(x['protein_id'])
                    else:
                        hitmRNAs = [x['mRNA'] for x in hit_pairedList]
                        if oldGenes[gene[2]]['type'] == 'mRNA':
                            hitCDSs = [x['CDS'] for x in hit_pairedList]
                    
                        #get new model information
                        pairedList = matchCDS2mRNA(gene[2], newGenes, newProteins, newTranscripts)
                        mRNAs = [x['mRNA'] for x in pairedList]
                        if newGenes[gene[2]]['type'] == 'mRNA':
                            CDSs = [x['CDS'] for x in pairedList]
                    
                        #calculate AEDs
                        exonAED = pairwiseAED(mRNAs, hitmRNAs)
                        if hitCDSs and CDSs:
                            cdsAED = pairwiseAED(CDSs, hitCDSs)
                        else:
                            cdsAED = '0.000'
                    
                        #determine what happened to gene model
                        msg = message(pairedList[0]['location'], hit_pairedList[0]['location'], cdsAED, exonAED)
                        no_change, UTR_added, yardSale, exonChange, modelChangeNotProt = outputCounter(msg, no_change, UTR_added, yardSale, exonChange, modelChangeNotProt)

                        #populate output dictionary with results
                        total_transcripts += len(pairedList)
                        for x in pairedList:
                            if not gene[2] in result:
                                result[gene[2]] = {'contig': x['contig'], 'old_num_transcripts': len(hit_pairedList), 'old_location': hit_pairedList[0]['location'], 'num_transcripts': len(pairedList), 'strand': x['strand'], 'mRNA': [x['mRNA']], 'location': x['location'], 'CDS': [x['CDS']], 'message': msg, 'cdsAED': cdsAED, 'exonAED': exonAED, 'transcript_id': [x['transcript_id']], 'protein_id':[x['protein_id']], 'seq': [x['seq']]}
                            else:
                                result[gene[2]]['mRNA'].append(x['mRNA'])
                                result[gene[2]]['CDS'].append(x['CDS'])
                                result[gene[2]]['seq'].append(x['seq'])
                                result[gene[2]]['transcript_id'].append(x['transcript_id'])
                                result[gene[2]]['protein_id'].append(x['protein_id'])             
    total_cdsAED = []
    total_exonAED = []               
    with open(output, 'w') as out:
        out.write('Locus_tag\tOrig_Location\tOrig_Num_Transcripts\tContig:start-end\tStrand\tGene_Length\tNum_Transcripts\tmRNA_AED\tCDS_AED\tDescription\n')
        for k,v in natsorted(result.items()):
            start = str(v['location'][0])
            end = str(v['location'][1])
            GeneLength = int(end) - int(start)
            total_cdsAED.append(float(v['cdsAED']))
            total_exonAED.append(float(v['exonAED']))
            out.write('{:}\t{:}:{:}-{:}\t{:}\t{:}:{:}-{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(k, v['contig'],v['old_location'][0],v['old_location'][1],v['old_num_transcripts'],v['contig'],start,end,v['strand'],GeneLength,v['num_transcripts'],v['exonAED'],v['cdsAED'], v['message']))
    Avg_cdsAED = sum(total_cdsAED) / float(len(total_cdsAED))
    Avg_exonAED = sum(total_exonAED) / float(len(total_exonAED))
    #output some simple stats to cmd line
    lib.log.info("Updated annotation complete:\n\
-------------------------------------------------------\n\
Total Gene Models:\t{:,}\n\
Total transcripts:\t{:,}\n\
New Gene Models:\t{:,}\n\
No Change:\t\t{:,}\n\
Update UTRs:\t\t{:,}\n\
Exons Changed:\t\t{:,}\n\
Exons/CDS Changed:\t{:,}\n\
Exon Changed/Prot same:\t{:,}\n\
Dropped Models:\t\t{:,}\n\
CDS AED:\t\t{:.3f}\n\
mRNA AED:\t\t{:.3f}\n\
-------------------------------------------------------".format(total_genes, total_transcripts, added, no_change, UTR_added, exonChange, yardSale, modelChangeNotProt, dropped, Avg_cdsAED, Avg_exonAED))


def pairwiseAED(query, reference):
    import itertools
    '''
    takes a multiple transcripts and sums AED from lowest pairwise comparison
    '''
    AEDsum = []
    pAED = [float(getAED(a,b)) for a, b in itertools.product(query, reference)]
    #split into parts to get lowest AED
    splitAED = [pAED[i:i+len(query)] for i  in range(0, len(pAED), len(query))]
    for pair in splitAED:
        AEDsum.append(min(pair))
    AEDavg = sum(AEDsum) / len(query)
    return '{:.3f}'.format(AEDavg)
     
def getAED(query, reference):
    '''
    function to calcuate annotation edit distance between two mRNA transcript coordinates
    AED = 1 - (SN + SP / 2)
    SN = fraction of ref predicted
    SP = fraction prediction overlapping the ref
    '''
    def _length(listTup):
        len = 0
        for i in listTup:
            l = abs(i[0] - i[1])
            len += l
        return len
    #make sure sorted
    s_query = sorted(query, key=lambda tup: tup[0])
    s_ref = sorted(reference, key=lambda tup: tup[0])
    rLen = _length(s_ref)
    refInterlap = InterLap(s_ref)
    RefPred = 0
    QueryOverlap = 0
    qLen = 0
    for exon in s_query:
        qLen += abs(exon[0] - exon[1])
        if exon in refInterlap: #exon overlaps at least partially with reference
            hit = list(refInterlap.find(exon))
            for h in hit:
                diff = np.subtract(exon,h) #will return array of exon minus hit at each pos
                if diff[0] < 1 and diff[1] > 0: #then query exon covers ref exon
                    QueryOverlap += abs(h[0] - h[1])
                elif diff[0] < 1 and diff[1] < 1: # means query partial covers ref
                    QueryOverlap = abs(h[0] - exon[1])
                elif diff[0] > 0 and diff[1] > 0: #means query partial covers ref
                    QueryOverlap = abs(exon[0] - h[1])
                elif diff[0] > 0 and diff[1] < 1: 
                    QueryOverlap = abs(exon[0] - exon[1])
    #calculate AED
    SP = QueryOverlap / float(qLen)
    SN = QueryOverlap / float(rLen)
    AED = 1 - ((SN + SP) / 2)
    return '{:.3f}'.format(AED)
    
#create folder structure
if os.path.isdir(args.input): #then funannoate folder is passed
    args.out = args.input

if not args.out:
    lib.log.error("No output folder specified, -o, --out.")
    sys.exit(1)
if not os.path.isdir(args.out):
    os.makedirs(args.out)
    os.makedirs(os.path.join(args.out, 'update_misc'))
    os.makedirs(os.path.join(args.out, 'update_results'))
    os.makedirs(os.path.join(args.out, 'logfiles'))
else:
    #make sure subdirectories exist
    dirs = [os.path.join(args.out, 'update_misc'), os.path.join(args.out, 'logfiles'), os.path.join(args.out, 'update_results')]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)

#assign temp directory
tmpdir = os.path.join(args.out, 'update_misc')

#create log file
log_name = os.path.join(args.out, 'logfiles', 'funannotate-update.log')
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print("-------------------------------------------------------")
lib.SystemInfo()

#get version of funannotate
version = lib.get_version()
lib.log.info("Running %s" % version)
   
#do some checks and balances
if not args.PASAHOME:
    try:
        PASA = os.environ["PASAHOME"]
    except KeyError:
        lib.log.error("$PASAHOME environmental variable not found, PASA is not properly configured.  You can use the --PASAHOME argument to specifiy a path at runtime")
        sys.exit(1)
else:
    PASA = args.PASAHOME

if not args.TRINITYHOME:
    try:
        TRINITY = os.environ["TRINITYHOME"]
    except KeyError:
        lib.log.error("$TRINITYHOME environmental variable not found, TRINITY is not properly configured. You can use the --TRINITYHOME argument to specify a path at runtime.")
        sys.exit(1)
else:
    TRINITY = args.TRINITYHOME
        
programs = ['fasta', 'mysql', 'gmap', 'blat', 'tbl2asn', 'hisat2', 'hisat2-build', 'kallisto', 'Trinity', 'bedtools', 'java']
lib.CheckDependencies(programs)

#take care of some preliminary checks
if args.sbt == 'SBT':
    SBT = os.path.join(parentdir, 'lib', 'test.sbt')
    lib.log.info("No NCBI SBT file given, will use default, for NCBI submissions pass one here '--sbt'")
else:
    SBT = args.sbt

#check input, allow for passing the output directory of funannotate, otherwise must be gbk or gbff files
#set read inputs to None, populate as you go
s_reads, l_reads, r_reads, trim_left, trim_right, trim_single, left_norm, right_norm, single_norm, all_reads, trim_reads, norm_reads, GBK, trinity_results, pasaConfigFile = (None,)*15
if os.path.isdir(args.input):
    if os.path.isdir(os.path.join(args.input, 'predict_results')):
        for file in os.listdir(os.path.join(args.input, 'predict_results')):
            if file.endswith('.gbk'):
                GBK = os.path.join(args.input, 'predict_results', file)
    #now lets also check if training folder/files are present, as then can pull all the data you need for update directly
    if os.path.isdir(os.path.join(args.input, 'training')): #then funannotate train has been run, try to get reads, trinity, PASA
        inputDir = os.path.join(args.input, 'training')
        if os.path.isfile(os.path.join(inputDir, 'left.fq.gz')):
            l_reads = os.path.join(inputDir, 'left.fq.gz')
        if os.path.isfile(os.path.join(inputDir, 'right.fq.gz')):
            r_reads = os.path.join(inputDir, 'right.fq.gz')
        if os.path.isfile(os.path.join(inputDir, 'single.fq.gz')):
            s_reads = os.path.join(inputDir, 'single.fq.gz')
        if os.path.isfile(os.path.join(inputDir, 'trimmomatic', 'trimmed_left.fastq.gz')):
            trim_left = os.path.join(inputDir, 'trimmomatic', 'trimmed_left.fastq.gz')
        if os.path.isfile(os.path.join(inputDir, 'trimmomatic', 'trimmed_right.fastq.gz')):
            trim_right = os.path.join(inputDir, 'trimmomatic', 'trimmed_right.fastq.gz')
        if os.path.isfile(os.path.join(inputDir, 'trimmomatic', 'trimmed_single.fastq.gz')):
            trim_single = os.path.join(inputDir, 'trimmomatic', 'trimmed_single.fastq.gz')
        if os.path.isfile(os.path.join(inputDir, 'normalize', 'left.norm.fq')):
            left_norm = os.path.join(inputDir, 'normalize', 'left.norm.fq')
        if os.path.isfile(os.path.join(inputDir, 'normalize', 'right.norm.fq')):
            right_norm = os.path.join(inputDir, 'normalize', 'right.norm.fq')
        if os.path.isfile(os.path.join(inputDir, 'normalize', 'single.norm.fq')):
            single_norm = os.path.join(inputDir, 'normalize', 'single.norm.fq')
        if l_reads or s_reads:
            all_reads = (l_reads, r_reads, s_reads)
        if trim_left or trim_single:
            trim_reads = (trim_left, trim_right, trim_single)
        if left_norm or single_norm:
            norm_reads = (left_norm, right_norm, single_norm)
        if os.path.isfile(os.path.join(inputDir, 'trinity.fasta')):
            trinity_results = os.path.join(inputDir, 'trinity.fasta')
        if os.path.isfile(os.path.join(inputDir, 'pasa', 'alignAssembly.txt')):
            pasaConfigFile = os.path.join(inputDir, 'pasa', 'alignAssembly.txt')
else:
    GBK = args.input

if not GBK:
    lib.log.error("Error in input (-i, --input): pass either funannotate directory or GenBank file")
    sys.exit(1)

#check if RefSeq --> NCBI does not want you to reannotate RefSeq genomes
if lib.checkRefSeq(GBK):
    lib.log.error('%s is a NCBI RefSeq genome, to reannotate please use original submission.' % GBK)
    sys.exit(1)

#split GenBank into parts
#setup output files
gffout = os.path.join(tmpdir, 'genome.gff3')
proteinsout = os.path.join(tmpdir, 'genome.proteins.fa')
trnaout = os.path.join(tmpdir, 'genome.trna.gff3')
fastaout = os.path.join(tmpdir, 'genome.fa')
spliceout = os.path.join(tmpdir, 'genome.ss')
exonout = os.path.join(tmpdir, 'genome.exons')
locustag, genenumber, justify = gbk2pasa(GBK, gffout, trnaout, fastaout, spliceout, exonout, proteinsout)
organism, strain, isolate, accession, WGS_accession, gb_gi, version = lib.getGBKinfo(GBK)
lib.log.info("Reannotating %s, NCBI accession: %s" % (organism, WGS_accession))
lib.log.info("Previous annotation consists of: {:,} protein coding gene models and {:,} non-coding gene models".format(lib.countGFFgenes(gffout), lib.countGFFgenes(trnaout)))

#check if organism/species/isolate passed at command line, if so, overwrite what you detected.
if args.species:
    organism = args.species
if args.strain:
    strain = args.strain
if args.isolate:
    isolate = args.isolate
if strain:
    organism_name = organism+'_'+strain
elif isolate:
    organism_name = organism+'_'+isolate
else:
    organism_name = organism
organism_name = organism_name.replace(' ', '_')

#check input reads
#get absolute paths for reads and concate if there are multiple
if not all_reads:
    if not os.path.isfile(os.path.join(args.out, 'update_misc', 'single.fq.gz')):
        if args.single:
            single_reads = []
            for y in args.single:
                single_reads.append(os.path.abspath(y))
            if single_reads[0].endswith('.gz'):
                ending = '.fq.gz'
            else:
                ending = '.fq'
            s_reads = os.path.join(tmpdir, 'single'+ending)
            if len(single_reads) > 1:
                lib.log.info("Multiple inputs for --single detected, concatenating SE reads")
                concatenateReads(single_reads, s_reads)
            else:
                s_reads = single_reads[0]
    else:
        s_reads = os.path.join(args.out, 'update_misc', 'single.fq.gz')
    
    if not os.path.isfile(os.path.join(args.out, 'update_misc', 'left.fq.gz')) or not os.path.isfile(os.path.join(args.out, 'update_misc', 'right.fq.gz')):
        if args.left and args.right:
            left_reads = []
            for i in args.left:
                left_reads.append(os.path.abspath(i))
            right_reads = []
            for x in args.right:
                right_reads.append(os.path.abspath(x))
            #since I can't get the comma separated input to work through subprocess, lets concatenate reads
            if left_reads[0].endswith('.gz'):
                ending = '.fq.gz'
            else:
                ending = '.fq'
            l_reads = os.path.join(tmpdir, 'left'+ending)
            r_reads = os.path.join(tmpdir, 'right'+ending)
            if len(left_reads) > 1:
                lib.log.info("Multiple inputs for --left and --right detected, concatenating PE reads")
                concatenateReads(left_reads, l_reads)
                concatenateReads(right_reads, r_reads)
            else:
                l_reads = left_reads[0]
                r_reads = right_reads[0]
    else:
        l_reads = os.path.join(args.out, 'update_misc', 'left.fq.gz')
        r_reads = os.path.join(args.out, 'update_misc', 'right.fq.gz')
    
    #get tuple of input reads so you can parse them in downstream tools
    all_reads = (l_reads, r_reads, s_reads)

lib.log.debug(all_reads)
#trimmomatic on reads, first run PE
if not trim_reads:
    if args.no_trimmomatic or args.trinity or args.left_norm or args.single_norm:
        lib.log.info("Trimmomatic will be skipped")
        trim_left = l_reads
        trim_right = r_reads
        trim_single = s_reads
    else:
        #check if they exist already in folder
        if not os.path.isfile(os.path.join(args.out, 'update_misc', 'trimmomatic', 'trimmed_left.fastq.gz')) or not os.path.isfile(os.path.join(args.out, 'update_misc', 'trimmomatic', 'trimmed_right.fastq.gz')):
            if all_reads[0] and all_reads[1]:
                trim_left, trim_right = runTrimmomaticPE(l_reads, r_reads)
            else:
                trim_left, trim_right = (None,)*2
        else:
            trim_left, trim_right = os.path.join(args.out, 'update_misc', 'trimmomatic', 'trimmed_left.fastq.gz'), os.path.join(args.out, 'update_misc', 'trimmomatic', 'trimmed_right.fastq.gz')
        if not os.path.isfile(os.path.join(args.out, 'update_misc', 'trimmomatic', 'trimmed_single.fastq.gz')) and s_reads:
            if all_reads[2]:
                trim_single = runTrimmomaticSE(s_reads)
            else:
                trim_single = None
        else:
            if s_reads:
                trim_single = os.path.join(args.out, 'update_misc', 'trimmomatic', 'trimmed_single.fastq.gz')
            else:
                trim_single = None
    #get tuple of trimmed reads
    trim_reads = (trim_left, trim_right, trim_single)
lib.log.debug(trim_reads)
#check that reads are present and make sure they follow trinity naming conventions, i.e. either illumina default or /1 and /2 to PE reads
for read in trim_reads:
    if read:
        if not os.path.isfile(read):
            lib.log.error("Trimmomatic failed, %s does not exist." % read)
            sys.exit(1)
if trim_reads[0] and trim_reads[1]: #PE reads are passed, lets make sure they have proper naming
    lib.CheckFASTQandFix(trim_reads[0], trim_reads[1]) #if needed to fix they will be fixed in place  

#normalize reads
if not norm_reads:
    if args.no_normalize_reads or args.trinity or args.left_norm or args.single_norm:
        lib.log.info("Read normalization will be skipped")
        if args.left_norm:
            left_norm = args.left_norm
            right_norm = args.right_norm
        else:   
            left_norm = trim_left
            right_norm = trim_right
        if args.single_norm:
            single_norm = args.single_norm
        else:
            single_norm = trim_single
    else:
        #check if exists
        if not os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'left.norm.fq')) or not os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'right.norm.fq')):
            lib.log.info("Running read normalization with Trinity")
            left_norm, right_norm, single_norm = runNormalization(trim_reads, args.memory)
        else:
            left_norm, right_norm = os.path.join(args.out, 'update_misc', 'normalize', 'left.norm.fq'), os.path.join(args.out, 'update_misc', 'normalize', 'right.norm.fq')
            if os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'single.norm.fq')):
                single_norm = os.path.join(args.out, 'update_misc', 'normalize', 'single.norm.fq')
        if not os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'single.norm.fq')) and not trim_left and not trim_right and trim_single:
            lib.log.info("Running read normalization with Trinity")
            left_norm, right_norm, single_norm = runNormalization(trim_reads, args.memory)
        else:
            if os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'single.norm.fq')):
                single_norm = os.path.join(args.out, 'update_misc', 'normalize', 'single.norm.fq')
            else:
                single_norm = None

    norm_reads = (left_norm, right_norm, single_norm)
lib.log.debug(norm_reads)

#now run Trinity with trimmomatic and read normalization
PASA_gff = os.path.join(tmpdir, 'pasa_final.gff3')
if not trinity_results:
    trinity_transcripts = os.path.join(tmpdir, 'trinity.fasta')
    trinity_tmp = os.path.join(tmpdir, 'trinity.tmp')
    if not lib.checkannotations(trinity_tmp):
        if args.trinity:
            shutil.copyfile(os.path.abspath(args.trinity), trinity_tmp)
        else:
            #run trinity genome guided
            runTrinityGG(fastaout, norm_reads, exonout, spliceout, trinity_tmp)
        
    #clip polyA tails
    polyAclip(trinity_tmp, trinity_transcripts)

    if not lib.checkannotations(trinity_tmp):
        lib.log.error("TRINITY step failed, check logfile, exiting")
        sys.exit(1)

    #if RNA is stranded, remove anti-sense transcripts by mapping back reads to transcripts and investigating strandeness
    if args.stranded != 'no' and not args.no_antisense_filter and not args.pasa_gff and not args.single and not lib.checkannotations(PASA_gff):
        trinity_transcripts_backup = trinity_transcripts+'.bak'
        os.rename(trinity_transcripts, trinity_transcripts_backup)  
        removeAntiSense(trinity_transcripts_backup, norm_reads, trinity_transcripts)
else:
    trinity_transcripts = trinity_results
    
#now run PASA steps
if not pasaConfigFile:
    if args.pasa_gff:
        lib.log.info("You passed a --pasa_gff file; are you sure this is a good idea?")
        shutil.copyfile(args.pasa_gff, PASA_gff)
    if not lib.checkannotations(PASA_gff):
        runPASA(fastaout, trinity_transcripts, args.stranded, args.max_intronlen, args.cpus, gffout, organism_name, PASA_gff, args.pasa_config)
else:
    if not lib.checkannotations(PASA_gff):
        runPASA(fastaout, trinity_transcripts, args.stranded, args.max_intronlen, args.cpus, gffout, organism_name, PASA_gff, pasaConfigFile)
    else:
        lib.log.info('Skipping PASA, found existing output: %s' % PASA_gff)
        
#now run Kallisto steps, if mixed PE and SE reads, only PE reads will be used for Kallisto as there isn't a reasonable way to combine them
KallistoAbundance = os.path.join(tmpdir, 'kallisto.tsv')
if args.kallisto:
    lib.log.info("You passed a --kallisto file; are you sure this is a good idea?")
    shutil.copyfile(args.kallisto, KallistoAbundance)
if not lib.checkannotations(KallistoAbundance):
    runKallisto(PASA_gff, fastaout, trim_reads, args.stranded, args.cpus, KallistoAbundance)

#parse Kallisto results with PASA GFF
global ExpressionValues
ExpressionValues = {}
BestModelGFF = os.path.join(tmpdir, 'bestmodels.gff3')
getBestModels(PASA_gff, fastaout, KallistoAbundance, args.alt_transcripts, BestModelGFF)

#make sure tRNA models don't overlap new gene models
cleanTRNA = os.path.join(tmpdir, 'trna.no-overlaps.gff')
cmd = ['bedtools', 'intersect', '-v', '-a', trnaout, '-b', BestModelGFF]
lib.runSubprocess2(cmd, '.', lib.log, cleanTRNA)

#generate tbl file
gagdir = os.path.join(tmpdir, 'tbl2asn')
TBLFile = os.path.join(gagdir, 'genome.tbl')
if not lib.checkannotations(TBLFile):
    if os.path.isdir(gagdir):
        shutil.rmtree(gagdir)
    os.makedirs(gagdir)
    GFF2tblCombined(BestModelGFF, fastaout, cleanTRNA, locustag, genenumber, justify, args.SeqCenter, args.SeqAccession, TBLFile)

#need a function here to clean up the ncbi tbl file if this is a reannotation
#a reannotation would have a WGS_accession, if none, then it is a first pass and not from genbank
if WGS_accession:
    os.rename(os.path.join(tmpdir, 'tbl2asn', 'genome.tbl'), os.path.join(tmpdir, 'tbl2asn', 'genome.tbl.bak'))
    p2g = {}
    if args.p2g: #load into dictionary
        shutil.copyfile(args.p2g, os.path.join(args.out, 'update_results', 'ncbi.p2g'))
        with open(args.p2g, 'rU') as input:
            for line in input:
                cols = line.split('\t')
                if not cols[0] in p2g:
                    p2g[cols[0]] = cols[1]
    with open(os.path.join(tmpdir, 'tbl2asn', 'genome.tbl'), 'w') as outfile:
        with open(os.path.join(tmpdir, 'tbl2asn', 'genome.tbl.bak'), 'rU') as infile:
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


#run tbl2asn in new directory directory
shutil.copyfile(fastaout, os.path.join(gagdir, 'genome.fsa'))
lib.log.info("Converting to Genbank format")
discrep = os.path.join(args.out, 'update_results', organism_name + '.discrepency.report.txt')
if version and WGS_accession: #this would mean it is a GenBank reannotation, so update accordingly. else it is just 1st version.
    rev_version = int(version) + 1
else:
    rev_version = 1
tbl2asn_cmd = lib.runtbl2asn(gagdir, SBT, discrep, organism, isolate, strain, args.tbl2asn, rev_version)

#grab results, populate results output directory
final_fasta = os.path.join(args.out, 'update_results', organism_name + '.scaffolds.fa')
final_gff = os.path.join(args.out, 'update_results', organism_name + '.gff3')
final_gbk = os.path.join(args.out, 'update_results', organism_name + '.gbk')
final_tbl = os.path.join(args.out, 'update_results', organism_name + '.tbl')
final_proteins = os.path.join(args.out, 'update_results', organism_name + '.proteins.fa')
final_transcripts = os.path.join(args.out, 'update_results', organism_name + '.transcripts.fa')
final_validation = os.path.join(args.out, 'update_results', organism_name+'.validation.txt')
final_error = os.path.join(args.out, 'update_results', organism_name+'.error.summary.txt')
final_fixes = os.path.join(args.out, 'update_results', organism_name+'.models-need-fixing.txt')

#retrieve files/reorganize
shutil.copyfile(os.path.join(gagdir, 'genome.gbf'), final_gbk)
shutil.copyfile(os.path.join(gagdir, 'genome.tbl'), final_tbl)
shutil.copyfile(os.path.join(gagdir, 'genome.val'), final_validation)
shutil.copyfile(os.path.join(gagdir, 'errorsummary.val'), final_error)
lib.log.info("Collecting final annotation files")
lib.gb2allout(final_gbk, final_gff, final_proteins, final_transcripts, final_fasta)

#since no place to write the WGS accession to, output to a file for reading by funannotate annotate
with open(os.path.join(args.out, 'update_results', 'WGS_accession.txt'), 'w') as out:
    out.write('%s\n' % WGS_accession)

#last step will be to compare old gene models versus updated ones, outputting a tsv file describing the changes
Changes = os.path.join(args.out, 'update_results', organism_name + '.pasa-reannotation.changes.txt')
compareAnnotations2(GBK, final_gbk, Changes)

lib.log.info("Funannotate update is finished, output files are in the %s/update_results folder" % (args.out))
errors = lib.ncbiCheckErrors(final_error, final_validation, locustag, final_fixes)
if errors > 0:
    print('-------------------------------------------------------')
    lib.log.info("Manually edit the tbl file %s, then run:\n\nfunannotate fix -i %s -t %s\n" % (final_tbl, final_gbk, final_tbl))
    lib.log.info("After the problematic gene models are fixed, you can proceed with functional annotation.")

lib.log.info("Your next step might be functional annotation, suggested commands:\n\
-------------------------------------------------------\n\
Run InterProScan (Docker required): \nfunannotate iprscan -i {:} -m docker -c {:}\n\n\
Run antiSMASH: \nfunannotate remote -i {:} -m antismash -e youremail@server.edu\n\n\
Annotate Genome: \nfunannotate annotate -i {:} --cpus {:} --sbt yourSBTfile.txt\n\
-------------------------------------------------------\n\
            ".format(args.out, \
            args.cpus, \
            args.out, \
            args.out, \
            args.cpus))
sys.exit(1)
