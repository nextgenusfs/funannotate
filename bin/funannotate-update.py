#!/usr/bin/env python

import sys, os, subprocess, inspect, shutil, argparse, fnmatch
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib
from natsort import natsorted

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
parser.add_argument('--name', help='Shortname for genes, perhaps assigned by NCBI, eg. VC83')
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
                                    current_phase = 0
                                    gff.write("%s\tGenBank\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n" % (chr, start, end, strand, ID))
                                    gff.write("%s\tGenBank\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, ID, product))
                                    if num_exons < 2: #only a single exon
                                        ex_start = str(f.location.nofuzzy_start + 1)
                                        ex_end = str(f.location.nofuzzy_end)
                                        gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon1;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ID))
                                        gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, current_phase, ID, ID))
                                    else: #more than 1 exon, so parts sub_features
                                        splices = []
                                        for i in range(0,num_exons):
                                            ex_start = str(f.location.parts[i].nofuzzy_start + 1)
                                            ex_end = str(f.location.parts[i].nofuzzy_end)
                                            ex_num = i + 1
                                            gff.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon%i;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, ID, ex_num, ID))
                                            gff.write("%s\tGenBank\tCDS\t%s\t%s\t.\t%s\t%i\tID=%s-T1.cds;Parent=%s-T1\n" % (chr, ex_start, ex_end, strand, current_phase, ID, ID))
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
                                        trna.write("%s\tGenBank\ttRNA\t%s\t%s\t.\t%s\t.\tID=%s-T1;Parent=%s;product=%s\n" % (chr, start, end, strand, ID, ID, product))
                                        trna.write("%s\tGenBank\texon\t%s\t%s\t.\t%s\t.\tID=%s-T1.exon1;Parent=%s-T1\n" % (chr, start, end, strand, ID, ID))
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
    tag, count = lastTag.split('_')
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
            cmd = [os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl'), '-c', os.path.abspath(alignConfig), '-r', '-C', '-R', '-g', os.path.abspath(genome), '--ALIGNERS', 'blat,gmap', '-t', os.path.abspath(transcripts), '--stringent_alignment_overlap', args.pasa_alignment_overlap, '--TRANSDECODER', '--MAX_INTRON_LENGTH', str(intronlen), '--CPU', str(cpus)]
        else:
            cmd = [os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl'), '-c', os.path.abspath(alignConfig), '-r', '-C', '-R', '-g', os.path.abspath(genome), '--ALIGNERS', 'blat,gmap', '-t', os.path.abspath(transcripts), '--transcribed_is_aligned_orient', '--stringent_alignment_overlap', args.pasa_alignment_overlap, '--TRANSDECODER', '--MAX_INTRON_LENGTH', str(intronlen), '--CPU', str(cpus)]
        lib.runSubprocess(cmd, folder, lib.log)
        
    #generate comparison template file
    with open(annotConfig, 'w') as config2:
        with open(os.path.join(PASA, 'pasa_conf', 'pasa.annotationCompare.Template.txt'), 'rU') as template2:
            for line in template2:
                line = line.replace('<__MYSQLDB__>', DataBaseName)
                config2.write(line)
    
    #now run Annotation comparisons
    lib.log.info("Running PASA annotation comparison step 1")
    cmd = [os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl'), '-c', os.path.abspath(annotConfig), '-g', os.path.abspath(genome), '-t', os.path.abspath(transcripts), '-A', '-L', '--annots_gff3', os.path.abspath(previousGFF), '--CPU', str(cpus)]
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
    cmd = [os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl'), '-c', os.path.abspath(annotConfig), '-g', os.path.abspath(genome), '-t', os.path.abspath(transcripts), '-A', '-L', '--annots_gff3', os.path.abspath(round1GFF), '--CPU', str(cpus)]
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
    #since mRNA is unique, parse the transcript file which has mRNAID genID in header
    with open(PASAtranscripts, 'rU') as transin:
        for line in transin:
            if line.startswith('>'):
                line = line.replace('>', '')
                cols = line.split(' ')
                mRNAID = cols[0]
                geneID = cols[1]
                if not mRNAID in mRNADict:
                    mRNADict[mRNAID] = geneID
    #now make new tsv file with #mRNAID geneID  TPM
    with open(output, 'w') as outfile:
        outfile.write("#mRNA-ID\tgene-ID\tTPM\n")
        with open(os.path.join(folder, 'kallisto', 'abundance.tsv'), 'rU') as infile:
            for line in infile:
                if line.startswith('targed_id'):
                    continue
                line = line.replace('\n', '')
                cols = line.split('\t')
                geneID = mRNADict.get(cols[0])
                outfile.write('%s\t%s\t%s\n' % (cols[0], geneID, cols[4]))

    
def getBestModel(input, fasta, abundances, outfile):
    #enforce the minimum protein length, generate list of IDs to ignore
    lib.log.info("Parsing Kallisto results and selecting best gene model")
    ignore = []
    allproteins = os.path.join(tmpdir, 'all_proteins.fasta')
    cmd = [os.path.join(PASA, 'misc_utilities', 'gff3_file_to_proteins.pl'), input, fasta, 'prot']
    lib.runSubprocess2(cmd, '.', lib.log, allproteins)
    with open(allproteins, 'rU') as protfile:
        for record in SeqIO.parse(protfile, 'fasta'):
            Seq = str(record.seq)
            Seq = Seq[:-1] # remove stop codon
            if len(record.seq) < args.min_protlen:
                ignore.append(record.id)
    if len(ignore) > 0:
        lib.log.debug("%i model(s) less than minimum protein length (%i), dropping" % (len(ignore), args.min_protlen))
    else:
        lib.log.debug("0 models less than minimum protein length")           
    #now parse the output, get list of bestModels
    bestModels = {}
    with open(abundances, 'rU') as tpms:
        for line in tpms:
            line = line.replace('\n', '')
            if line.startswith('#') or line.startswith('target_id'):
                continue
            cols = line.split('\t')
            #check ignore list
            if cols[0] in ignore:
                continue
            if not cols[1] in bestModels:
                bestModels[cols[1]] = (cols[0], cols[2])
            else:
                oldTPM = bestModels.get(cols[1])[1]
                if float(cols[2]) > float(oldTPM):
                    bestModels[cols[1]] = (cols[0], cols[2])
    bestModelFile = os.path.join(tmpdir, 'best_hits.txt')
    extractList = []
    with open(bestModelFile, 'w') as bmout:
        for k,v in natsorted(bestModels.items()):
            extractList.append(v[0])
            bmout.write("%s\t%s\t%s\n" % (v[0], k, v[1]))
    extractList = set(extractList)
    with open(outfile, 'w') as output:
        output.write("##gff-version 3\n")
        with open(input, 'rU') as gff:
            for line in gff:
                if line.startswith("#") or line.startswith('\n'):
                    continue
                line = line.replace('\n', '')
                cols = line.split('\t')
                gffID = cols[8].split(';Parent')[0].replace('ID=', '')
                if 'gene' in cols[2]:
                    continue
                elif 'mRNA' in cols[2]:
                    if gffID in extractList:
                        geneID = cols[8].split(';Name=')[0]
                        geneID = 'ID='+ geneID.split(';Parent=')[-1]
                        mRNAID = cols[8].split(';Name=')[0]
                        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cols[0], 'PASA', 'gene', cols[3], cols[4], cols[5], cols[6], cols[7], geneID))
                        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], mRNAID))
                elif '_prime_UTR' in cols[2]:
                    utrID = gffID.split('.utr')[0]
                    if utrID in extractList:
                        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
                elif 'exon' in cols[2]:
                    exonID = gffID.split('.exon')[0]
                    if exonID in extractList:
                        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
                elif 'CDS' in cols[2]:
                    cdsID = gffID.split('cds.')[-1]
                    if cdsID in extractList:
                        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
    count = lib.countGFFgenes(outfile)
    lib.log.info("Wrote {:,} gene models to {:}".format(count, outfile))

def GFF2tblCombined(evm, genome, trnascan, proteins, prefix, genenumber, justify, SeqCenter, SeqRefNum, tblout):
    from collections import OrderedDict
    '''
    function to take GFF3 annotation to produce a GBK tbl file.
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['start'])
    genenumber = int(genenumber)
    #generate genome length dictionary used for feature tbl generation
    scaffLen = {}
    with open(genome, 'rU') as fastain:
        for record in SeqIO.parse(fastain, 'fasta'):
            if not record.id in scaffLen:
                scaffLen[record.id] = len(record.seq)
    #generate protein start/stop data
    Proteins = {}
    with open(proteins, 'rU') as prots:
        for record in SeqIO.parse(prots, 'fasta'):
            ID = record.description.split(' ')[1]
            start, stop = (True,)*2
            Seq = str(record.seq)
            if not Seq.endswith('*'):
                stop = False
            if not Seq.startswith('M'):
                start = False
            if not ID in Proteins:
                Proteins[ID] = {'start': start, 'stop': stop}    
    Genes = {}
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
                    Genes[ID] = {'type': 'mRNA', 'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'mRNA': [], 'CDS': [], 'proper_start': Proteins[ID]['start'], 'proper_stop': Proteins[ID]['stop'], 'phase': [], 'product': 'hypothetical protein', 'note': '' }
                else:
                    print("Duplicate Gene IDs found, %s" % ID)
            else: #meaning needs to append to a gene ID
                info = attributes.split(';')
                ID,Parent = (None,)*2
                for x in info:
                    if x.startswith('ID='):
                        ID = x.replace('ID=', '')
                    if x.startswith('Parent='):
                        Parent = x.replace('Parent=', '')
                if not ID or not Parent:
                    print("Error, can't find ID or Parent. Malformed GFF file.")
                    sys.exit(1)
                if feature == 'mRNA':
                    if not ID in idParent:
                        idParent[ID] = Parent
                elif feature == 'exon':
                    if Parent in idParent:
                        geneName = idParent.get(Parent)
                    if geneName in Genes:
                        Genes[geneName]['mRNA'].append((start, end))
                elif feature == 'CDS':
                    if Parent in idParent:
                        geneName = idParent.get(Parent)
                    if geneName in Genes:
                        Genes[geneName]['CDS'].append((start, end))
                        Genes[geneName]['phase'].append(phase)
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
                    Genes[ID] = {'type':'tRNA', 'contig': contig, 'source': source, 'start': start, 'end': end, 'strand': strand, 'mRNA': [], 'CDS': [], 'proper_start': False, 'proper_stop': False, 'phase': [], 'product': '', 'note': ''  }
            elif feature == 'tRNA':
                info = attributes.split(';')
                ID = info[0].replace('ID=', '')
                Parent = info[1].replace('Parent=', '')
                Product = info[2].replace('product=', '')
                Genes[Parent]['product'] = Product
            elif feature == 'exon':
                ID = info[0].replace('ID=', '')
                Parent = info[1].replace('Parent=', '')
                if Parent in Genes:
                    Genes[Parent]['mRNA'].append((start, end))
                    Genes[Parent]['CDS'].append((start, end))
                    Genes[Parent]['phase'].append('0')

    #now sort dictionary by contig and location, rename using prefix
    sGenes = sorted(Genes.iteritems(), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    scaff2genes = {}
    for k,v in sortedGenes.items():
        #renaming scheme, leave if startswith locustag
        if k.startswith('novel_gene') or k.startswith('temp_gene') or k.startswith('split_gene'):
            genenumber += 1
            locusTag = prefix+'_'+str(genenumber).zfill(justify)            
        elif k.startswith(prefix) and '_'+prefix in k: #means models were merged
            locusTag = k.split('_'+prefix)[0]
        else:
            locusTag = k
        #get orientation and order exons/CDS
        if v['strand'] == '+':
            sortedExons = sorted(v['mRNA'], key=lambda tup: tup[0])
            sortedCDS = sorted(v['CDS'], key=lambda tup: tup[0])
        else:
            sortedExons = sorted(v['mRNA'], key=lambda tup: tup[0], reverse=True)
            sortedCDS = sorted(v['CDS'], key=lambda tup: tup[0], reverse=True)
        #get the codon_start by getting first CDS phase + 1
        indexStart = [x for x, y in enumerate(v['CDS']) if y[0] == sortedCDS[0][0]]
        codon_start = int(v['phase'][indexStart[0]]) + 1
        v['mRNA'] = sortedExons
        v['CDS'] = sortedCDS
        renamedGenes[locusTag] = v
        renamedGenes[locusTag]['codon_start'] = codon_start
        if not v['contig'] in scaff2genes:
            scaff2genes[v['contig']] = [locusTag]
        else:
            scaff2genes[v['contig']].append(locusTag)
    #now have scaffolds dict and gene dict, loop through scaff dict printing tbl
    with open(tblout, 'w') as tbl:
        for k,v in natsorted(scaff2genes.items()):
            tbl.write('>Feature %s\n' % k)
            tbl.write('1\t%s\tREFERENCE\n' % scaffLen.get(k))
            tbl.write('\t\t\t%s\t%s\n' % (SeqCenter, SeqRefNum))
            for genes in v: #now loop through each gene on the scaffold
                geneInfo = renamedGenes.get(genes)
                if geneInfo['type'] == 'mRNA':
                    #check for partial models
                    if geneInfo['proper_start'] and geneInfo['codon_start'] == 1:
                        partialStart = ''
                    else:
                        partialStart = '<'
                    if geneInfo['proper_stop']:
                        partialStop = ''
                    else:
                        partialStop = '>'
                    if geneInfo['strand'] == '+':
                        tbl.write('%s%i\t%s%i\tgene\n' % (partialStart, geneInfo['start'], partialStop, geneInfo['end']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0 and num == len(geneInfo['mRNA']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tmRNA\n' % (partialStart, exon[0], partialStop, exon[1]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tmRNA\n' % (partialStart, exon[0], exon[1]))
                            elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (exon[0], partialStop, exon[1]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T1_mrna\n' % genes)
                        for num, cds in enumerate(geneInfo['CDS']):
                            if num == 0 and num == len(geneInfo['CDS']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tCDS\n' % (partialStart, cds[0], partialStop, cds[1]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tCDS\n' % (partialStart, cds[0], cds[1]))
                            elif num == len(geneInfo['CDS']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (cds[0], partialStop, cds[1]))
                            else:
                                tbl.write('%s\t%s\n' % (cds[0], cds[1]))
                        tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'])
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T1\n' % genes)                                      
                    else:
                        tbl.write('%s%i\t%s%i\tgene\n' % (partialStart, geneInfo['end'], partialStop, geneInfo['start']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0 and num == len(geneInfo['mRNA']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tmRNA\n' % (partialStart, exon[1], partialStop, exon[0]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tmRNA\n' % (partialStart, exon[1], exon[0]))
                            elif num == len(geneInfo['mRNA']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (exon[1], partialStop, exon[0]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[1], exon[0]))                 
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T1_mrna\n' % genes)
                        for num, cds in enumerate(geneInfo['CDS']):
                            if num == 0 and num == len(geneInfo['CDS']) - 1: #single exon, so slightly differnt method
                                tbl.write('%s%s\t%s%s\tCDS\n' % (partialStart, cds[1], partialStop, cds[0]))
                            elif num == 0:
                                tbl.write('%s%s\t%s\tCDS\n' % (partialStart, cds[1], cds[0]))
                            elif num == len(geneInfo['CDS']) - 1: #this is last one
                                tbl.write('%s\t%s%s\n' % (cds[1], partialStop, cds[0]))
                            else:
                                tbl.write('%s\t%s\n' % (cds[1], cds[0]))
                        tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'])
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T1\n' % genes)
                elif geneInfo['type'] == 'tRNA':
                    if geneInfo['strand'] == '+':
                        tbl.write('<%i\t>%i\tgene\n' % (geneInfo['start'], geneInfo['end']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0:
                                tbl.write('<%s\t>%s\ttRNA\n' % (exon[0], exon[1]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        if geneInfo['product'] == 'tRNA-Xxx':
                            tbl.write('\t\t\tpseudo\n')
                        tbl.write('\t\t\tnote\t%s\n' % geneInfo['note'])                                    
                    else:
                        tbl.write('<%i\t>%i\tgene\n' % (geneInfo['end'], geneInfo['start']))
                        tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                        for num, exon in enumerate(geneInfo['mRNA']):
                            if num == 0:
                                tbl.write('<%s\t>%s\ttRNA\n' % (exon[1], exon[0]))
                            else:
                                tbl.write('%s\t%s\n' % (exon[1], exon[0]))
                        tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'])
                        if geneInfo['product'] == 'tRNA-Xxx':
                            tbl.write('\t\t\tpseudo\n')
                        tbl.write('\t\t\tnote\t%s\n' % geneInfo['note']) 

def getGBKmodels(input):
    '''
    function returns a dictionary of all gene models from a genbank file
    dictionary structure looks like:
    { 'gene_id': { 'contig': string,
                 'location': (start, end),
                   'strand': plus or minus,    
                     'mRNA': [(exon1_start, exon1_end), (exon2_start, exon2_end)...],
                     'CDS': [(cds1_start, cds1_end), (cds2_start, cds2_end)...],
                 'protein' : string(translation)}}
    '''
    result = {}
    with open(input, 'rU') as filein:
        for record in SeqIO.parse(filein, 'genbank'):
            Contig = record.id
            for f in record.features:
                #reset for every feature
                ID,start,end,strand,num_parts,exons,cds,protID,transcriptID,protSeq = (None,)*10
                if not f.type in ['gene', 'mRNA', 'tRNA', 'CDS']:
                    continue
                #get info from features
                ID = f.qualifiers['locus_tag'][0]
                start = f.location.nofuzzy_start
                end = f.location.nofuzzy_end
                strand = f.location.strand
                if strand == 1:
                    str = '+'
                else:
                    str = '-'
                num_parts = len(f.location.parts)           
                if f.type == "gene":
                    if not ID in result:
                        result[ID] = { 'contig': Contig, 'location': (int(start),int(end)), 'strand': str, 'mRNA': [], 'CDS': [], 'protein': protSeq, 'protein_id': protID, 'transcript_id': transcriptID }
                    else:
                        result[ID]['location'] = (int(start),int(end))
                        result[ID]['strand'] = strand
                
                elif f.type == "mRNA" or f.type == 'tRNA':
                    try:
                        transcriptID = f.qualifiers['transcript_id'][0]
                    except KeyError:
                        transcriptID = None
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
                    #update dictionary
                    if not ID in result:
                        result[ID] = { 'contig': Contig, 'location': '', 'strand': '', 'mRNA': sortedExons, 'CDS': [], 'protein': protSeq, 'protein_id': protID, 'transcript_id': transcriptID }
                    else:
                        result[ID]['mRNA'] = sortedExons
                        result[ID]['transcript_id'] = transcriptID
                
                elif f.type == 'CDS':
                    protSeq = f.qualifiers['translation'][0]
                    try:
                        protID = f.qualifiers['protein_id'][0]
                    except KeyError:
                        protID = None
                    cdsTuples = []
                    if num_parts < 2: #only single CDS
                        cdsTuples.append((int(start),int(end)))
                    else:
                        for i in range(0, num_parts):
                            ex_start = f.location.parts[i].nofuzzy_start
                            ex_end = f.location.parts[i].nofuzzy_end
                            cdsTuples.append((int(ex_start),int(ex_end)))
                    sortedCDS = sorted(cdsTuples, key=lambda tup: tup[0])
                    #update dictionary
                    if not ID in result:
                        result[ID] = { 'contig': Contig, 'location': '', 'strand': '', 'mRNA': [], 'CDS': sortedCDS, 'protein': protSeq, 'protein_id': protID, 'transcript_id': transcriptID }
                    else:
                        result[ID]['CDS'] = sortedCDS
                        result[ID]['protein'] = protSeq
                        result[ID]['protein_id'] = protID
    #return the dictionary
    return result

def compareAnnotations(old, new, output):
    '''
    function takes two GenBank annotated genomes and compares gene models
    output is a tsv file for each locus and a description of what is different
    '''
    result = {}
    dropped = 0
    added = 0
    no_change = 0
    UTR_added = 0
    yardSale = 0
    exonChange = 0
    modelChangeNotProt = 0
    lib.log.info("Parsing GenBank files...comparing annotation")
    oldAnnotation = getGBKmodels(old)
    newAnnotation = getGBKmodels(new)
    #do the simple stuff first, find models that were deleted
    for key,value in oldAnnotation.items():
        if not key in newAnnotation:
            result[key] = [value['contig'], value['location'], value['strand'], value['mRNA'], value['CDS'], value['protein'], value['protein_id'], value['transcript_id'], 'gene model removed']
            dropped += 1
    #now look for new ones as well as do comparison of gene models
    for key,newModel in newAnnotation.items():
        if not key in oldAnnotation:
            result[key] = [newModel['contig'], newModel['location'], newModel['strand'], newModel['mRNA'], newModel['CDS'], newModel['protein'], '', '', 'new gene model']
            added += 1
        else: #now look to see how its changed
            oldModel = oldAnnotation.get(key) #returns dictionary
            #run comparisons of the data
            contigs = oldModel['contig'] == newModel['contig']
            location = oldModel['location'] == newModel['location']
            strand = oldModel['strand'] == newModel['strand']
            mRNA = oldModel['mRNA'] == newModel['mRNA']
            CDS = oldModel['CDS'] == newModel['CDS']
            protein = oldModel['protein'] == newModel['protein']
            compare = [contigs, location, strand, mRNA, CDS, protein]
            if compare == [True,True,True,True,True,True]:
                result[key] = [newModel['contig'], newModel['location'], newModel['strand'], newModel['mRNA'], newModel['CDS'], newModel['protein'], oldModel['protein_id'], oldModel['transcript_id'], 'gene model has not changed']
                no_change += 1
            elif compare == [True,False,True,False,True,True]:
                result[key] = [newModel['contig'], newModel['location'], newModel['strand'], newModel['mRNA'], newModel['CDS'], newModel['protein'], oldModel['protein_id'], oldModel['transcript_id'], 'gene coordinates changed; UTR has been added to gene model; CDS unchanged; translation unchanged']
                UTR_added += 1
            elif compare == [True,False,True,False,False,False]:
                result[key] = [newModel['contig'], newModel['location'], newModel['strand'], newModel['mRNA'], newModel['CDS'], newModel['protein'], oldModel['protein_id'], oldModel['transcript_id'], 'gene coordinates changed; mRNA and CDS has changed; translation changed']
                yardSale += 1
            elif compare == [True,True,True,False,False,False]:
                result[key] = [newModel['contig'], newModel['location'], newModel['strand'], newModel['mRNA'], newModel['CDS'], newModel['protein'], oldModel['protein_id'], oldModel['transcript_id'], 'gene coordinates unchanged; mRNA and CDS has changed; translation changed']
                exonChange += 1
            elif compare == [True, False, True, False, False, True]:
                result[key] = [newModel['contig'], newModel['location'], newModel['strand'], newModel['mRNA'], newModel['CDS'], newModel['protein'], oldModel['protein_id'], oldModel['transcript_id'], 'gene coordinates changed; mRNA and CDS has changed; translation unchanged']
                modelChangeNotProt += 1
                
    with open(output, 'w') as out:
        out.write('Locus_tag\tContig:start-end\tStrand\tTranscript_ID\tLength\tNum_Exons\tProtein_ID\tProt_Length\tDescription\n')
        for k,v in natsorted(result.items()):
            contig, loc, strand, exons, cds, protSeq, protID, tranID, description = v
            try:
                ProtLength = len(protSeq)
            except TypeError:
                ProtLength = 0
            GeneLength = int(loc[1]) - int(loc[0])
            NumExons = len(exons)
            start = str(loc[0])
            end = str(loc[1])
            out.write('%s\t%s:%s-%s\t%s\t%s\t%i\t%i\t%s\t%i\t%s\n' % (k, contig, start, end, strand, tranID, GeneLength, NumExons, protID, ProtLength, description))
    #output some simple stats to cmd line
    lib.log.info("Updated annotation complete:\n\
-------------------------------------------------------\n\
Total Gene Models:\t{:,}\n\
New Gene Models:\t{:,}\n\
No Change in Model:\t{:,}\n\
Update UTRs:\t\t{:,}\n\
Exons Changed:\t\t{:,}\n\
Exons/CDS Changed:\t{:,}\n\
Exon Changed/Prot same:\t{:,}\n\
Dropped Models:\t\t{:,}\n\
-------------------------------------------------------".format(len(newAnnotation), added, no_change, UTR_added, exonChange, yardSale, modelChangeNotProt, dropped))

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
        
programs = ['fasta', 'mysql', 'gmap', 'blat', 'tbl2asn', 'hisat2', 'hisat2-build', 'kallisto', 'Trinity', 'bedtools']
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
        shutil.copyfile(args.pasa_gff, PASA_gff)
    if not lib.checkannotations(PASA_gff):
        runPASA(fastaout, trinity_transcripts, args.stranded, args.max_intronlen, args.cpus, gffout, organism_name, PASA_gff, args.pasa_config)
else:
    if not lib.checkannotations(PASA_gff):
        runPASA(fastaout, trinity_transcripts, args.stranded, args.max_intronlen, args.cpus, gffout, organism_name, PASA_gff, pasaConfigFile)
    
#now run Kallisto steps, if mixed PE and SE reads, only PE reads will be used for Kallisto as there isn't a reasonable way to combine them
KallistoAbundance = os.path.join(tmpdir, 'kallisto.tsv')
if args.kallisto:
    shutil.copyfile(args.kallisto, KallistoAbundance)
if not lib.checkannotations(KallistoAbundance):
    runKallisto(PASA_gff, fastaout, trim_reads, args.stranded, args.cpus, KallistoAbundance)

#parse Kallisto results with PASA GFF
BestModelGFF = os.path.join(tmpdir, 'bestmodels.gff3')
getBestModel(PASA_gff, fastaout, KallistoAbundance, BestModelGFF)

#generate proteins 
finalProts = os.path.join(tmpdir, 'final.proteins.fasta')
cmd = [os.path.join(PASA, 'misc_utilities', 'gff3_file_to_proteins.pl'), BestModelGFF, fastaout, 'prot']
lib.runSubprocess2(cmd, '.', lib.log, finalProts)

#make sure tRNA models don't overlap new gene models
cleanTRNA = os.path.join(tmpdir, 'trna.no-overlaps.gff')
cmd = ['bedtools', 'intersect', '-v', '-a', trnaout, '-b', BestModelGFF]
lib.runSubprocess2(cmd, '.', lib.log, cleanTRNA)

#generate tbl file
gagdir = os.path.join(tmpdir, 'tbl2asn')
if os.path.isdir(gagdir):
    shutil.rmtree(gagdir)
os.makedirs(gagdir)
TBLFile = os.path.join(gagdir, 'genome.tbl')
GFF2tblCombined(BestModelGFF, fastaout, cleanTRNA, finalProts, locustag, genenumber, justify, args.SeqCenter, args.SeqAccession, TBLFile)

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

#retrieve files/reorganize
#shutil.copyfile(os.path.join(gagdir, 'genome.gff'), final_gff)
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
compareAnnotations(GBK, final_gbk, Changes)

lib.log.info("Funannotate update is finished, output files are in the %s/update_results folder" % (args.out))
lib.log.info("Your next step might be functional annotation, suggested commands:\n\
-------------------------------------------------------\n\
Run InterProScan (Docker required): \n{:} -i={:} -c={:}\n\n\
Run antiSMASH: \nfunannotate remote -i {:} -m antismash -e youremail@server.edu\n\n\
Annotate Genome: \nfunannotate annotate -i {:} --iprscan {:} --cpus {:} --sbt yourSBTfile.txt\n\
-------------------------------------------------------\n\
            ".format(os.path.join(parentdir, 'util', 'interproscan_docker.sh'), \
            os.path.join(args.out, 'update_results', organism_name+'.proteins.fa'), \
            args.cpus, \
            args.out, \
            args.out, \
            os.path.join(args.out, 'update_results', organism_name+'.proteins.fa.xml'), \
            args.cpus))
sys.exit(1)
