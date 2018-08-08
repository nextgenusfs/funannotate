#!/usr/bin/env python

import sys, os, subprocess, inspect, shutil, argparse, fnmatch
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib
from natsort import natsorted
from lib.interlap import InterLap
from collections import defaultdict

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='funannotate-train.py', usage="%(prog)s [options] -i genome.fasta",
    description = '''Script is a wrapper for automated Trinity/PASA generation of training data.''',
    epilog = """Written by Jon Palmer (2017-2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--input', required=True, help='Genome in FASTA format')
parser.add_argument('-l', '--left', nargs='+', help='Left (R1) FASTQ Reads')
parser.add_argument('--left_norm', help='Left (R1) FASTQ Reads')
parser.add_argument('--right_norm', help='Right (R2) normalized FASTQ Reads')
parser.add_argument('--single_norm', help='single normalized FASTQ Reads')
parser.add_argument('-r', '--right', nargs='+', help='Right (R2) FASTQ Reads')
parser.add_argument('-s', '--single', nargs='+', help='Single ended FASTQ Reads')
parser.add_argument('--pacbio_isoseq', help='PacBio Iso-seq data')
parser.add_argument('--nanopore_cdna', help='Nanopore 2d cDNA data')
parser.add_argument('--nanopore_mrna', help='Nanopore direct mRNA data')
parser.add_argument('-o', '--out', required=True, help='Basename folder of output files')
parser.add_argument('-c', '--coverage', default=50, type=int, help='Depth to normalize reads to')
parser.add_argument('--trinity', help='Trinity genome guided FASTA results')
parser.add_argument('--memory', default='50G', help='RAM to use for Jellyfish/Trinity')
parser.add_argument('--no_normalize_reads', action='store_true', help='skip normalization')
parser.add_argument('--no_trimmomatic', action='store_true', help='skip quality trimming via trimmomatic')
parser.add_argument('--jaccard_clip', action='store_true', help='Turn on jaccard_clip for dense genomes')
parser.add_argument('--pasa_alignment_overlap', default='30.0', help='PASA --stringent_alingment_overlap')
parser.add_argument('--pasa_db', default='sqlite', choices=['mysql', 'sqlite'], help='PASA SQL database to use')
parser.add_argument('--max_intronlen', default=3000, help='Maximum intron length for gene models')
parser.add_argument('--stranded', default = 'no', choices=['RF','FR','F','R','no'], help='RNA seq strandedness')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to use')
parser.add_argument('--header_length', default=16, type=int, help='Max length for fasta headers')
parser.add_argument('--species', help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('--isolate', help='Isolate name (e.g. Af293)')
parser.add_argument('--strain', help='Strain name (e.g. CEA10)')
parser.add_argument('--PASAHOME', help='Path to PASA home directory, $PASAHOME')
parser.add_argument('--TRINITYHOME', help='Path to Trinity config directory, $TRINITYHOME')
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
        lib.Fzip_inplace(x, args.cpus)
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
    lib.Fzip_inplace(output, args.cpus)
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
    lib.log.info("Running read normalization with Trinity")
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

def long2fasta(readTuple, cpus, tmpdir, combined, combinedClean):
    '''
    Run SeqClean on long reads, return cleaned tuple and combined output
    tuple is (pb_iso, nano_cdna, nano_mrna)
    '''
    def _convert2fasta(file, output):
        messy = []
        with open(output, 'w') as outfile:
            if file.endswith('.gz'):
                newfile = file.replace('.gz', '')
                messy.append(newfile)
                lib.Funzip(file, newfile, cpus)
                file = newfile                  
            if file.endswith('.fa') or file.endswith('.fasta'):
                with open(file, 'rU') as infile:
                    for record in SeqIO.parse(infile, 'fasta'):
                        SeqIO.write(record, outfile, 'fasta')
            elif file.endswith('.fq') or file.endswith('.fastq'):
                with open(file, 'rU') as infile:
                    for record in SeqIO.parse(infile, 'fastq'):
                        SeqIO.write(record, outfile, 'fasta')
        #clean up
        for x in messy:
            lib.SafeRemove(x)
    if os.path.islink(combined) or os.path.isfile(combined):
        results = []
        originals = []
        for i in [PBiso+'.clean', nanocdna+'.clean', nanomrna+'.clean']:
            if lib.checkannotations(i):
                results.append(i)
                originals.append(i.replace('.clean', ''))
            else:
                results.append(None)
                originals.append(None)
        return tuple(originals), tuple(results)
    else:
        lib.log.info('Processing long reads: converting to fasta and running SeqClean')
        results = []
        if readTuple[0] and not lib.checkannotations(PBiso+'.clean'):
            _convert2fasta(readTuple[0], PBiso)
            runSeqClean(PBiso, tmpdir)
        if readTuple[1] and not lib.checkannotations(nanocdna+'.clean'):
            _convert2fasta(readTuple[1], nanocdna)
            runSeqClean(nanocdna, tmpdir)
        if readTuple[2] and not lib.checkannotations(nanomrna+'.clean'):
            _convert2fasta(readTuple[2], nanomrna)
            runSeqClean(nanomrna, tmpdir)
        for i in [PBiso+'.clean', nanocdna+'.clean', nanomrna+'.clean']:
            if lib.checkannotations(i):
                results.append(i)
            else:
                results.append(None)
        validResults = [x for x in results if x is not None]
        validOriginal = [x.replace('.clean', '') for x in validResults]
        validCln = [x.replace('.clean', '.cln') for x in validResults]
        ClnOut = combined+'.cln'
        if len(validResults) > 1:
            lib.catFiles(*validResult, output=combinedClean)
            lib.catFiles(*validOriginal, output=combined)
            lib.catFiles(*validCln, output=ClnOut)
        else:
            if not lib.checkannotations(combinedClean):
                os.symlink(os.path.abspath(validResults[0]), combinedClean)
            if not lib.checkannotations(combined):  
                os.symlink(os.path.abspath(validOriginal[0]), combined)
            if not lib.checkannotations(ClnOut):
                os.symlink(os.path.abspath(validCln[0]), ClnOut)
        return tuple(validOriginal), tuple(results)

def runSeqClean(input, folder):
    '''
    wrapper to run PASA seqclean on Trinity transcripts
    '''
    if args.cpus > 16:
        cpus = 16
    else:
        cpus = args.cpus
    if os.path.isfile(input + ".clean"):
        lib.log.info('Existing SeqClean output found: {:}'.format(os.path.join(folder, input + ".clean")))
    else:
        cmd = [os.path.join(PASA, 'bin', 'seqclean'), os.path.basename(input), '-c', str(cpus)]
        lib.runSubprocess(cmd, folder, lib.log)
    for f in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, f)):
            if f.startswith('cleaning'):
                lib.SafeRemove(os.path.join(folder, f))

def bam2fasta(input, output):
    cmd = ['samtools', 'fasta', '-@', str(args.cpus), '-F', '0x4', input]
    lib.runSubprocess2(cmd, '.', lib.log, output)

def bam2fasta_unmapped(input, output):
    cmd = ['samtools', 'fasta', '-@', str(args.cpus), '-f', '0x4', input]
    lib.runSubprocess2(cmd, '.', lib.log, output)

def mapTranscripts(genome, longTuple, assembled, tmpdir, trinityBAM, allBAM):
    '''
    function will map long reads and trinity to genome, return sorted BAM
    '''
    isoBAM = os.path.join(tmpdir, 'isoseq.coordSorted.bam')
    isoSeqs = os.path.join(tmpdir, 'isoseq.coordSorted.fasta')
    nano_cdnaBAM = os.path.join(tmpdir, 'nano_cDNA.coordSorted.bam')
    nano_cdnaSeqs = os.path.join(tmpdir, 'nano_cDNA.coordSorted.fasta')
    nano_mrnaBAM = os.path.join(tmpdir, 'nano_mRNA.coordSorted.bam')
    nano_mrnaSeqs = os.path.join(tmpdir, 'nano_mRNA.coordSorted.fasta')
    mappedSeqs = []
    mappedLong = os.path.join(tmpdir, 'long-reads.mapped.fasta')
    if not all(v is None for v in longTuple): #tuple is (iso-seq, nanopore_cDNA, nanopore_mRNA)
        #run minimap2 alignment
        lib.log.info('Aligning long reads to genome with minimap2')
        if longTuple[0]: #run iso-seq method
            lib.iso_seq_minimap2(longTuple[0], genome, args.cpus, args.max_intronlen, isoBAM)
            bam2fasta(isoBAM, isoSeqs)
        if longTuple[1]: #run nano cDNA
           lib.nanopore_cDNA_minimap2(longTuple[1], genome, args.cpus, args.max_intronlen, nano_cdnaBAM)
           bam2fasta(nano_cdnaBAM, nano_cdnaSeqs)
        if longTuple[2]: #run nano mRNA
            lib.nanopore_mRNA_minimap2(longTuple[2], genome, args.cpus, args.max_intronlen, nano_mrnaBAM)
            bam2fasta(nano_mrnaBAM, nano_mrnaSeqs)
        for x in [isoSeqs, nano_cdnaSeqs, nano_mrnaSeqs]:
            if lib.checkannotations(x):
                mappedSeqs.append(x)
        if len(mappedSeqs) > 0:
            lib.catFiles(*mappedSeqs, output=mappedLong)
            lib.SafeRemove(isoSeqs)
            lib.SafeRemove(nano_cdnaSeqs)
            lib.SafeRemove(nano_mrnaSeqs)
    if lib.checkannotations(assembled): #Trinity transcripts
        #want to recover any long-reads that don't map to Trinity transcripts but do map to genome
        crosscheckBAM = os.path.join(tmpdir, 'trinity.vs.long-reads.bam')
        unmappedLong = os.path.join(tmpdir, 'long-reads.trinity.unique.fasta')
        if lib.checkannotations(mappedLong):
            lib.log.info('Finding long-reads not represented in Trinity assemblies')
            minimap_cmd = ['minimap2', '-ax', 'map-ont', '-t', str(args.cpus), '--secondary=no', assembled, mappedLong]
            cmd = [os.path.join(parentdir, 'util', 'sam2bam.sh'), " ".join(minimap_cmd), str(args.cpus // 2), crosscheckBAM]
            if not lib.checkannotations(crosscheckBAM):
                lib.runSubprocess(cmd, '.', lib.log)
            bam2fasta_unmapped(crosscheckBAM, unmappedLong)
            lib.log.info('Adding {:,} unique long-reads to Trinity assemblies'.format(lib.countfasta(unmappedLong)))
            lib.SafeRemove(crosscheckBAM)
        if lib.checkannotations(unmappedLong):
            trinityCombined = os.path.join(tmpdir, 'trinity.long-reads.fasta')
            trinityCombinedClean = trinityCombined+'.clean'
            lib.catFiles(*[assembled,unmappedLong], output=trinityCombined)
            runSeqClean(trinityCombined, tmpdir)
        else:
            trinityCombinedClean = assembled
            trinityCombined = assembled.replace('.clean', '')
        #finally run trinity mapping
        lib.minimap2Align(trinityCombinedClean, genome, args.cpus, args.max_intronlen, trinityBAM)
        
    bamResults = [isoBAM, nano_cdnaBAM, nano_mrnaBAM, trinityBAM]
    foundResults = []
    for r in bamResults:
        if lib.checkannotations(r):
            foundResults.append(r)
    if len(foundResults) > 1:
        lib.log.info('Merging BAM files: {:}'.format(', '.join(foundResults)))
        lib.mergeBAMs(*foundResults, cpus=args.cpus, output=allBAM)
    elif len(foundResults) == 0:
        lib.log.error('Alignment failed, BAM files empty. Please check logfile')
        sys.exit(1)
    else:
        os.symlink(os.path.abspath(foundResults[0]), os.path.abspath(allBAM))
    return trinityCombined, trinityCombinedClean


def runTrinityGG(genome, readTuple, longReads, shortBAM, output):
    '''
    function will run genome guided Trinity. First step will be to run hisat2 to align reads
    to the genome, then pass that BAM file to Trinity to generate assemblies
    '''
    #build hisat2 index, using exons and splice sites
    lib.log.info("Building Hisat2 genome index")
    cmd = ['hisat2-build', genome, os.path.join(tmpdir, 'hisat2.genome')]
    lib.runSubprocess4(cmd, '.', lib.log)
    #align reads using hisat2
    lib.log.info("Aligning reads to genome using Hisat2")
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
    cmd = [os.path.join(parentdir, 'util', 'sam2bam.sh'), " ".join(hisat2cmd), str(bamthreads), shortBAM]
    lib.runSubprocess(cmd, '.', lib.log)

    #now launch Trinity genome guided
    TrinityLog = os.path.join(tmpdir, 'Trinity-gg.log')
    lib.log.info("Running genome-guided Trinity, logfile: %s" % TrinityLog)
    lib.log.info("Clustering of reads from BAM and preparing assembly commands")
    jaccard_clip = []
    if args.jaccard_clip:
        jaccard_clip = ['--jaccard_clip']
    if args.stranded != 'no':
        cmd = ['Trinity', '--SS_lib_type', args.stranded, '--no_distributed_trinity_exec', '--genome_guided_bam', shortBAM, '--genome_guided_max_intron', str(args.max_intronlen), '--CPU', str(args.cpus), '--max_memory', args.memory, '--output', os.path.join(tmpdir, 'trinity_gg')]
    else:
        cmd = ['Trinity', '--no_distributed_trinity_exec', '--genome_guided_bam', shortBAM, '--genome_guided_max_intron', str(args.max_intronlen), '--CPU', str(args.cpus), '--max_memory', args.memory, '--output', os.path.join(tmpdir, 'trinity_gg')]
    cmd = cmd + jaccard_clip
    if lib.checkannotations(longReads):
        cmd = cmd + ['--long_reads', os.path.realpath(longReads)]
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
    
def runPASAtrain(genome, transcripts, cleaned_transcripts, gff3_alignments, stringtie_gtf, stranded, intronlen, cpus, dbname, output):
    '''
    function will run PASA align assembly and then choose best gene models for training
    '''
    if cpus > 2:
        pasa_cpus = cpus
    else:
        pasa_cpus = 2
    #create tmpdir
    folder = os.path.join(tmpdir, 'pasa')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    
    #get config files and edit
    alignConfig = os.path.join(folder, 'alignAssembly.txt')
    pasaDBname = dbname.replace('-', '_')
    if args.pasa_db == 'sqlite':
        pasaDBname_path = os.path.abspath(os.path.join(folder, pasaDBname))
    else:
        pasaDBname_path = pasaDBname
    with open(alignConfig, 'w') as config1:
        with open(os.path.join(PASA, 'pasa_conf', 'pasa.alignAssembly.Template.txt'), 'rU') as template1:
            for line in template1:
                line = line.replace('<__DATABASE__>', pasaDBname_path)
                line = line.replace('<__MYSQLDB__>', pasaDBname_path)
                config1.write(line)
    if not os.path.isfile(os.path.join(folder, pasaDBname+'.assemblies.fasta')):
        #now run first PASA step, note this will dump any database with same name 
        lib.log.info("Running PASA alignment step using {:,} transcripts".format(lib.countfasta(cleaned_transcripts)))
        cmd = [LAUNCHPASA, '-c', os.path.abspath(alignConfig), '-r', '-C', '-R', '-g', os.path.abspath(genome), '--ALIGNERS', 'blat', '--IMPORT_CUSTOM_ALIGNMENTS', gff3_alignments, '-T','-t', os.path.abspath(cleaned_transcripts), '-u', os.path.abspath(transcripts), '--stringent_alignment_overlap', args.pasa_alignment_overlap, '--TRANSDECODER', '--ALT_SPLICE', '--MAX_INTRON_LENGTH', str(intronlen), '--CPU', str(pasa_cpus)]
        if stranded != 'no':
            cmd = cmd + ['--transcribed_is_aligned_orient']
        if lib.checkannotations(stringtie_gtf):
            cmd = cmd + ['--trans_gtf', os.path.abspath(stringtie_gtf)]
        lib.runSubprocess(cmd, folder, lib.log)
    else:
        lib.log.info('Existing PASA assemblies found: {:}'.format(os.path.join(folder, pasaDBname+'.assemblies.fasta')))
    #generate TSV gene-transcripts
    Loci = []
    numTranscripts = 0
    with open(os.path.join(folder, 'pasa.gene2transcripts.tsv'), 'w') as gene2transcripts:
        with open(os.path.join(folder, pasaDBname+'.pasa_assemblies_described.txt'), 'rU') as description:
            for line in description:
                if not line.startswith('#'):
                    cols = line.split('\t')
                    gene2transcripts.write('g_%s\t%s\n' % (cols[1], cols[2]))
                    numTranscripts += 1
                    if not cols[1] in Loci:
                        Loci.append(cols[1])
    lib.log.info("PASA assigned {:,} transcipts to {:,} loci (genes)".format(numTranscripts, len(Loci)))
    lib.log.info("Getting PASA models for training with TransDecoder")
    pasa_training_gff = os.path.join(folder, pasaDBname+'.assemblies.fasta.transdecoder.genome.gff3') 
    cmd = [os.path.join(PASA, 'scripts', 'pasa_asmbls_to_training_set.dbi'), '--pasa_transcripts_fasta', pasaDBname+'.assemblies.fasta', '--pasa_transcripts_gff3', pasaDBname+'.pasa_assemblies.gff3']
    lib.runSubprocess(cmd, folder, lib.log)
    #grab final result
    shutil.copyfile(pasa_training_gff, output)
    lib.log.info('PASA finished. PASAweb accessible via: localhost:port/cgi-bin/index.cgi?db=%s' % pasaDBname_path)

def pasa_transcript2gene(input):
    #modify kallisto ouput to map gene names to each mRNA ID so you know what locus they have come from
    mRNADict = {}
    #since mRNA is unique, parse the transcript file which has mRNAID geneID in header
    with open(input, 'rU') as transin:
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
    return mRNADict

def runKallisto(input, fasta, readTuple, stranded, cpus, folder, output):
    '''
    function takes GFF3 output from PASA compare, extracts transcripts, and then calculates TPM
    using Kallisto to idenitfy the best scoring gene model for each locus, the left and right
    these should be the adapter cleaned non-normalized Illumina reads
    '''
    lib.log.info("Using Kallisto TPM data to determine which PASA gene models to select at each locus")
    #convert GFF to transcripts
    if not os.path.exists(folder):
        os.makedirs(folder) # handle already existing folder okay? could also delete it
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
    mRNADict = pasa_transcript2gene(PASAtranscripts)
    
    #some PASA models can have incomplete CDS and are wrong, get list of incompletes to ignore list
    ignore = []
    with open(input, 'rU') as infile:
        for line in infile:
            if line.startswith('#PROT'):
                if line.endswith('\t\n'):
                    ID = line.split(' ')[1]
                    ignore.append(ID)
    if len(ignore) > 0:
        lib.log.debug("Ignoring %i incomplete PASA models: %s" % (len(ignore), ','.join(ignore)))
                    
    #now make new tsv file with #mRNAID geneID location TPM
    with open(output, 'w') as outfile:
        outfile.write("#mRNA-ID\tgene-ID\tLocation\tTPM\n")
        with open(os.path.join(folder, 'kallisto', 'abundance.tsv'), 'rU') as infile:
            for line in infile:
                if line.startswith('targed_id'):
                    continue
                line = line.rstrip()
                cols = line.split('\t')
                if cols[0] in ignore:
                    continue
                if cols[0] in mRNADict:
                    geneHit = mRNADict.get(cols[0])
                    geneID = geneHit[0]
                    location = geneHit[1]
                    outfile.write('%s\t%s\t%s\t%s\n' % (cols[0], geneID, location, cols[4]))
   

def getPASAtranscripts2genes(input, output):
    '''
    function to parse PASA assemblies GFF3 to generate TSV file for transdecoder
    GFF format is non-standard and looks like this, the transcript IDs are in the Target field
    CM002236    assembler-Neurospora_crassa_train2  cDNA_match  933 2973    .   -   .   ID=align_64070;Target=asmbl_1 1 2041 +
    CM002236    assembler-Neurospora_crassa_train2  cDNA_match  933 1449    .   -   .   ID=align_64071;Target=asmbl_2 1447 1963 +
    CM002236    assembler-Neurospora_crassa_train2  cDNA_match  1528    2973    .   -   .   ID=align_64071;Target=asmbl_2 1 1446 +
    '''
    Genes = {}
    with open(input, 'rU') as infile:
        for line in infile:
            line = line.rstrip()
            contig, source, feature, start, end, score, strand, phase, attributes = line.split('\t')
            ID,Target = (None,)*2
            info = attributes.split(';')
            for x in info:
                if x.startswith('ID='):
                    ID = x.replace('ID=', '')
                elif x.startswith('Target='):
                    tmp = x.replace('Target=', '')
                    Target = tmp.split(' ')[0]
            if ID and Target:
                if not ID in Genes:
                    Genes[ID] = {'contig': contig, 'ids': Target, 'mRNA': [(int(start), int(end))]}
                else:
                    Genes[ID]['mRNA'].append((int(start), int(end)))
    #after all positions added, now create interlap on start stop positions
    inter = defaultdict(InterLap)
    for k,v in natsorted(Genes.items()):
        sortedExons = sorted(v['mRNA'], key=lambda tup: tup[0])
        inter[v['contig']].add((sortedExons[0][0], sortedExons[-1][1], k, v['ids']))
    #now loop through interlap object and create a gene2transcript dictionary
    Transcript2Gene = {}
    counter = 1
    for scaffold in inter:
        for x in inter[scaffold]:
            loc = [x[0],x[1]]
            hits = list(inter[scaffold].find(loc))
            Overlap = []
            for y in hits:
                percentOverlap = pOverlap(loc,[y[0],y[1]])
                if percentOverlap >= (float(args.pasa_alignment_overlap) / 100) and not y[3] in Transcript2Gene:
                    Overlap.append(y[3])
            if len(Overlap) > 0:
                for transcript in Overlap:
                    if not transcript in Transcript2Gene:
                        Transcript2Gene[transcript] = 'g_'+str(counter)
            counter += 1
    #finally print out TSV file
    unique = []
    with open(output, 'w') as outfile:
        for k,v in natsorted(Transcript2Gene.items()):
            if not v in unique:
                unique.append(v)
            outfile.write('{:}\t{:}\n'.format(v,k))
    return len(unique)
    
def pOverlap(one, two):
    rone = set(range(one[0],one[1]))
    rtwo = set(range(two[0],two[1]))
    overlap = rone & rtwo
    return len(overlap) / float(len(rone))                      

def getBestModel(input, fasta, abundances, outfile):
    #function to parse PASA results and generate GFF3; supports multiple transcripts
    lib.log.info("Parsing expression value results. Keeping best transcript at each locus.")
    Expression = {}
    with open(abundances, 'rU') as tpms:
        for line in tpms:
            line = line.rstrip()
            if line.startswith('#') or line.startswith('target_id'):
                continue
            transcriptID, geneID, Loc, TPM = line.split('\t')
            if not transcriptID in Expression:
                Expression[geneID] = float(TPM)
    
    #load GFF3 output into annotation and interlap dictionaries.
    inter_gene, Genes = lib.gff2interlap(input, fasta)
    bestHits = []
    overlap = []
    for scaffold in inter_gene:
        for x in inter_gene[scaffold]:
            loc = [x[0],x[1]]
            hits = list(inter_gene[scaffold].find(loc))
            ExpHits = []
            for y in hits:
                percentOverlap = pOverlap(loc,[y[0],y[1]])
                if percentOverlap >= (float(args.pasa_alignment_overlap) / 100): #overlap more than args.pasa_alignment_overlap
                    if y[2] in Expression:
                        ExpHits.append((y[2], Expression[y[2]], percentOverlap))
                    else:
                        ExpHits.append((y[2], 0.00, percentOverlap))
            sortedExpHits = sorted(ExpHits, key=lambda x: x[1], reverse=True)
            for i in range(0, len(sortedExpHits)):
                if i == 0:
                    bestHits.append(sortedExpHits[i][0])
                else:
                    overlap.append(sortedExpHits[i][0])
    bestModels = {}
    for k,v in natsorted(Genes.items()):
        if k in bestHits:
            bestModels[k] = v
    lib.dict2gff3(bestModels, outfile)
    lib.log.info('Wrote {:,} PASA gene models'.format(len(bestModels)))     


#create folder structure
if not os.path.isdir(args.out):
    os.makedirs(args.out)
    os.makedirs(os.path.join(args.out, 'training'))
    os.makedirs(os.path.join(args.out, 'logfiles'))
else:
    #make sure subdirectories exist
    dirs = [os.path.join(args.out, 'training'), os.path.join(args.out, 'logfiles')]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)

tmpdir = os.path.join(args.out, 'training')

#create log file
log_name = os.path.join(args.out, 'logfiles', 'funannotate-train.log')
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
   
#do some checks and balances
if not args.PASAHOME:
    try:
        PASA = os.environ["PASAHOME"]
    except KeyError:
        lib.log.error("$PASAHOME environmental variable not found, PASA is not properly configured.  You can use the --PASAHOME argument to specifiy a path at runtime")
        sys.exit(1)
else:
    PASA = args.PASAHOME

#try to autodetect different PASA distributions
if os.path.isfile(os.path.join(PASA,'Launch_PASA_pipeline.pl')): #then v2.3.0 or newer
    LAUNCHPASA = os.path.join(PASA,'Launch_PASA_pipeline.pl')
    PASAVERSION = '2.3.0'
elif os.path.isfile(os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl')): #older version
    LAUNCHPASA = os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl')
    args.pasa_db = 'mysql' #sqlite not available
    PASAVERSION = '2.2.0'

if not args.TRINITYHOME:
    try:
        TRINITY = os.environ["TRINITYHOME"]
    except KeyError:
        lib.log.error("$TRINITYHOME environmental variable not found, TRINITY is not properly configured. You can use the --TRINITYHOME argument to specify a path at runtime.")
        sys.exit(1)
else:
    TRINITY = args.TRINITYHOME
        
programs = ['fasta', 'minimap2', 'blat', 'hisat2', 'hisat2-build', 'Trinity', 'java', 'kallisto', LAUNCHPASA, os.path.join(PASA, 'bin', 'seqclean')]
lib.CheckDependencies(programs)

#see if organism/species/isolate was passed at command line, build PASA naming scheme
organism, strain, isolate = (None,)*3
if args.species:
    organism = args.species
else:
    organism = os.path.basename(args.input).split('.fa')[0]
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

#check input, make sure fasta headers are compatible
header_test = lib.checkFastaHeaders(args.input, args.header_length)
#move into tmpfolder
genome = os.path.join(tmpdir, 'genome.fasta')
shutil.copyfile(args.input, genome)

if args.left and args.right and args.single:
    lib.log.info("Combining PE and SE reads supported, but you will lose stranded information, setting --stranded no")
    args.stranded = 'no'
    
#check input reads
#get absolute paths for reads and concate if there are multiple
s_reads, l_reads, r_reads = (None,)*3
if not lib.checkannotations(os.path.join(tmpdir, 'single.fq.gz')):
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
            lib.concatenateReads(single_reads, s_reads)
        else:
            s_reads = single_reads[0]
        if s_reads.endswith('.fq'):
            lib.Fzip_inplace(s_reads, args.cpus)
            s_reads = s_reads+'.gz'
        if not lib.checkannotations(os.path.join(tmpdir, 'single.fq.gz')):   
            if os.path.dirname(os.path.abspath(tmpdir)) != os.path.dirname(os.path.abspath(s_reads)):
                try:
                    os.symlink(os.path.realpath(s_reads), os.path.join(tmpdir, 'single.fq.gz'))
                except OSError:
                    pass
else:
    s_reads = os.path.join(tmpdir, 'single.fq.gz')
    
if not lib.checkannotations(os.path.join(tmpdir, 'left.fq.gz')) or not lib.checkannotations(os.path.join(tmpdir, 'right.fq.gz')):
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
            lib.concatenateReads(left_reads, l_reads)
            lib.concatenateReads(right_reads, r_reads)
        else:
            l_reads = left_reads[0]
            r_reads = right_reads[0]
        if l_reads.endswith('.fq'):
            lib.Fzip_inplace(l_reads, args.cpus)
            l_reads = l_reads+'.gz'
        if r_reads.endswith('.fq'):
            lib.Fzip_inplace(r_reads, args.cpus)
            r_reads = r_reads+'.gz'
        if not lib.checkannotations(os.path.join(tmpdir, 'left.fq.gz')):
            if os.path.dirname(os.path.abspath(tmpdir)) != os.path.dirname(os.path.abspath(l_reads)):
                try:
                    os.symlink(os.path.realpath(l_reads), os.path.join(tmpdir, 'left.fq.gz'))
                except OSError:
                    pass
        if not lib.checkannotations(os.path.join(tmpdir, 'right.fq.gz')):
            if os.path.dirname(os.path.abspath(tmpdir)) != os.path.dirname(os.path.abspath(r_reads)):
                try:
                    os.symlink(os.path.realpath(r_reads), os.path.join(tmpdir, 'right.fq.gz'))
                except OSError:
                    pass       
else:
    l_reads = os.path.join(tmpdir, 'left.fq.gz')
    r_reads = os.path.join(tmpdir, 'right.fq.gz')
    
#get tuple of input reads so you can parse them in downstream tools
all_reads = (l_reads, r_reads, s_reads)
lib.log.debug('Input reads: {:}'.format(all_reads))

#trimmomatic on reads, first run PE
if args.no_trimmomatic or args.trinity or args.left_norm or args.single_norm:
    lib.log.info("Trimmomatic will be skipped")
    trim_left = l_reads
    trim_right = r_reads
    trim_single = s_reads
else:
    #check if they exist already in folder
    if not os.path.isfile(os.path.join(tmpdir, 'trimmomatic', 'trimmed_left.fastq.gz')) or not os.path.isfile(os.path.join(tmpdir, 'trimmomatic', 'trimmed_right.fastq.gz')):
        if all_reads[0] and all_reads[1]:
            trim_left, trim_right = runTrimmomaticPE(l_reads, r_reads)
        else:
            trim_left, trim_right = (None,)*2
    else:
        trim_left, trim_right = os.path.join(tmpdir, 'trimmomatic', 'trimmed_left.fastq.gz'), os.path.join(tmpdir, 'trimmomatic', 'trimmed_right.fastq.gz')
    if not os.path.isfile(os.path.join(tmpdir, 'trimmomatic', 'trimmed_single.fastq.gz')) and s_reads:
        if all_reads[2]:
            trim_single = runTrimmomaticSE(s_reads)
        else:
            trim_single = None
    else:
        if s_reads:
            trim_single = os.path.join(tmpdir, 'trimmomatic', 'trimmed_single.fastq.gz')
        else:
            trim_single = None
#get tuple of trimmed reads
trim_reads = (trim_left, trim_right, trim_single)
lib.log.debug('Quality trimmed reads: {:}'.format(trim_reads))

#check that reads are present and make sure they follow trinity naming conventions, i.e. either illumina default or /1 and /2 to PE reads
for read in trim_reads:
    if read:
        if not os.path.isfile(read):
            lib.log.error("Trimmomatic failed, %s does not exist." % read)
            sys.exit(1)
if trim_reads[0] and trim_reads[1]: #PE reads are passed, lets make sure they have proper naming
    lib.CheckFASTQandFix(trim_reads[0], trim_reads[1]) #if needed to fix they will be fixed in place  

#normalize reads
left_norm, right_norm, single_norm = (None,)*3
if not os.path.isdir(os.path.join(tmpdir, 'normalize')):
    os.makedirs(os.path.join(tmpdir, 'normalize'))
if args.no_normalize_reads or args.trinity or args.left_norm or args.single_norm:
    lib.log.info("Read normalization will be skipped")
    if args.left_norm:
        left_norm = args.left_norm
        right_norm = args.right_norm 
        lib.SafeRemove(os.path.join(tmpdir, 'normalize', 'left.norm.fq')) 
        lib.SafeRemove(os.path.join(tmpdir, 'normalize', 'right.norm.fq'))
        if os.path.dirname(os.path.abspath(tmpdir)) != os.path.dirname(os.path.abspath(args.left_norm)):
            os.symlink(os.path.realpath(args.left_norm), os.path.join(tmpdir, 'normalize', 'left.norm.fq'))
        if os.path.dirname(os.path.abspath(tmpdir)) != os.path.dirname(os.path.abspath(args.right_norm)):
            os.symlink(os.path.realpath(args.right_norm), os.path.join(tmpdir, 'normalize', 'right.norm.fq'))  
    else:   
        left_norm = trim_left
        right_norm = trim_right
    if args.single_norm:
        single_norm = args.single_norm
        lib.SafeRemove(os.path.join(tmpdir, 'normalize', 'single.norm.fq'))
        if os.path.dirname(os.path.abspath(tmpdir)) != os.path.dirname(os.path.abspath(args.single_norm)):
            os.symlink(os.path.realpath(args.single_norm), os.path.join(tmpdir, 'normalize', 'single.norm.fq')) 
    else:
        single_norm = trim_single
else:
    #check if exists
    if trim_left and trim_right:
        if not os.path.islink(os.path.join(tmpdir, 'normalize', 'left.norm.fq')) or not os.path.islink(os.path.join(tmpdir, 'normalize', 'right.norm.fq')):
            if not all(v is None for v in trim_reads):
                left_norm, right_norm, single_norm = runNormalization(trim_reads, args.memory)
        else:
            left_norm, right_norm = os.path.join(tmpdir, 'normalize', 'left.norm.fq'), os.path.join(tmpdir, 'normalize', 'right.norm.fq')
            if os.path.islink(os.path.join(tmpdir, 'normalize', 'single.norm.fq')):
                single_norm = os.path.join(tmpdir, 'normalize', 'single.norm.fq')
    if trim_single:
        if not os.path.islink(os.path.join(tmpdir, 'normalize', 'single.norm.fq')) and not trim_left and not trim_right and trim_single:
            if not all(v is None for v in trim_reads):
                left_norm, right_norm, single_norm = runNormalization(trim_reads, args.memory)
        else:
            if os.path.islink(os.path.join(tmpdir, 'normalize', 'single.norm.fq')):
                single_norm = os.path.join(tmpdir, 'normalize', 'single.norm.fq')
            else:
                single_norm = None

#setup reads and check if normalization worked
norm_reads = (left_norm, right_norm, single_norm)
lib.log.debug('Normalized reads: {:}'.format(norm_reads))
if all(v is None for v in norm_reads):
    lib.log.error('No short reads detected, Trinity will be skipped.')
for read in norm_reads:
    if read:
        if not os.path.isfile(read):
            lib.log.error("Read normalization failed, %s does not exist." % read)
            sys.exit(1)

#check if long reads are passed, get full path
pb_iso, nano_cdna, nano_mrna = (None,)*3
if args.pacbio_isoseq:
    pb_iso = os.path.abspath(args.pacbio_isoseq)
if args.nanopore_cdna:
    nano_cdna = os.path.abspath(args.nanopore_cdna)
if args.nanopore_mrna:
    nano_mrna = os.path.abspath(args.nanopore_mrna)
long_reads = (pb_iso, nano_cdna, nano_mrna)
lib.log.debug('Long reads: {:}'.format(long_reads))

#get long read FASTA file
longReadFA = os.path.join(tmpdir, 'long-reads.fasta')
longReadClean = os.path.join(tmpdir, 'long-reads.fasta.clean')
PBiso = os.path.join(tmpdir, 'iso-seq.fasta')
nanocdna = os.path.join(tmpdir, 'nano-cdna.fasta')
nanomrna = os.path.join(tmpdir, 'nano-mrna.fasta')
long_clean = (None, None, None)
if not all(v is None for v in long_reads):
    if not lib.checkannotations(longReadFA):
        long_readsFA, long_clean = long2fasta(long_reads, args.cpus, tmpdir, os.path.abspath(longReadFA), os.path.abspath(longReadClean))
    else:
        found_clean = []
        for x in [PBiso, nanocdna, nanomrna]:
            if lib.checkannotations(x):
                found_clean.append(x)
            else:
                found_clean.append(None)
        long_clean = tuple(found_clean)
if not lib.checkannotations(longReadFA):
    longReadFA = None
lib.log.debug('Long reads FASTA format: {:}'.format(long_reads))
lib.log.debug('Long SeqCleaned reads: {:}'.format(long_clean))

#now run Trinity with trimmomatic and read normalization 
shortBAM = os.path.join(tmpdir, 'hisat2.coordSorted.bam')
trinity_transcripts = os.path.join(tmpdir, 'trinity.fasta')
if not lib.checkannotations(trinity_transcripts):
    if args.trinity:
        lib.log.info("Parsing assembled trinity data : {:}".format(args.trinity))
        shutil.copyfile(os.path.abspath(args.trinity), trinity_transcripts)
    else:
        if not all(v is None for v in norm_reads):
            #run trinity genome guided
            runTrinityGG(genome, norm_reads, longReadClean, shortBAM, trinity_transcripts)
else:
    lib.log.info("Existing Trinity results found: {:}".format(trinity_transcripts))

#if stringtie installed, run on shortBAM incorporate into PASA later on
stringtieGTF = os.path.join(tmpdir, 'funannotate_train.stringtie.gtf')
if lib.which('stringtie') and lib.checkannotations(shortBAM):
    lib.log.info('StringTie installed, running StringTie on Hisat2 coordsorted BAM')
    cmd = ['stringtie', '-p', str(args.cpus)]
    if args.stranded != 'no':
        if args.stranded.startswith('R'):
            cmd = cmd + ['--rf']
        else:
            cmd = cmd + ['--fr']
    cmd = cmd + [shortBAM]
    lib.runSubprocess8(cmd, '.', lib.log, stringtieGTF)
    
#run SeqClean to clip polyA tails and remove low quality seqs.
cleanTranscripts = os.path.join(tmpdir, 'trinity.fasta.clean')
if lib.checkannotations(trinity_transcripts):
    runSeqClean(trinity_transcripts, tmpdir)

#map long reads and Trinity transcripts to genome for PASA
allBAM = os.path.join(tmpdir, 'transcript.alignments.bam')
trinityBAM = os.path.join(tmpdir, 'trinity.alignments.bam')
if not lib.checkannotations(allBAM):
    trinity_transcripts, cleanTranscripts = mapTranscripts(genome, long_clean, cleanTranscripts, tmpdir, trinityBAM, allBAM)
else:
    if lib.checkannotations(trinityBAM):
        lib.log.info("Existing BAM alignments found: {:}, {:}".format(trinityBAM, allBAM))
    else:
        lib.log.info("Existing BAM alignments found: {:}".format(allBAM))
    
#convert BAM to GFF3
allGFF3 = os.path.join(tmpdir, 'transcript.alignments.gff3')
trinityGFF3 = os.path.join(tmpdir, 'trinity.alignments.gff3')
if not lib.checkannotations(allGFF3) and lib.checkannotations(allBAM):
    lib.log.info('Converting transcript alignments to GFF3 format')
    lib.bam2gff3(allBAM, allGFF3)
if not lib.checkannotations(trinityGFF3) and lib.checkannotations(trinityBAM):
    lib.log.info('Converting Trinity transcript alignments to GFF3 format')
    lib.bam2gff3(trinityBAM, trinityGFF3)

#now run PASA steps
PASA_gff = os.path.join(tmpdir, 'funannotate_train.pasa.gff3')
PASA_tmp = os.path.join(tmpdir, 'pasa.step1.gff3')
if not lib.checkannotations(PASA_tmp):
    if lib.checkannotations(trinityBAM):
        runPASAtrain(genome, trinity_transcripts, cleanTranscripts, os.path.abspath(trinityGFF3), stringtieGTF, args.stranded, args.max_intronlen, args.cpus, organism_name, PASA_tmp)
    elif lib.checkannotations(longReadFA): #no trinity seqs, so running PASA with only long reads
        runPASAtrain(genome, os.path.abspath(longReadFA), os.path.abspath(longReadClean), os.path.abspath(allGFF3), stringtieGTF, args.stranded, args.max_intronlen, args.cpus, organism_name, PASA_tmp)
else:
    lib.log.info("Existing PASA output found: {:}".format(PASA_tmp))

#Refine PASA models (there are many overlapping transcripts run kallisto and choose best model at each location)
KallistoAbundance = os.path.join(tmpdir, 'kallisto.tsv')
if all(v is None for v in trim_reads):
    kallistoreads = norm_reads
else:
    kallistoreads = trim_reads
if all(v is None for v in kallistoreads) and lib.checkannotations(longReadClean):
    #use minimap to count reads
    lib.log.info('Generating relative expression values to PASA transcripts')
    PASAtranscripts = os.path.join(tmpdir, 'transcripts.fa')
    cmd = [os.path.join(PASA, 'misc_utilities', 'gff3_file_to_proteins.pl'), PASA_tmp, genome, 'cDNA']
    if not lib.checkannotations(PASAtranscripts):
        lib.runSubprocess2(cmd, '.', lib.log, PASAtranscripts)
    PASAdict = pasa_transcript2gene(PASAtranscripts)
    minimapBAM = os.path.join(tmpdir, 'long-reads_transcripts.bam')
    minimap_cmd = ['minimap2', '-ax' 'map-ont', '-t', str(args.cpus), '--secondary=no',PASAtranscripts, longReadClean]
    cmd = [os.path.join(parentdir, 'util', 'sam2bam.sh'), " ".join(minimap_cmd), str(args.cpus // 2), minimapBAM]
    if not lib.checkannotations(minimapBAM):
        lib.runSubprocess(cmd, '.', lib.log)
    if not lib.checkannotations(KallistoAbundance):
        lib.mapCount(minimapBAM, PASAdict, KallistoAbundance)
else:
    if not lib.checkannotations(KallistoAbundance):
        runKallisto(PASA_tmp, genome, kallistoreads, args.stranded, args.cpus, os.path.join(tmpdir, 'getBestModel'), KallistoAbundance)
    else:
        lib.log.info("Existing Kallisto output found: {:}".format(KallistoAbundance))
    
#parse Kallisto results with PASA GFF
getBestModel(PASA_tmp, genome, KallistoAbundance, PASA_gff)

#collect final output files
BAMfinal = os.path.join(tmpdir, 'funannotate_train.coordSorted.bam')
TranscriptFinal = os.path.join(tmpdir, 'funannotate_train.trinity-GG.fasta')
LongFinal = os.path.join(tmpdir, 'funannotate_long-reads.fasta')
TranscriptAlignments = os.path.join(tmpdir, 'funannotate_train.transcripts.gff3')
#remove symlinks if from previous run
for x in [BAMfinal, TranscriptFinal, LongFinal, TranscriptAlignments]:
    lib.SafeRemove(x)
if lib.checkannotations(allBAM):
    os.symlink(os.path.realpath(allBAM), os.path.abspath(BAMfinal))
if longReadFA:
    os.symlink(os.path.realpath(longReadClean), os.path.abspath(LongFinal))
if lib.checkannotations(allGFF3):
    os.symlink(os.path.realpath(allGFF3), os.path.abspath(TranscriptAlignments))
if lib.checkannotations(trinity_transcripts):
    os.symlink(os.path.realpath(trinity_transcripts), os.path.abspath(TranscriptFinal))
lib.log.info('PASA database name: {:}'.format(organism_name.replace('-', '_')))
if args.strain:
    lib.log.info('Trinity/PASA has completed, you are now ready to run funanotate predict, for example:\n\n\
funannotate predict -i {:} \\\n\
        -o {:} -s "{:}" --strain {:} --cpus {:}\n'.format(args.input, args.out, organism, args.strain, args.cpus))
elif args.isolate:
    lib.log.info('Trinity/PASA has completed, you are now ready to run funanotate predict, for example:\n\n\
funannotate predict -i {:} \\\n\
        -o {:} -s "{:}" --isolate {:} --cpus {:}\n'.format(args.input, args.out, organism, args.isolate, args.cpus))
else:
    lib.log.info('Trinity/PASA has completed, you are now ready to run funanotate predict, for example:\n\n\
funannotate predict -i {:} \\\n\
        -o {:} -s "{:}" --cpus {:}\n'.format(args.input, args.out, organism, args.cpus))
print("-------------------------------------------------------")
sys.exit(1)
