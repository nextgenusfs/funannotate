#!/usr/bin/env python
from __future__ import division
import sys, os, subprocess, inspect, shutil, argparse, fnmatch, itertools
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
parser.add_argument('-i', '--input', help='Genome in GBK format or funannotate folder')
parser.add_argument('-g', '--gff', help='Genome annotation in GFF3 format')
parser.add_argument('-f', '--fasta', help='Genome sequence in FASTA format')
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

def validateCDSmRNAPairs(gene, cds, mrna, strand):
    '''
    function to make sure CDS is contained inside mRNA exons
    input is a list of lists of tuples for each
    return the correctly phased lists, only move the mRNA as CDS is tied to the id
    '''
    def _order(lst):
        return lst == sorted(lst) or lst == sorted(lst)[::-1]
    #load first mRNA exon into InterLap
    combined = []
    num = len(cds)
    warning = False
    for i in range(0,num):
        if strand == '+':
            sortCDS = sorted(cds[i], key=lambda tup: tup[0])
        else:
            sortCDS = sorted(cds[i], key=lambda tup: tup[0], reverse=True)
        compatible = []
        for x in range(0,num):
            if strand == '+':
                sortExon = sorted(mrna[x], key=lambda tup: tup[0])
            else:
                sortExon = sorted(mrna[x], key=lambda tup: tup[0], reverse=True)
            #simple first, if more cds than exons it is not compatible
            if len(sortCDS) > len(sortExon):
                compatible.append(False)
                continue
            result = True
            inter = InterLap(mrna[x])
            for i,coord in enumerate(sortCDS):
                if coord in inter:
                    hit = list(inter.find(coord))[0]
                    diff = np.subtract(coord,hit)
                    if diff[0] < 0 or diff[1] > 0: #then it cannot contain the cds so has to be wrong
                        result = False
                    if len(sortCDS) > 1:  
                        if i != 0 or (i+1) != len(sortCDS): #if an internal CDS, then must match perfectly or its wrong
                            if diff[0] != 0 and diff[1] != 0:
                                result = False
                        elif i == 0:
                            if strand == '+':
                                if diff[1] != 0:
                                    result = False
                            else:
                                if diff[0] != 0:
                                    result = False
                        elif (i+1) == len(sortCDS):
                            if strand == '+':    
                                if diff[0] != 0:
                                    return False
                            else:
                                if diff[1] != 0:
                                    return False             
            compatible.append(result) 
        combined.append(compatible)
    valid_orders = []
    for test in list(itertools.permutations(range(0,len(combined)),len(combined))):
        #test is a tuple, slice list to see if all True
        tester = []
        for num,x in enumerate(test):
            tester.append(combined[num][x])
        if all(tester):
            valid_orders.append(list(test))
    mRNA_order = valid_orders[0]      
    if not _order(mRNA_order):
        lib.log.debug('%s CDS/mRNA features out of phase, trying to fix. %s' % (gene, mRNA_order))
        if len(valid_orders) > 1:
            lib.log.debug('%s had %i possible solutions for CDS/mRNA, expect errors...' % (gene, len(valid_orders)))
            warning = True
    mRNAout = []
    for y in mRNA_order:
        mRNAout.append(mrna[y])
    return cds, mRNAout, warning

def gbk2pasaNEW(input, gff, trnaout, fastaout, spliceout, exonout, proteinsout):
    '''
    function is to parse a genbank flat file and move protein coding genes into GFF3
    and then parse out splice sites for hisat2. also filter tRNA gene models and move to
    new GFF3 file for adding back after PASA
    '''
    LocusTags = []
    multiExon = {}
    genes = {}
    with open(fastaout, 'w') as fasta:
        with open(input, 'rU') as gbk:
            for record in SeqIO.parse(gbk, 'genbank'):
                fasta.write(">%s\n%s\n" % (record.id, record.seq))
                for f in record.features:
                    lib.gb_feature_add2dict(f, record, genes)
    #out of order mRNA/CDS in genbank files can break this... so try to validate those with multiple transcripts
    warn = False
    for k,v in natsorted(genes.items()):
        if v['type'] == 'mRNA' and len(v['ids']) > 1:
            confirmedCDS,confirmedExons, warning = validateCDSmRNAPairs(k, v['CDS'], v['mRNA'], v['strand'])
            if warning:
                warn = True
            genes[k]['CDS'] = confirmedCDS
            genes[k]['mRNA'] = confirmedExons
    if warn:
        lib.log.info("GenBank file has multiple transcripts per locus, I tried my hardest to match them up but can't gaurantee there aren't errors. You can blame NCBI. You may want to try to pass a GFF3 + FASTA files instead of GBK.")
    with open(gff, 'w') as gffout:
        gffout.write('##gff-version 3\n')
        with open(trnaout, 'w') as trna:
            with open(proteinsout, 'w') as protout:
                for k,v in natsorted(genes.items()):
                    if not k in LocusTags:
                        LocusTags.append(k)
                    if v['type'] == 'mRNA':
                        #write GFF gene feature
                        if v['name']:
                            gffout.write("{:}\tGenBank\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(v['contig'], v['location'][0], v['location'][1], v['strand'], k, v['name']))
                        else:
                            gffout.write("{:}\tGenBank\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(v['contig'], v['location'][0], v['location'][1], v['strand'], k))
                        for i in range(0,len(v['ids'])):                 
                            #now write mRNA feature
                            gffout.write("{:}\tGenBank\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(v['contig'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i]))
                            protout.write('>%s %s\n%s\n' % (v['ids'][i], k, v['protein'][i]))
                            #write the exons and CDS features
                            num_exons = len(v['mRNA'][i])
                            for x in range(0,num_exons):
                                ex_num = x + 1
                                gffout.write("{:}\tGenBank\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(v['contig'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))                          
                                if num_exons > 1:
                                    #ss and exons are 0-based position, so 1 less than GFF
                                    exons_start = int(v['mRNA'][i][x][0]) - 1
                                    exons_end = int(v['mRNA'][i][x][1]) -1                                            
                                    #add to exon dictionary
                                    if not v['ids'][i] in multiExon:
                                        multiExon[v['ids'][i]] = [v['contig'], v['strand'], [(exons_start,exons_end)]]
                                    else:
                                        multiExon[v['ids'][i]][2].append((exons_start,exons_end))                               
                            num_cds = len(v['CDS'][i])
                            current_phase = v['codon_start'][i] - 1 #GFF3 phase is 1 less than flat file
                            for y in range(0,num_cds):
                                gffout.write("{:}\tGenBank\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(v['contig'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                                current_phase = (current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                                if current_phase == 3:
                                    current_phase = 0
                    elif v['type'] == 'tRNA' or v['type'] == 'rRNA':
                        #check length of tRNA gene should be between 50 and 150
                        if v['type'] == 'tRNA':
                            if v['strand'] == '+':
                                length = abs(int(v['location'][1]) - int(v['location'][0]))
                            else:
                                length = abs(int(v['location'][0]) - int(v['location'][1]))
                        else:
                            length = 100 #just a placeholder for rRNA features --> not sure if they have length requirements?
                        if length < 50 or length > 150:
                            continue
                        trna.write("{:}\tGenBank\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(v['contig'], v['location'][0], v['location'][1], v['strand'], k))
                        for i in range(0, len(v['ids'])):
                            trna.write("{:}\tGenBank\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(v['contig'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i]))
                            if v['type'] == 'tRNA':
                                num_exons = len(v['mRNA'][i])
                                for x in range(0,num_exons):
                                    ex_num = x + 1
                                    trna.write("{:}\tGenBank\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(v['contig'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))                   
    
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
    
def gff2pasa(gff_in, fasta, gff_out, trnaout, spliceout, exonout):
    '''
    function to parse GFF3 input file and split protein coding models from tRNA and/or rRNA
    models. Generate Hisat2 splice and exon files for mapping.
    '''
    #load into funanotate structured dictionary
    LocusTags = []
    multiExon = {}  
    genes = {}
    genes = lib.gff2dict(gff_in, fasta, genes)
    #now loop through dictionary and output desired files
    with open(gff_out, 'w') as gffout:
        gffout.write('##gff-version 3\n')
        with open(trnaout, 'w') as trna:
            for k,v in natsorted(genes.items()):
                if not k in LocusTags:
                    LocusTags.append(k)
                if v['type'] == 'mRNA':
                    #write GFF gene feature
                    if v['name']:
                        gffout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k, v['name']))
                    else:
                        gffout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
                    for i in range(0,len(v['ids'])):                 
                        #now write mRNA feature
                        gffout.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i]))
                        #write the exons and CDS features
                        num_exons = len(v['mRNA'][i])
                        for x in range(0,num_exons):
                            ex_num = x + 1
                            gffout.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))                          
                            if num_exons > 1:
                                #ss and exons are 0-based position, so 1 less than GFF
                                exons_start = int(v['mRNA'][i][x][0]) - 1
                                exons_end = int(v['mRNA'][i][x][1]) -1                                            
                                #add to exon dictionary
                                if not v['ids'][i] in multiExon:
                                    multiExon[v['ids'][i]] = [v['contig'], v['strand'], [(exons_start,exons_end)]]
                                else:
                                    multiExon[v['ids'][i]][2].append((exons_start,exons_end))                               
                        num_cds = len(v['CDS'][i])
                        current_phase = v['codon_start'][i] - 1 #GFF3 phase is 1 less than flat file
                        for y in range(0,num_cds):
                            gffout.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(v['contig'], v['source'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                            current_phase = (current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                            if current_phase == 3:
                                current_phase = 0
                elif v['type'] == 'tRNA' or v['type'] == 'rRNA':
                    #check length of tRNA gene should be between 50 and 150
                    if v['type'] == 'tRNA':
                        if v['strand'] == '+':
                            length = abs(int(v['location'][1]) - int(v['location'][0]))
                        else:
                            length = abs(int(v['location'][0]) - int(v['location'][1]))
                    else:
                        length = 100 #just a placeholder for rRNA features --> not sure if they have length requirements?
                    if length < 50 or length > 150:
                        continue
                    trna.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
                    for i in range(0, len(v['ids'])):
                        trna.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i]))
                        if v['type'] == 'tRNA':
                            num_exons = len(v['mRNA'][i])
                            for x in range(0,num_exons):
                                ex_num = x + 1
                                trna.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
    
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
    if bamthreads > 4:
        bamthreads = 4
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
    if bamthreads > 4:
        bamthreads = 4
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
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};Note=TPM:{:0.2f};\n'.format(cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], mRNAID, expression))
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
                         
def GFF2tblCombinedNEW(evm, genome, trnascan, prefix, genenumber, justify, SeqCenter, SeqRefNum, tblout):
    from collections import OrderedDict
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    '''
    function to take GFF3 annotation to produce a GBK tbl file, support multiple transcripts per locus.
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    #make sure genenumber is integer 
    genenumber = int(genenumber)
    
    #generate genome length dictionary used for feature tbl generation
    scaffLen = {}
    with open(genome, 'rU') as fastain:
        for record in SeqIO.parse(fastain, 'fasta'):
            if not record.id in scaffLen:
                scaffLen[record.id] = len(record.seq)
    #setup interlap database for genes on each chromosome and load EVM models into dictionary
    gene_inter = defaultdict(InterLap)
    Genes = {}
    gene_inter, Genes = lib.gff2interlapDict(evm, gene_inter, Genes)
    #now load tRNA predictions
    gene_inter, Genes = lib.gff2interlapDict(trnascan, gene_inter, Genes)
    #now sort dictionary by contig and location, rename using prefix
    sGenes = sorted(Genes.iteritems(), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    scaff2genes = {}
    SeqRecords = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
    inter = defaultdict(InterLap)
    skipList = []
    dropped = 0
    keeper = 0
    tooShort = 0
    internalStop = 0
    lib.log.info("Validating gene models (renaming, checking translations, filtering, etc)")
    for k,v in sortedGenes.items():
        GoodModel = True
        #check if gene model completely contained inside another one on same strand
        if args.alt_transcripts == '1':
            loc = sorted([v['location'][0],v['location'][1]])
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
        for i in range(0,len(v['ids'])):
            protSeq = None
            if v['type'] == 'mRNA': #get transcript for valid models
                cdsSeq = lib.getSeqRegions(SeqRecords, v['contig'], v['CDS'][i])
                protSeq = lib.translate(cdsSeq, v['strand'], v['codon_start'][i]-1)
                v['protein'].append(protSeq)
            if protSeq and len(protSeq) - 1 < 50:
                tooShort += 1
                continue
            if protSeq and '*' in protSeq[:-1]:
                internalStop += 1
                continue
            if protSeq:
                if protSeq.endswith('*'):
                    v['partialStop'][i] = False
                else: #try to extend the CDS up to 20 codons to see if you can find valid stop codon
                    '''
                    result, newCoords = lib.extend2stop(SeqRecords, v['contig'], v['CDS'][i], v['strand'], v['codon_start'][i]-1, len(protSeq))
                    if result == True:
                        print k, v
                        v['partialStop'][i] = False
                        v['CDS'][i] = newCoords
                        #update the location and mRNA if necessary
                        if v['strand'] == '+':
                            if newCoords[-1][1] > v['mRNA'][i][-1][1]:
                                v['mRNA'][i][-1] = (v['mRNA'][i][-1][0], newCoords[-1][1])
                            if newCoords[-1][1] > v['location'][1]:
                                v['location'] = (v['location'][0], newCoords[-1][1])
                        else:
                            if newCoords[-1][0] < v['mRNA'][i][-1][0]:
                                v['mRNA'][i][-1] = (newCoords[-1][0], v['mRNA'][i][-1][1])
                            if newCoords[-1][0] < v['location'][0]:
                                v['location'] = (newCoords[-1][0], v['location'][1])
                        print v
                    else:
                        v['partialStop'][i] = True
                    '''
                    v['partialStop'][i] = True
                if v['codon_start'][i] == 1 and v['protein'][i].startswith('M'):
                    v['partialStart'][i] = False
                else:
                    v['partialStart'][i] = True
        if not locusTag in renamedGenes:
            renamedGenes[locusTag] = v
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
                geneInfo = renamedGenes.get(genes) #single funannotate standard dictionary
                #check for partial models
                if True in geneInfo['partialStart']:
                    ps = '<'
                else:
                    ps = ''
                if True in geneInfo['partialStop']:
                    pss = '>'
                else:
                    pss = ''
                if geneInfo['type'] == 'rRNA':
                    ps = '<'
                    pss = '>'              
                #now write gene model
                if geneInfo['strand'] == '+':
                    tbl.write('%s%i\t%s%i\tgene\n' % (ps, geneInfo['location'][0], pss, geneInfo['location'][1]))
                    tbl.write('\t\t\tlocus_tag\t%s\n' % genes)
                else:
                    tbl.write('%s%i\t%s%i\tgene\n' % (ps, geneInfo['location'][1], pss, geneInfo['location'][0]))
                    tbl.write('\t\t\tlocus_tag\t%s\n' % genes)                                 
                #now will output the gene models with -T1, -T2, -T3 annotations based on expression values
                #means need to get the order
                order = []
                if len(geneInfo['ids']) > 1: #multiple transcripts, so get order of highest TPM
                    tpms = []
                    for num,tpm in enumerate(geneInfo['note']):
                        for item in tpm:
                            if item.startswith('TPM:'):
                                value = float(item.split(':')[-1])
                                tpms.append((value,num))
                    for x in sorted(tpms, reverse=True):
                        order.append(x[1])  
                else:
                    order.append(0)
                for i in order: #now write mRNA and CDS features
                    if geneInfo['type'] == 'mRNA':
                        if geneInfo['strand'] == '+':
                            if geneInfo['partialStart'][i] == False:
                                ps = ''
                            else:
                                ps = '<'
                            if geneInfo['partialStop'][i] == False:
                                pss = ''
                            else:
                                pss = '>'
                            for num, exon in enumerate(geneInfo['mRNA'][i]):
                                if num == 0 and num == len(geneInfo['mRNA'][i]) - 1: #single exon, so slightly differnt method
                                    tbl.write('%s%s\t%s%s\tmRNA\n' % (ps, exon[0], pss, exon[1]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tmRNA\n' % (ps, exon[0], exon[1]))
                                elif num == len(geneInfo['mRNA'][i]) - 1: #this is last one
                                    tbl.write('%s\t%s%s\n' % (exon[0], pss, exon[1]))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,i+1))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,i+1))
                            for num, cds in enumerate(geneInfo['CDS'][i]):
                                if num == 0 and num == len(geneInfo['CDS'][i]) - 1: #single exon, so slightly differnt method
                                    tbl.write('%s%s\t%s%s\tCDS\n' % (ps, cds[0], pss, cds[1]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tCDS\n' % (ps, cds[0], cds[1]))
                                elif num == len(geneInfo['CDS'][i]) - 1: #this is last one
                                    tbl.write('%s\t%s%s\n' % (cds[0], pss, cds[1]))
                                else:
                                    tbl.write('%s\t%s\n' % (cds[0], cds[1]))
                            tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'][i])
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,i+1))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,i+1))                                 
                        else: #means this is on crick strand
                            if geneInfo['partialStart'][i] == False:
                                ps = ''
                            else:
                                ps = '<'
                            if geneInfo['partialStop'][i] == False:
                                pss = ''
                            else:
                                pss = '>'               
                            for num, exon in enumerate(geneInfo['mRNA'][i]):
                                if num == 0 and num == len(geneInfo['mRNA'][i]) - 1: #single exon, so slightly differnt method
                                    tbl.write('%s%s\t%s%s\tmRNA\n' % (ps, exon[1], pss, exon[0]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tmRNA\n' % (ps, exon[1], exon[0]))
                                elif num == len(geneInfo['mRNA'][i]) - 1: #this is last one
                                    tbl.write('%s\t%s%s\n' % (exon[1], pss, exon[0]))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[1], exon[0]))                 
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,i+1))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,i+1))
                            for num, cds in enumerate(geneInfo['CDS'][i]):
                                if num == 0 and num == len(geneInfo['CDS'][i]) - 1: #single exon, so slightly differnt method
                                    tbl.write('%s%s\t%s%s\tCDS\n' % (ps, cds[1], pss, cds[0]))
                                elif num == 0:
                                    tbl.write('%s%s\t%s\tCDS\n' % (ps, cds[1], cds[0]))
                                elif num == (len(geneInfo['CDS'][i]) - 1): #this is last one
                                    tbl.write('%s\t%s%s\n' % (cds[1], pss, cds[0]))
                                else:
                                    tbl.write('%s\t%s\n' % (cds[1], cds[0]))
                            tbl.write('\t\t\tcodon_start\t%i\n' % geneInfo['codon_start'][i])
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                            tbl.write('\t\t\ttranscript_id\tgnl|ncbi|%s-T%i_mrna\n' % (genes,i+1))
                            tbl.write('\t\t\tprotein_id\tgnl|ncbi|%s-T%i\n' % (genes,i+1))
                    elif geneInfo['type'] == 'tRNA':
                        if geneInfo['strand'] == '+':
                            for num, exon in enumerate(geneInfo['mRNA'][i]):
                                if num == 0:
                                    tbl.write('<%s\t>%s\t%s\n' % (exon[0], exon[1], geneInfo['type']))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[0], exon[1]))
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tbl.write('\t\t\tpseudo\n')        
                        else:
                            for num, exon in enumerate(geneInfo['mRNA'][i]):
                                if num == 0:
                                    tbl.write('<%s\t>%s\t%s\n' % (exon[1], exon[0], geneInfo['type']))
                                else:
                                    tbl.write('%s\t%s\n' % (exon[1], exon[0]))
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])
                            if geneInfo['product'] == 'tRNA-Xxx':
                                tbl.write('\t\t\tpseudo\n')
                    elif geneInfo['type'] == 'rRNA':
                        if geneInfo['strand'] == '+':
                            tbl.write('<%s\t>%s\t%s\n' % (geneInfo['location'][0],geneInfo['location'][1], geneInfo['type']))
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])   
                        else:
                            tbl.write('<%s\t>%s\t%s\n' % (geneInfo['location'][1],geneInfo['location'][0], geneInfo['type']))
                            tbl.write('\t\t\tproduct\t%s\n' % geneInfo['product'][i])

def gbk2interlap(input):
    '''
    function to parse GBK file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    '''
    inter = defaultdict(InterLap)
    Genes = {}
    with open(input, 'rU') as filein:
        for record in SeqIO.parse(filein, 'genbank'):
            for f in record.features:
                if f.type == 'gene':
                    locusTag, ID, Parent = lib.getID(f, f.type)
                    start = int(f.location.nofuzzy_start)
                    end = int(f.location.nofuzzy_end)
                    inter[record.id].add((start,end,locusTag))      
                lib.gb_feature_add2dict(f, record, Genes)
    return inter, Genes

def gff2interlap(input, fasta):
    '''
    function to parse GFF3 file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    '''
    inter = defaultdict(InterLap)
    Genes = {}
    Genes = lib.gff2dict(input, fasta, Genes)
    for k,v in natsorted(Genes.items()):
        inter[v['contig']].add((v['location'][0],v['location'][1],k))
    return inter, Genes
        
   
def merge_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z
    
def message(loc1, loc2, cdsAED, mrnaAED, protMatches, UTRs, no_change, UTR_added, yardSale, exonChange):
    msg = []
    if not cdsAED or cdsAED == '':
        cds = 0
    else:
        cds = float(cdsAED)
    mrna = float(mrnaAED)
    pos = loc1 == loc2
    #structured message, coordinates, 5prime, 3prime, exon, cds, pident
    if not pos: #coordinates changed
        msg.append('gene coordinates updated')
        for u in UTRs:
            if u[0]:
                if not any('5prime' in x for x in msg):
                    msg.append('5prime UTR added')
            if u[1]:
                if not any('3prime' in x for x in msg):
                    msg.append('3prime UTR added')
    if mrna > 0:
        msg.append('mRNA updated')
    if cds > 0:
        pidentmsg = []
        for x in protMatches:
            pidentmsg.append('{0:.0f}%'.format(x))
        msg.append('CDS update [translation pident: %s]' % ', '.join(pidentmsg))
    #now work on the counter
    if len(msg) < 1:
        msg = ['no change']
        no_change += 1
    elif any('UTR' in x for x in msg):
        UTR_added += 1
    elif any('mRNA' in x for x in msg):
        exonChange += 1
    else:
        yardSale += 1
    final_message = ';'.join(msg)
    return final_message, no_change, UTR_added, yardSale, exonChange
  
def pairwiseAlign(query, ref):
    from Bio import pairwise2
    '''
    do global alignment and return pident
    '''
    if query == ref:
        return 100.0
    align = pairwise2.align.globalxx(query, ref)
    length = max(len(query), len(ref))
    pident = (align[0][2] / float(length)) * 100
    return pident

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
    if args.gff and args.fasta:
        oldInter, oldGenes = gff2interlap(old, args.fasta)
    else:
        oldInter, oldGenes = gbk2interlap(old)
    newInter, newGenes = gbk2interlap(new)
    #do the simple stuff first, find models that were deleted
    for contig in oldInter:
        for gene in oldInter[contig]:
            if not gene in newInter[contig]: #these models are removed
                dropped += 1
                if not gene[2] in oldGenes:
                    continue
                #populate output dictionary with results
                if not gene[2] in result:
                    #dropped model has AED of 1.000
                    cdsAED = '1.000'
                    exonAED = '1.000'
                    result[gene[2]] = {'contig': oldGenes[gene[2]]['contig'], 'old_num_transcripts': len(oldGenes[gene[2]]['ids']), 
                                'old_location': oldGenes[gene[2]]['location'], 'num_transcripts': len(oldGenes[gene[2]]['ids']), 'strand': oldGenes[gene[2]]['strand'], 
                                'mRNA': oldGenes[gene[2]]['mRNA'], 'location': oldGenes[gene[2]]['location'], 'CDS': oldGenes[gene[2]]['CDS'], 'message': 'gene model removed', 
                                'cdsAED': cdsAED, 'exonAED': exonAED, 'transcript_id': oldGenes[gene[2]]['ids'], 'pident': [],
                                'protein_id':oldGenes[gene[2]]['ids'], 'seq': oldGenes[gene[2]]['protein']}

    #now go through the updated annotation, comparing to old annot
    for contig in newInter:
        for gene in newInter[contig]:
            if not gene in oldInter[contig]: #means this is a new model, so add it
                added += 1
                total_genes += 1
                if not gene[2] in newGenes:
                    continue
                total_transcripts += len(newGenes[gene[2]]['ids'])
                if not gene[2] in result:
                    result[gene[2]] = {'contig': newGenes[gene[2]]['contig'], 'old_num_transcripts': 0, 
                                'old_location': newGenes[gene[2]]['location'], 'num_transcripts': len(newGenes[gene[2]]['ids']), 'strand': newGenes[gene[2]]['strand'], 
                                'mRNA': newGenes[gene[2]]['mRNA'], 'location': newGenes[gene[2]]['location'], 'CDS': newGenes[gene[2]]['CDS'], 'message': 'new gene model', 
                                'cdsAED': '0.000', 'exonAED': '0.000', 'transcript_id': newGenes[gene[2]]['ids'], 
                                'protein_id':newGenes[gene[2]]['ids'], 'seq': newGenes[gene[2]]['protein'], 'pident': []}
            else: #means this is existing model, and need to do some comparisons
                hitList = list(oldInter[contig].find(gene))
                #there might be some overlapping transcripts, so enforce locus name
                hit = None
                for z in hitList:
                    if gene[2] == z[2]:
                        hit = z
                if not hit:
                    #there is no real hit, so this a new gene
                    total_transcripts += len(newGenes[gene[2]]['ids'])
                    added += 1
                    total_genes += 1
                    if not gene[2] in result:
                        result[gene[2]] = {'contig': newGenes[gene[2]]['contig'], 'old_num_transcripts': 0, 
                                    'old_location': newGenes[gene[2]]['location'], 'num_transcripts': len(newGenes[gene[2]]['ids']), 'strand': newGenes[gene[2]]['strand'], 
                                    'mRNA': newGenes[gene[2]]['mRNA'], 'location': newGenes[gene[2]]['location'], 'CDS': newGenes[gene[2]]['CDS'], 'message': 'new gene model', 
                                    'cdsAED': '0.000', 'exonAED': '0.000', 'transcript_id': newGenes[gene[2]]['ids'], 
                                    'protein_id':newGenes[gene[2]]['ids'], 'seq': newGenes[gene[2]]['protein'], 'pident': []}
                else:
                    #since we may have multiple transcripts from hit as well as new annotation we need to be aware of that
                    #also, tRNA annotations do not exist in Proteins dictionary, so process them differently
                    #get the reference hits, pull out CDS and mRNA for pairwiseAED calculation
                    total_genes += 1
                    total_transcripts += len(newGenes[gene[2]]['ids'])
                    
                    #get the old annotation
                    hitInfo = oldGenes.get(gene[2])
                    
                    #calculate AED
                    exonAED = pairwiseAED(newGenes[gene[2]]['mRNA'], hitInfo['mRNA'])
                    if newGenes[gene[2]]['type'] == 'mRNA' and hitInfo['type'] == 'mRNA':
                        cdsAED = pairwiseAED(newGenes[gene[2]]['CDS'], hitInfo['CDS'])
                    else:
                        cdsAED = '0.000'
                    
                    #check translation, to deal with multiple transcripts, lets loop through new
                    protMatches = []
                    if newGenes[gene[2]]['type'] == 'mRNA' and hitInfo['type'] == 'mRNA':
                        for i in range(0,len(newGenes[gene[2]]['ids'])):
                            protMatch = None
                            for y in range(0,len(oldGenes[gene[2]]['ids'])):
                                pident = pairwiseAlign(newGenes[gene[2]]['protein'][i], oldGenes[gene[2]]['protein'][y])
                                if not protMatch:
                                    protMatch = pident
                                else:
                                    if pident > protMatch:
                                        protMatch = pident
                            protMatches.append(protMatch)
                    #summarize UTRs
                    UTRs = findUTRs(newGenes[gene[2]]['CDS'], newGenes[gene[2]]['mRNA'], newGenes[gene[2]]['strand'])
                    
                    #structured comments/counts for gene models
                    msg, no_change, UTR_added, yardSale, exonChange = message(newGenes[gene[2]]['location'], oldGenes[gene[2]]['location'], cdsAED, exonAED, protMatches, UTRs, no_change, UTR_added, yardSale, exonChange)
                    
                    if not gene[2] in result:
                        result[gene[2]] = {'contig': newGenes[gene[2]]['contig'], 'old_num_transcripts': len(oldGenes[gene[2]]['ids']), 
                                    'old_location': oldGenes[gene[2]]['location'], 'num_transcripts': len(newGenes[gene[2]]['ids']), 'strand': newGenes[gene[2]]['strand'], 
                                    'mRNA': newGenes[gene[2]]['mRNA'], 'location': newGenes[gene[2]]['location'], 'CDS': newGenes[gene[2]]['CDS'], 'message': msg, 
                                    'cdsAED': cdsAED, 'exonAED': exonAED, 'transcript_id': newGenes[gene[2]]['ids'], 
                                    'protein_id':newGenes[gene[2]]['ids'], 'seq': newGenes[gene[2]]['protein'], 'pident': protMatches}
         
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
Dropped Models:\t\t{:,}\n\
CDS AED:\t\t{:.3f}\n\
mRNA AED:\t\t{:.3f}\n\
-------------------------------------------------------".format(total_genes, total_transcripts, added, no_change, UTR_added, exonChange, yardSale, dropped, Avg_cdsAED, Avg_exonAED))

def findUTRs(cds, mrna, strand):
    '''
    take list of list of CDS coordiantes and compare to list of list of mRNA coordinates to
    determine if 5 prime or 3 prime UTR exist
    '''
    #supporting multiple transcripts, however, they are already matched up and sorted
    UTRs = []
    for i in range(0, len(cds)):
        Fiveprime = False
        Threeprime = False
        refInterlap = InterLap(mrna[i])
        if strand == '+': #look at first CDS for 5 prime and last CDS for 3 prime
            if cds[i][0] in refInterlap: #means it overlaps with mrNA (which it obviously should)
                hit = list(refInterlap.find(cds[i][0]))[0]
                loc = mrna[i].index(hit) #if first exon, then compare, if not first then there is 5prime UTR
                if loc == 0:            
                    diff = np.subtract(cds[i][0],hit) #will return array of exon minus hit at each pos
                    if diff[0] > 0:
                        Fiveprime = True
                else:
                    Fiveprime = True
            #check for 3 prime UTR
            if cds[i][-1] in refInterlap:
                hit = list(refInterlap.find(cds[i][-1]))[0]
                loc = mrna[i].index(hit)
                if len(mrna[i]) == loc+1:
                    diff = np.subtract(cds[i][-1],hit) #will return array of exon minus hit at each pos
                    if diff[1] < 0:
                        Threeprime = True
                else:
                    Threeprime = True
        else:
            if cds[i][0] in refInterlap: #means it overlaps with mrNA (which it obviously should)
                hit = list(refInterlap.find(cds[i][0]))[0]
                loc = mrna[i].index(hit) #if first exon, then compare, if not first then there is 5prime UTR
                if loc == 0:       
                    diff = np.subtract(cds[i][0],hit) #will return array of exon minus hit at each pos
                    if diff[1] < 0:
                        Fiveprime = True
                else:
                    Fiveprime = True
            #check for 3 prime UTR
            if cds[i][-1] in refInterlap:
                hit = list(refInterlap.find(cds[i][-1]))[0]
                loc = mrna[i].index(hit)
                if len(mrna[i]) == loc+1:
                    diff = np.subtract(cds[i][-1],hit) #will return array of exon minus hit at each pos
                    if diff[0] > 0:
                        Threeprime = True
                else:
                    Threeprime = True     
        UTRs.append((Fiveprime,Threeprime))
    return UTRs

def pairwiseAED(query, reference):
    '''
    takes a multiple transcripts and sums AED from lowest pairwise comparison and then calculates
    the average based on number of transcripts in the query
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
    #check if identical
    if query == reference:
       return '0.000'
    #make sure sorted
    rLen = _length(reference)
    refInterlap = InterLap(reference)
    RefPred = 0
    QueryOverlap = 0
    qLen = 0
    for exon in query:
        qLen += abs(exon[0] - exon[1])
        if exon in refInterlap: #exon overlaps at least partially with reference
            hit = list(refInterlap.find(exon))
            for h in hit:
                diff = np.subtract(exon,h) #will return array of exon minus hit at each pos
                if diff[0] <= 0 and diff[1] >= 0: #then query exon covers ref exon
                    cov = abs(h[0] - h[1])
                    QueryOverlap += cov
                elif diff[0] <= 0 and diff[1] < 0: # means query partial covers ref
                    cov = abs(h[0] - exon[1])
                    QueryOverlap += cov
                elif diff[0] > 0 and diff[1] >= 0: #means query partial covers ref
                    cov = abs(exon[0] - h[1])
                    QueryOverlap += cov        
                elif diff[0] > 0 and diff[1] < 1:
                    cov = abs(exon[0] - exon[1])
                    QueryOverlap += cov
    #calculate AED
    SP = QueryOverlap / float(qLen)
    SN = QueryOverlap / float(rLen)
    AED = 1 - ((SN + SP) / 2)
    return '{:.3f}'.format(AED)
    
#create folder structure
if args.input:
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
        
programs = ['fasta', 'gmap', 'blat', 'tbl2asn', 'hisat2', 'hisat2-build', 'kallisto', 'Trinity', 'bedtools', 'java']
lib.CheckDependencies(programs)

#take care of some preliminary checks
if args.sbt == 'SBT':
    SBT = os.path.join(parentdir, 'lib', 'test.sbt')
    lib.log.info("No NCBI SBT file given, will use default, for NCBI submissions pass one here '--sbt'")
else:
    SBT = args.sbt

#setup output files
gffout = os.path.join(tmpdir, 'genome.gff3')
proteinsout = os.path.join(tmpdir, 'genome.proteins.fa')
trnaout = os.path.join(tmpdir, 'genome.trna.gff3')
fastaout = os.path.join(tmpdir, 'genome.fa')
spliceout = os.path.join(tmpdir, 'genome.ss')
exonout = os.path.join(tmpdir, 'genome.exons')

#check input, allow for passing the output directory of funannotate, otherwise must be gbk or gbff files
#set read inputs to None, populate as you go
s_reads, l_reads, r_reads, trim_left, trim_right, trim_single, left_norm, right_norm, single_norm, all_reads, trim_reads, norm_reads, GBK, trinity_results, pasaConfigFile = (None,)*15
if args.input:
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
        #let user know which files are being re-used
        look4files = [s_reads, l_reads, r_reads, trim_left, trim_right, trim_single, left_norm, right_norm, single_norm, trinity_results, pasaConfigFile]
        look4files_names = ['Single reads', 'Forward reads', 'Reverse reads', 'Forward Q-trimmed reads', 'Reverse Q-trimmed reads', 'Single Q-trimmed reads', 'Forward normalized reads', 'Reverse normalized reads', 'Single normalized reads', 'Trinity results', 'PASA config file']
        files_used = []
        for i,x in enumerate(look4files):
            if x is not None:
                files_used.append('\t'+look4files_names[i]+': '+x)
        #files_used = [x for x in [s_reads, l_reads, r_reads, trim_left, trim_right, trim_single, left_norm, right_norm, single_norm, trinity_results, pasaConfigFile] if x is not None]
        if len(files_used) > 0:
            lib.log.info('Found relevent files in %s, will re-use them:\n%s' % (inputDir, '\n'.join(files_used)))
    else:
        GBK = args.input
    #check if RefSeq --> NCBI does not want you to reannotate RefSeq genomes
    if lib.checkRefSeq(GBK):
        lib.log.error('%s is a NCBI RefSeq genome, to reannotate please use original submission.' % GBK)
        sys.exit(1)
    #split GenBank into parts
    locustag, genenumber, justify = gbk2pasaNEW(GBK, gffout, trnaout, fastaout, spliceout, exonout, proteinsout)
    organism, strain, isolate, accession, WGS_accession, gb_gi, version = lib.getGBKinfo(GBK)
    lib.log.info("Reannotating %s, NCBI accession: %s" % (organism, WGS_accession))

else:
    if args.gff and args.fasta:
        if not args.species:
            lib.log.error("Input error: please enter a name for -s,--species")
            sys.exit(1)
        shutil.copyfile(args.fasta, fastaout)
        locustag, genenumber, justify = gff2pasa(args.gff, fastaout, gffout, trnaout, spliceout, exonout)
        organism,strain,isolate,accession,WGS_accession,gb_gi,version = (None,)*7
    else:
        lib.log.error("Error in input: pass either funannotate directory or GenBank file to -i,--input; or GFF3 to -g,--gff and genome FASTA to -f,--fasta.")
        sys.exit(1)
        
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
        #check if files already exist
        #check if exists
        if not os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'left.norm.fq')) or not os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'right.norm.fq')):
            if trim_reads[0] and trim_reads[1]:
                lib.log.info("Running read normalization with Trinity")
                left_norm, right_norm, single_norm = runNormalization(trim_reads, args.memory)
            else:
                left_norm, right_norm = (None,)*2
        else:
            left_norm, right_norm = os.path.join(args.out, 'update_misc', 'normalize', 'left.norm.fq'), os.path.join(args.out, 'update_misc', 'normalize', 'right.norm.fq')
            if os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'single.norm.fq')):
                single_norm = os.path.join(args.out, 'update_misc', 'normalize', 'single.norm.fq')
        if not os.path.isfile(os.path.join(args.out, 'update_misc', 'normalize', 'single.norm.fq')) and not trim_reads[0] and not trim_reads[1] and trim_single:
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
    GFF2tblCombinedNEW(BestModelGFF, fastaout, cleanTRNA, locustag, genenumber, justify, args.SeqCenter, args.SeqAccession, TBLFile)

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
if args.gff and args.fasta:
    compareAnnotations2(args.gff, final_gbk, Changes)
else:
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
