#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os
import subprocess
import shutil
import argparse
import itertools
from Bio import SeqIO
import funannotate.library as lib
from funannotate.interlap import InterLap
from collections import defaultdict
from natsort import natsorted
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def validateCDSmRNAPairs(gene, cds, mrna, strand):
    '''
    function to make sure CDS is contained inside mRNA exons
    input is a list of lists of tuples for each
    return the correctly phased lists, only move the mRNA as CDS is tied to the id
    '''
    def _order(lst):
        return lst == sorted(lst) or lst == sorted(lst)[::-1]
    # load first mRNA exon into InterLap
    combined = []
    num = len(cds)
    warning = False
    for i in range(0, num):
        if strand == '+':
            sortCDS = sorted(cds[i], key=lambda tup: tup[0])
        else:
            sortCDS = sorted(cds[i], key=lambda tup: tup[0], reverse=True)
        compatible = []
        for x in range(0, num):
            if strand == '+':
                sortExon = sorted(mrna[x], key=lambda tup: tup[0])
            else:
                sortExon = sorted(
                    mrna[x], key=lambda tup: tup[0], reverse=True)
            # simple first, if more cds than exons it is not compatible
            if len(sortCDS) > len(sortExon):
                compatible.append(False)
                continue
            result = True
            inter = InterLap(mrna[x])
            for i, coord in enumerate(sortCDS):
                if coord in inter:
                    hit = list(inter.find(coord))[0]
                    diff = np.subtract(coord, hit)
                    # then it cannot contain the cds so has to be wrong
                    if diff[0] < 0 or diff[1] > 0:
                        result = False
                    if len(sortCDS) > 1:
                        # if an internal CDS, then must match perfectly or its wrong
                        if i != 0 or (i+1) != len(sortCDS):
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
    for test in list(itertools.permutations(list(range(0, len(combined))), len(combined))):
        # test is a tuple, slice list to see if all True
        tester = []
        for num, x in enumerate(test):
            tester.append(combined[num][x])
        if all(tester):
            valid_orders.append(list(test))
    mRNA_order = valid_orders[0]
    if not _order(mRNA_order):
        lib.log.debug(
            '%s CDS/mRNA features out of phase, trying to fix. %s' % (gene, mRNA_order))
        if len(valid_orders) > 1:
            lib.log.debug(
                '%s had %i possible solutions for CDS/mRNA, expect errors...' % (gene, len(valid_orders)))
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
        with open(input, 'r') as gbk:
            for record in SeqIO.parse(gbk, 'genbank'):
                fasta.write(">%s\n%s\n" % (record.id, record.seq))
                for f in record.features:
                    lib.gb_feature_add2dict(f, record, genes)
    # out of order mRNA/CDS in genbank files can break this... so try to validate those with multiple transcripts
    warn = False
    for k, v in natsorted(list(genes.items())):
        if v['type'] == 'mRNA' and len(v['ids']) > 1:
            confirmedCDS, confirmedExons, warning = validateCDSmRNAPairs(
                k, v['CDS'], v['mRNA'], v['strand'])
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
                for k, v in natsorted(list(genes.items())):
                    if not k in LocusTags:
                        LocusTags.append(k)
                    if v['type'] == 'mRNA':
                        # write GFF gene feature
                        if v['name']:
                            gffout.write("{:}\tGenBank\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(
                                v['contig'], v['location'][0], v['location'][1], v['strand'], k, v['name']))
                        else:
                            gffout.write("{:}\tGenBank\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                                v['contig'], v['location'][0], v['location'][1], v['strand'], k))
                        for i in range(0, len(v['ids'])):
                            # now write mRNA feature
                            gffout.write("{:}\tGenBank\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(
                                v['contig'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i]))
                            protout.write('>%s %s\n%s\n' %
                                          (v['ids'][i], k, v['protein'][i]))
                            # write the exons and CDS features
                            num_exons = len(v['mRNA'][i])
                            for x in range(0, num_exons):
                                ex_num = x + 1
                                gffout.write("{:}\tGenBank\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                                    v['contig'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
                                if num_exons > 1:
                                    # ss and exons are 0-based position, so 1 less than GFF
                                    exons_start = int(v['mRNA'][i][x][0]) - 1
                                    exons_end = int(v['mRNA'][i][x][1]) - 1
                                    # add to exon dictionary
                                    if not v['ids'][i] in multiExon:
                                        multiExon[v['ids'][i]] = [
                                            v['contig'], v['strand'], [(exons_start, exons_end)]]
                                    else:
                                        multiExon[v['ids'][i]][2].append(
                                            (exons_start, exons_end))
                            num_cds = len(v['CDS'][i])
                            # GFF3 phase is 1 less than flat file
                            current_phase = v['codon_start'][i] - 1
                            for y in range(0, num_cds):
                                gffout.write("{:}\tGenBank\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(
                                    v['contig'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                                current_phase = (
                                    current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                                if current_phase == 3:
                                    current_phase = 0
                    elif v['type'] in ['tRNA', 'rRNA', 'ncRNA']:
                        # check length of tRNA gene should be between 50 and 150
                        if v['type'] == 'tRNA':
                            if v['strand'] == '+':
                                length = abs(
                                    int(v['location'][1]) - int(v['location'][0]))
                            else:
                                length = abs(
                                    int(v['location'][0]) - int(v['location'][1]))
                        else:
                            length = 100  # just a placeholder for rRNA features --> not sure if they have length requirements?
                        if length < 50 or length > 150:
                            continue
                        trna.write("{:}\tGenBank\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                            v['contig'], v['location'][0], v['location'][1], v['strand'], k))
                        for i in range(0, len(v['ids'])):
                            trna.write("{:}\tGenBank\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(
                                v['contig'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i]))
                            if v['type'] == 'tRNA':
                                num_exons = len(v['mRNA'][i])
                                for x in range(0, num_exons):
                                    ex_num = x + 1
                                    trna.write("{:}\tGenBank\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                                        v['contig'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))

    # parse splice sites and write to file
    with open(exonout, 'w') as exon:
        with open(spliceout, 'w') as splicer:
            for k, v in natsorted(list(multiExon.items())):
                sortedList = sorted(v[2], key=lambda tup: tup[0])
                for y in sortedList:
                    exon.write("%s\t%i\t%i\t%s\n" % (v[0], y[0], y[1], v[1]))
                splices = []
                for i in range(1, len(sortedList)):
                    splices.append((sortedList[i-1][1], sortedList[i][0]))
                for x in splices:
                    splicer.write('%s\t%i\t%i\t%s\n' %
                                  (v[0], x[0], x[1], v[1]))
    # finally lets return the base locus tag name and the last number
    lastTag = natsorted(LocusTags)[-1]
    if '_' in lastTag:
        tag, count = lastTag.split('_')
        tag = tag+'_'
    else:
        for i, c in enumerate(lastTag):
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
    # load into funanotate structured dictionary
    LocusTags = []
    multiExon = {}
    genes = {}
    genes = lib.gff2dict(gff_in, fasta, genes)
    # now loop through dictionary and output desired files
    with open(gff_out, 'w') as gffout:
        gffout.write('##gff-version 3\n')
        with open(trnaout, 'w') as trna:
            for k, v in natsorted(list(genes.items())):
                if not k in LocusTags:
                    LocusTags.append(k)
                if v['type'] == 'mRNA':
                    # write GFF gene feature
                    if v['name']:
                        gffout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(
                            v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k, v['name']))
                    else:
                        gffout.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                            v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
                    for i in range(0, len(v['ids'])):
                        # now write mRNA feature
                        gffout.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(
                            v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i]))
                        # write the exons and CDS features
                        num_exons = len(v['mRNA'][i])
                        for x in range(0, num_exons):
                            ex_num = x + 1
                            gffout.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                                v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))
                            if num_exons > 1:
                                # ss and exons are 0-based position, so 1 less than GFF
                                exons_start = int(v['mRNA'][i][x][0]) - 1
                                exons_end = int(v['mRNA'][i][x][1]) - 1
                                # add to exon dictionary
                                if not v['ids'][i] in multiExon:
                                    multiExon[v['ids'][i]] = [
                                        v['contig'], v['strand'], [(exons_start, exons_end)]]
                                else:
                                    multiExon[v['ids'][i]][2].append(
                                        (exons_start, exons_end))
                        num_cds = len(v['CDS'][i])
                        # GFF3 phase is 1 less than flat file
                        current_phase = v['codon_start'][i] - 1
                        for y in range(0, num_cds):
                            gffout.write("{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(
                                v['contig'], v['source'], v['CDS'][i][y][0], v['CDS'][i][y][1], v['strand'], current_phase, v['ids'][i], v['ids'][i]))
                            current_phase = (
                                current_phase - (int(v['CDS'][i][y][1]) - int(v['CDS'][i][y][0]) + 1)) % 3
                            if current_phase == 3:
                                current_phase = 0
                elif v['type'] in ['tRNA', 'rRNA', 'ncRNA']:
                    # check length of tRNA gene should be between 50 and 150
                    if v['type'] == 'tRNA':
                        if v['strand'] == '+':
                            length = abs(
                                int(v['location'][1]) - int(v['location'][0]))
                        else:
                            length = abs(
                                int(v['location'][0]) - int(v['location'][1]))
                    else:
                        length = 100  # just a placeholder for rRNA features --> not sure if they have length requirements?
                    if length < 50 or length > 150:
                        continue
                    trna.write("{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                        v['contig'], v['source'], v['location'][0], v['location'][1], v['strand'], k))
                    for i in range(0, len(v['ids'])):
                        trna.write("{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};\n".format(
                            v['contig'], v['source'], v['type'], v['location'][0], v['location'][1], v['strand'], v['ids'][i], k, v['product'][i]))
                        if v['type'] == 'tRNA':
                            num_exons = len(v['mRNA'][i])
                            for x in range(0, num_exons):
                                ex_num = x + 1
                                trna.write("{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                                    v['contig'], v['source'], v['mRNA'][i][x][0], v['mRNA'][i][x][1], v['strand'], v['ids'][i], ex_num, v['ids'][i]))

    # parse splice sites and write to file
    with open(exonout, 'w') as exon:
        with open(spliceout, 'w') as splicer:
            for k, v in natsorted(list(multiExon.items())):
                sortedList = sorted(v[2], key=lambda tup: tup[0])
                for y in sortedList:
                    exon.write("%s\t%i\t%i\t%s\n" % (v[0], y[0], y[1], v[1]))
                splices = []
                for i in range(1, len(sortedList)):
                    splices.append((sortedList[i-1][1], sortedList[i][0]))
                for x in splices:
                    splicer.write('%s\t%i\t%i\t%s\n' %
                                  (v[0], x[0], x[1], v[1]))
    # finally lets return the base locus tag name and the last number
    lastTag = natsorted(LocusTags)[-1]
    if '_' in lastTag:
        tag, count = lastTag.split('_')
        tag = tag+'_'
        try:
            count = int(count)
        except ValueError: #means it is not a number, so then count gens
            count = len(LocusTags) + 1
    else:
        for i, c in enumerate(lastTag):
            if c.isdigit():
                tag = lastTag[:i]
                count = lastTag[i:]
                break
    count = str(count)
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


def runTrimmomaticPE(left, right, cpus=1):
    '''
    function is wrapper for Trinity trimmomatic
    '''
    # create tmpdir
    folder = os.path.join(tmpdir, 'trimmomatic')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    lib.log.info("Adapter and Quality trimming PE reads with Trimmomatic")
    left_paired = os.path.join(folder, 'trimmed_left.fastq')
    left_single = os.path.join(folder, 'trimmed_left.unpaired.fastq')
    right_paired = os.path.join(folder, 'trimmed_right.fastq')
    right_single = os.path.join(folder, 'trimmed_right.unpaired.fastq')
    cmd = ['trimmomatic', 'PE', '-threads', str(cpus), '-phred33',
           left, right, left_paired, left_single, right_paired, right_single,
           'ILLUMINACLIP:' +
           os.path.join(parentdir, 'config', 'TruSeq3-PE.fa')+':2:30:10',
           'SLIDINGWINDOW:4:5', 'LEADING:5', 'TRAILING:5', 'MINLEN:25']
    lib.runSubprocess(cmd, '.', lib.log)
    for x in [left_paired, left_single, right_paired, right_single]:
        lib.Fzip_inplace(x, cpus)
    trim_left = os.path.join(folder, 'trimmed_left.fastq.gz')
    trim_right = os.path.join(folder, 'trimmed_right.fastq.gz')
    return trim_left, trim_right


def runTrimmomaticSE(reads, cpus=1):
    '''
    function is wrapper for Trinity trimmomatic
    '''
    # create tmpdir
    folder = os.path.join(tmpdir, 'trimmomatic')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    lib.log.info("Adapter and Quality trimming SE reads with Trimmomatic")
    output = os.path.join(folder, 'trimmed_single.fastq')
    cmd = ['trimmomatic', 'SE', '-threads', str(cpus), '-phred33',
           reads, output, 'ILLUMINACLIP:' +
           os.path.join(parentdir, 'config', 'TruSeq3-SE.fa')+':2:30:10',
           'SLIDINGWINDOW:4:5', 'LEADING:5', 'TRAILING:5', 'MINLEN:25']
    lib.runSubprocess(cmd, '.', lib.log)
    lib.Fzip_inplace(output, cpus)
    trim_single = os.path.join(folder, 'trimmed_single.fastq.gz')
    return trim_single


def runNormalization(readTuple, memory, min_coverage=5, coverage=50, cpus=1, stranded='no'):
    '''
    function is wrapper for Trinity read normalization
    have to run normalization separately for PE versus single
    '''
    left_norm, right_norm, single_norm = (None,)*3
    SENormalLog = os.path.join(tmpdir, 'trinity_normalization.SE.log')
    PENormalLog = os.path.join(tmpdir, 'trinity_normalization.PE.log')
    lib.log.info("Running read normalization with Trinity")
    if stranded != 'no':
        cmd = [os.path.join(TRINITY, 'util', 'insilico_read_normalization.pl'), '--PARALLEL_STATS',
               '--JM', memory, '--min_cov', str(
                   min_coverage), '--max_cov', str(coverage),
               '--seqType', 'fq', '--output', os.path.join(
                   tmpdir, 'normalize'), '--CPU', str(cpus),
               '--SS_lib_type', stranded]
    else:
        cmd = [os.path.join(TRINITY, 'util', 'insilico_read_normalization.pl'), '--PARALLEL_STATS',
               '--JM', memory, '--min_cov', str(
                   min_coverage), '--max_cov', str(coverage),
               '--seqType', 'fq', '--output', os.path.join(tmpdir, 'normalize'), '--CPU', str(cpus)]
    if readTuple[2]:  # single reads present, so run normalization just on those reads
        cmd = cmd + ['--single', readTuple[2]]
        lib.runSubprocess2(cmd, '.', lib.log, SENormalLog)
        single_norm = os.path.join(tmpdir, 'normalize', 'single.norm.fq')
    if readTuple[0] and readTuple[1]:
        cmd = cmd + ['--pairs_together', '--left',
                     readTuple[0], '--right', readTuple[1]]
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


def getPASAinformation(configFile, DBname, folder, genome):
    '''
    function to dump GFF from existing PASA database, compare genome headers to what is in PASA
    DB to make sure genome is same, return True if good to go, else error out
    '''
    # run some checks of the data to make sure it is same assembly
    mysqlDB, mysqlUser, mysqlPass = (None,)*3
    pasaconf_file = os.path.join(PASA, 'pasa_conf', 'conf.txt')
    if os.environ.get('PASACONF'):
        pasaconf_file = os.environ.get('PASACONF').strip()
    with open(pasaconf_file, 'r') as pasaconf:
        for line in pasaconf:
            line = line.replace('\n', '')
            if line.startswith('MYSQLSERVER='):
                mysqlDB = line.split('=')[-1]
            if line.startswith('MYSQL_RW_USER='):
                mysqlUser = line.split('=')[-1]
            if line.startswith('MYSQL_RW_PASSWORD='):
                mysqlPass = line.split('=')[-1]
    pasaExistingGFF = os.path.join(folder, 'existing_pasa.gff3')
    cmd = [os.path.join(PASA, 'scripts', 'pasa_asmbl_genes_to_GFF3.dbi'),
           '-M', DBname+':'+mysqlDB, '-p', mysqlUser+':'+mysqlPass]
    lib.runSubprocess2(cmd, folder, lib.log, pasaExistingGFF)
    if not lib.checkannotations(pasaExistingGFF):
        return False
    # now get number of genes and list of contigs
    pasaContigs = []
    geneCount = 0
    with open(pasaExistingGFF, 'r') as infile:
        for line in infile:
            if line.startswith('\n'):
                continue
            cols = line.split('\t')
            if not cols[0] in pasaContigs:
                pasaContigs.append(cols[0])
            if cols[2] == 'gene':
                geneCount += 1
    # now get fasta headers from genome
    genomeContigs = []
    with open(genome, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                line = line.replace('\n', '')
                line = line.replace('>', '')
                if not line in genomeContigs:
                    genomeContigs.append(line)
    # now make sure PASA headers in genome
    genomeContigs = set(genomeContigs)
    for contig in pasaContigs:
        if not contig in genomeContigs:
            return False
    lib.log.info(
        "Existing PASA database contains {:,} gene models, validated FASTA headers match".format(geneCount))
    return True


def runPASA(genome, transcripts, cleanTranscripts, gff3_alignments,
            stringtie_gtf, stranded, intronlen, cpus, previousGFF, dbname,
            output, configFile, pasa_db='sqlite',
            pasa_alignment_overlap=30, aligners=['blat'], min_pct_aligned=90,
            min_avg_id=95, num_bp_perfect=3):
    '''
    function will run PASA align assembly, followed by 2 rounds of comparison to update
    annotations for preexisting gene models
    '''
    pasa_cpus = int(cpus)
    # create tmpdir
    folder = os.path.join(tmpdir, 'pasa')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    pasaLOG = os.path.join(folder, 'pasa-assembly.log')
    pasaLOG1 = os.path.join(folder, 'pasa-comparison1.log')
    pasaLOG2 = os.path.join(folder, 'pasa-comparison2.log')

    # get config files and edit
    alignConfig = os.path.join(folder, 'alignAssembly.txt')
    annotConfig = os.path.join(folder, 'annotCompare.txt')

    # check if config file is passed, if so, get databasename and copy to assembly config file
    # dashes will get stripped in MySQL
    DataBaseName = dbname.replace('-', '_')
    DataBaseName += '_pasa'
    if pasa_db == 'sqlite':
        DataBaseName = os.path.abspath(os.path.join(folder, DataBaseName))
    if configFile:
        with open(configFile, 'r') as infile:
            for line in infile:
                line = line.replace('\n', '')
                if line.startswith('DATABASE=') or line.startswith('MYSQLDB='):
                    DataBaseName = line.split('=')[-1]
        shutil.copyfile(configFile, alignConfig)
        if pasa_db == 'mysql':
            # check existing database
            if not getPASAinformation(configFile, DataBaseName, folder, genome):
                lib.log.error(
                    "MySQL database not found or headers in PASA database, do not match those in FASTA.")
                # now run PASA alignment step
                lib.log.info("Running PASA alignment step using " +
                             "{0:,}".format(lib.countfasta(cleanTranscripts))+" transcripts")
                cmd = [LAUNCHPASA, '-c', os.path.abspath(alignConfig), '-r', '-C', '-R', '-g', os.path.abspath(genome),
                       '--IMPORT_CUSTOM_ALIGNMENTS', gff3_alignments, '-T',
                       '-t', os.path.abspath(
                           cleanTranscripts), '-u', os.path.abspath(transcripts),
                       '--stringent_alignment_overlap', pasa_alignment_overlap, '--TRANSDECODER',
                       '--MAX_INTRON_LENGTH', str(intronlen), '--CPU', str(pasa_cpus)]
                if 'minimap2' in aligners:
                    aligners.remove('minimap2')
                if aligners:
                    cmd.append('--ALIGNERS')
                    cmd.append(','.join(aligners))
                if stranded != 'no':
                    cmd = cmd + ['--transcribed_is_aligned_orient']
                if lib.checkannotations(stringtie_gtf):
                    cmd = cmd + ['--trans_gtf', stringtie_gtf]
                lib.runSubprocess6(cmd, folder, lib.log, pasaLOG)
        else:
            lib.log.info('PASA database is SQLite: {:}'.format(DataBaseName))
        # finally need to index the genome using cdbfasta so lookups can be done
        CDBFASTA = lib.which_path('cdbfasta')
        if not CDBFASTA:
            CDBFASTA = os.path.join(PASA, 'bin', 'cdbfasta')
        cmd = [CDBFASTA, genome]
        lib.runSubprocess(cmd, '.', lib.log)
    else:
        # create new config file from template
        with open(alignConfig, 'w') as config1:
            with open(os.path.join(PASA, 'pasa_conf', 'pasa.alignAssembly.Template.txt'), 'r') as template1:
                for line in template1:
                    if '<__DATABASE__>' in line:
                        line = line.replace('<__DATABASE__>', DataBaseName)
                    elif '<__MYSQLDB__>' in line:
                        line = line.replace('<__MYSQLDB__>', DataBaseName)
                    elif line.startswith('#script validate_alignments_in_db.dbi'):
                        line = line + '\n' + 'validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY={}\n'.format(num_bp_perfect)
                    elif '<__MIN_PERCENT_ALIGNED__>' in line:
                        line = line.replace('<__MIN_PERCENT_ALIGNED__>', str(min_pct_aligned))
                    elif '<__MIN_AVG_PER_ID__>' in line:
                        line = line.replace('<__MIN_AVG_PER_ID__>', str(min_avg_id))
                    config1.write(line)
        # align transcripts using minimap2
        # now run PASA alignment step
        lib.log.info("Running PASA alignment step using " +
                     "{0:,}".format(lib.countfasta(cleanTranscripts))+" transcripts")
        cmd = [LAUNCHPASA, '-c', os.path.abspath(alignConfig), '-r', '-C',
               '-R', '-g', os.path.abspath(genome),
               '--IMPORT_CUSTOM_ALIGNMENTS', gff3_alignments, '-T',
               '-t', os.path.abspath(cleanTranscripts),
               '-u', os.path.abspath(transcripts),
               '--stringent_alignment_overlap', pasa_alignment_overlap,
               '--TRANSDECODER',
               '--MAX_INTRON_LENGTH', str(intronlen), '--CPU', str(pasa_cpus)]
        cmd += ['--ALIGNERS']
        filtaligners = []
        for x in aligners:
            if x != 'minimap2':
                filtaligners.append(x)
        cmd.append(','.join(filtaligners))
        if stranded != 'no':
            cmd = cmd + ['--transcribed_is_aligned_orient']
        if lib.checkannotations(stringtie_gtf):
            cmd = cmd + ['--trans_gtf', stringtie_gtf]
        lib.runSubprocess6(cmd, folder, lib.log, pasaLOG)

    # generate comparison template file
    with open(annotConfig, 'w') as config2:
        with open(os.path.join(PASA, 'pasa_conf', 'pasa.annotationCompare.Template.txt'), 'r') as template2:
            for line in template2:
                line = line.replace('<__MYSQLDB__>', DataBaseName)
                line = line.replace('<__DATABASE__>', DataBaseName)
                config2.write(line)

    # now run Annotation comparisons
    lib.log.info("Running PASA annotation comparison step 1")
    cmd = [LAUNCHPASA, '-c', os.path.abspath(annotConfig),
           '-g', os.path.abspath(genome),
           '-t', os.path.abspath(cleanTranscripts),
           '-A', '-L', '--CPU', str(pasa_cpus)]
    if lib.versionCheck(PASAVERSION, '2.3.0'):
        cmd = cmd + ['--annots', os.path.abspath(previousGFF)]
    else:
        cmd = cmd + ['--annots_gff3', os.path.abspath(previousGFF)]
    lib.runSubprocess6(cmd, folder, lib.log, pasaLOG1)
    round1GFF = None
    for file in os.listdir(folder):
        if not file.endswith('.gff3'):
            continue
        if 'gene_structures_post_PASA_updates' in file:
            round1GFF = os.path.join(folder, file)
    if not round1GFF:
        lib.log.error("PASA failed, check log, exiting")
        sys.exit(1)
    # run round 2 comparison
    lib.log.info("Running PASA annotation comparison step 2")
    cmd = [LAUNCHPASA, '-c', os.path.abspath(annotConfig),
           '-g', os.path.abspath(genome),
           '-t', os.path.abspath(cleanTranscripts),
           '-A', '-L', '--CPU', str(pasa_cpus)]
    if lib.versionCheck(PASAVERSION, '2.3.0'):
        cmd = cmd + ['--annots', os.path.abspath(round1GFF)]
    else:
        cmd = cmd + ['--annots_gff3', os.path.abspath(round1GFF)]
    lib.runSubprocess6(cmd, folder, lib.log, pasaLOG2)
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
    # grab final result
    shutil.copyfile(round2GFF, output)


def pasa_transcript2gene(input):
    # modify kallisto ouput to map gene names to each mRNA ID so you know what locus they have come from
    mRNADict = {}
    # since mRNA is unique, parse the transcript file which has mRNAID geneID in header
    with open(input, 'r') as transin:
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
                with open(file, 'r') as infile:
                    for title, seq in SimpleFastaParser(infile):
                        if '/' in title:
                            title = title.replace('/', '_')
                        outfile.write('>{:}\n{:}\n'.format(title, lib.softwrap(seq)))
            elif file.endswith('.fq') or file.endswith('.fastq'):
                with open(file, 'r') as infile:
                    for title, seq, qual in FastqGeneralIterator(infile):
                        if '/' in title:
                            title = title.replace('/', '_')
                        outfile.write('>{:}\n{:}\n'.format(title, lib.softwrap(seq)))
        # clean up
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
        lib.log.info(
            'Processing long reads: converting to fasta and running SeqClean')
        results = []
        if readTuple[0] and not lib.checkannotations(PBiso+'.clean'):
            _convert2fasta(readTuple[0], PBiso)
            runSeqClean(PBiso, tmpdir, cpus=cpus)
        if readTuple[1] and not lib.checkannotations(nanocdna+'.clean'):
            _convert2fasta(readTuple[1], nanocdna)
            runSeqClean(nanocdna, tmpdir, cpus=cpus)
        if readTuple[2] and not lib.checkannotations(nanomrna+'.clean'):
            _convert2fasta(readTuple[2], nanomrna)
            runSeqClean(nanomrna, tmpdir, cpus=cpus)
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
            lib.catFiles(*validResults, output=combinedClean)
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


def runSeqClean(input, folder, cpus=1):
    '''
    wrapper to run PASA seqclean on Trinity transcripts
    '''
    if cpus > 16:
        cpus = 16

    if os.path.isfile(input + ".clean"):
        lib.log.info('Existing SeqClean output found: {:}'.format(
            os.path.join(folder, input + ".clean")))
    else:
        cmd = [os.path.join(PASA, 'bin', 'seqclean'),
               os.path.basename(input), '-c', str(cpus)]
        lib.runSubprocess(cmd, folder, lib.log)
    for f in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, f)):
            if f.startswith('cleaning'):
                lib.SafeRemove(os.path.join(folder, f))


def bam2fasta(input, output, cpus=1):
    tmpout = output+'.tmp'
    cmd = ['samtools', 'fasta', '-@', str(cpus), '-F', '0x4', input]
    lib.runSubprocess2(cmd, '.', lib.log, tmpout)
    # make sure no empty sequences
    with open(output, 'w') as outfile:
        with open(tmpout, 'r') as infile:
            for header, seq in SimpleFastaParser(infile):
                if len(seq) > 0:
                    outfile.write('>{:}\n{:}\n'.format(header, lib.softwrap(seq)))
    os.remove(tmpout)


def bam2fasta_unmapped(input, output, cpus=1):
    tmpout = output+'.tmp'
    cmd = ['samtools', 'fasta', '-@', str(cpus), '-f', '0x4', input]
    lib.runSubprocess2(cmd, '.', lib.log, tmpout)
    # make sure no empty sequences
    with open(output, 'w') as outfile:
        with open(tmpout, 'r') as infile:
            for header, seq in SimpleFastaParser(infile):
                if len(seq) > 0:
                    outfile.write('>{:}\n{:}\n'.format(header, lib.softwrap(seq)))
    os.remove(tmpout)


def longReadFilter(fastx, reference, output, cpus=8, method='map-ont',
                   min_pident=80, secondary=False, options=[]):
    cmd = ['minimap2', '-x', method, '-c', '-t', str(cpus), '--paf-no-hit']
    if not secondary:
        cmd.append('--secondary=no')
    if options:
        cmd += options
    cmd += [reference, fastx]
    keep = set()
    for line in lib.execute(cmd):
        line = line.rstrip()
        query, qlen, qstart, qend, strand, target, tlen, tstart, tend, matches, alnlen, qual = line.split('\t')[:12]
        if query in keep:
            continue
        if target == '*':
            keep.add(query)
        extras = line.split('\t')[12:]
        gap_exl_pident = None
        for x in extras:
            if x.startswith('de:f:'):
                gap_exl_pident = 100 - (float(x.replace('de:f:', ''))*100)
        if gap_exl_pident and gap_exl_pident < min_pident:
            keep.add(query)
        elif not gap_exl_pident:
            keep.add(query)
    with open(output, 'w') as outfile:
        with open(fastx, 'r') as infile:
            for title, seq in SimpleFastaParser(infile):
                if title.split()[0] in keep:
                    outfile.write('>{}\n{}\n'.format(title, lib.softwrap(seq)))
    return len(keep)


def longReadMap(fastx, reference, output, maxintronlen=3000,
                cpus=8, method='splice', minqual=2, options=[]):
    cmd = ['minimap2', '-x', method, '-t', str(cpus), '-G', str(maxintronlen)]
    if options:
        cmd += options
    cmd += [reference, fastx]
    keep = set()
    for line in lib.execute(cmd):
        line = line.rstrip()
        query, qlen, qstart, qend, strand, target, tlen, tstart, tend, matches, alnlen, qual = line.split('\t')[:12]
        if query in keep:
            continue
        if int(qual) > minqual:
            keep.add(query)
    with open(output, 'w') as outfile:
        with open(fastx, 'r') as infile:
            for title, seq in SimpleFastaParser(infile):
                if title.split()[0] in keep:
                    outfile.write('>{}\n{}\n'.format(title, lib.softwrap(seq)))
    return len(keep)


def mapTranscripts(genome, longTuple, assembled, tmpdir, trinityBAM, allBAM, cpus=1, max_intronlen=3000):
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
    # tuple is (iso-seq, nanopore_cDNA, nanopore_mRNA)
    if not all(v is None for v in longTuple):
        # run minimap2 alignment
        lib.log.info('Aligning long reads to genome with minimap2')
        if longTuple[0]:  # run iso-seq method
            isoMap = os.path.join(tmpdir, 'isoseq.mapped.fasta')
            mapped = longReadMap(longTuple[0], genome, isoMap, cpus=cpus,
                                 maxintronlen=max_intronlen,
                                 options=['-uf', '-C5'])
            lib.log.debug('{:,} IsoSeq reads mapped to genome'.format(mapped))
            if lib.checkannotations(assembled):
                unmapped_count = longReadFilter(isoMap, assembled, isoSeqs,
                                                cpus=cpus)
                lib.log.debug('{:,} IsoSeq reads not found in Trinity'.format(unmapped_count))
            else:
                isoSeqs = isoMap
            lib.iso_seq_minimap2(isoSeqs, genome, cpus, max_intronlen, isoBAM)
        if longTuple[1]:  # run nano cDNA
            nano_cDNA_Map = os.path.join(tmpdir, 'nano_cDNA.mapped.fasta')
            mapped = longReadMap(longTuple[1], genome, nano_cDNA_Map, cpus=cpus,
                                 maxintronlen=max_intronlen)
            lib.log.debug('{:,} ONT cDNA reads mapped to genome'.format(mapped))
            if lib.checkannotations(assembled):
                unmapped_count = longReadFilter(nano_cDNA_Map, assembled,
                                                nano_cdnaSeqs, cpus=cpus)
                lib.log.debug('{:,} ONT cDNA reads not found in Trinity'.format(unmapped_count))
            else:
                nano_cdnaSeqs = nano_cDNA_Map
            lib.nanopore_cDNA_minimap2(nano_cdnaSeqs, genome, cpus,
                                       max_intronlen, nano_cdnaBAM)
        if longTuple[2]:  # run nano mRNA
            nano_mrna_Map = os.path.join(tmpdir, 'nano_mRNA.mapped.fasta')
            mapped = longReadMap(longTuple[2], genome, nano_mrna_Map, cpus=cpus,
                                 maxintronlen=max_intronlen,
                                 options=['-uf', '-k14'])
            lib.log.debug('{:,} ONT mRNA reads mapped to genome'.format(mapped))
            if lib.checkannotations(assembled):
                unmapped_count = longReadFilter(nano_mrna_Map, assembled,
                                                nano_mrnaSeqs, cpus=cpus)
                lib.log.debug('{:,} ONT mRNA reads not found in Trinity'.format(unmapped_count))
            else:
                nano_mrnaSeqs = nano_mrna_Map
            lib.nanopore_mRNA_minimap2(nano_mrnaSeqs, genome, cpus,
                                       max_intronlen, nano_mrnaBAM)
        for x in [isoSeqs, nano_cdnaSeqs, nano_mrnaSeqs]:
            if lib.checkannotations(x):
                mappedSeqs.append(x)
        if len(mappedSeqs) > 0:
            lib.catFiles(*mappedSeqs, output=mappedLong)
            lib.SafeRemove(isoSeqs)
            lib.SafeRemove(nano_cdnaSeqs)
            lib.SafeRemove(nano_mrnaSeqs)
            lib.log.info('Adding {:,} unique long-reads to Trinity assemblies'.format(
                lib.countfasta(mappedLong)))

    if lib.checkannotations(assembled):  # Trinity transcripts
        if lib.checkannotations(mappedLong):
            trinityCombined = os.path.join(tmpdir, 'trinity.long-reads.fasta')
            trinityCombinedClean = trinityCombined+'.clean'
            lib.catFiles(*[assembled, mappedLong], output=trinityCombined)
            runSeqClean(trinityCombined, tmpdir, cpus=cpus)
        else:
            trinityCombinedClean = assembled
            trinityCombined = assembled.replace('.clean', '')
        # finally run trinity mapping
        lib.minimap2Align(trinityCombinedClean, genome,
                          cpus, max_intronlen, trinityBAM)
    else:
        trinityCombined = mappedLong
        trinityCombinedClean = trinityCombined+'.clean'
        runSeqClean(trinityCombined, tmpdir, cpus=cpus)

    bamResults = [isoBAM, nano_cdnaBAM, nano_mrnaBAM, trinityBAM]
    foundResults = []
    for r in bamResults:
        if lib.checkannotations(r):
            foundResults.append(r)
    if len(foundResults) > 1:
        lib.log.info('Merging BAM files: {:}'.format(', '.join(foundResults)))
        lib.mergeBAMs(*foundResults, cpus=cpus, output=allBAM)
    elif len(foundResults) == 0:
        lib.log.error(
            'Alignment failed, BAM files empty. Please check logfile')
        sys.exit(1)
    else:
        os.symlink(os.path.abspath(foundResults[0]), os.path.abspath(allBAM))
    return trinityCombined, trinityCombinedClean


def runKallisto(input, fasta, readTuple, stranded, cpus, output):
    '''
    function takes GFF3 output from PASA compare, extracts transcripts, and then calculates TPM
    using Kallisto to idenitfy the best scoring gene model for each locus, the left and right
    these should be the adapter cleaned non-normalized Illumina reads
    '''
    lib.log.info(
        "Using Kallisto TPM data to determine which PASA gene models to select at each locus")
    # convert GFF to transcripts
    folder = os.path.join(tmpdir, 'getBestModel')
    os.makedirs(folder)
    PASAtranscripts = os.path.join(folder, 'transcripts.fa')
    cmd = [os.path.join(PASA, 'misc_utilities',
                        'gff3_file_to_proteins.pl'), input, fasta, 'cDNA']
    lib.log.info("Building Kallisto index")
    lib.runSubprocess2(cmd, '.', lib.log, PASAtranscripts)
    # generate kallisto index
    cmd = ['kallisto', 'index', '-i',
           os.path.join(folder, 'bestModel'), PASAtranscripts]
    lib.runSubprocess(cmd, '.', lib.log)
    # use kallisto to map reads to index
    # base command
    cmd = ['kallisto', 'quant', '-i', os.path.join(folder, 'bestModel'), '-o', os.path.join(
        folder, 'kallisto'), '--plaintext', '-t', str(cpus)]
    # parse the strand information
    if stranded == 'RF':
        strandcmd = ['--rf-stranded']
    elif stranded == 'FR':
        strandcmd = ['--fr-stranded']
    else:
        strandcmd = []
    # adapt command for input, i.e. single or PE ends -> what do you do if you have both?
    # single, not just using estimated lengths and SD, I think this is okay? can make this an option otherwise
    if readTuple[2] and not readTuple[0] and not readTuple[1]:
        cmd = cmd + ['--single', '-l', '200', '-s', '20', readTuple[2]]
    elif readTuple[0] and readTuple[1]:
        cmd = cmd + strandcmd + [readTuple[0], readTuple[1]]
    lib.log.info("Mapping reads using pseudoalignment in Kallisto")
    lib.runSubprocess(cmd, '.', lib.log)

    # modify kallisto ouput to map gene names to each mRNA ID so you know what locus they have come from
    mRNADict = pasa_transcript2gene(PASAtranscripts)

    # some PASA models can have incomplete CDS and are wrong, get list of incompletes to ignore list
    ignore = []
    with open(input, 'r') as infile:
        for line in infile:
            if line.startswith('#PROT'):
                if line.endswith('\t\n'):
                    ID = line.split(' ')[1]
                    ignore.append(ID)
    if len(ignore) > 0:
        lib.log.debug("Ignoring %i incomplete PASA models: %s" %
                      (len(ignore), ','.join(ignore)))

    # now make new tsv file with #mRNAID geneID location TPM
    with open(output, 'w') as outfile:
        outfile.write("#mRNA-ID\tgene-ID\tLocation\tTPM\n")
        with open(os.path.join(folder, 'kallisto', 'abundance.tsv'), 'r') as infile:
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
                    outfile.write('%s\t%s\t%s\t%s\n' %
                                  (cols[0], geneID, location, cols[4]))


def getBestModels(input, fasta, abundances, alt_transcripts, outfile):
    # function to parse PASA results and generate GFF3; supports multiple transcripts
    if float(alt_transcripts) == 0:
        lib.log.info(
            "Parsing Kallisto results. Keeping all alt-splicing transcripts at each locus.")
    elif float(alt_transcripts) < 1:
        lib.log.info(
            "Parsing Kallisto results. Keeping alt-splicing transcripts if expressed at least {0:.1f}% of highest transcript per locus.".format(float(alt_transcripts)*100))
    else:
        lib.log.info(
            "Parsing Kallisto results. Keeping best transcript at each locus.")
    bestModels = {}
    locations = {}
    with open(abundances, 'r') as tpms:
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
                bestModels[geneLocus].append((transcriptID, float(TPM)))
    # now we have geneID dictionary containing list of tuples of of transcripts
    # loop through each locus grabbing transcript IDs of those models to keep
    # use best expression value * alt_transcript threshold to filter models
    # alt_transcript == 1 would be then only keeping best hit
    extractList = []
    ExpValues = {}
    for k, v in natsorted(list(bestModels.items())):
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
    # now go through the PASA GFF file and generate filtered GFF3 file composed of extractList
    extractList = set(extractList)
    with open(outfile, 'w') as output:
        output.write("##gff-version 3\n")
        with open(input, 'r') as gff:
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
                        geneID = 'ID=' + geneID.split(';Parent=')[-1]
                        mRNAID = cols[8].split(';Name=')[0]
                        expression = ExpValues.get(gffID)
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};\n'.format(
                            cols[0], 'PASA', 'gene', cols[3], cols[4], cols[5], cols[6], cols[7], geneID))
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};Note=TPM:{:0.2f};\n'.format(
                            cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], mRNAID, expression))
                elif '_prime_UTR' in cols[2]:
                    utrID = gffID.split('.utr')[0]
                    if utrID in extractList:
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};\n'.format(
                            cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
                elif 'exon' in cols[2]:
                    exonID = gffID.split('.exon')[0]
                    if exonID in extractList:
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};\n'.format(
                            cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
                elif 'CDS' in cols[2]:
                    cdsID = gffID.split('cds.')[-1]
                    if cdsID in extractList:
                        output.write('{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:};\n'.format(
                            cols[0], 'PASA', cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]))
    lib.log.info('Wrote {:,} transcripts derived from {:,} protein coding loci.'.format(
        len(extractList), len(bestModels)))


def GFF2tblCombinedNEW(evm, genome, trnascan, prefix, genenumber, justify, SeqCenter, SeqRefNum, tblout, alt_transcripts='1'):
    from collections import OrderedDict
    '''
    function to take GFF3 annotation to produce a GBK tbl file, support multiple transcripts per locus.
    '''
    def _sortDict(d):
        return (d[1]['contig'], d[1]['location'][0])
    # make sure genenumber is integer
    genenumber = int(genenumber)
    # generate genome length dictionary used for feature tbl generation
    scaffLen = {}
    with open(genome, 'r') as fastain:
        for record in SeqIO.parse(fastain, 'fasta'):
            if not record.id in scaffLen:
                scaffLen[record.id] = len(record.seq)
    # setup interlap database for genes on each chromosome and load EVM models into dictionary
    gene_inter = defaultdict(InterLap)
    Genes = {}
    gene_inter, Genes = lib.gff2interlapDict(
        evm, genome, gene_inter, Genes)
    # now load tRNA predictions
    gene_inter, Genes = lib.gff2interlapDict(
        trnascan, genome, gene_inter, Genes)
    # now sort dictionary by contig and location, rename using prefix
    sGenes = sorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    scaff2genes = {}
    SeqRecords = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
    skipList = []
    dropped = 0
    keeper = 0
    tooShort = 0
    internalStop = 0
    lib.log.info(
        "Validating gene models (renaming, checking translations, filtering, etc)")
    for k, v in list(sortedGenes.items()):
        GoodModel = True
        # check if gene model completely contained inside another one on same strand
        if alt_transcripts == '1':
            loc = sorted([v['location'][0], v['location'][1]])
            if loc in gene_inter[v['contig']]:
                for hit in list(gene_inter[v['contig']].find(loc)):
                    if hit[3] != k and hit[2] == v['strand']:  # same strand but diff gene
                        sortedhit = sorted([hit[0], hit[1]])
                        # then this gene is fully contained, skip it
                        if loc[0] >= sortedhit[0] and loc[1] <= sortedhit[1]:
                            # if two gene models have exact same start stop they will both be removed, not really what I want, so run check
                            # exact same, then choose which has higher TPM
                            if loc[0] == sortedhit[0] and loc[1] == sortedhit[1]:
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
        # rename gene locus here
        keeper += 1
        # renaming scheme, leave if startswith locustag
        if k.startswith('novel_gene') or k.startswith('temp_gene') or k.startswith('split_gene'):
            genenumber += 1
            locusTag = prefix+str(genenumber).zfill(justify)
        elif k.startswith(prefix) and '_'+prefix in k:  # means models were merged
            locusTag = k.split('_'+prefix)[0]
        else:
            locusTag = k
        # translate to protein space, drop if less than minimum
        # translate to protein sequence, construct cDNA Seq object, translate
        # if passes then add to final output dictionary
        for i in range(0, len(v['ids'])):
            protSeq = None
            if v['type'] == 'mRNA':  # get transcript for valid models
                cdsSeq = lib.getSeqRegions(
                    SeqRecords, v['contig'], v['CDS'][i])
                protSeq = lib.translate(
                    cdsSeq, v['strand'], v['codon_start'][i]-1)
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
                else:  # try to extend the CDS up to 20 codons to see if you can find valid stop codon
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
    lib.log.info('Writing {:,} loci to TBL format: dropped {:,} overlapping, {:,} too short, and {:,} frameshift gene models'.format(
        len(renamedGenes), dropped, tooShort, internalStop))
    lib.dicts2tbl(renamedGenes, scaff2genes, scaffLen,
                  SeqCenter, SeqRefNum, skipList, tblout)


def gbk2interlap(input):
    '''
    function to parse GBK file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    '''
    inter = defaultdict(InterLap)
    Genes = {}
    with open(input, 'r') as filein:
        for record in SeqIO.parse(filein, 'genbank'):
            for f in record.features:
                if f.type == 'gene':
                    locusTag, ID, Parent = lib.getID(f, f.type)
                    start = int(f.location.nofuzzy_start)
                    end = int(f.location.nofuzzy_end)
                    inter[record.id].add((start, end, locusTag))
                lib.gb_feature_add2dict(f, record, Genes)
    return inter, Genes


def gff2interlap(input, fasta):
    '''
    function to parse GFF3 file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    '''
    inter = defaultdict(InterLap)
    Genes = {}
    Genes = lib.gff2dict(input, fasta, Genes)
    for k, v in natsorted(list(Genes.items())):
        inter[v['contig']].add((v['location'][0], v['location'][1], k))
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
    # structured message, coordinates, 5prime, 3prime, exon, cds, pident
    if not pos:  # coordinates changed
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
        msg.append('CDS update [translation pident: %s]' %
                   ', '.join(pidentmsg))
    # now work on the counter
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


def compareAnnotations2(old, new, output, args={}):
    '''
    function takes two GenBank annotated genomes and compares gene models
    output is a tsv file for each locus and a description of what is different
    can handle multiple transcripts per locus
    '''
    result = {}
    global no_change, UTR_added, yardSale, exonChange, modelChangeNotProt, dropped, added, total_transcripts, total_genes
    no_change, UTR_added, yardSale, exonChange, modelChangeNotProt, dropped, added, total_transcripts, total_genes = (
        0,)*9
    lib.log.info("Parsing GenBank files...comparing annotation")
    if args.gff and args.fasta:
        oldInter, oldGenes = gff2interlap(old, args.fasta)
    else:
        oldInter, oldGenes = gbk2interlap(old)
    newInter, newGenes = gbk2interlap(new)
    # do the simple stuff first, find models that were deleted
    for contig in oldInter:
        for gene in oldInter[contig]:
            if not gene in newInter[contig]:  # these models are removed
                dropped += 1
                if not gene[2] in oldGenes:
                    continue
                # populate output dictionary with results
                if not gene[2] in result:
                    # dropped model has AED of 1.000
                    cdsAED = '1.000'
                    exonAED = '1.000'
                    result[gene[2]] = {'contig': oldGenes[gene[2]]['contig'], 'old_num_transcripts': len(oldGenes[gene[2]]['ids']),
                                       'old_location': oldGenes[gene[2]]['location'], 'num_transcripts': len(oldGenes[gene[2]]['ids']), 'strand': oldGenes[gene[2]]['strand'],
                                       'mRNA': oldGenes[gene[2]]['mRNA'], 'location': oldGenes[gene[2]]['location'], 'CDS': oldGenes[gene[2]]['CDS'], 'message': 'gene model removed',
                                       'cdsAED': cdsAED, 'exonAED': exonAED, 'transcript_id': oldGenes[gene[2]]['ids'], 'pident': [],
                                       'protein_id': oldGenes[gene[2]]['ids'], 'seq': oldGenes[gene[2]]['protein']}

    # now go through the updated annotation, comparing to old annot
    for contig in newInter:
        for gene in newInter[contig]:
            # means this is a new model, so add it
            if not gene in oldInter[contig]:
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
                                       'protein_id': newGenes[gene[2]]['ids'], 'seq': newGenes[gene[2]]['protein'], 'pident': []}
            else:  # means this is existing model, and need to do some comparisons
                hitList = list(oldInter[contig].find(gene))
                # there might be some overlapping transcripts, so enforce locus name
                hit = None
                for z in hitList:
                    if gene[2] == z[2]:
                        hit = z
                if not hit:
                    # there is no real hit, so this a new gene
                    total_transcripts += len(newGenes[gene[2]]['ids'])
                    added += 1
                    total_genes += 1
                    if not gene[2] in result:
                        result[gene[2]] = {'contig': newGenes[gene[2]]['contig'], 'old_num_transcripts': 0,
                                           'old_location': newGenes[gene[2]]['location'], 'num_transcripts': len(newGenes[gene[2]]['ids']), 'strand': newGenes[gene[2]]['strand'],
                                           'mRNA': newGenes[gene[2]]['mRNA'], 'location': newGenes[gene[2]]['location'], 'CDS': newGenes[gene[2]]['CDS'], 'message': 'new gene model',
                                           'cdsAED': '0.000', 'exonAED': '0.000', 'transcript_id': newGenes[gene[2]]['ids'],
                                           'protein_id': newGenes[gene[2]]['ids'], 'seq': newGenes[gene[2]]['protein'], 'pident': []}
                else:
                    # since we may have multiple transcripts from hit as well as new annotation we need to be aware of that
                    # also, tRNA annotations do not exist in Proteins dictionary, so process them differently
                    # get the reference hits, pull out CDS and mRNA for pairwiseAED calculation
                    total_genes += 1
                    total_transcripts += len(newGenes[gene[2]]['ids'])

                    # get the old annotation
                    hitInfo = oldGenes.get(gene[2])

                    # calculate AED
                    exonAED = pairwiseAED(
                        newGenes[gene[2]]['mRNA'], hitInfo['mRNA'])
                    if newGenes[gene[2]]['type'] == 'mRNA' and hitInfo['type'] == 'mRNA':
                        cdsAED = pairwiseAED(
                            newGenes[gene[2]]['CDS'], hitInfo['CDS'])
                    else:
                        cdsAED = '0.000'

                    # check translation, to deal with multiple transcripts, lets loop through new
                    protMatches = []
                    if newGenes[gene[2]]['type'] == 'mRNA' and hitInfo['type'] == 'mRNA':
                        for i in range(0, len(newGenes[gene[2]]['ids'])):
                            protMatch = None
                            for y in range(0, len(oldGenes[gene[2]]['ids'])):
                                pident = pairwiseAlign(
                                    newGenes[gene[2]]['protein'][i], oldGenes[gene[2]]['protein'][y])
                                if not protMatch:
                                    protMatch = pident
                                else:
                                    if pident > protMatch:
                                        protMatch = pident
                            protMatches.append(protMatch)
                    # summarize UTRs
                    UTRs = findUTRs(
                        newGenes[gene[2]]['CDS'], newGenes[gene[2]]['mRNA'], newGenes[gene[2]]['strand'])

                    # structured comments/counts for gene models
                    msg, no_change, UTR_added, yardSale, exonChange = message(
                        newGenes[gene[2]]['location'], oldGenes[gene[2]]['location'], cdsAED, exonAED, protMatches, UTRs, no_change, UTR_added, yardSale, exonChange)

                    if not gene[2] in result:
                        result[gene[2]] = {'contig': newGenes[gene[2]]['contig'], 'old_num_transcripts': len(oldGenes[gene[2]]['ids']),
                                           'old_location': oldGenes[gene[2]]['location'], 'num_transcripts': len(newGenes[gene[2]]['ids']), 'strand': newGenes[gene[2]]['strand'],
                                           'mRNA': newGenes[gene[2]]['mRNA'], 'location': newGenes[gene[2]]['location'], 'CDS': newGenes[gene[2]]['CDS'], 'message': msg,
                                           'cdsAED': cdsAED, 'exonAED': exonAED, 'transcript_id': newGenes[gene[2]]['ids'],
                                           'protein_id': newGenes[gene[2]]['ids'], 'seq': newGenes[gene[2]]['protein'], 'pident': protMatches}

    total_cdsAED = []
    total_exonAED = []
    with open(output, 'w') as out:
        out.write('Locus_tag\tOrig_Location\tOrig_Num_Transcripts\tContig:start-end\tStrand\tGene_Length\tNum_Transcripts\tmRNA_AED\tCDS_AED\tDescription\n')
        for k, v in natsorted(list(result.items())):
            start = str(v['location'][0])
            end = str(v['location'][1])
            GeneLength = int(end) - int(start)
            total_cdsAED.append(float(v['cdsAED']))
            total_exonAED.append(float(v['exonAED']))
            out.write('{:}\t{:}:{:}-{:}\t{:}\t{:}:{:}-{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format(k, v['contig'], v['old_location'][0], v['old_location'][
                      1], v['old_num_transcripts'], v['contig'], start, end, v['strand'], GeneLength, v['num_transcripts'], v['exonAED'], v['cdsAED'], v['message']))
    Avg_cdsAED = sum(total_cdsAED) / float(len(total_cdsAED))
    Avg_exonAED = sum(total_exonAED) / float(len(total_exonAED))
    # output some simple stats to cmd line
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
    # supporting multiple transcripts, however, they are already matched up and sorted
    UTRs = []
    for i in range(0, len(cds)):
        Fiveprime = False
        Threeprime = False
        refInterlap = InterLap(mrna[i])
        if strand == '+':  # look at first CDS for 5 prime and last CDS for 3 prime
            # means it overlaps with mrNA (which it obviously should)
            if cds[i][0] in refInterlap:
                hit = list(refInterlap.find(cds[i][0]))[0]
                # if first exon, then compare, if not first then there is 5prime UTR
                loc = mrna[i].index(hit)
                if loc == 0:
                    # will return array of exon minus hit at each pos
                    diff = np.subtract(cds[i][0], hit)
                    if diff[0] > 0:
                        Fiveprime = True
                else:
                    Fiveprime = True
            # check for 3 prime UTR
            if cds[i][-1] in refInterlap:
                hit = list(refInterlap.find(cds[i][-1]))[0]
                loc = mrna[i].index(hit)
                if len(mrna[i]) == loc+1:
                    # will return array of exon minus hit at each pos
                    diff = np.subtract(cds[i][-1], hit)
                    if diff[1] < 0:
                        Threeprime = True
                else:
                    Threeprime = True
        else:
            # means it overlaps with mrNA (which it obviously should)
            if cds[i][0] in refInterlap:
                hit = list(refInterlap.find(cds[i][0]))[0]
                # if first exon, then compare, if not first then there is 5prime UTR
                loc = mrna[i].index(hit)
                if loc == 0:
                    # will return array of exon minus hit at each pos
                    diff = np.subtract(cds[i][0], hit)
                    if diff[1] < 0:
                        Fiveprime = True
                else:
                    Fiveprime = True
            # check for 3 prime UTR
            if cds[i][-1] in refInterlap:
                hit = list(refInterlap.find(cds[i][-1]))[0]
                loc = mrna[i].index(hit)
                if len(mrna[i]) == loc+1:
                    # will return array of exon minus hit at each pos
                    diff = np.subtract(cds[i][-1], hit)
                    if diff[0] > 0:
                        Threeprime = True
                else:
                    Threeprime = True
        UTRs.append((Fiveprime, Threeprime))
    return UTRs


def pairwiseAED(query, reference):
    '''
    takes a multiple transcripts and sums AED from lowest pairwise comparison and then calculates
    the average based on number of transcripts in the query
    '''
    AEDsum = []
    pAED = [float(getAED(a, b))
            for a, b in itertools.product(query, reference)]
    # split into parts to get lowest AED
    splitAED = [pAED[i:i+len(query)] for i in range(0, len(pAED), len(query))]
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
    # check if identical
    if query == reference:
        return '0.000'
    # make sure sorted
    rLen = _length(reference)
    refInterlap = InterLap(reference)
    QueryOverlap = 0
    qLen = 0
    for exon in query:
        qLen += abs(exon[0] - exon[1])
        if exon in refInterlap:  # exon overlaps at least partially with reference
            hit = list(refInterlap.find(exon))
            for h in hit:
                # will return array of exon minus hit at each pos
                diff = np.subtract(exon, h)
                if diff[0] <= 0 and diff[1] >= 0:  # then query exon covers ref exon
                    cov = abs(h[0] - h[1])
                    QueryOverlap += cov
                elif diff[0] <= 0 and diff[1] < 0:  # means query partial covers ref
                    cov = abs(h[0] - exon[1])
                    QueryOverlap += cov
                elif diff[0] > 0 and diff[1] >= 0:  # means query partial covers ref
                    cov = abs(exon[0] - h[1])
                    QueryOverlap += cov
                elif diff[0] > 0 and diff[1] < 1:
                    cov = abs(exon[0] - exon[1])
                    QueryOverlap += cov
    # calculate AED
    SP = QueryOverlap / float(qLen)
    SN = QueryOverlap / float(rLen)
    AED = 1 - ((SN + SP) / 2)
    return '{:.3f}'.format(AED)


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='funannotate-update.py', usage="%(prog)s [options] -i genome.gbk -l left.fq.gz -r right.fg.gz",
                                     description='''Script is a wrapper for automated Trinity/PASA reannotation.''',
                                     epilog="""Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--input',
                        help='Genome in GBK format or funannotate folder')
    parser.add_argument('-g', '--gff', help='Genome annotation in GFF3 format')
    parser.add_argument('-f', '--fasta',
                        help='Genome sequence in FASTA format')
    parser.add_argument('-l', '--left', nargs='+',
                        help='Left (R1) FASTQ Reads')
    parser.add_argument('--left_norm', help='Left (R1) FASTQ Reads')
    parser.add_argument('--right_norm', help='Right (R2) normalized FASTQ Reads')
    parser.add_argument('--single_norm', help='single normalized FASTQ Reads')
    parser.add_argument('-r', '--right', nargs='+',
                        help='Right (R2) FASTQ Reads')
    parser.add_argument('-s', '--single', nargs='+',
                        help='Single ended FASTQ Reads')
    parser.add_argument('--pacbio_isoseq', help='PacBio Iso-seq data')
    parser.add_argument('--nanopore_cdna', help='Nanopore 2d cDNA data')
    parser.add_argument('--nanopore_mrna', help='Nanopore direct mRNA data')
    parser.add_argument('-o', '--out', help='Basename of output files')
    parser.add_argument('--species',
                        help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
    parser.add_argument('-c', '--coverage', default=50,
                        type=int, help='Depth to normalize reads to')
    parser.add_argument('-m', '--min_coverage', default=5, type=int,
                        help='Minimum depth to pass to Trinity during normalization')
    parser.add_argument('--isolate', help='Isolate name (e.g. Af293)')
    parser.add_argument('--strain', help='Strain name (e.g. CEA10)')
    parser.add_argument('--trinity',
                        help='Trinity genome guided FASTA results')
    parser.add_argument('--pasa_gff', help='PASA GFF')
    parser.add_argument('--pasa_alignment_overlap', default='30.0',
                        help='PASA --stringent_alingment_overlap')
    parser.add_argument('--pasa_min_pct_aligned', default='90',
                        help='PASA --MIN_PERCENT_ALIGNED')
    parser.add_argument('--pasa_min_avg_per_id', default='95',
                        help='PASA --MIN_AVG_PER_ID')
    parser.add_argument('--pasa_num_bp_splice', default='3',
                        help='PASA --NUM_BP_PERFECT_SPLICE_BOUNDARY')
    parser.add_argument('--pasa_config',
                        help='PASA assembly configuration file')
    parser.add_argument('--pasa_db', default='sqlite',
                        choices=['mysql', 'sqlite'], help='PASA SQL database to use')
    parser.add_argument('--memory', default='50G',
                        help='RAM to use for Jellyfish/Trinity')
    parser.add_argument('--no_normalize_reads',
                        action='store_true', help='skip normalization')
    parser.add_argument('--no_trimmomatic', action='store_true',
                        help='skip quality trimming via trimmomatic')
    parser.add_argument('--jaccard_clip', action='store_true',
                        help='Turn on jaccard_clip for dense genomes')
    parser.add_argument('--kallisto', help='Kallisto abundances table')
    parser.add_argument('--name',
                        help='Shortname for genes, perhaps assigned by NCBI, eg. VC83_')
    parser.add_argument('--max_intronlen', default=3000,
                        help='Maximum intron length for gene models')
    parser.add_argument('--min_protlen', default=50, type=int,
                        help='Minimum amino acid length for valid gene model')
    parser.add_argument('--stranded', default='no',
                        choices=['RF', 'FR', 'F', 'R', 'no'], help='RNA seq strandedness')
    parser.add_argument('--cpus', default=2, type=int,
                        help='Number of CPUs to use')
    parser.add_argument('-t', '--tbl2asn', default='-l paired-ends',
                        help='Parameters for tbl2asn, linkage and gap info')
    parser.add_argument('--sbt', default='SBT',
                        help='Basename of output files')
    parser.add_argument('--p2g', help='NCBI p2g file from previous annotation')
    parser.add_argument('--aligners', default=['minimap2', 'blat'], nargs='+', choices=[
                        'minimap2', 'gmap', 'blat'], help='transcript alignment programs')
    parser.add_argument('--PASAHOME',
                        help='Path to PASA home directory, $PASAHOME')
    parser.add_argument('--TRINITYHOME',
                        help='Path to Trinity config directory, $TRINITYHOME')
    parser.add_argument('--SeqCenter', default='CFMR',
                        help='Sequencing center for GenBank tbl file')
    parser.add_argument('--SeqAccession', default='12345',
                        help='Sequencing accession number')
    parser.add_argument('--alt_transcripts', default='0.10',
                        help='Threshold to keep alt-transcripts, percent highest expression')
    parser.add_argument('--no-progress', dest='progress', action='store_false',
                        help='no progress on multiprocessing')
    args = parser.parse_args(args)

    global FNULL
    FNULL = open(os.devnull, 'w')

    # create folder structure
    if args.input:
        if os.path.isdir(args.input):  # then funannoate folder is passed
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
        # make sure subdirectories exist
        dirs = [os.path.join(args.out, 'update_misc'), os.path.join(
            args.out, 'logfiles'), os.path.join(args.out, 'update_results')]
        for d in dirs:
            if not os.path.isdir(d):
                os.makedirs(d)

    # assign temp directory
    global tmpdir, PASA, LAUNCHPASA, PASAVERSION, TRINITY,  PBiso, nanocdna, nanomrna, parentdir
    parentdir = os.path.join(os.path.dirname(__file__))
    tmpdir = os.path.join(args.out, 'update_misc')

    # create log file
    log_name = os.path.join(args.out, 'logfiles', 'funannotate-update.log')
    if os.path.isfile(log_name):
        os.remove(log_name)

    # initialize script, log system info and cmd issue at runtime
    lib.setupLogging(log_name)
    cmd_args = " ".join(sys.argv)+'\n'
    lib.log.debug(cmd_args)
    print("-------------------------------------------------------")
    lib.SystemInfo()

    # get version of funannotate
    version = lib.get_version()
    lib.log.info("Running %s" % version)

    # do some checks and balances
    if not args.PASAHOME:
        try:
            PASA = os.environ["PASAHOME"].strip()
        except KeyError:
            lib.log.error(
                "$PASAHOME environmental variable not found, PASA is not properly configured.  You can use the --PASAHOME argument to specifiy a path at runtime")
            sys.exit(1)
    else:
        PASA = args.PASAHOME.strip()

    # try to autodetect different PASA distributions
    if os.path.isfile(os.path.join(PASA, 'Launch_PASA_pipeline.pl')):  # then v2.3.0 or newer
        LAUNCHPASA = os.path.join(PASA, 'Launch_PASA_pipeline.pl')
        PASAVERSION = '2.3.0'
    elif os.path.isfile(os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl')):  # older version
        LAUNCHPASA = os.path.join(PASA, 'scripts', 'Launch_PASA_pipeline.pl')
        args.pasa_db = 'mysql'  # sqlite not available
        PASAVERSION = '2.2.0'

    if not args.TRINITYHOME:
        try:
            TRINITY = os.environ["TRINITYHOME"].strip()
        except KeyError:
            try:
                TRINITY = os.environ["TRINITY_HOME"].strip()
            except KeyError:
                lib.log.error(
                    "$TRINITYHOME nor $TRINITY_HOME environmental variable not found, TRINITY is not properly configured. You can use the --TRINITYHOME argument to specify a path at runtime.")
                sys.exit(1)
    else:
        TRINITY = args.TRINITYHOME.strip()

    programs = ['fasta', 'minimap2', 'tbl2asn', 'hisat2', 'hisat2-build', 'kallisto',
                'Trinity', 'bedtools', 'java', LAUNCHPASA, os.path.join(PASA, 'bin', 'seqclean')]
    if not args.no_trimmomatic:
        programs.append('trimmomatic')
    programs += args.aligners
    lib.CheckDependencies(programs)

    # take care of some preliminary checks
    if args.sbt == 'SBT':
        SBT = os.path.join(parentdir, 'config', 'test.sbt')
        lib.log.info(
            "No NCBI SBT file given, will use default, for NCBI submissions pass one here '--sbt'")
    else:
        SBT = args.sbt

    # setup output files
    gffout = os.path.join(tmpdir, 'genome.gff3')
    proteinsout = os.path.join(tmpdir, 'genome.proteins.fa')
    trnaout = os.path.join(tmpdir, 'genome.trna.gff3')
    fastaout = os.path.join(tmpdir, 'genome.fa')
    spliceout = os.path.join(tmpdir, 'genome.ss')
    exonout = os.path.join(tmpdir, 'genome.exons')

    # check input, allow for passing the output directory of funannotate, otherwise must be gbk or gbff files
    # set read inputs to None, populate as you go
    existingStats = False
    s_reads, l_reads, r_reads, trim_left, trim_right, trim_single, left_norm, right_norm, single_norm, all_reads, trim_reads, norm_reads, GBK, trinity_results, pasaConfigFile, PBiso, nanocdna, nanomrna, long_clean, pb_iso, nano_cdna, nano_mrna, stringtieGTF, longReadClean, shortBAM = (
        None,)*25
    if args.input:
        if os.path.isdir(args.input):
            if os.path.isdir(os.path.join(args.input, 'predict_results')):
                for file in os.listdir(os.path.join(args.input, 'predict_results')):
                    if file.endswith('.gbk'):
                        GBK = os.path.join(args.input, 'predict_results', file)
                    if file.endswith('.stats.json'):
                        existingStats = os.path.join(args.input, 'predict_results', file)
            # now lets also check if training folder/files are present, as then can pull all the data you need for update directly
            # then funannotate train has been run, try to get reads, trinity, PASA
            if os.path.isdir(os.path.join(args.input, 'training')):
                inputDir = os.path.join(args.input, 'training')
                if lib.checkannotations(os.path.join(inputDir, 'left.fq.gz')):
                    l_reads = os.path.join(inputDir, 'left.fq.gz')
                if lib.checkannotations(os.path.join(inputDir, 'right.fq.gz')):
                    r_reads = os.path.join(inputDir, 'right.fq.gz')
                if lib.checkannotations(os.path.join(inputDir, 'single.fq.gz')):
                    s_reads = os.path.join(inputDir, 'single.fq.gz')
                if lib.checkannotations(os.path.join(inputDir, 'trimmomatic', 'trimmed_left.fastq.gz')):
                    trim_left = os.path.join(
                        inputDir, 'trimmomatic', 'trimmed_left.fastq.gz')
                if lib.checkannotations(os.path.join(inputDir, 'trimmomatic', 'trimmed_right.fastq.gz')):
                    trim_right = os.path.join(
                        inputDir, 'trimmomatic', 'trimmed_right.fastq.gz')
                if lib.checkannotations(os.path.join(inputDir, 'trimmomatic', 'trimmed_single.fastq.gz')):
                    trim_single = os.path.join(
                        inputDir, 'trimmomatic', 'trimmed_single.fastq.gz')
                if lib.checkannotations(os.path.join(inputDir, 'normalize', 'left.norm.fq')):
                    left_norm = os.path.join(
                        inputDir, 'normalize', 'left.norm.fq')
                if lib.checkannotations(os.path.join(inputDir, 'normalize', 'right.norm.fq')):
                    right_norm = os.path.join(
                        inputDir, 'normalize', 'right.norm.fq')
                if lib.checkannotations(os.path.join(inputDir, 'normalize', 'single.norm.fq')):
                    single_norm = os.path.join(
                        inputDir, 'normalize', 'single.norm.fq')
                if lib.checkannotations(os.path.join(inputDir, 'nano-mrna.fasta')):
                    nanomrna = os.path.join(inputDir, 'nano-mrna.fasta')
                if lib.checkannotations(os.path.join(inputDir, 'nano-cdna.fasta')):
                    nanocdna = os.path.join(inputDir, 'nano-cdna.fasta')
                if lib.checkannotations(os.path.join(inputDir, 'iso-seq.fasta')):
                    PBiso = os.path.join(inputDir, 'iso-seq.fasta')
                long_clean = (PBiso, nanocdna, nanomrna)
                if l_reads or s_reads:
                    all_reads = (l_reads, r_reads, s_reads)
                if trim_left or trim_single:
                    trim_reads = (trim_left, trim_right, trim_single)
                if left_norm or single_norm:
                    norm_reads = (left_norm, right_norm, single_norm)
                if lib.checkannotations(os.path.join(inputDir, 'funannotate_train.stringtie.gtf')):
                    stringtieGTF = os.path.join(
                        inputDir, 'funannotate_train.stringtie.gtf')
                if lib.checkannotations(os.path.join(inputDir, 'funannotate_long-reads.fasta')):
                    longReadClean = os.path.join(
                        inputDir, 'funannotate_long-reads.fasta')
                if lib.checkannotations(os.path.join(inputDir, 'funannotate_train.trinity-GG.fasta')):
                    trinity_results = os.path.join(
                        inputDir, 'funannotate_train.trinity-GG.fasta')
                if lib.checkannotations(os.path.join(inputDir, 'funannotate_train.coordSorted.bam')):
                    shortBAM = os.path.join(
                        inputDir, 'funannotate_train.coordSorted.bam')
                if args.pasa_config:
                    pasaConfigFile = args.pasa_config
                elif os.path.isfile(os.path.join(inputDir, 'pasa', 'alignAssembly.txt')):
                    pasaConfigFile = os.path.join(
                        inputDir, 'pasa', 'alignAssembly.txt')
            # let user know which files are being re-used
            look4files = [s_reads, l_reads, r_reads, trim_left, trim_right, trim_single, left_norm,
                          right_norm, single_norm, trinity_results, longReadClean, pasaConfigFile, shortBAM, stringtieGTF]
            look4files_names = ['Single reads', 'Forward reads', 'Reverse reads', 'Forward Q-trimmed reads', 'Reverse Q-trimmed reads', 'Single Q-trimmed reads', 'Forward normalized reads',
                                'Reverse normalized reads', 'Single normalized reads', 'Trinity results', 'Long-read results', 'PASA config file', 'BAM alignments', 'StringTie GTF']
            files_used = []
            for i, x in enumerate(look4files):
                if x is not None:
                    files_used.append('\t'+look4files_names[i]+': '+x)
            if len(files_used) > 0:
                lib.log.info('Found relevant files in %s, will re-use them:\n%s' %
                             (inputDir, '\n'.join(files_used)))
        else:
            GBK = args.input
        # check if RefSeq --> NCBI does not want you to reannotate RefSeq genomes
        if GBK is None:
            print("No GBK file found")
            sys.exit(1)
        elif lib.checkRefSeq(GBK):
            lib.log.error(
                '%s is a NCBI RefSeq genome, to reannotate please use original submission.' % GBK)
            sys.exit(1)
        # split GenBank into parts
        locustag, genenumber, justify = gbk2pasaNEW(
            GBK, gffout, trnaout, fastaout, spliceout, exonout, proteinsout)
        organism, strain, isolate, accession, WGS_accession, gb_gi, version = lib.getGBKinfo(
            GBK)
        lib.log.info("Reannotating %s, NCBI accession: %s" %
                     (organism, WGS_accession))
    else:
        if args.gff and args.fasta:
            if not args.species:
                lib.log.error(
                    "Input error: please enter a name for -s,--species")
                sys.exit(1)
            shutil.copyfile(args.fasta, fastaout)
            locustag, genenumber, justify = gff2pasa(
                args.gff, fastaout, gffout, trnaout, spliceout, exonout)
            organism, strain, isolate, accession, WGS_accession, gb_gi, version = (
                None,)*7
        else:
            lib.log.error(
                "Error in input: pass either funannotate directory or GenBank file to -i,--input; or GFF3 to -g,--gff and genome FASTA to -f,--fasta.")
            sys.exit(1)

    lib.log.info("Previous annotation consists of: {:,} protein coding gene models and {:,} non-coding gene models".format(
        lib.countGFFgenes(gffout), lib.countGFFgenes(trnaout)))

    # check if organism/species/isolate passed at command line, if so, overwrite what you detected.
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

    # check input reads
    # get absolute paths for reads and concate if there are multiple
    if not all_reads:
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
                    lib.log.info(
                        "Multiple inputs for --single detected, concatenating SE reads")
                    lib.concatenateReads(single_reads, s_reads)
                else:
                    s_reads = single_reads[0]
                if s_reads.endswith('.fq'):
                    lib.Fzip_inplace(s_reads, args.cpus)
                    s_reads = s_reads+'.gz'
                if not lib.checkannotations(os.path.join(tmpdir, 'single.fq.gz')):
                    if os.path.dirname(os.path.abspath(tmpdir)) != os.path.dirname(os.path.abspath(s_reads)):
                        lib.SafeRemove(os.path.join(tmpdir, 'single.fq.gz'))
                        os.symlink(os.path.realpath(s_reads),
                                   os.path.join(tmpdir, 'single.fq.gz'))
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
                # since I can't get the comma separated input to work through subprocess, lets concatenate reads
                if left_reads[0].endswith('.gz'):
                    ending = '.fq.gz'
                else:
                    ending = '.fq'
                l_reads = os.path.join(tmpdir, 'left'+ending)
                r_reads = os.path.join(tmpdir, 'right'+ending)
                if len(left_reads) > 1:
                    lib.log.info(
                        "Multiple inputs for --left and --right detected, concatenating PE reads")
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
                        lib.SafeRemove(os.path.join(tmpdir, 'left.fq.gz'))
                        os.symlink(os.path.realpath(l_reads),
                                   os.path.join(tmpdir, 'left.fq.gz'))
                if not lib.checkannotations(os.path.join(tmpdir, 'right.fq.gz')):
                    if os.path.dirname(os.path.abspath(tmpdir)) != os.path.dirname(os.path.abspath(r_reads)):
                        lib.SafeRemove(os.path.join(tmpdir, 'right.fq.gz'))
                        os.symlink(os.path.realpath(r_reads),
                                   os.path.join(tmpdir, 'right.fq.gz'))
        else:
            l_reads = os.path.join(tmpdir, 'left.fq.gz')
            r_reads = os.path.join(tmpdir, 'right.fq.gz')

        # get tuple of input reads so you can parse them in downstream tools
        all_reads = (l_reads, r_reads, s_reads)

    lib.log.debug('Input reads: {:}'.format(all_reads))
    # trimmomatic on reads, first run PE
    if not trim_reads:
        if args.no_trimmomatic or args.trinity or left_norm or single_norm or args.left_norm or args.single_norm:
            lib.log.info("Trimmomatic will be skipped")
            trim_left = l_reads
            trim_right = r_reads
            trim_single = s_reads
        else:
            # check if they exist already in folder
            if not os.path.isfile(os.path.join(tmpdir, 'trimmomatic', 'trimmed_left.fastq.gz')) or not os.path.isfile(os.path.join(tmpdir, 'trimmomatic', 'trimmed_right.fastq.gz')):
                if all_reads[0] and all_reads[1]:
                    trim_left, trim_right = runTrimmomaticPE(
                        l_reads, r_reads, cpus=args.cpus)
                else:
                    trim_left, trim_right = (None,)*2
            else:
                trim_left, trim_right = os.path.join(tmpdir, 'trimmomatic', 'trimmed_left.fastq.gz'), os.path.join(
                    tmpdir, 'trimmomatic', 'trimmed_right.fastq.gz')
            if not os.path.isfile(os.path.join(tmpdir, 'trimmomatic', 'trimmed_single.fastq.gz')) and s_reads:
                if all_reads[2]:
                    trim_single = runTrimmomaticSE(s_reads, cpus=args.cpus)
                else:
                    trim_single = None
            else:
                if s_reads:
                    trim_single = os.path.join(
                        tmpdir, 'trimmomatic', 'trimmed_single.fastq.gz')
                else:
                    trim_single = None
        # get tuple of trimmed reads
        trim_reads = (trim_left, trim_right, trim_single)

    lib.log.debug('Quality trimmed reads: {:}'.format(trim_reads))
    # check that reads are present and make sure they follow trinity naming conventions, i.e. either illumina default or /1 and /2 to PE reads
    for read in trim_reads:
        if read:
            if not os.path.isfile(read):
                lib.log.error("Trimmomatic failed, %s does not exist." % read)
                sys.exit(1)
    # PE reads are passed, lets make sure they have proper naming
    if trim_reads[0] and trim_reads[1]:
        # if needed to fix they will be fixed in place
        lib.CheckFASTQandFix(trim_reads[0], trim_reads[1])

    # normalize reads
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
            # check if exists
            if trim_left and trim_right:
                if not os.path.islink(os.path.join(tmpdir, 'normalize', 'left.norm.fq')) or not os.path.islink(os.path.join(tmpdir, 'normalize', 'right.norm.fq')):
                    if not all(v is None for v in trim_reads):
                        left_norm, right_norm, single_norm = runNormalization(trim_reads, args.memory, cpus=args.cpus,
                                                                              stranded=args.stranded, min_coverage=args.min_coverage, coverage=args.coverage)
                else:
                    left_norm, right_norm = os.path.join(tmpdir, 'normalize', 'left.norm.fq'), os.path.join(
                        tmpdir, 'normalize', 'right.norm.fq')
                    if os.path.islink(os.path.join(tmpdir, 'normalize', 'single.norm.fq')):
                        single_norm = os.path.join(
                            tmpdir, 'normalize', 'single.norm.fq')
            if trim_single:
                if not os.path.islink(os.path.join(tmpdir, 'normalize', 'single.norm.fq')) and not trim_left and not trim_right and trim_single:
                    if not all(v is None for v in trim_reads):
                        left_norm, right_norm, single_norm = runNormalization(trim_reads, args.memory, cpus=args.cpus,
                                                                              stranded=args.stranded, min_coverage=args.min_coverage, coverage=args.coverage)
                else:
                    if os.path.islink(os.path.join(tmpdir, 'normalize', 'single.norm.fq')):
                        single_norm = os.path.join(
                            tmpdir, 'normalize', 'single.norm.fq')
                    else:
                        single_norm = None
        norm_reads = (left_norm, right_norm, single_norm)

    lib.log.debug('Normalized reads: {:}'.format(norm_reads))

    # check if long reads are passed, get full path
    if args.pacbio_isoseq:
        pb_iso = os.path.abspath(args.pacbio_isoseq)
    if args.nanopore_cdna:
        nano_cdna = os.path.abspath(args.nanopore_cdna)
    if args.nanopore_mrna:
        nano_mrna = os.path.abspath(args.nanopore_mrna)
    long_reads = (pb_iso, nano_cdna, nano_mrna)
    lib.log.debug('Long reads: {:}'.format(long_reads))

    if not trinity_results and all(v is None for v in norm_reads) and all(v is None for v in long_reads):
        lib.log.error('No reads to generate transcriptome assemblies, exiting')
        sys.exit(1)
    if trinity_results and all(v is None for v in norm_reads) and all(v is None for v in long_reads):
        lib.log.error(
            'Trinity results detected, but no RNA-seq reads detected, exiting')
        sys.exit(1)
    if not long_clean:
        long_clean = (None, None, None)
    if not longReadClean:
        if not all(v is None for v in long_reads):
            # get long read FASTA file
            longReadFA = os.path.join(tmpdir, 'long-reads.fasta')
            longReadClean = os.path.join(tmpdir, 'long-reads.fasta.clean')
            PBiso = os.path.join(tmpdir, 'iso-seq.fasta')
            nanocdna = os.path.join(tmpdir, 'nano-cdna.fasta')
            nanomrna = os.path.join(tmpdir, 'nano-mrna.fasta')
            if not all(v is None for v in long_reads):
                if not lib.checkannotations(longReadFA):
                    long_readsFA, long_clean = long2fasta(long_reads, args.cpus, tmpdir, os.path.abspath(
                        longReadFA), os.path.abspath(longReadClean))
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

    trinity_transcripts = os.path.join(tmpdir, 'trinity.fasta')
    if not trinity_results:
        # now run Trinity with trimmomatic and read normalization
        shortBAM = os.path.join(tmpdir, 'hisat2.coordSorted.bam')
        if not lib.checkannotations(trinity_transcripts):
            if args.trinity:
                lib.log.info(
                    "Parsing assembled trinity data: {:}".format(args.trinity))
                shutil.copyfile(os.path.abspath(
                    args.trinity), trinity_transcripts)
            else:
                if not all(v is None for v in norm_reads):
                    # run trinity genome guided
                    # runTrinityGG(genome, norm_reads, longReadClean, shortBAM, trinity_transcripts)
                    cmd = [sys.executable, os.path.join(parentdir, 'aux_scripts', 'trinity.py'),
                           '-f', fastaout, '-o', trinity_transcripts, '-b', shortBAM, '-t', tmpdir,
                           '--stranded', args.stranded, '--max_intronlen', str(
                               args.max_intronlen),
                           '--cpus', str(args.cpus), '--TRINITYHOME', TRINITY, '--memory', args.memory,
                           '--logfile', os.path.join(args.out, 'logfiles', 'funannotate-trinity.log')]
                    if args.jaccard_clip:
                        cmd.append('--jaccard_clip')
                    if norm_reads[2]:  # single
                        cmd += ['-s', norm_reads[2]]
                    else:
                        cmd += ['-l', norm_reads[0], '-r', norm_reads[1]]
                    if lib.checkannotations(longReadClean):
                        cmd += ['--long', longReadClean]
                    if not args.progress:
                        cmd.append('--no-progress')
                    # run trinity
                    subprocess.call(cmd)
                    if not lib.checkannotations(trinity_transcripts):
                        lib.log.info('ERROR: Trinity de novo assembly failed')
                        sys.exit(1)
        else:
            lib.log.info("Existing Trinity results found: {:}".format(
                trinity_transcripts))
    else:
        shutil.copyfile(trinity_results, trinity_transcripts)

    if not stringtieGTF:
        # if stringtie installed, run on shortBAM incorporate into PASA later on
        stringtieGTF = os.path.join(tmpdir, 'funannotate_train.stringtie.gtf')
        if not lib.checkannotations(stringtieGTF):
            if lib.which('stringtie') and lib.checkannotations(shortBAM):
                lib.log.info(
                    'StringTie installed, running StringTie on Hisat2 coordsorted BAM')
                cmd = ['stringtie', '-p', str(args.cpus)]
                if args.stranded != 'no':
                    if args.stranded.startswith('R'):
                        cmd = cmd + ['--rf']
                    else:
                        cmd = cmd + ['--fr']
                cmd = cmd + [shortBAM]
                lib.runSubprocess8(cmd, '.', lib.log, stringtieGTF)

    # run SeqClean to clip polyA tails and remove low quality seqs.
    cleanTranscripts = trinity_transcripts+'.clean'
    if not lib.checkannotations(cleanTranscripts) and lib.checkannotations(trinity_transcripts):
        runSeqClean(trinity_transcripts, tmpdir, cpus=args.cpus)

    # map long reads and Trinity transcripts to genome for PASA
    allBAM = os.path.join(tmpdir, 'transcript.alignments.bam')
    trinityBAM = os.path.join(tmpdir, 'trinity.alignments.bam')
    if not lib.checkannotations(allBAM):
        trinity_transcripts, cleanTranscripts = mapTranscripts(
            fastaout, long_clean, cleanTranscripts, tmpdir, trinityBAM, allBAM, cpus=args.cpus, max_intronlen=args.max_intronlen)
    else:
        if lib.checkannotations(trinityBAM):
            lib.log.info("Existing BAM alignments found: {:}, {:}".format(
                trinityBAM, allBAM))
        else:
            lib.log.info("Existing BAM alignments found: {:}".format(allBAM))

    # convert BAM to GFF3
    allGFF3 = os.path.join(tmpdir, 'transcript.alignments.gff3')
    trinityGFF3 = os.path.join(tmpdir, 'trinity.alignments.gff3')
    if not lib.checkannotations(allGFF3) and lib.checkannotations(allBAM):
        lib.log.info('Converting transcript alignments to GFF3 format')
        lib.bam2gff3(allBAM, allGFF3)
    if not lib.checkannotations(trinityGFF3) and lib.checkannotations(trinityBAM):
        lib.log.info('Converting Trinity transcript alignments to GFF3 format')
        lib.bam2gff3(trinityBAM, trinityGFF3)

    # now run PASA steps
    PASA_gff = os.path.join(tmpdir, 'pasa_final.gff3')
    if not pasaConfigFile:
        if args.pasa_gff:
            lib.log.info(
                "You passed a --pasa_gff file; are you sure this is a good idea?")
            shutil.copyfile(args.pasa_gff, PASA_gff)
        if not lib.checkannotations(PASA_gff):
            if lib.checkannotations(trinityBAM):
                runPASA(fastaout, trinity_transcripts, cleanTranscripts, os.path.abspath(trinityGFF3),
                        os.path.abspath(stringtieGTF), args.stranded, args.max_intronlen, args.cpus, gffout, organism_name,
                        PASA_gff, args.pasa_config, pasa_db=args.pasa_db, pasa_alignment_overlap=args.pasa_alignment_overlap,
                        aligners=args.aligners)
            else:
                runPASA(fastaout, os.path.abspath(longReadFA), os.path.abspath(longReadClean),
                        os.path.abspath(allGFF3), os.path.abspath(stringtieGTF), args.stranded, args.max_intronlen, args.cpus,
                        gffout, organism_name, PASA_gff, args.pasa_config, pasa_db=args.pasa_db,
                        pasa_alignment_overlap=args.pasa_alignment_overlap, aligners=args.aligners)
    else:
        if not lib.checkannotations(PASA_gff):
            if lib.checkannotations(trinityBAM):
                runPASA(fastaout, trinity_transcripts, cleanTranscripts, os.path.abspath(trinityGFF3),
                        os.path.abspath(stringtieGTF), args.stranded, args.max_intronlen, args.cpus, gffout, organism_name,
                        PASA_gff, pasaConfigFile, pasa_db=args.pasa_db, pasa_alignment_overlap=args.pasa_alignment_overlap,
                        aligners=args.aligners)
            else:
                runPASA(fastaout, os.path.abspath(longReadFA), os.path.abspath(longReadClean),
                        os.path.abspath(allGFF3), os.path.abspath(stringtieGTF), args.stranded, args.max_intronlen,
                        args.cpus, gffout, organism_name, PASA_gff, pasaConfigFile, pasa_db=args.pasa_db,
                        pasa_alignment_overlap=args.pasa_alignment_overlap, aligners=args.aligners)
        else:
            lib.log.info('Skipping PASA, found existing output: %s' % PASA_gff)

    # now run Kallisto steps, if mixed PE and SE reads, only PE reads will be used for Kallisto as there isn't a reasonable way to combine them
    KallistoAbundance = os.path.join(tmpdir, 'kallisto.tsv')
    if args.kallisto:
        lib.log.info("You passed a --kallisto file; are you sure this is a good idea?")
        shutil.copyfile(args.kallisto, KallistoAbundance)

    if all(v is None for v in trim_reads):
        kallistoreads = norm_reads
    else:
        kallistoreads = trim_reads
    if all(v is None for v in kallistoreads) and lib.checkannotations(longReadClean):
        # use minimap to count reads
        lib.log.info(
            'Generating relative expression values to PASA transcripts')
        PASAtranscripts = os.path.join(tmpdir, 'transcripts.fa')
        cmd = [os.path.join(
            PASA, 'misc_utilities', 'gff3_file_to_proteins.pl'), PASA_gff, fastaout, 'cDNA']
        if not lib.checkannotations(PASAtranscripts):
            lib.runSubprocess2(cmd, '.', lib.log, PASAtranscripts)
        PASAdict = pasa_transcript2gene(PASAtranscripts)
        minimapBAM = os.path.join(tmpdir, 'long-reads_transcripts.bam')
        minimap2_cmd = ['minimap2', '-ax' 'map-ont', '-t',
                        str(args.cpus), '--secondary=no', PASAtranscripts, longReadClean]
        samtools_cmd = ['samtools', 'sort', '--reference', PASAtranscripts,
                        '-@', '2', '-o', minimapBAM, '-']
        if not lib.checkannotations(minimapBAM):
            lib.log.debug('{} | {}'.format(' '.join(minimap2_cmd), ' '. join(samtools_cmd)))
            p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=FNULL)
            p2 = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=FNULL, stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()
        if not lib.checkannotations(KallistoAbundance):
            lib.mapCount(minimapBAM, PASAdict, KallistoAbundance)
    else:
        if not lib.checkannotations(KallistoAbundance):
            runKallisto(PASA_gff, fastaout, kallistoreads,
                        args.stranded, args.cpus, KallistoAbundance)
        else:
            lib.log.info(
                "Existing Kallisto output found: {:}".format(KallistoAbundance))

    # parse Kallisto results with PASA GFF
    global ExpressionValues
    ExpressionValues = {}
    BestModelGFF = os.path.join(tmpdir, 'bestmodels.gff3')
    getBestModels(PASA_gff, fastaout, KallistoAbundance,
                  args.alt_transcripts, BestModelGFF)

    # make sure tRNA models don't overlap new gene models
    cleanTRNA = os.path.join(tmpdir, 'trna.no-overlaps.gff')
    lib.validate_tRNA(trnaout, BestModelGFF, False, cleanTRNA)

    # generate tbl file
    gagdir = os.path.join(tmpdir, 'tbl2asn')
    TBLFile = os.path.join(gagdir, 'genome.tbl')
    if os.path.isdir(gagdir):
        shutil.rmtree(gagdir)
    os.makedirs(gagdir)
    GFF2tblCombinedNEW(BestModelGFF, fastaout, cleanTRNA, locustag, genenumber,
                       justify, args.SeqCenter, args.SeqAccession, TBLFile,
                       alt_transcripts=args.alt_transcripts)

    # need a function here to clean up the ncbi tbl file if this is a reannotation
    # a reannotation would have a WGS_accession, if none, then it is a first pass and not from genbank
    if WGS_accession:
        os.rename(os.path.join(tmpdir, 'tbl2asn', 'genome.tbl'),
                  os.path.join(tmpdir, 'tbl2asn', 'genome.tbl.bak'))
        p2g = {}
        if args.p2g:  # load into dictionary
            shutil.copyfile(args.p2g, os.path.join(
                args.out, 'update_results', 'ncbi.p2g'))
            with open(args.p2g, 'r') as input:
                for line in input:
                    cols = line.split('\t')
                    if not cols[0] in p2g:
                        p2g[cols[0]] = cols[1]
        with open(os.path.join(tmpdir, 'tbl2asn', 'genome.tbl'), 'w') as outfile:
            with open(os.path.join(tmpdir, 'tbl2asn', 'genome.tbl.bak'), 'r') as infile:
                for line in infile:
                    line = line.replace('\n', '')
                    if line.startswith('\t\t\tprotein_id') or line.startswith('\t\t\ttranscript_id'):
                        ID = line.rsplit('|', 1)[-1].replace('_mrna', '')
                        type = 'prot'
                        if 'transcript_id' in line:
                            type = 'transcript'
                        if not ID in p2g:
                            if type == 'prot':
                                outfile.write(
                                    '\t\t\tprotein_id\tgnl|%s|%s\n' % (WGS_accession, ID))
                            elif type == 'transcript':
                                outfile.write(
                                    '\t\t\ttranscript_id\tgnl|%s|%s_mrna\n' % (WGS_accession, ID))
                        else:
                            p2gID = p2g.get(ID)
                            if type == 'prot':
                                outfile.write('\t\t\tprotein_id\tgnl|%s|%s|gb|%s\n' % (
                                    WGS_accession, ID, p2gID))
                            elif type == 'transcript':
                                outfile.write(
                                    '\t\t\ttranscript_id\tgnl|%s|%s_mrna\n' % (WGS_accession, ID))
                    else:
                        outfile.write('%s\n' % line)

    # run tbl2asn in new directory directory
    lib.log.info("Converting to Genbank format")
    discrep = os.path.join(args.out, 'update_results',
                           organism_name + '.discrepency.report.txt')
    # this would mean it is a GenBank reannotation, so update accordingly. else it is just 1st version.
    if version and WGS_accession:
        rev_version = int(version) + 1
    else:
        rev_version = 1

    # have to run as subprocess because of multiprocessing issues
    cmd = [sys.executable, os.path.join(parentdir, 'aux_scripts', 'tbl2asn_parallel.py'),
           '-i', TBLFile, '-f', fastaout, '-o', gagdir, '--sbt', SBT, '-d', discrep,
           '-s', organism, '-t', args.tbl2asn, '-v', str(rev_version), '-c', str(args.cpus)]
    if isolate:
        cmd += ['--isolate', isolate]
    if strain:
        cmd += ['--strain', strain]
    lib.log.debug(' '.join(cmd))
    subprocess.call(cmd)

    # grab results, populate results output directory
    final_fasta = os.path.join(
        args.out, 'update_results', organism_name + '.scaffolds.fa')
    final_gff = os.path.join(
        args.out, 'update_results', organism_name + '.gff3')
    final_gbk = os.path.join(
        args.out, 'update_results', organism_name + '.gbk')
    final_tbl = os.path.join(
        args.out, 'update_results', organism_name + '.tbl')
    final_proteins = os.path.join(
        args.out, 'update_results', organism_name + '.proteins.fa')
    final_transcripts = os.path.join(
        args.out, 'update_results', organism_name + '.mrna-transcripts.fa')
    final_cds_transcripts = os.path.join(
        args.out, 'update_results', organism_name + '.cds-transcripts.fa')
    final_validation = os.path.join(
        args.out, 'update_results', organism_name+'.validation.txt')
    final_error = os.path.join(
        args.out, 'update_results', organism_name+'.error.summary.txt')
    final_fixes = os.path.join(
        args.out, 'update_results', organism_name+'.models-need-fixing.txt')
    final_stats = os.path.join(args.out, 'update_results', organism_name+'.stats.json')

    # retrieve files/reorganize
    shutil.copyfile(os.path.join(gagdir, 'genome.gbf'), final_gbk)
    shutil.copyfile(os.path.join(gagdir, 'genome.tbl'), final_tbl)
    shutil.copyfile(os.path.join(gagdir, 'genome.val'), final_validation)
    shutil.copyfile(os.path.join(gagdir, 'errorsummary.val'), final_error)

    lib.log.info("Collecting final annotation files")
    lib.tbl2allout(final_tbl, fastaout, final_gff, final_proteins,
                   final_transcripts, final_cds_transcripts, final_fasta)
    lib.annotation_summary(fastaout, final_stats, tbl=final_tbl,
                           transcripts=allGFF3, previous=existingStats)
    # since no place to write the WGS accession to, output to a file for reading by funannotate annotate
    with open(os.path.join(args.out, 'update_results', 'WGS_accession.txt'), 'w') as out:
        out.write('%s\n' % WGS_accession)

    # last step will be to compare old gene models versus updated ones, outputting a tsv file describing the changes
    Changes = os.path.join(args.out, 'update_results',
                           organism_name + '.pasa-reannotation.changes.txt')
    if args.gff and args.fasta:
        compareAnnotations2(args.gff, final_gbk, Changes, args=args)
    else:
        compareAnnotations2(GBK, final_gbk, Changes, args=args)

    lib.log.info(
        "Funannotate update is finished, output files are in the %s/update_results folder" % (args.out))
    errors = lib.ncbiCheckErrors(
        final_error, final_validation, locustag, final_fixes)
    if errors > 0:
        print('-------------------------------------------------------')
        lib.log.info("Manually edit the tbl file %s, then run:\n\nfunannotate fix -i %s -t %s\n" %
                     (final_tbl, final_gbk, final_tbl))
        lib.log.info(
            "After the problematic gene models are fixed, you can proceed with functional annotation.")

    lib.log.info("Your next step might be functional annotation, suggested commands:\n\
-------------------------------------------------------\n\
Run InterProScan (Docker required): \nfunannotate iprscan -i {:} -m docker -c {:}\n\n\
Run antiSMASH: \nfunannotate remote -i {:} -m antismash -e youremail@server.edu\n\n\
Annotate Genome: \nfunannotate annotate -i {:} --cpus {:} --sbt yourSBTfile.txt\n\
-------------------------------------------------------\n\
                ".format(args.out,
                         args.cpus,
                         args.out,
                         args.out,
                         args.cpus))


if __name__ == "__main__":
    main(sys.argv[1:])
