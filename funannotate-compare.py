#!/usr/bin/env python
from __future__ import division

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time, csv, glob, warnings
from Bio import SeqIO
from natsort import natsorted
import pandas as pd
import numpy as np

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-compare.py', usage="%(prog)s [options] genome1.gbk genome2.gbk",
    description='''Funannotate comparative genomics.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', nargs='+', help='List of annotated genomes in GenBank format')
parser.add_argument('-o','--out', default='funannotate_compare', help='Name of output folder')
args=parser.parse_args()

#create log file
log_name = 'funnannotate-compare.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.log.info("Operating system: %s, %i cores, ~ %i GB RAM" % (sys.platform, multiprocessing.cpu_count(), lib.MemoryCheck()))

#make output folder
if not os.path.isdir(args.out):
    os.makedirs(args.out)
go_folder = os.path.join(args.out, 'go_terms')
if os.path.isdir(go_folder):
    shutil.rmtree(go_folder)
    os.makedirs(go_folder)
else:
    os.makedirs(go_folder)

#loop through each genome
stats = []
merops = []
ipr = []
cazy = []
num_input = len(args.input)
lib.log.info("Now parsing %i genomes" % num_input)
for i in range(0,num_input):
    if not args.input[i].endswith('.gbk'):
        lib.log.error("Error, one of input files is not in GenBank format")
        os._exit(1)
    else: #split arguments into genomes and run a bunch of stats/comparisons
        stats.append(lib.genomeStats(args.input[i]))
        merops.append(lib.getStatsfromNote(args.input[i], 'MEROPS'))
        ipr.append(lib.getStatsfromDbxref(args.input[i], 'InterPro'))
        cazy.append(lib.getStatsfromNote(args.input[i], 'CAZy'))
        lib.parseGOterms(args.input[i], go_folder, stats[i][0].replace(' ', '_'))

#add species names to pandas table
names = []
for i in stats:
    sci_name = i[0]
    genus = sci_name.split(' ')[0]
    species = ' '.join(sci_name.split(' ')[1:])
    abbrev = genus[:1] + '.'
    final_name = abbrev + ' ' + species
    names.append(final_name)
    
####InterProScan##################################
lib.log.info("Summarizing InterProScan results")
#convert to counts
IPRdf = lib.convert2counts(ipr)
IPRdf.fillna(0, inplace=True) #fill in zeros for missing data
IPRdf['species'] = names
IPRdf.set_index('species', inplace=True)

#transpose and write to csv file
ipr2 = IPRdf.transpose()
ipr2.to_csv(os.path.join(args.out, 'interproscan.results.csv'))

#NMDS analysis of InterPro Domains
lib.distance2mds(IPRdf, 'braycurtis', 'InterProScan', os.path.join(args.out, 'InterProScan.nmds.pdf'))

    
####MEROPS################################
lib.log.info("Summarizing MEROPS protease results")
MEROPS = {'A': 'Aspartic Peptidase', 'C': 'Cysteine Peptidase', 'G': 'Glutamic Peptidase', 'M': 'Metallo Peptidase', 'N': 'Asparagine Peptide Lyase', 'P': 'Mixed Peptidase','S': 'Serine Peptidase', 'T': 'Threonine Peptidase', 'U': 'Unknown Peptidase'}
#convert to counts
meropsdf = lib.convert2counts(merops)
meropsdf.fillna(0, inplace=True)
meropsdf['species'] = names
meropsdf.set_index('species', inplace=True)

#make a simple table with just these numbers
meropsA = meropsdf.filter(regex='A').sum(numeric_only=True, axis=1)
meropsC = meropsdf.filter(regex='C').sum(numeric_only=True, axis=1)
meropsG = meropsdf.filter(regex='G').sum(numeric_only=True, axis=1)
meropsM = meropsdf.filter(regex='M').sum(numeric_only=True, axis=1)
meropsN = meropsdf.filter(regex='N').sum(numeric_only=True, axis=1)
meropsP = meropsdf.filter(regex='P').sum(numeric_only=True, axis=1)
meropsS = meropsdf.filter(regex='S').sum(numeric_only=True, axis=1)
meropsT = meropsdf.filter(regex='T').sum(numeric_only=True, axis=1)
meropsU = meropsdf.filter(regex='U').sum(numeric_only=True, axis=1)
#get totals for determining height of y-axis
totals = meropsdf.sum(numeric_only=True, axis=1)
max_num = max(totals)
round_max = int(lib.roundup(max_num))
diff = round_max - int(max_num)
if diff < 100:
    ymax = round_max + 100
else:
    ymax = round_max
#recombine sums
enzymes = ['A', 'C', 'G', 'M', 'N', 'P', 'S', 'T', 'U']
meropsShort = pd.concat([meropsA, meropsC, meropsG, meropsM, meropsN, meropsP, meropsS, meropsT, meropsU], axis=1, keys=enzymes)
meropsShort['species'] = names
meropsShort.set_index('species', inplace=True)
#remove any columns with no hits
meropsShort = meropsShort.loc[:, (meropsShort != 0).any(axis=0)]
meropsall = meropsdf.transpose()

#write to file
meropsdf.transpose().to_csv(os.path.join(args.out, 'MEROPS.all.results.csv'))
meropsShort.transpose().to_csv(os.path.join(args.out, 'MEROPS.summary.results.csv'))

#draw plots for merops data
#stackedbar graph
lib.drawStackedBar(meropsShort, 'MEROPS', MEROPS, ymax, os.path.join(args.out, 'MEROPS.graph.pdf'))

#drawheatmap of all merops families
lib.drawHeatmap(meropsall, 'BuPu', os.path.join(args.out, 'MEROPS.heatmap.pdf'), False)
#######################################################

#####run CAZy routine#################################
lib.log.info("Summarizing CAZyme results")
#convert to counts
CAZydf = lib.convert2counts(cazy)

#with CAZy there are 7 possible families
CAZY = {'CBM': 'Carbohydrate-binding module', 'CE': 'Carbohydrate esterase','GH': 'Glycoside hydrolase', 'GT': 'Glycosyltransferase', 'PL': 'Polysaccharide lyase', 'AA': 'Auxillary activities'}
#make a simple table with just these numbers
cazyAA = CAZydf.filter(regex='AA').sum(numeric_only=True, axis=1)
cazyGT = CAZydf.filter(regex='GT').sum(numeric_only=True, axis=1)
cazyPL = CAZydf.filter(regex='PL').sum(numeric_only=True, axis=1)
cazyCE = CAZydf.filter(regex='CE').sum(numeric_only=True, axis=1)
cazyCBM = CAZydf.filter(regex='CBM').sum(numeric_only=True, axis=1)
cazyGH = CAZydf.filter(regex='GH').sum(numeric_only=True, axis=1)
#get totals for determining height of y-axis
totals = CAZydf.sum(numeric_only=True, axis=1)
max_num = max(totals)
round_max = int(lib.roundup(max_num))
diff = round_max - int(max_num)
if diff < 100:
    ymax = round_max + 100
else:
    ymax = round_max
#print max_num, round_max, diff, ymax
enzymes = ['AA', 'CBM', 'CE', 'GH', 'GT', 'PL']
CAZyShort = pd.concat([cazyAA, cazyCBM, cazyCE, cazyGH, cazyGT, cazyPL], axis=1, keys=enzymes)
CAZydf['species'] = names
CAZyShort['species'] = names
CAZydf.set_index('species', inplace=True)
CAZyShort.set_index('species', inplace=True)
cazyall = CAZydf.transpose()

#write to file
CAZydf.transpose().to_csv(os.path.join(args.out, 'CAZyme.all.results.csv'))
CAZyShort.transpose().to_csv(os.path.join(args.out, 'CAZyme.summary.results.csv'))

#draw stacked bar graph for CAZY's
lib.drawStackedBar(CAZyShort, 'CAZyme', CAZY, ymax, os.path.join(args.out, 'CAZy.graph.pdf'))

#drawheatmap of all CAZys
lib.drawHeatmap(cazyall, 'YlOrRd', os.path.join(args.out, 'CAZy.heatmap.pdf'), False)
########################################################

####GO Terms, GO enrichment
#concatenate all genomes into a population file
lib.log.info("Running GO enrichment for each genome")
with open(os.path.join(go_folder, 'population.txt'), 'w') as pop:
    for file in os.listdir(go_folder):
        if not file.startswith('associations'):
            file = os.path.join(go_folder, file)
            with open(file) as input:
                pop.write(input.read())

#now loop through each genome comparing to population
for f in os.listdir(go_folder):
    if f.startswith('associations'):
        continue
    if f.startswith('population'):
        continue
    file = os.path.join(go_folder, f)
    base = f.replace('.txt', '')
    proc = subprocess.Popen(['find_enrichment.py', '--obo', os.path.join(currentdir, 'DB', 'go.obo'), '--pval', '0.001', '--alpha', '0.001', '--fdr', file, os.path.join(go_folder, 'population.txt'), os.path.join(go_folder, 'associations.txt')], stderr=FNULL, stdout=subprocess.PIPE)
    result = proc.communicate()[0].split('\n')
    #print result
    under = []
    over = []
    #now parse result
    for line in result:
        if line.startswith('#'):
            continue
        if line.startswith('id'):
            header = line
        if line.startswith('GO'):
            cols = line.split('\t')
            #print cols
            if float(cols[9]) <= 0.001: #if fdr is significant, then save line
                if cols[1] == 'p': #this is under-represented
                    under.append(line)
                if cols[1] == 'e': #over-represented
                    over.append(line)
    lib.log.info("Found %i over-represented and %i under-represented GO terms in %s" % (len(over), len(under), base))   
    if over:
        with open(os.path.join(args.out, base+'.goterms.over_represented.txt'), 'w') as output:
            output.write("%s\n" % (header))
            for line in over:
                output.write("%s\n" % (line))
    if under:
        with open(os.path.join(args.out, base+'.goterms.under_represented.txt'), 'w') as output:
            output.write("%s\n" % (header))
            for line in under:
                output.write("%s\n" % (line))
#################################################### 

#test script up to here
os._exit(1)

