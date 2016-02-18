#!/usr/bin/env python
from __future__ import division

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time, csv, glob, warnings, re
from goatools import obo_parser
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
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to utilize')
parser.add_argument('--go_fdr', default=0.05, type=float, help='P-value for FDR GO-enrichment')
parser.add_argument('--heatmap_stdev', default=1.0, type=float, help='Standard Deviation threshold for heatmap retention')
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

#check dependencies
PROTORTHO = os.path.join(currentdir, 'util', 'proteinortho_v5.11', 'proteinortho5.pl')
programs = ['find_enrichment.py']
lib.CheckDependencies(programs)

#make output folder
if not os.path.isdir(args.out):
    os.makedirs(args.out)
go_folder = os.path.join(args.out, 'go_terms')
protortho = os.path.join(args.out, 'protortho')
if os.path.isdir(go_folder):
    shutil.rmtree(go_folder)
    os.makedirs(go_folder)
else:
    os.makedirs(go_folder)
if not os.path.isdir(protortho):
    os.makedirs(protortho)

#copy over html file
if not os.path.isdir(os.path.join(args.out,'css')):
    lib.copyDirectory(os.path.join(currentdir, 'html_template', 'css'), os.path.join(args.out, 'css'))
if not os.path.isdir(os.path.join(args.out, 'js')):
    lib.copyDirectory(os.path.join(currentdir, 'html_template', 'js'), os.path.join(args.out, 'js'))
    
#loop through each genome
stats = []
merops = []
ipr = []
cazy = []
pfam = []
eggnog = []
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
        pfam.append(lib.getStatsfromDbxref(args.input[i], 'PFAM'))
        cazy.append(lib.getStatsfromNote(args.input[i], 'CAZy'))
        lib.parseGOterms(args.input[i], go_folder, stats[i][0].replace(' ', '_'))
        lib.gb2proteinortho(args.input[i], protortho, stats[i][0].replace(' ', '_'))
        eggnog.append(lib.getEggNogfromNote(args.input[i]))

#convert eggnog to a single dictionary
EGGNOG = { k: v for d in eggnog for k, v in d.items() }

#add species names to pandas table
names = []
for i in stats:
    sci_name = i[0]
    genus = sci_name.split(' ')[0]
    species = ' '.join(sci_name.split(' ')[1:])
    abbrev = genus[:1] + '.'
    final_name = abbrev + ' ' + species
    names.append(final_name)


#PFAM#############################################
lib.log.info("Summarizing PFAM domain results")
if not os.path.isdir(os.path.join(args.out, 'pfam')):
    os.makedirs(os.path.join(args.out, 'pfam'))

#convert to counts
pfamdf = lib.convert2counts(pfam)
pfamdf.fillna(0, inplace=True)
pfamdf['species'] = names
pfamdf.set_index('species', inplace=True)

#make an nmds
lib.distance2mds(pfamdf, 'braycurtis', 'PFAM', os.path.join(args.out, 'pfam','PFAM.nmds.pdf'))

#get the PFAM descriptions
pfamdf2 = pfamdf.transpose()
PFAM = lib.pfam2dict(os.path.join(currentdir, 'DB', 'Pfam-A.clans.tsv'))
pfam_desc = []
for i in pfamdf2.index.values:
    pfam_desc.append(PFAM.get(i))
pfamdf2['descriptions'] = pfam_desc
#write to file
pfamdf2.to_csv(os.path.join(args.out, 'pfam', 'pfam.results.csv'))
pfamdf2.reset_index(inplace=True)
pfamdf2.rename(columns = {'index':'PFAM'}, inplace=True)
pfamdf2['PFAM'] = '<a href="http://pfam.xfam.org/family/'+ pfamdf2['PFAM'].astype(str)+'">'+pfamdf2['PFAM']+'</a>'
#create html output
with open(os.path.join(args.out, 'pfam.html'), 'w') as output:
    pd.set_option('display.max_colwidth', -1)
    output.write(lib.HEADER)
    output.write(lib.PFAM)
    output.write(pfamdf2.to_html(index=False, escape=False, classes='table table-hover'))
    output.write(lib.FOOTER)

##################################################
  
####InterProScan##################################
lib.log.info("Summarizing InterProScan results")
if not os.path.isdir(os.path.join(args.out, 'interpro')):
    os.makedirs(os.path.join(args.out, 'interpro'))

#convert to counts
IPRdf = lib.convert2counts(ipr)
IPRdf.fillna(0, inplace=True) #fill in zeros for missing data
IPRdf['species'] = names
IPRdf.set_index('species', inplace=True)

#NMDS analysis of InterPro Domains
lib.distance2mds(IPRdf, 'braycurtis', 'InterProScan', os.path.join(args.out, 'interpro', 'InterProScan.nmds.pdf'))

#write to csv file
ipr2 = IPRdf.transpose()
#get IPR descriptions
INTERPRO = lib.iprxml2dict(os.path.join(currentdir, 'DB', 'interpro.xml'))
ipr_desc = []
for i in ipr2.index.values:
    ipr_desc.append(INTERPRO.get(i))
ipr2['descriptions'] = ipr_desc
ipr2.to_csv(os.path.join(args.out, 'interpro','interproscan.results.csv'))
ipr2.reset_index(inplace=True)
ipr2.rename(columns = {'index':'InterPro'}, inplace=True)
ipr2['InterPro'] = '<a href="http://www.ebi.ac.uk/interpro/entry/'+ ipr2['InterPro'].astype(str)+'">'+ipr2['InterPro']+'</a>'

#create html output
with open(os.path.join(args.out, 'interpro.html'), 'w') as output:
    pd.set_option('display.max_colwidth', -1)
    output.write(lib.HEADER)
    output.write(lib.INTERPRO)
    output.write(ipr2.to_html(index=False, escape=False, classes='table table-hover'))
    output.write(lib.FOOTER)

##############################################

    
####MEROPS################################
lib.log.info("Summarizing MEROPS protease results")
if not os.path.isdir(os.path.join(args.out, 'merops')):
    os.makedirs(os.path.join(args.out, 'merops'))

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
meropsdf.transpose().to_csv(os.path.join(args.out, 'merops', 'MEROPS.all.results.csv'))
meropsShort.transpose().to_csv(os.path.join(args.out, 'merops', 'MEROPS.summary.results.csv'))

#draw plots for merops data
#stackedbar graph
lib.drawStackedBar(meropsShort, 'MEROPS', MEROPS, ymax, os.path.join(args.out, 'merops', 'MEROPS.graph.pdf'))

#drawheatmap of all merops families where there are any differences 
stdev = meropsall.std(axis=1)
meropsall['stdev'] = stdev
df2 = meropsall[meropsall.stdev >= args.heatmap_stdev ]
meropsplot = df2.drop('stdev', axis=1)
lib.log.info("found %i/%i MEROPS familes with stdev >= %f" % (len(meropsplot), len(meropsall), args.heatmap_stdev))
lib.drawHeatmap(meropsplot, 'BuPu', os.path.join(args.out, 'merops', 'MEROPS.heatmap.pdf'), False)

meropsall.reset_index(inplace=True)
meropsall.rename(columns = {'index':'MEROPS'}, inplace=True)
meropsall['MEROPS'] = '<a href="https://merops.sanger.ac.uk/cgi-bin/famsum?family='+ meropsall['MEROPS'].astype(str)+'">'+meropsall['MEROPS']+'</a>'

#create html output
with open(os.path.join(args.out, 'merops.html'), 'w') as output:
    pd.set_option('display.max_colwidth', -1)
    output.write(lib.HEADER)
    output.write(lib.MEROPS)
    output.write(meropsall.to_html(escape=False, index=False, classes='table table-hover'))
    output.write(lib.FOOTER)

#######################################################

#####run CAZy routine#################################
lib.log.info("Summarizing CAZyme results")
if not os.path.isdir(os.path.join(args.out, 'cazy')):
    os.makedirs(os.path.join(args.out, 'cazy'))
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
CAZydf.transpose().to_csv(os.path.join(args.out, 'cazy', 'CAZyme.all.results.csv'))
CAZyShort.transpose().to_csv(os.path.join(args.out, 'cazy', 'CAZyme.summary.results.csv'))

#draw stacked bar graph for CAZY's
lib.drawStackedBar(CAZyShort, 'CAZyme', CAZY, ymax, os.path.join(args.out, 'cazy', 'CAZy.graph.pdf'))

#drawheatmap of all CAZys that have standard deviation > 0.5
stdev = cazyall.std(axis=1)
cazyall['stdev'] = stdev
df2 = cazyall[cazyall.stdev >= args.heatmap_stdev ]
cazyplot = df2.drop('stdev', axis=1)
lib.log.info("found %i/%i CAZy familes with stdev >= %f" % (len(cazyplot), len(cazyall), args.heatmap_stdev))
lib.drawHeatmap(cazyplot, 'YlOrRd', os.path.join(args.out, 'cazy', 'CAZy.heatmap.pdf'), False)

cazyall.reset_index(inplace=True)
cazyall.rename(columns = {'index':'CAZy'}, inplace=True)
cazyall['CAZy'] = '<a href="http://www.cazy.org/'+ cazyall['CAZy'].astype(str)+'.html">'+cazyall['CAZy']+'</a>'

#create html output
with open(os.path.join(args.out, 'cazy.html'), 'w') as output:
    pd.set_option('display.max_colwidth', -1)
    output.write(lib.HEADER)
    output.write(lib.CAZY)
    output.write(cazyall.to_html(escape=False, index=False, classes='table table-hover'))
    output.write(lib.FOOTER)
########################################################


####GO Terms, GO enrichment############################
if not os.path.isdir(os.path.join(args.out, 'go_enrichment')):
    os.makedirs(os.path.join(args.out, 'go_enrichment'))

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
    goa_out = os.path.join(args.out, 'go_enrichment', base+'.go.enrichment.txt')
    with open(goa_out, 'w') as output:
        subprocess.call(['find_enrichment.py', '--obo', os.path.join(currentdir, 'DB', 'go.obo'), '--pval', '0.001', '--alpha', '0.001', '--fdr', file, os.path.join(go_folder, 'population.txt'), os.path.join(go_folder, 'associations.txt')], stderr=FNULL, stdout=output)

#load into pandas and write to html
with open(os.path.join(args.out, 'go.html'), 'w') as output:
    pd.set_option('display.max_colwidth', -1)
    pd.options.mode.chained_assignment = None #turn off warning
    output.write(lib.HEADER)
    output.write(lib.GO)
    for f in os.listdir(os.path.join(args.out, 'go_enrichment')):
        if f.endswith('go.enrichment.txt'):
            file = os.path.join(args.out, 'go_enrichment', f)
            base = file.strip('.go_enrichment.txt')
            name = base.split('/')[-1]
            #check output, > 3 lines means there is some data, otherwise nothing.
            num_lines = sum(1 for line in open(file))
            output.write('<h4 class="sub-header" align="left">GO Enrichment: '+name+'</h4>')
            if num_lines > 3:
                df = pd.read_csv(file, sep='\t', skiprows=2)
                df['enrichment'].replace('p', 'under', inplace=True)
                df['enrichment'].replace('e', 'over', inplace=True)
                df2 = df.loc[df['p_fdr'] < args.go_fdr]
                df2.sort_values(by='enrichment', inplace=True)
                if len(df2) > 0:
                    df2.to_csv(base+'.fdr_enriched.csv', index=False)
                    df2['id'] = '<a href="http://amigo.geneontology.org/amigo/search/ontology?q='+ df2['id'].astype(str)+'">'+df2['id']+'</a>'
                    output.write(df2.to_html(escape=False, index=False, classes='table table-hover'))
                else:
                    output.write('<table border="1" class="dataframe table table-hover">\n<th>No enrichment found</th></table>')
            else:
                output.write('<table border="1" class="dataframe table table-hover">\n<th>No enrichment found</th></table>')
    output.write(lib.FOOTER)
        
#################################################### 

##ProteinOrtho################################
if not os.path.isdir(os.path.join(args.out, 'annotations')):
    os.makedirs(os.path.join(args.out, 'annotations'))
lib.log.info("Running orthologous clustering tool, ProteinOrtho5.  This may take awhile...")
#setup protein ortho inputs, some are a bit strange in the sense that they use equals signs
log = os.path.join(protortho, 'proteinortho.log')
#get list of files in folder
filelist = []
for file in os.listdir(protortho):
    if file.endswith('.faa'):
        filelist.append(file)
fileinput = ' '.join(filelist)
cmd = [PROTORTHO, '-project=funannotate', '-synteny', '-cpus='+str(args.cpus), '-singles', '-selfblast']
cmd2 = cmd + filelist
if not os.path.isfile(os.path.join(args.out, 'protortho', 'funannotate.poff')):
    with open(log, 'w') as logfile:
        subprocess.call(cmd2, cwd = protortho, stderr = logfile, stdout = logfile)

#now process the output, get # of singletons per genome, total orthologs, single-copy orthologs and append to stats, output text file with groups
orthologs = os.path.join(args.out, 'annotations','orthology_groups.txt')
with open(orthologs, 'w') as output:
    with open(os.path.join(args.out, 'protortho', 'funannotate.poff'), 'rU') as input:
        count = 0
        scoCount = 0
        for line in input:
            line = line.replace('\n', '') #strip line ending
            if line.startswith('#'):
                header = line
                species = header.split('\t')[3:]
                num_species = header.count('\t') - 2
                continue
            col = re.split(r'[,\t]', line)
            if col[0] != '1':
                count +=1
                ID = 'orth'+str(count)
                prots = col[3:]
                prots = [x for x in prots if x != '*']
                eggs = []
                for i in prots:
                    hit = EGGNOG.get(i)
                    if not hit in eggs:
                        eggs.append(hit)
                eggs = [x for x in eggs if x is not None]
                if col[0] == str(num_species) and col[1] == str(num_species):
                    scoCount += 1
                if len(eggs) > 0:
                    output.write("%s\t%s\t%s\n" % (ID, ', '.join(str(v) for v in eggs), ', '.join(prots)))
                else:
                    output.write("%s\t%s\t%s\n" % (ID, 'None', ', '.join(prots)))

if not os.path.isdir(os.path.join(args.out, 'stats')):
    os.makedirs(os.path.join(args.out, 'stats'))
summary = []
for i in stats:
    singles = lib.singletons(os.path.join(args.out, 'protortho', 'funannotate.poff'), i[0])
    i.append("{0:,}".format(singles))
    orthos = lib.orthologs(os.path.join(args.out, 'protortho', 'funannotate.poff'), i[0])
    i.append("{0:,}".format(orthos))
    i.append("{0:,}".format(scoCount))
    summary.append(i)
#convert to dataframe for easy output
header = ['species', 'isolate', 'Assembly Size', 'Largest Scaffold', 'Average Scaffold', 'Num Scaffolds', 'Scaffold N50', 'Percent GC', 'Num Genes', 'Num Proteins', 'Num tRNA', 'Unique Proteins', 'Prots atleast 1 ortholog', 'Single-copy orthologs']
df = pd.DataFrame(summary, columns=header)
df.set_index('species', inplace=True)
df.transpose().to_csv(os.path.join(args.out, 'stats','genome.stats.summary.csv'))
with open(os.path.join(args.out, 'index.html'), 'w') as output:
    pd.set_option('display.max_colwidth', -1)
    output.write(lib.HEADER)
    output.write(lib.SUMMARY)
    output.write(df.transpose().to_html(classes='table table-condensed'))
    output.write(lib.FOOTER)
############################################
######summarize all annotation for each gene in a table
lib.log.info("Compiling all annotations for each genome")

#get orthology into dictionary
orthoDict = {}
with open(orthologs, 'rU') as input:
    for line in input:
        line = line.replace('\n', '')
        col = line.split('\t')
        genes = col[1].split(',')
        for i in genes:
            orthoDict[i] = col[0]
            
#get GO associations into dictionary as well
with lib.suppress_stdout_stderr():
    goLookup = obo_parser.GODag(os.path.join(currentdir, 'DB', 'go.obo'))
goDict = {}
with open(os.path.join(go_folder, 'associations.txt'), 'rU') as input:
    for line in input:
        line = line.replace('\n', '')
        col = line.split('\t')
        gos = col[1].split(';')
        goList = []
        for i in gos:
            description = i+' '+goLookup[i].name
            goList.append(description)
        goDict[col[0]] = goList

EggNog = lib.eggnog2dict()
iprDict = lib.dictFlipLookup(ipr, INTERPRO)
pfamDict = lib.dictFlipLookup(pfam, PFAM)
meropsDict = lib.dictFlip(merops)  
cazyDict = lib.dictFlip(cazy)

table = []
header = ['GeneID','length','description', 'Ortho Group', 'EggNog', 'Protease family', 'CAZyme family', 'InterPro Domains', 'PFAM Domains', 'GO terms', 'SecMet Cluster', 'SMCOG']
for i in range(0,num_input):
    outputname = os.path.join(args.out, 'annotations', stats[i][0].replace(' ', '_')+'.all.annotations.tsv')
    with open(outputname, 'w') as output:
        output.write("%s\n" % ('\t'.join(header)))
        with open(args.input[i], 'rU') as input:
            SeqRecords = SeqIO.parse(input, 'genbank')
            for record in SeqRecords:
                for f in record.features:
                    if f.type == 'CDS':
                        egg = ''
                        cluster = ''
                        smcog = ''
                        ID = f.qualifiers['locus_tag'][0]
                        length = len(f.qualifiers['translation'][0])
                        description = f.qualifiers['product'][0]
                        if ID in iprDict:
                            IPRdomains = "; ".join(iprDict.get(ID))
                        else:
                            IPRdomains = ''
                        if ID in pfamDict:
                            pfamdomains = "; ".join(pfamDict.get(ID))
                        else:
                            pfamdomains = ''
                        if ID in meropsDict:
                            meropsdomains = "; ".join(meropsDict.get(ID))
                        else:
                            meropsdomains = ''
                        if ID in cazyDict:
                            cazydomains = "; ".join(cazyDict.get(ID))
                        else:
                            cazydomains = ''
                        if ID in goDict:
                            goTerms = "; ".join(goDict.get(ID))
                        else:
                            goTerms = ''
                        if ID in orthoDict:
                            orthogroup = orthoDict.get(ID)
                        else:
                            orthogroup = ''
                        for k,v in f.qualifiers.items():
                            if k == 'note':
                                notes = v[0].split('; ')
                                for i in notes:
                                    if i.startswith('EggNog:'):
                                        hit = i.replace('EggNog:', '')
                                        egg = hit+': '+EggNog.get(hit)
                                    if i.startswith('antiSMASH:'):
                                        cluster = i.replace('antiSMASH:', '')
                                    if i.startswith('SMCOG:'):
                                        smcog = i

                        final_result = [ID, str(length), description, orthogroup, egg, meropsdomains, cazydomains, IPRdomains, pfamdomains, goTerms, cluster, smcog]
                        output.write("%s\n" % ('\t'.join(final_result)))        
############################################
#building remaining HTML output

with open(os.path.join(args.out, 'orthologs.html'), 'w') as output:  
    df = pd.read_csv(orthologs, sep='\t', header=None)
    df.columns = ['Orthology Group', 'EggNog Ref', 'Gene Names']
    pd.set_option('display.max_colwidth', -1)
    output.write(lib.HEADER)
    output.write(lib.ORTHOLOGS)
    output.write(df.to_html(index=False, classes='table table-hover'))
    output.write(lib.FOOTER)
    
with open(os.path.join(args.out, 'citation.html'), 'w') as output:
    output.write(lib.HEADER)
    output.write(lib.CITATION)
    output.write(lib.FOOTER)
                                     
lib.log.info("Finished!")
os._exit(1)


############################################



