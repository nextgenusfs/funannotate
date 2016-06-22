#!/usr/bin/env python

#BUSCO - Benchmarking sets of Universal Single-Copy Orthologs.

#Copyright (C) 2015 E. Zdobnov lab: F. Simao Neto
#<felipe.simao@unige.ch> based on code by R. Waterhouse.

#BUSCO is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#BUSCO is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

########################################################################

# Version 1.2 -- May 3rd 2016 15:37
# several bug fixes and changes June 6th, 2016 via Jon Palmer
#   - sed -i is not valid on all POSIX systems, exchanged for python version
#   - added smooth exit instead of error if type script with no options
#   - OGS is broken due to parsing hmmer results and unbound variable, fixed by setting busco_query to file input instead of parsing from results

#-------------------------------------------------------------------------------#

import os 
import subprocess
import argparse
import time
import threading
import sys
import heapq         

from collections import deque

#check python version
p_version = sys.version_info[0]
if p_version != 3:
    import Queue as queue
else:
    import queue

start_time = time.time()

#------------------------------------ Argument parser START ----------------------------------------#
parser=argparse.ArgumentParser(description = 'Welcome to the Benchmarking set of Universal Single Copy Orthologs (BUSCO).\n\n For further usage information, please check the README file provided with this distrubution.',
                               usage = 'BUSCO_v1.2.py -in [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] [OTHER OPTIONS]')
parser.add_argument('-g','--genome','-in', metavar = 'FASTA FILE' , 
                    help = 'Input file in fasta format.\nCan be a genome, proteome or transcriptome. Default analysis is run on the genome mode, for other files please specify the mode with (-m [MODE])\n')	#Input file in fasta
parser.add_argument('-c', '--cpu', metavar = 'N', help = 'Number of threads/cores to use.')	#Number of threads (CPUs)
parser.add_argument('-a', '--abrev','-o', metavar = 'output', help = 'How to name output and temporary files.')	#output folder
parser.add_argument('--ev', '-e', '-ev', metavar = 'N', type = float, help = 'E-value cutoff for BLAST searches. (Default: 0.01)')	#evalue option
parser.add_argument('-m', '--mode', metavar = 'mode', help = 'which module to run the analysis to run, valid modes are \'all\'(genome assembly), \'OGS\' (gene set / proteome) and \'Trans\' (transcriptome).\n Defaults to \'all\'')
parser.add_argument('-l', '--clade', '--lineage', metavar = 'lineage',  help = 'Which BUSCO lineage to be used.')	#lineage
parser.add_argument('-f', action = 'store_true', default = False, dest = 'force', help = 'Force rewrting of existing files. Must be used when output files with the provided name already exist.')
parser.add_argument('-sp', '--species', default = 'generic', metavar = 'species', help = 'Name of existing Augustus species gene finding metaparameters. (Default: generic)')	
parser.add_argument('-flank', '--flank', '-F', metavar = 'flanks', type = int, help = 'Flanking sequence size for candidate regions. If not provided, flank size is calculated based on genome size with a range from 5 to 20 Kbp.')
parser.add_argument('-Z', '--size', metavar = 'dbsize', type = int, help = 'HMM library total size (Z). Important if using external datasets')
parser.add_argument('-t', '--tmp', metavar = 'temp', default = './',  help = 'Where to store temporary files (Default: Current directory).')	#lineage
parser.add_argument('--limit', '--lim', metavar = 'region_limit', default = 3, type = int, help = 'How many candidate regions to consider (default: 3)')
parser.add_argument('--long', action = 'store_true', default = False, dest = 'long', help = 'Optimization mode Augustus self-training (Default: Off) adds ~20h extra run time, but can improve results for some non-model organisms')
args = vars(parser.parse_args()) #parse arguments 

#print(args) #DBG 

if len(sys.argv) < 2:
    print "Error on input: %s -h" % sys.argv[0]
    sys.exit(1)

#------------------------------------ Argument parser END ----------------------------------------#

#------------------------------------ Set-up the chosen parameters START  -------------------------------#


#Use an e-value cutoff of 0.01 unless user has supplied a custom value using "-ev float" option
ev_cutoff = 0.01	#default e-value cuttof
try:
  if args['ev'] != ev_cutoff and args['ev']!=None:
    print('WARNING: You are using a custom e-value cutoff')
    ev_cutoff = args['ev']
except:
  pass

maxflank = 20000
try:
  if args['clade'] != None:
      clade = args['clade']
      clade_name = clade.strip('/').split('/')[-1].lower()
except:
  print('Please indicate the full path to a BUSCO clade data: Eukaryota, Metazoa, Arthropoda, Vertebrata or Fungi\nExample: -l /path/to/clade')
  raise SystemExit

#Use "generic" as the Augustus species unless user has specified the desired species metaparameters using the "-sp species" option
if args['species'] == 'generic':
    if args['clade'].startswith('arthrop'):
        target_species = 'fly'
    elif args['clade'].startswith('vertebr'):
        target_species = 'human'
    elif args['clade'].startswith('fung'):
        target_species = 'aspergillus_nidulans'
    elif args['clade'].startswith('metazoa'):
        target_species = 'fly' #caenorhabditis
    elif args['clade'].startswith('bacteri'):
        target_species = 'E_coli_K12'
    elif args['clade'].startswith('plant'):
        target_species = 'maize'
    elif args['clade'].startswith('eukary'):
        target_species = 'fly'
    else:
        target_species = 'fly'
else:
     target_species = args['species']
     
#Set up the number of cores to be used
#Augustus uses the python 'threading' library to be run in parallel, blast and HMMer allow this by default
cpus = 1	#1 core default
try:
  if args['cpu'] != None:
    cpus = args['cpu']
except:
  pass


#BUSCO mode (valid modes are genome, transcriptome and ogs)
#Genome is run by default unless overriden by user (-m mode)
modes = ['all','ogs','OGS','transcriptome','trans','ogs','genome'] #valid modes
mode = 'genome'	#unless otherwise specified, run on all (mode for genome assembly)
try:
  if args['mode'] != None and args['mode'] in modes:
    mode = args['mode']
    if mode == 'ogs':
      mode = 'OGS'
    elif mode == 'all' or mode == 'genome':
      mode = 'genome'
    elif mode == 'transcriptome':
      mode = 'trans'
except:
  print('Error: Unknown mode specified * %s *, please check the documentation for valid modes.' % args['mode'])
  raise SystemExit

if mode == 'genome':
    if args['limit'] == 0 or args['limit'] > 20:
        print('ERROR: Limit must be an integer between 1 and 20 (you have used: %s).' % args['limit'])
        raise SystemExit
    else:
        region_limit = args['limit']
        print('Maximum number of regions limited to: %s' % region_limit)

#Get the flank size
#Minimum 5 Kbp and maximum 20 Kbp 
#Scaled as GenomeSize/50 
if mode == 'genome':   #scalled flanks
  f = open(args['genome'])
  size=0
  for line in f:
    if line.startswith('>'):
      pass
    else:
      size += len(line.strip())
  size = size/1000 #size in mb
  flank = int(size/50)	#proportional flank size
  if flank < 5000:
    flank = 5000
  elif flank > maxflank:
    flank = maxflank

#------------------------------------ Set-up the chosen parameters END  -------------------------------#

#------------------------------------ Check dependencies START --------------------#


#Check if command exists and is accessible from the command-line
def cmd_exists(cmd):
    return subprocess.call('type %s' % cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) == 0

#If genome or transcriptome mode, BLAST is required.
#Check if blast is acessible from command-line (tblastn)
if mode in ('genome','trans') and cmd_exists('tblastn') == False:
  print('Error: Blast is not accessible from the command-line, please add it to the environment')
  raise SystemExit

#HMMer 3.1 is always required, check if it is acessible from command-line (as 'hmmsearch')
#Also check if HMMer is the correct version (3.1+)
if mode in ('genome','trans','OGS') and cmd_exists('hmmsearch') == False:
  print('Error: HMMer is not accessible from the command-line, please add it to the environment')
  raise SystemExit
elif cmd_exists('hmmsearch') == True:
    hmmer_check = subprocess.check_output('hmmsearch -h', shell=True)
    hmmer_check = hmmer_check.decode('utf-8')
    hmmer_check = hmmer_check.split('\n')[1].split()[2]
    hmmer_check = float(hmmer_check[:3])
    if hmmer_check >= 3.1:
        pass
    else:
        print('Error: HMMer version detected is not unsupported, please use HMMer 3.1+')
        raise SystemExit

#If genome  mode, Augustus is required.
#Check if Augustus is acessible from command-line (as 'augustus')
if mode == 'genome' and cmd_exists('augustus') == False:
  print('Error: Augustus is not accessible from the command-line, please add it to the environment')
  raise SystemExit

#If genome  mode, Augustus training requires WRITE access to the Augustus directory
#Check if WRITE permissions are enabled for the Augustus installation directory
if mode == 'genome' and os.access(os.environ.get('AUGUSTUS_CONFIG_PATH'),os.W_OK) == False:
  print('Error: Cannot write to Augustus directory, please make sure you have write permissions to %s' % os.environ.get('AUGUSTUS_CONFIG_PATH'))


#------------------------------------ Check dependencies END --------------------#

def pythonSED(input):
    with open(input, 'r+') as f:
        d = f.readlines()
        f.seek(0)
        for i in d[3:]:
            f.write(i)
        f.truncate()
        
def getTrainResults(input):  
    with open(input, 'rU') as train:
        for line in train:
            if line.startswith('nucleotide level'):
                line = line.replace(' ', '')
                values1 = line.split('|') #get [1] and [2]
            if line.startswith('exon level'):
                line = line.replace(' ', '') #get [6] and [7]
                values2 = line.split('|')
            if line.startswith('gene level'):
                line = line.replace(' ', '')
                values3 = line.split('|') #get [6] and [7]
        return (values1[1], values1[2], values2[6], values2[7], values3[6], values3[7])

#------------------------------------- Run START  ---------------------------------------#**********************#


#create the run directory
mainout = './run_%s/' % args['abrev'] #final output directory
if os.path.exists(mainout)==False and args['abrev'] != None:
  subprocess.call(['mkdir', mainout])
else:
  if args['force'] != True:
    print('A run with that name already exists!\nIf are sure you wish to rewrite existing files please use the -f option')
    raise SystemExit
  elif args['force'] == True:
    subprocess.call('rm -rf %s/*' % (mainout), shell = True)

#create the tmp directory
if args['tmp']!='./':
    if os.path.exists(args['tmp'])==False:
        subprocess.call(['mkdir', args['tmp']])
    if args['tmp'][-1]!='/':
        args['tmp']+='/'
print(args['tmp'])


#------------------------------------ Necessary functions START  -------------------------------#

def extract(path,group):
  count = 0
  group_name = group.split('.')[0]
  
  try: 
      group_index = int(group[-1])
  except:
      group_index = '1'
      group = group + '.out.1'
  group_index = str(group_index)
  
  f = open('%saugustus/%s' % (path,group))
  global no_predictions
  written_check = 0
  check = 0;
  while True:
      line = f.readline()
      if not line: 
        break
      if line.startswith('# start gene'):
        line = f.readline(); line = line.split(); places = [line[0], line[3], line[4]]
      elif line.startswith('# protein'):
        line = line.strip().split('[')
        count += 1
        if written_check == 0:
            out = open('%saugustus_proteins/%s.fas.%s' % (path, group_name, group_index), 'w')
            written_check = 1
        out.write('>p%s[%s:%s-%s]\n' % (count, places[0], places[1], places[2]))
        if line[1][-1] == ']':
          line[1] = line[1][:-1]        
        out.write(line[1])
        check = 1
      else:
        if line.startswith(('# end','# sequence')):
          check = 0
          out.write('\n')
        elif check == 1:
          line = line.split()[1]
          if line[-1] == ']':
            line = line[:-1]
          out.write(line)
  if written_check == 1:
    out.close() 
  else:
    no_predictions.append('%s.fas.%s' % (group_name, group_index))  

def checkOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def defineBoundary(a, b):
    temp_start = a[0]; temp_end = a[1]
    current_start = b[0]; current_end = b[1]
    if temp_start < current_start and temp_end < current_start: 
        #i.e. entry is fully before
        #append left, IF entry is the first one; otherwise put into the proper position
        boundary = deque([a, b])
    elif temp_start < current_start and temp_end <= current_end and temp_end >= current_start: 
        #i.e. overlap starts before, but ends inside
        boundary = deque([temp_start, current_end])    
    elif temp_start >= current_start and temp_start <= current_end and temp_end > current_end:
        #overlap starts inside, but ends outside
        boundary = deque([current_start, temp_end])
    elif temp_start > current_end and temp_end > current_end:
        #i.e. query is fully after
        #append right; otherwise put into the proper position
        boundary = deque([b, a])
    elif temp_start >= current_start and temp_end <= current_end and temp_start <= current_end:
        #i.e. query is fully inside, no further operations needed
        boundary = deque(b)
    elif temp_start == current_start and temp_end == current_end:
        boundary = deque(a)
    elif temp_start <= current_start and temp_end >= current_end:
        #i.e. query is longer and contains all coordinates
        #replace by the query
        boundary = deque(a)
    return boundary

          
def gargantua(deck):
  total = 0
  for entry in deck:
    total += entry[1] - entry[0]
  return(total)  

#Compacts numbers into 2 digits.
def shrink (number):
  number = number / totalbuscos
  number = number*100
  if number >= 10:
    if number != 100:
        number = str(number)[:2]
    else:
        number = 100
  elif number < 10 and number > 0:
    number = str(number)[:3]
  return(number)

#------------------------------------ Necessary functions END  -------------------------------#

codons = {'TTT':'F','TTC':'F',
       'TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
       'ATT':'I','ATC':'I','ATA':'I',
       'ATG':'M',
       'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
       'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
       'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
       'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
       'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
       'TAT':'Y','TAC':'Y',
       'TAA':'X','TAG':'X','TGA':'X',
       'CAT':'H','CAC':'H',
       'CAA':'Q','CAG':'Q',
       'AAT':'N','AAC':'N',
       'AAA':'K','AAG':'K',
       'GAT':'D','GAC':'D',
       'GAA':'E','GAG':'E',
       'TGT':'C','TGC':'C',
       'TGG':'W',
       'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
       'AGT':'S','AGC':'S',
       'AGA':'R','AGG':'R',
       'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

comp = {'A':'T','T':'A','C':'G','G':'C','N':'N'}

#gets the sixframe translation
def sixpack(seq):
    s1 = seq
    s2 = seq[1:]
    s3 = seq[2:]
    rev = ''
    for letter in seq[::-1]: 
        rev += comp[letter]
    r1 = rev
    r2 = rev[1:]
    r3 = rev[2:]
    transc = []
    frames = [s1,s2,s3,r1,r3,r2]
    for sequence in frames:
        part = ''
        new = ''
        for letter in sequence:
            if len(part)==3:
                try:
                    new += codons[part]
                except:
                    new += 'NNN'
                part = ''
                part += letter
            else:
                part += letter
        if len(part)==3:
            new += codons[part]
        transc.append(new)
    return(transc)
        

#---------------------------BLAST START -------------------------------------------# <<<<<<<<<<<<<<<<<<<<<<<<<<<

#Make a blast database and run tblastn
if mode == 'genome' or mode == 'blast' or mode == 'trans':
  print('    Phase One')  
  print('*** Running tBlastN ***')
  subprocess.call('makeblastdb -in %s -dbtype nucl -out %s' % 
                  (args['genome'], args['abrev']), shell = True)
  subprocess.call('tblastn -num_threads %s -query %s/ancestral -db %s -out %s_tblastn -outfmt 7' % 
                  (cpus, clade, args['abrev'], args['abrev']), shell = True)


#Get coordinates for a genome analysis
if mode == 'genome' or mode == 'blast':  
  print('*** Getting coordinates for candidate regions! ***')
  
  blast_file = open('%s_tblastn' % args['abrev'])	#open input file
  
  dic = {}; coords = {}
  for line in blast_file:
    if line.startswith('#'):
        pass
    else:
        line = line.strip().split()
        
        busco_name = line[0]; contig = line[1] #busco_og and contig name, respectively 
        busco_start = int(line[6]); busco_end = int(line[7]) #busco hit positions
        contig_start = int(line[8]); contig_end = int(line[9]) #contig postions
        ev = float(line[10]); aln_len = int(line[3]) #e_value and alignment length
        
        #check evalue cutoff 
        if ev <= ev_cutoff:
            if contig_end < contig_start:	#for minus-strand genes, invert coordinates for convenience
                temp = contig_end; contig_end = contig_start; contig_start = temp         
            if busco_name not in dic.keys(): #create new entry in dictionary for current BUSCO
                dic[busco_name] = [contig]; 
                coords[busco_name] = {}; 
                coords[busco_name][contig] = [contig_start,contig_end,deque([[busco_start,busco_end]]),aln_len]  
                
                
            elif contig not in dic[busco_name] and len(dic[busco_name]) < region_limit:	#get just the top3 scoring regions 
                dic[busco_name].append(contig); coords[busco_name][contig] = [contig_start,contig_end,deque([[busco_start,busco_end]]),aln_len]
           
            elif contig in dic[busco_name] and ev < ev_cutoff:	#contigold already checked, now update coordinates
                if contig_start < coords[busco_name][contig][0] and coords[busco_name][contig][0] - contig_start <= 50000:	#starts before, and withing 50kb of current position
                    coords[busco_name][contig][0] = contig_start; coords[busco_name][contig][2].append([busco_start,busco_end]);
                if contig_end > coords[busco_name][contig][1] and contig_end - coords[busco_name][contig][1] <= 50000:	#ends after and within 50 kbs
                    coords[busco_name][contig][1] = contig_end; coords[busco_name][contig][3] = busco_end; coords[busco_name][contig][2].append([busco_start,busco_end]);  
                elif contig_start > coords[busco_name][contig][0] and contig_start < coords[busco_name][contig][1]:#starts inside current coordinates 
                    if contig_end < coords[busco_name][contig][1]:
                        coords[busco_name][contig][2].append([busco_start,busco_end])  #if ending inside, just add alignemnt positions to deque
                    elif contig_end > coords[busco_name][contig][1]: #if ending after current coordinates, extend
                        coords[busco_name][contig][2][1] = contig_end; coords[busco_name][contig][2].append([busco_start,busco_end])             

if mode == 'genome' or mode == 'blast':  
  final_locations = {}
  out = open('coordinates_%s' % args['abrev'],'w')          	#open Coordinates output file
  
  for busco_group in coords:
      final_locations[busco_group] = []
      candidate_contigs = list(coords[busco_group].keys()) #list of candidate contigs 
      size_lists = []
      
      for contig in candidate_contigs:
          potential_locations = coords[busco_group][contig][2]
          max_iterations = len(potential_locations); iter_count = 0 
          
          final_regions = [] #nested list of regions 
          used_pieces = []; non_used = []
          
          while iter_count < max_iterations:
            currently = potential_locations[iter_count]
            if final_regions == []:
                final_regions.append(currently)
            else:
                for region in final_regions:
                    if checkOverlap(currently, region) != 0:
                        gg = defineBoundary(currently, region)
                        region_index = final_regions.index(region)
                        final_regions[region_index] = gg
                        used_pieces.append(iter_count)
                    else:
                        non_used.append(iter_count)
            iter_count += 1
            
          ##done for this contig, now consolidate
          for entry_index in non_used:
              entry = potential_locations[entry_index]
              if entry in used_pieces:
                  pass #already used
              else:
                  ok = []
                  for region in final_regions: 
                      checking = checkOverlap(entry, region)
                      if checking == 0:
                          #i.e. no overlap
                          pass
                      else:
                          ok.append([region,entry])
                  if ok == []:
                      #no overlaps at all (i.e. unique boundary)
                      final_regions.append(entry)
                  else:
                      region = ok[0][0]
                      currently = ok[0][1]
                      gg = defineBoundary(currently, region)
                      final_regions[final_regions.index(region)] = gg
                      
          size_lists.append(gargantua(final_regions))
          
      max_size = max(size_lists)
      size_cutoff = int(0.7 * max_size)
      index_passed_cutoffs = []
      index_passed_cutoffs = heapq.nlargest(region_limit, range(len(size_lists)), size_lists.__getitem__)
      
      for candidate in index_passed_cutoffs:
          if size_lists[candidate] >= size_cutoff: 
              seq_name = candidate_contigs[candidate]
              seq_start = int(coords[busco_group][candidate_contigs[candidate]][0]) - flank
              if seq_start < 0:
                seq_start = 0
              seq_end = int(coords[busco_group][candidate_contigs[candidate]][1]) + flank
              final_locations[busco_group].append([seq_name,seq_start,seq_end])
              out.write('%s\t%s\t%s\t%s\n' % (busco_group, seq_name, seq_start, seq_end))
  out.close()





#Get coordinates, candidate regions and translate sequences (transcriptome analysis)
if mode == 'transcriptome' or mode == 'trans':
  print('*** Getting coordinates for candidate transcripts! ***')
  f = open('%s_tblastn' % args['abrev'])	#open input file
  dic = {}; transdic = {}
  for i in f:				#get a dictionary of BUSCO matches vs candidate scaffolds
    if i.startswith('#'):
      pass
    else:
      line = i.strip().split()
      name = line[0]; scaff = line[1]; e_val = float(line[10]); leng = int(line[3])
      if name not in dic.keys() and e_val <= ev_cutoff:
        dic[name] = [scaff]; maxi = leng
        transdic[scaff] = name
      elif e_val <= ev_cutoff and scaff not in dic[name] and len(dic[name]) < 3 and leng >= 0.7*maxi:
        dic[name].append(scaff); transdic[scaff] = name

  scaff_list = [] #list of unique scaffolds with buscos matches
  for busco in dic:
    for scaff in dic[busco]:
       if scaff not in scaff_list:
         scaff_list.append(scaff)
  print('*** Extracting candidate transcripts! ***')
  f = open(args['genome']);check=0
  for i in f:
    if i.startswith('>'):
      i = i.strip().split(); i = i[0][1:];
      if i in scaff_list:
        out = open('%s%s%s_.temp' % (args['tmp'],i,args['abrev']),'w')
        out.write('>%s\n' % (i))
        check = 1
      else:
        check = 0
    elif check == 1:
      out.write(i)
  out.close()
  if os.path.exists('%stranslated_proteins' % mainout) == False:
    subprocess.call(['mkdir', '%stranslated_proteins' % mainout])
  files = os.listdir('.'); lista = []
  for entry in files:
    if entry.endswith(args['abrev']+'_.temp'):
      lista.append(entry)
   
if mode == 'transcriptome' or mode == 'trans':
  print('Translating candidate transcripts !')    
  for entry in lista:
    raw_seq = open(args['tmp'] + entry)
    trans_seq = open(mainout + 'translated_proteins/' + entry.split(args['abrev'])[0] + '.fas','w')
    
    nucl_seq = ''
    name = ''
    for line in raw_seq:
        if line.startswith('>'):
            name=line.strip()[1:]
        else:
            nucl_seq+=line.strip()
    seq_count = 0
    for translation in sixpack(nucl_seq):
        seq_count += 1
        trans_seq.write('>transcript_%s\n%s\n' % (seq_count, translation))
    trans_seq.close()
    
if mode == 'transcriptome' or mode == 'trans':
  f2 = open('%s/scores_cutoff' % (clade))	#open target scores file
  #Load dictionary of HMM expected scores and full list of groups
  score_dic = {};
  for i in f2:
    i = i.strip().split()
    try:
      score_dic[i[0]] = float(i[1]); 	#float [1] = mean value; [2] = minimum value
    except:
      pass
  totalbuscos = len(list(score_dic.keys())) 


#---------------------------BLAST END -------------------------------------------# <<<<<<<<<<<<<<<<<<<<<<<<<<<

#---------------------------AUGUSTUS steps START -------------------------------------------#
#Run Augustus on all candidate regions
#1- Get the temporary sequence files (no multi-fasta support in Augustus)
#2- Build a list with the running commands (for threading)
#3- Launch Augustus in paralell using Threading
#4- Prepare the sequence files to be analysed with HMMer 3.1


exitFlag=0

#Threading class
class augustusThreads (threading.Thread):
    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
    def run(self):
        process_data(self.name, self.q)

#Threading function
total = 0
state = 0
slate = [100,75,50,25,10]
def process_data(threadName, q):
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            data = q.get()
            queueLock.release()
            subprocess.call('%s' % (data), shell = True)
            check = len([name for name in os.listdir('%saugustus'% mainout) if os.path.isfile(os.path.join('%saugustus'% mainout, name))])
            state = 100 * check / total
            if state > slate[-1]:
                print('=>\t%s%% of predictions performed' % slate.pop())
        else:
            queueLock.release()
        #time.sleep(1)

#Extract candidate contigs/scaffolds from genome assembly 
## Augustus can't handle multi-fasta files, each sequence has to be present in its own file
### Write the temporary sequence files 
if mode == 'genome':
  print('*** pre-Augustus scaffold extraction ***')
  coord = open('coordinates_%s' % args['abrev'])
  dic = {}; scaff_list = []
  for i in coord:
    i = i.strip().split()
    if len(i) != 2:
      dic[i[0]] = [i[1],i[2],i[3]]
      if i[1] not in scaff_list:
        scaff_list.append(i[1])
  f = open(args['genome']); check = 0
  for i in f:
    if i.startswith('>'):
      i = i.split(); i = i[0][1:]
      if i in scaff_list:
        out = open('%s%s%s_.temp' % (args['tmp'],i,args['abrev']),'w')
        out.write('>%s\n' % (i))
        check = 1
      else:
        check = 0
    elif check == 1:
      out.write(i)
  out.close()

#Now run Augustus on each candidate region with its respective Block-profile

  print('*** Running Augustus prediction ***')
  if os.path.exists('%saugustus' % mainout) == False:
    subprocess.call(['mkdir','%saugustus' % mainout])
    

#coordinates of hits by BUSCO
  location_dic = {}
  f = open('coordinates_%s' % args['abrev'])
  for line in f:
        
        line = line.strip().split('\t')
        scaff_id = line[1]
        scaff_start = line[2]
        scaff_end = line[3]
        group_name = line[0]
        
        if group_name not in location_dic:  
            location_dic[group_name] = [ [scaff_id, scaff_start, scaff_end] ] #scaffold,start and end
        elif group_name in location_dic:
            location_dic[group_name].append([scaff_id, scaff_start, scaff_end])
      

  #Make a list containing the commands to be executed in parallel with threading.
  augustus_first_run_strings = []

  for entry in location_dic:
      
      for location in location_dic[entry]:
          
          scaff = location[0] + args['abrev'] + '_.temp'
          scaff_start = location[1]
          scaff_end = location[2]
          output_index = location_dic[entry].index(location)+1
          
          out_name = '%s/augustus/%s.out.%s' % (mainout, entry, output_index)
          
          augustus_call = 'augustus --stopCodonExcludedFromCDS=False --proteinprofile=%(clade)s/prfl/%(busco_group)s.prfl --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' % \
              {'clade' : clade, 'species' : target_species, 'busco_group' : entry, 
               'start_coord' : scaff_start, 'end_coord' : scaff_end, 'scaffold' : scaff, 
               'output' : out_name}
          augustus_first_run_strings.append(augustus_call) #list of call strings
    
  #Create X number of threads    
  threadList = []
  for i in range(int(cpus)):
    threadList.append("Thread-%s" % str(i + 1))
  nameList = list(dic.keys())
  queueLock = threading.Lock()
  workQueue = queue.Queue(len(augustus_first_run_strings))
  threads = []
  threadID = 1
  total = int(len(augustus_first_run_strings))

  #Generate the new threads
  for tName in threadList:
      thread = augustusThreads(threadID, tName, workQueue)
      thread.start()
      threads.append(thread)
      threadID += 1

  #Fill the queue with the Augustus commands
  queueLock.acquire()
  for word in augustus_first_run_strings:
      workQueue.put(word)  
  queueLock.release()

  #Wait for all Augustus jobs to finish (i.e. queue being empty)
  while not workQueue.empty():
      pass
  #Send exit signal
  exitFlag = 1

  #Wait for all threads to finish
  for t in threads:
      t.join()
  exitFlag = 0 #reset the exit flag for next threading step
  print('=>\t100%% of predictions performed')
  slate = [100.0,75.0,50.0,25.0,10.0]


#Preparation of sequences for use with HMMer

#Parse Augustus output files ('run_XXXX/augustus') and extract protein sequences to a FASTA file ('run_XXXX/augustus_proteins').
if mode == 'genome': 
  print('*** Extracting predicted proteins ***')
  output_dir = mainout + 'augustus'
  files = os.listdir(output_dir)
  count = 0; check = 0
  for i in files:
    FileIn = output_dir + '/' + i
    pythonSED(FileIn)
  if os.path.exists(mainout+'augustus_proteins') == False:
    os.makedirs(mainout+'augustus_proteins')

  
  for entry in files:
    f = open(mainout+'augustus/'+entry)
    
    group_name = entry.split('.')[0]
    group_index = entry[-1]
    out = open('%saugustus_proteins/%s.fas.%s' % (mainout, group_name, group_index), 'w')
    
      
    count = 0; tr = 0
    for line in f:
      if line.startswith('# start gene'):
        tr = 1;
      elif tr == 1:
        line = line.split(); places = [line[0], line[3], line[4]]; tr = 0
      elif line.startswith('# protein'):
        line = line.strip().split('[')
        count += 1
        out.write('>g%s[%s:%s-%s]\n' % (count, places[0], places[1], places[2]))
        if line[1][-1] == ']':
          line[1] = line[1][:-1]        
        out.write(line[1])
        check = 1
      else:
        if line.startswith(('# end','# sequence')):
          check = 0
          out.write('\n')
        elif check == 1:
          line = line.split()[1]
          if line[-1] == ']':
            line = line[:-1]
          out.write(line)
  out.close() 

#---------------------------AUGUSTUS steps END -------------------------------------------#

if mode == 'genome':
    subprocess.call('find %saugustus_proteins -size 0 -delete   ' % mainout, shell = True)

#---------------------------HMMER steps START -------------------------------------------#

#Just run HMMer 3.1, slightly different approach for Genome, Transcriptome and Gene Set (OGS)
#Genome mode HMMer
if mode == 'genome':  
  print('*** Running HMMER to confirm orthology of predicted proteins ***')

  files = os.listdir(mainout + 'augustus_proteins/')
  if os.path.exists(mainout + 'hmmer_output') == False:
    subprocess.call(['mkdir', '%shmmer_output' % mainout])
  
  for entry in files:
    if entry.startswith(('BUSCO','EOG','COG')):
        group_name = entry.split('.')[0]
        group_index = entry[-1]
        
        subprocess.call('hmmsearch --domtblout %(folder)s/hmmer_output/%(group_file)s.out.%(index)s -o temp --cpu %(cpu)s %(clade)s/hmms/%(group_file)s.hmm %(folder)s/augustus_proteins/\'%(input_file)s\'' %  {'folder' : mainout, 'index' : group_index, 'input_file' : entry,
      'cpu' : cpus, 'group_file' : group_name, 'clade' : clade}, shell = True)
                

#Transcriptome mode hmmer
if mode == 'trans' or mode == 'transcriptome':  
  print('*** Running HMMER to confirm transcript orthology ***')
  files = os.listdir('%stranslated_proteins/' % mainout)
  if os.path.exists('%shmmer_output' % mainout) == False:
    subprocess.call(['mkdir', '%shmmer_output' % mainout])
  group = ''; grouplist = []
  for i in files:
    if i.endswith('.fas'):
      f = open('%stranslated_proteins/%s' % (mainout,i))
      name = i[:-4]; group = transdic[name]
      if group not in grouplist:
        grouplist.append(group)
        subprocess.call("hmmsearch --domtblout %(output_file)s.out.1 -o temp --cpu %(cpu)s %(group_file)s.hmm \'%(input_file)s\'" % 
	    {'input_file' : mainout + '/translated_proteins/' + i, 'cpu' : cpus, 'group_file' : clade + '/hmms/' + group, 'output_file' : mainout + 'hmmer_output/' + group}, shell = True)
      else:
        grouplist.append(group)
        subprocess.call("hmmsearch --domtblout %(output_file)s.out.%(count)s -o temp --cpu %(cpu)s %(group_file)s.hmm \'%(input_file)s\' " % 
	    {'input_file' : mainout + '/translated_proteins/' + i, 'cpu' : cpus, 'group_file' : clade + '/hmms/' + group, 'output_file' : mainout + 'hmmer_output/' + group, 'count' : str(grouplist.count(group))}, shell = True)

#OGS/Proteome module
if mode == 'OGS':
  if os.path.exists(mainout + 'hmmer_output') == False:
    subprocess.call(['mkdir', '%shmmer_output' % mainout])
  files = os.listdir(clade + '/hmms')
  f2 = open('%s/scores_cutoff' % clade)	#open target scores file
  #Load dictionary of HMM expected scores and full list of groups
  score_dic = {};
  for i in f2:
    i = i.strip().split()
    try:
      score_dic[i[0]] = float(i[1]); 	#[1] = mean value; [2] = minimum value
    except:
      pass
  totalbuscos = len(list(score_dic.keys()))
  for i in files:
    name = i[:-4]
    if name in score_dic:
      subprocess.call('hmmsearch --domtblout %(output_file)s.out.1 -o temp  --cpu %(cpu)s %(group_file)s.hmm \'%(input_file)s\' ' % 
            {'input_file' : args['genome'] , 'cpu' : cpus, 'group_file' : clade + '/hmms/' + name, 'output_file' : mainout + 'hmmer_output/' + name}, shell = True)

#---------------------------HMMER steps END -------------------------------------------#

#executable to test the hmmer parsing functionality


#load scores 
score_file = open('%s/scores_cutoff' % clade)	#open target scores file
cutoff_dictionary = {}
score_dic = {}
for entry in score_file:
    entry = entry.strip().split()
    try:
        score_dic[entry[0]] = float(entry[1]) #name : score
        cutoff_dictionary[entry[0]] = {'score':float(entry[1])}
    except:
        pass
    totalbuscos = len(list(score_dic.keys())) #legacy
    totalbuscos = len(list(cutoff_dictionary.keys()))

#load lengths
leng_dic = {}
sd_dic = {}
f = open('%s/lengths_cutoff' % clade)
for line in f:
    line = line.strip().split()
    
    leng_dic[line[0]] = float(line[3]) #legacy
    sd_dic[line[0]] = float(line[2]) #legacy
    
    cutoff_dictionary[line[0]]['sigma'] = float(line[2])
    #there is an arthropod profile with sigma 0 that causes a crash on divisions
    if float(line[2]) == 0.0:
        cutoff_dictionary[line[0]]['sigma'] = 1
        
    cutoff_dictionary[line[0]]['length'] = float(line[3])
    


def parse_hmmer(hmmer_results_files, mainout, mode, location_dic=False):
    def measuring (nested):
        if isinstance(nested,str):
            return('0')
        scaffolds = list(nested.keys())
        if len(nested) == 1:
            total_len = [0]
            for hit in nested[scaffolds[0]]:
                total_len[0] += hit[1]-hit[0]
        elif len(nested) > 1:
            total_len = [0]*len(nested)
            for entry in range(0, len(scaffolds)):
                for hit in nested[scaffolds[entry]]:
                    total_len[entry] += hit[1]-hit[0]
        try:
            return(total_len)
        except:
            pass
    
    env = []
    fragment = []
    complete_groups = []  #busco groups that are complete (BUSCOxxxx)
    complete_files = [] #output files for groups that are complete (BUSCOxxxx.out.y or BUSCOxxxx.out)
    everything = {} #all info from hit_dic + lengths
    filtered_everything = {}
        
    for file_name in hmmer_results_files:
        group_name = file_name.split('.')[0]
        group_lim_id = int(file_name[-1])-1
        busco_query = group_name
        f = open('%shmmer_output/%s'% (mainout, file_name))
        
        hit_dic = {}
        bit_score_list = []
        
        for line in f:
            if line.startswith('#'):
                pass
            else:
                line = line.strip().split()
                if mode == 'genome':
                    prot_id = line[0] + '-' + file_name
                else:
                    prot_id = line[0];

                tlen = int(line[2]); qlen = int(line[5])
                
                
                bit_score = float(line[7])
                bit_score_list.append(bit_score)
                hmm_start = int(line[15])
                hmm_end = int(line[16])
                
                #new protein that passes score cutoff
                if bit_score >= cutoff_dictionary[busco_query]['score']:
                    if prot_id not in hit_dic.keys():
                        hit_dic[prot_id] = [ [hmm_start, hmm_end, bit_score] ]
                    else:
                        hit_dic[prot_id].append([hmm_start, hmm_end, bit_score])
        f.close()
        
        length = measuring(hit_dic)  

        length_count = 0
        
        if busco_query not in everything:
            everything[busco_query] = hit_dic
        else:
            for part in hit_dic:
                everything[busco_query][part] = hit_dic[part]
        for hit in hit_dic:
            everything[busco_query][hit][0].append(length[length_count])
            length_count += 1
        #classify genes using sigmas
        for entry in everything[busco_query]:

            size = everything[busco_query][entry][0][3]
#            print('BUSCO: %s and Length things are : %s | %s | %s and sigma is: %s' % 
#                  (busco_query, cutoff_dictionary[busco_query]['length'], size, abs(cutoff_dictionary[busco_query]['length'] - size), cutoff_dictionary[busco_query]['sigma']))
            sigma = abs(cutoff_dictionary[busco_query]['length'] - size) / cutoff_dictionary[busco_query]['sigma']
            
            everything[busco_query][entry][0].append(sigma)
            everything[busco_query][entry][0].append(file_name)

        

        
    ##REFINE CLASSIFICATION 
    
    #separate complete into multi and single-copy
    is_complete = {}
    is_fragment = {}
    
    for thing in everything:
        for sequence in everything[thing]: 
                
            all_data = everything[thing][sequence][0]

            sigma = everything[thing][sequence][0][-2]
            seq_name = sequence 
            bit_score = everything[thing][sequence][0][2]
                
            if sigma <= 2:
                if thing not in is_complete:
                    is_complete[thing] = {}
                is_complete[thing][seq_name] = all_data
            elif sigma > 2:
                if thing not in is_fragment:
                    is_fragment[thing] = {}
                is_fragment[thing][seq_name] = all_data

    the_sc = {}
    the_mc = {}
    the_fg = {}
    
    sc_count = 0
    mc_count = 0
    fg_count = 0
    
    has_complete_match = []
    
    for entity in is_complete:
        #single copy
        if len(is_complete[entity]) == 1: #e.g. BUSCOaEOG7QCM97
            the_sc[entity] = is_complete[entity]
            sc_count += 1
            has_complete_match.append(entity)
        elif len(is_complete[entity]) >= 2:
            the_mc[entity] = is_complete[entity]
            mc_count += 1
            has_complete_match.append(entity)
            
    for entity in is_fragment:
        if entity not in has_complete_match:
            fg_count += 1
            the_fg[entity] = is_fragment[entity]
    
    env.append(sc_count)
    env.append(mc_count)
    env.append(fg_count)
    
    out = open('%s/full_table_%s' % (mainout,args['abrev']), 'w')
    
    not_missing = []
    complete_and_singlecopy = []
    csc = {}
    
    
    for entity in the_sc:
        
        for seq_id in the_sc[entity]:
            bit_score = the_sc[entity][seq_id][2]
            seq_len = the_sc[entity][seq_id][3]
            
            not_missing.append(entity)
            complete_and_singlecopy.append(entity)
            
            if mode == 'OGS' or mode == 'trans':
                out.write('%s\tComplete\t%s\t%s\t%s\n' % (entity, seq_id, bit_score, seq_len))
            elif mode == 'genome':
                #predicted_id = seq_id.strip(']').split('[')[0]
                scaff = seq_id.strip(']').split('[')[1].split(':')[0]
                scaff_start = seq_id.strip(']').split('[')[1].split(':')[1].split('-')[0]
                scaff_end = seq_id.strip(']').split('[')[1].split(':')[1].split('-')[1]
                
                out.write('%s\tComplete\t%s\t%s\t%s\t%s\t%s\n' % (entity, scaff, scaff_start, scaff_end, bit_score, seq_len))
                csc[entity] = seq_id
                
    for entity in the_mc:
        for seq_id in the_mc[entity]:
            bit_score = the_mc[entity][seq_id][2]
            seq_len = the_mc[entity][seq_id][3]
            
            not_missing.append(entity)
            
            if mode == 'OGS' or mode == 'trans':
                out.write('%s\tDuplicated\t%s\t%s\t%s\n' % (entity, seq_id, bit_score, seq_len))
            elif mode == 'genome':
                scaff = seq_id.strip(']').split('[')[1].split(':')[0]
                scaff_start = seq_id.strip(']').split('[')[1].split(':')[1].split('-')[0]
                scaff_end = seq_id.strip(']').split('[')[1].split(':')[1].split('-')[1]
                out.write('%s\tDuplicated\t%s\t%s\t%s\t%s\t%s\n' % (entity, scaff, scaff_start, scaff_end, bit_score, seq_len))
            
                
    for entity in the_fg:
        for seq_id in the_fg[entity]:
                bit_score = the_fg[entity][seq_id][2]
                seq_len = the_fg[entity][seq_id][3]
                
                not_missing.append(entity)
                
                if mode == 'OGS' or mode == 'trans':
                    out.write('%s\tFragmented\t%s\t%s\t%s\n' % (entity, seq_id, bit_score, seq_len))
                elif mode == 'genome':
                    scaff = seq_id.strip(']').split('[')[1].split(':')[0]
                    scaff_start = seq_id.strip(']').split('[')[1].split(':')[1].split('-')[0]
                    scaff_end = seq_id.strip(']').split('[')[1].split(':')[1].split('-')[1]
                    out.write('%s\tFragmented\t%s\t%s\t%s\t%s\t%s\n' % (entity, scaff, scaff_start, scaff_end, bit_score, seq_len))
                    
    missing = []
    miss_file = open('%s/missing_buscos_list_%s' % (mainout, args['abrev']), 'w')
    for busco_group in cutoff_dictionary:
        if busco_group in not_missing:
            pass
        else:
            out.write('%s\tMissing\n' % busco_group)
            miss_file.write('%s\n' % busco_group)
            missing.append(busco_group)
            
    env.append(missing)
    env.append(complete_and_singlecopy)
    env.append(csc)
    
    out.close()
    return(env) #DBG



hmmer_results = os.listdir('%shmmer_output' % mainout)
hmmer_results_files = []
for entry in hmmer_results:
    if entry.startswith(('BUSCO','EOG','COG')):
        hmmer_results_files.append(entry)

results_from_hmmer = parse_hmmer(hmmer_results_files, mainout, mode)

single_copy = results_from_hmmer[0] #int
multi_copy = results_from_hmmer[1] #int
only_fragments = results_from_hmmer[2] #int
missing_busco_list = results_from_hmmer[3] #list of BUSCO ids
complete_and_singlecopy = results_from_hmmer[4]
single_copy_files = results_from_hmmer[5]

summary_file = open('%s/short_summary_%s' % (mainout, args['abrev']), 'w')
summary_file.write('#Summarized BUSCO benchmarking for file: %s\n#BUSCO was run in mode: %s\n\n' % (args['genome'], mode))

#summary_file.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % ('CC','DD','FF','MM',totalbuscos))
summary_file.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % 
                   (shrink(single_copy + multi_copy), shrink(multi_copy), shrink(only_fragments), shrink(totalbuscos - single_copy - multi_copy - only_fragments),totalbuscos))  
  
summary_file.write('\t%s\tComplete BUSCOs\n' % str(single_copy + multi_copy) )
summary_file.write('\t%s\tComplete and single-copy BUSCOs\n' % single_copy)
summary_file.write('\t%s\tComplete and duplicated BUSCOs\n' % multi_copy)
summary_file.write('\t%s\tFragmented BUSCOs\n' % only_fragments)
summary_file.write('\t%s\tMissing BUSCOs\n' % str(totalbuscos - single_copy - multi_copy - only_fragments))
summary_file.write('\t%s\tTotal BUSCO groups searched\n' % (totalbuscos))
summary_file.close()

print('Complete: %s'  % str(single_copy + multi_copy) )
print('Single: %s' % single_copy)
print('Multi: %s' % multi_copy)
print('Fragment: %s'  % only_fragments)
print('Missing: %s'  % str(totalbuscos - single_copy - multi_copy - only_fragments))
print('Total BUSCO groups searched: %s' % (totalbuscos))

###############################################

# Retrain Augustus (Genome mode only)
exitFlag = 0 #reset the exitFlag  

#Create the necessary subdirectories
if mode == 'genome':
  if os.path.exists('%sgffs'% mainout) == False:
    subprocess.call(['mkdir','%sgffs' % mainout])
  if os.path.exists('%sgb' % mainout) == False:
    subprocess.call(['mkdir','%sgb' % mainout])
    
  print('    Phase Two')    
  print('*** Training Augustus using Single-Copy Complete BUSCOs ***')  
  try:
    AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
  except KeyError:
    print("$AUGUSTUS_CONFIG_PATH environmental variable not found, Augustus is not properly configured.")
    os._exit(1)
  if AUGUSTUS.endswith('config'):
    AUGUSTUS_BASE = AUGUSTUS.replace('config', '')
  elif AUGUSTUS.endswith('config'+os.sep):
    AUGUSTUS_BASE = AUGUSTUS.replace('config'+os.sep, '')  
  species = '--species='+args['abrev']
  aug_cpus = '--cpus='+str(cpus)
  trainingset = 'training_set_'+args['abrev']
  
  for entry in single_copy_files:
    check = 0
    file_name = single_copy_files[entry].split('-')[-1]
    target_seq_name = single_copy_files[entry].split('[')[0]
    group_name = file_name.split('.')[0]
    
    #create GFFs with only the "single-copy" BUSCO sequences
    pred_file = open('%s/augustus/%s' % (mainout, file_name))
    gff_file = open('%s/gffs/%s.gff' % (mainout, group_name), 'w')
    for line in pred_file:
        if line.startswith('# start gene'):
            pred_name = line.strip().split()[-1]
            if pred_name == target_seq_name:
                check = 1
        elif line.startswith('#'):
            check = 0
        elif check == 1:
            gff_file.write(line)
    gff_file.close()
        
    subprocess.call('$AUGUSTUS_CONFIG_PATH/../scripts/gff2gbSmallDNA.pl %s/gffs/%s.gff %s 1000 %s/gb/%s.raw.gb 1>/dev/null 2>/dev/null' %
                    (mainout, entry, args['genome'], mainout, entry), shell = True)
  #bacteria clade needs to be flagged as "prokaryotic"
  NEW_SPECIES = os.path.join(AUGUSTUS_BASE, 'scripts', 'new_species.pl')
  if clade.startswith('bacteria'):
    subprocess.call([NEW_SPECIES, '--prokaryotic', species])
  else:
    subprocess.call([NEW_SPECIES, species])
  #create training set from BUSCO models, concatentate each model, make sure it is not repeated
  locus_seen = {}
  keepers = []
  for file in os.listdir(os.path.join('run_'+args['abrev'], 'gb')):
    file = os.path.join('run_'+args['abrev'], 'gb', file)
    with open(file, 'rU') as input:
        line = input.readline().strip()
        if not line in locus_seen:
            keepers.append(file)
            locus_seen[line] = ''
  with open(trainingset, 'w') as output:
      for i in keepers:
        with open(i) as input2:
            output.write(input2.read())
  #subprocess.call('cat %sgb/*.gb > training_set_%s' % (mainout, args['abrev']), shell = True)
  subprocess.call(['etraining', species, trainingset])
  
  '''
  #split into training and test set
  RANDOMSPLIT = os.path.join(AUGUSTUS_BASE, 'scripts', 'randomSplit.pl')
  subprocess.call([RANDOMSPLIT, trainingset, '100'])
  #run initial training
  subprocess.call(['etraining', species, trainingset+'.train'])
  
  print ("----------------------------------------------------------")
  print ("***  Initial training set created, now measuring accuracy")
  #check accuracy of original training
  with open('initial_training_results.txt', 'w') as initialtraining:
    subprocess.call(['augustus', species, trainingset+'.test'], stdout=initialtraining)
  train_results = getTrainResults('initial_training_results.txt')
  print ("%f%% genes predicted exactly\n%f%% of exons predicted exactly\n%f%% of the predicted exons were exact in test set" % (float(train_results[4])*100, float(train_results[2])*100, float(train_results[3])*100))  
  print ("----------------------------------------------------------")         
  '''
  #long mode (--long) option - runs all the Augustus optimization scripts (adds ~1 day of runtime)
  if args['long'] == True:
    print('***  Optimizing augustus metaparameters, this may take around 20 hours')
    OPTIMIZE = os.path.join(AUGUSTUS_BASE, 'scripts', 'optimize_augustus.pl')
    with open('optimize_augustus.log', 'w') as logfile:
        subprocess.call([OPTIMIZE, species, aug_cpus, '--UTR=off', trainingset+'.train'], stderr = logfile, stdout = logfile)
    subprocess.call(['etraining', species, trainingset+'.train'])
    with open('optimized_training_results.txt', 'w') as finaltraining:
        subprocess.call(['augustus', species, trainingset+'.test'], stdout=finaltraining)
    train_results = getTrainResults('optimized_training_results.txt')
    print ("%f%% genes predicted exactly\n%f%% of exons predicted exactly\n%f%% of the predicted exons were exact in test set" % (float(train_results[4])*100, float(train_results[2])*100, float(train_results[3])*100))
    
  print('*** Re-running Augustus with the new metaparameters, number targets: %s ***' % len(missing_busco_list))
  
  augustus_rerun_strings = []
  augustus_rerun_seds = []
  hmmer_rerun_strings = []

  for entry in missing_busco_list:
      if entry in location_dic:
        for location in location_dic[entry]:
          
            scaff = location[0] + args['abrev'] + '_.temp'
            scaff_start = location[1]
            scaff_end = location[2]
            output_index = location_dic[entry].index(location)+1
          
            out_name = '%s/augustus/%s.out.%s' % (mainout, entry, output_index)
          
          
            augustus_call = 'augustus --proteinprofile=%(clade)s/prfl/%(busco_group)s.prfl --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' % \
              {'clade' : clade, 'species' : target_species, 'busco_group' : entry, 'start_coord' : scaff_start, 'end_coord' : scaff_end, 'scaffold' : scaff, 
               'output' : out_name}
            augustus_rerun_strings.append(augustus_call) #list of call strings

            sed_call = pythonSED(out_name)
            augustus_rerun_seds.append(sed_call)

            out_name = '%s/hmmer_output/%s.out.%s' % (mainout, entry, output_index)
            augustus_fasta = '%s/augustus_proteins/%s.fas.%s' % (mainout, entry, output_index)
            hmmer_call = 'hmmsearch --domtblout %(output_file)s -o temp --cpu 1 %(clade)s/hmms/%(busco_group)s.hmm \'%(input_file)s\'' % \
                         {'output_file' : out_name, 'busco_group' : entry, 'input_file' : augustus_fasta, 'clade' : clade}
            hmmer_rerun_strings.append(hmmer_call)    
      else:
         pass
          
##########

#pass the strings constructed above to the threading function
if mode == 'genome': 
  # augustus threading block  
  queueLock = threading.Lock()
  workQueue = queue.Queue(len(augustus_rerun_strings))
  threads = []
  threadID = 1
  
  needed = len(augustus_rerun_strings); mark = 0

  # Create new threads
  for tName in threadList:
      mark += 1
      thread = augustusThreads(threadID, tName, workQueue)
      thread.start()
      threads.append(thread)
      threadID += 1
      if mark >= needed:
          break

  # Fill the queue
  queueLock.acquire()
  for word in augustus_rerun_strings:
      workQueue.put(word)
  queueLock.release()

  # Wait for queue to empty
  while not workQueue.empty():
      pass
  # Notify threads it's time to exit
  exitFlag = 1

# Wait for all threads to complete
  for t in threads:
      t.join()

# augustus is done
# sed the results (remove weird nordic characters in augustus output)
  print('=>\t100%% of predictions performed')
  for sed_string in augustus_rerun_seds:
      subprocess.call('%s' % sed_string, shell=True)
  
# Extract fasta files from augustus output
  no_predictions = []
  for entry in missing_busco_list:      
     if entry in location_dic: 
        for location in location_dic[entry]:
          output_index = location_dic[entry].index(location)+1
          out_name = '%s/augustus/%s.out.%s' % (mainout, entry, output_index)
          plain_name = '%s.out.%s' % (entry, output_index) #when extract gets reworked to not need MAINOUT, change to OUT_NAME
          extract(mainout, plain_name)
     else:
         pass
  
# Run hmmer
  print('*** Running HMMER on the new predictions ***')
  exitFlag = 0
  
  queueLock = threading.Lock()
  workQueue = queue.Queue(len(hmmer_rerun_strings))
  threads = []
  threadID = 1
  
  needed = len(hmmer_rerun_strings); mark = 0

  # Create new threads
  for tName in threadList:
      mark += 1
      thread = augustusThreads(threadID, tName, workQueue)
      thread.start()
      threads.append(thread)
      threadID += 1
      if mark >= needed:
          break
  
  
  queueLock.acquire()
  for word in hmmer_rerun_strings:
      target_seq = word.split('/')[-1][:-1]
      if target_seq not in no_predictions:
        workQueue.put(word)
  queueLock.release()

  # Wait for queue to empty
  while not workQueue.empty():
      pass
  # Notify threads it's time to exit
  exitFlag = 1

# Wait for all threads to complete
  for t in threads:
      t.join()  

##########
if mode == 'genome':  #can remove this IF statement
  hmmer_results = os.listdir('%shmmer_output' % mainout)
  hmmer_results_files = []
  for entry in hmmer_results:
    if entry.startswith(('BUSCO','EOG','COG')):
        hmmer_results_files.append(entry)

  results_from_hmmer = parse_hmmer(hmmer_results_files, mainout, mode)

  single_copy = results_from_hmmer[0] #int
  multi_copy = results_from_hmmer[1] #int
  only_fragments = results_from_hmmer[2] #int
  missing_busco_list = results_from_hmmer[3] #list of BUSCO ids
  complete_and_singlecopy = results_from_hmmer[4]
  single_copy_files = results_from_hmmer[5]

  summary_file = open('%s/short_summary_%s' % (mainout, args['abrev']), 'w')
#    summary_file = open('%s/short_summary_%s' % (mainout, args['abrev']), 'w')

  summary_file.write('#Summarized BUSCO benchmarking for file: %s\n#BUSCO was run in mode: %s\n\n' % (args['genome'], mode))

  summary_file.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % 
                   (shrink(single_copy + multi_copy), shrink(multi_copy), shrink(only_fragments), shrink(totalbuscos - single_copy - multi_copy - only_fragments),totalbuscos))  
  
  summary_file.write('\t%s\tComplete BUSCOs\n' % str(single_copy + multi_copy) )
  summary_file.write('\t%s\tComplete and single-copy BUSCOs\n' % single_copy)
  summary_file.write('\t%s\tComplete and duplicated BUSCOs\n' % multi_copy)
  summary_file.write('\t%s\tFragmented BUSCOs\n' % only_fragments)
  summary_file.write('\t%s\tMissing BUSCOs\n' % str(totalbuscos - single_copy - multi_copy - only_fragments))
  summary_file.write('\t%s\tTotal BUSCO groups searched\n' % (totalbuscos))
  summary_file.close()

  print('Complete: %s'  % str(single_copy + multi_copy) )
  print('Single: %s' % single_copy)
  print('Multi: %s' % multi_copy)
  print('Fragment: %s'  % only_fragments)
  print('Missing: %s'  % str(totalbuscos - single_copy - multi_copy - only_fragments))
  print('Total BUSCO groups searched: %s' % (totalbuscos))

  ##get single-copy files as fasta
  if os.path.exists('%ssingle_copy' % mainout) == False:
    subprocess.call(['mkdir','%ssingle_copy' % mainout])
   
  #print('Getting single-copy files...') 
  for entry in single_copy_files:
    check = 0
    
    
    file_name = single_copy_files[entry].split('-')[-1].replace('out','fas')
    target_seq_name = single_copy_files[entry].split('[')[0]
    group_name = file_name.split('.')[0]
    seq_coord_start = single_copy_files[entry].split(']-')[0].split('[')[1]

    
    #create GFFs with only the "single-copy" BUSCO sequences
    pred_fasta_file = open('%s/augustus_proteins/%s' % (mainout, file_name))
    single_copy_outfile = open('%s/single_copy/%s.fas' % (mainout, group_name), 'w')
    
    for line in pred_fasta_file:
        if line.startswith('>%s' % target_seq_name):
            single_copy_outfile.write('>%s:%s:%s\n' % (group_name, args['genome'], seq_coord_start))
            check = 1
        elif line.startswith('>'):
            check = 0
        elif check == 1:
            single_copy_outfile.write(line)
    single_copy_outfile.close()

###############################################


#--------------------------- Cleaning up temporary files START -------------------------------------------#

#Clean up temporary files
if mode != 'OGS':
  subprocess.call('rm %s*%s_.temp' % (args['tmp'],args['abrev']), shell = True)
  subprocess.call('rm %s.nsq %s.nin %s.nhr'  % (args['abrev'], args['abrev'], args['abrev']), shell = True)
  subprocess.call('mv %s_tblastn run_%s' % (args['abrev'], args['abrev']), shell = True)
  if mode == 'trans':
    print("Total running time: %f seconds" % (time.time() - start_time))
  if mode != 'trans':
    subprocess.call('mv training_set_%s run_%s' % (args['abrev'], args['abrev']), shell = True)
    subprocess.call('mv coordinates_%s run_%s' % (args['abrev'], args['abrev']), shell = True)
else:
  print("Total running time: %f seconds" % (time.time() - start_time))





