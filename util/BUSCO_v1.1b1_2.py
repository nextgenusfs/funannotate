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


#Version 1.1c (may/15) - minor bug fixes

#                      - lineages may be specified using the full path (e.g. -l /path/to/lineage)
#                      - added threading support (signficant speed increase)
#                      - added checks for the necessary programs before running BUSCO 


#script modified by Jon Palmer 2016
#   - sed commands were causing error on OSX - I think they are fixed
#   - tried to fix the requirement of python3 - seems like nothing too 'fancy' is being used here, why won't 2.7 work?
#   - switch to mulitprocessing and not multi-threading - multi-threading does not use multiple CPUs

#-------------------------------------------------------------------------------#

import os 
import argparse
from collections import deque
import time
import multiprocessing
import subprocess

start_time = time.time()



#------------------------------------ Argument parser START ----------------------------------------#
parser=argparse.ArgumentParser(description='Welcome to the Benchmarking set of Universal Single Copy Orthologs (BUSCO).\n\n For further usage information, please check the README file provided with this distrubution.',
                               usage='BUSCO_v1.0.py -in [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] [OTHER OPTIONS]')
parser.add_argument('-g','--genome','-in',metavar='FASTA FILE',type=str,help='Input file in fasta format.\nCan be a genome, proteome or transcriptome. Default analysis is run on the genome mode, for other files please specify the mode with (-m [MODE])\n')	#genome assembly file
parser.add_argument('-c','--cpu',metavar='N',type=str,help='Number of threads/cores to use.')	#Number of available threads
parser.add_argument('-a','--abrev','-o',metavar='output',type=str,help='How to name output and temporary files.')	#Four letter abbreviation for use with genome assembly
parser.add_argument('--ev','-e','-ev',metavar='N', type=float, help='E-value cutoff for BLAST searches. (Default: 0.01)')	#evalue option
parser.add_argument('-m','--mode',metavar='mode', type=str, help='which module to run the analysis to run, valid modes are \'all\'(genome assembly), \'OGS\' (gene set / proteome) and \'Trans\' (transcriptome).\n Defaults to \'all\'')
parser.add_argument('-l','--clade','--lineage',metavar='lineage', type=str, help='Which BUSCO lineage to be used.')	#lineage
parser.add_argument('-f', action='store_true',default=False, dest='force',help='Force rewrting of existing files. Must be used when output files with the provided name already exist.')
parser.add_argument('-sp','--species',default='generic',metavar='species',type=str,help='Name of existing Augustus species gene finding metaparameters. (Default: generic)')	
parser.add_argument('-flank','--flank','-F',metavar='flanks', type=int, help='Flanking sequence size for candidate regions. If not provided, flank size is calculated based on genome size with a range from 5 to 20 Kbp.')
parser.add_argument('-Z','--size',metavar='dbsize', type=int, help='HMM library total size (Z). Important if using external datasets')
parser.add_argument('--long', action='store_true',default=False, dest='long',help='Optimization mode Augustus self-training (Default: Off) adds ~20h extra run time, but can improve results for some non-model organisms')


args=vars(parser.parse_args()) #parse the arguments
#print(args)
mainout='./run_%s/' % args['abrev'] #final output directory
if os.path.exists(mainout)==False and args['abrev']!=None:
  os.system('mkdir %s' % mainout)
else:
  if args['force']!=True:
    print('A run with that name already exists!\nIf are sure you wish to rewrite existing files please use the -f option')
    raise SystemExit

def CheckAugustusSpecies(input):
    #get the possible species from augustus
    augustus_list = []
    for i in os.listdir(os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], 'species')):
        if not i.startswith('.'):
            augustus_list.append(i)
    augustus_list = set(augustus_list)
    if input in augustus_list:
        return True
    else:
        return False

target_species='generic' 
if CheckAugustusSpecies(args['species']):
    target_species = args['species']

print(target_species)
  
 
ev_cut=0.01	#default e-value cuttof
try:
  if args['ev']!=ev_cut and args['ev']!=None:
    print('WARNING: You are using a custom e-value cutoff')
    ev_cut=args['ev']
except:
  pass



valid_clade_info={'arthropoda':102785,'metazoa':91897,'vertebrata':143785,'fungi':174195,'example':102785,'bacteria':107114,'eukaryota':41317}
maxflank=20000
#print(args['clade'])
try:
  if args['clade']!=None:
      clade=args['clade']
      clade_name=clade.strip('/').split('/')[-1].lower()
      print(clade_name)
      if clade_name in valid_clade_info:
          Z=valid_clade_info[clade_name]
      else:
        print('Using custom lineage data...')
        try:
          Z=args['dbsize'];  
        except:
          print('Please indicate the size of the custom HMM database using the (-Z integer)')
          raise SystemExit

except:
  print('Please indicate the full path to a BUSCO clade: Eukaryota, Metazoa, Arthropoda, Vertebrata or Fungi\nExample: -l /path/to/Arthropoda')
  raise SystemExit



cpus=1	#1 core default
try:
  if args['cpu']!=cpus and args['cpu']!=None:
    cpus=args['cpu']
except:
  pass

modes=['all','blast','hmmer','augustus','parser','hmmer+','OGS','transcriptome','trans','ogs','genome'] #valid modes
mode='genome'	#unless otherwise specified, run on all (mode for genome assembly)
try:
  if args['mode']!=None and args['mode'] in modes:
    mode=args['mode']
    if mode=='ogs':
      mode='OGS'
    elif mode=='all' or mode=='genome':
      mode='genome'
    elif mode=='transcriptome':
      mode='trans'
except:
  print('Error: Unknown mode specified * %s *, please check the documentation for valid modes.' % args['mode'])
  raise SystemExit

#------------------------------------ Check if necessary programs are acessible --------------------#
import subprocess
def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

if mode in ('genome','trans') and cmd_exists('tblastn')==False:
  print('Error: Blast is not accessible from the command-line, please add it to the environment')
  raise SystemExit


if mode in ('genome','trans','OGS') and cmd_exists('hmmsearch')==False:
  print('Error: HMMer is not accessible from the command-line, please add it to the environment')
  raise SystemExit
elif cmd_exists('hmmsearch')==True:
    hmmer_check=subprocess.check_output('hmmsearch -h', shell=True)
    hmmer_check=hmmer_check.decode('utf-8')
    hmmer_check=hmmer_check.split('\n')[1].split()[2]
    hmmer_check=float(hmmer_check[:3])
    if hmmer_check>=3.1:
        pass
    else:
        print('Error: HMMer version detected is not unsupported, please use HMMer 3.1+')
        raise SystemExit
        
if mode=='genome' and cmd_exists('augustus')==False:
  print('Error: Augustus is not accessible from the command-line, please add it to the environment')
  raise SystemExit

if mode=='genome' and os.access(os.environ.get('AUGUSTUS_CONFIG_PATH'),os.W_OK)==False:
  print('Error: Cannot write to Augustus directory, please make sure you have write permissions to %s' % os.environ.get('AUGUSTUS_CONFIG_PATH'))
    
if mode=='trans' and cmd_exists('transeq')==False:
  print('Error: EMBOSS transeq is not accessible from the commandline, please add it to the environment')
  raise SystemExit


#------------------------------------------------checks over ---------------------------------------#

#get clade specific parameters.
if mode=='genome' or mode=='blast':   #scalled flanks
  f=open(args['genome'])
  size=0
  for line in f:
    if line.startswith('>'):
      pass
    else:
      size+=len(line.strip())
  size=size/1000 #size in mb
  flank=int(size/50)	#proportional flank size
  if flank<5000:
    flank=5000
  elif flank>maxflank:
    flank=maxflank
 
#------------------------------------ Argument parser END ----------------------------------------#

def measuring (nested):
  if isinstance(nested,str):
    return('0')
  scaffolds=list(nested.keys())
  if len(nested)==1:
    total_len=[0]
    for hit in nested[scaffolds[0]]:
      total_len[0]+=hit[1]-hit[0]
  elif len(nested)>1:
    total_len=[0]*len(nested)
    for entry in range(0,len(scaffolds)):
      for hit in nested[scaffolds[entry]]:
       total_len[entry]+=hit[1]-hit[0]
  try:
    return(total_len)
  except:
    pass

def extract(path,group):
  count=0
  if group.endswith(('.1','.2','.3')):
    f=open('%saugustus/%s' % (path,group))
    out=open('%saugustus_proteins/%s.fas.%s' % (path,group[:-6],group[-1]),'w')
  else:
    f=open('%saugustus/%s.out' % (path,group))
    out=open('%saugustus_proteins/%s.fas' % (path,group),'w')
  check=0;
  while True:
      line = f.readline()
      if not line: 
        break
      if line.startswith('# start gene'):
        line=f.readline();line=line.split();places=[line[0],line[3],line[4]]
      elif line.startswith('# protein'):
        line=line.strip().split('[')
        count+=1
        out.write('>p%s[%s:%s-%s]\n' % (count,places[0],places[1],places[2]))
        if line[1][-1]==']':
          line[1]=line[1][:-1]        
        out.write(line[1])
        check=1
      else:
        if line.startswith('# end'):
          check=0
          out.write('\n')
        elif check==1:
          line=line.split()[1]
          if line[-1]==']':
            line=line[:-1]
          out.write(line)
  out.close() 

def disentangle(deck):
  structure=deque([deck.popleft()])
  try:
    while 1:
      temp=deck.popleft()
      start=temp[0];end=temp[1]
      for i in range(0,len(structure)):
        ds=structure[i][0];de=structure[i][1]
        if start<ds and end<ds: #fully before
          if i==0:  #first entry, just appendleft
            structure.appendleft(temp)
            break
          else:
            new=structure[0:i];new.append(temp)
            for z in range(i,len(structure)):
              new.append(structure[z])          
            break
        elif start<ds and end<de and end>ds: #end overlaps inside, but the start is before
          structure[i][0]=start
          break
        elif start>ds and start<de and end>de: #start overlaps inside, but the end is after
          structure[i][1]=end
          break
        elif start>de and end>de: #fully after
          if i==len(structure)-1: #only if its the last entry can it be safely added to structure
            structure.append(temp)
        elif start<ds and end>de: #current structure is found fully inside the current entry
         structure[i]=temp
  except:
    return(structure)
          
def gargantua(deck):
  total=0
  for entry in deck:
    total+=entry[1]-entry[0]
  return(total)  

def shrink (number):
  number=number*100
  if number>=10:
    number=str(number)[:2]
  elif number<10 and number>0:
    number=str(number)[:3]
  return(number)


def runAugustus(input):
    FNULL = open(os.devnull, 'w')
    protprofile = '--proteinprofile=' + input[3] + '/' + input[0]
    startcoords = '--predictionStart=' + input[1]
    endcoords = '--predictionEnd=' + input[2]
    Species = '--species=' + input[4]
    Scaffold = input[5]
    with open(input[6], 'w') as output:
        subprocess.call(['augustus', protprofile, startcoords, endcoords, Species, Scaffold], stdout = output, stderr = FNULL)

def runHMMerSearch(input):
    FNULL = open(os.devnull, 'w')
    HMM = input[2] + '.hmm'
    subprocess.call(['hmmsearch', '--domtblout', input[3], '-Z', input[1], '-o', 'temp', '--cpu', '1', HMM, input[0]], stdout = FNULL, stderr = FNULL)

def runSEDS(input):
    subprocess.call(['sed', '-i', '"1,3d"', input])

#---------------------------BLAST steps START -------------------------------------------# 

#Make a blast database and run tblastn
if mode=='genome' or mode=='blast' or mode=='trans':
  print('*** Running tBlastN ***')
  os.system('makeblastdb -in %s -dbtype nucl -out %s' % (args['genome'],args['abrev']))
  os.system('tblastn -num_threads %s -query %s/ancestral -db %s -out %s_tblastn -outfmt 7' % 	  (cpus,clade,args['abrev'],args['abrev']))

#Get coordinates for a genome analysis
if mode=='genome' or mode=='blast':  
  print('*** Getting coordinates for candidate regions! ***')
  f=open('%s_tblastn' % args['abrev'])	#open input file
  out=open('coordinates_%s' % args['abrev'],'w')          	#open Coordinates output file
  dic={};coords={}
  for i in f:
    if i.startswith('#'):
      pass
    else:
      line=i.strip().split()
      name=line[0];scaff=line[1];hitstart=int(line[6]);hitend=int(line[7]);postart=int(line[8]);posend=int(line[9])
      e_val=float(line[10]);sizer=int(line[3])
      if posend<postart:	#for minus-strand genes, invert coordinates for convenience
        temp=posend;posend=postart;postart=temp         
      if name not in dic.keys(): #create new entry in dictionary for current BUSCO
        dic[name]=[scaff];coords[name]={};coords[name][scaff]=[postart,posend,deque([[hitstart,hitend]]),sizer]  
      elif scaff not in dic[name] and len(dic[name])<3:	#get just the top3 scoring regions
        dic[name].append(scaff);coords[name][scaff]=[postart,posend,deque([[hitstart,hitend]]),sizer]
      elif scaff in dic[name] and e_val<ev_cut:	#scaffold already checked, now update coordinates
        if postart<coords[name][scaff][0] and coords[name][scaff][0]-postart<=50000:	#starts before, and withing 50kb of current position
            coords[name][scaff][0]=postart;coords[name][scaff][2].append([hitstart,hitend]);
        if posend>coords[name][scaff][1] and posend-coords[name][scaff][1]<=50000:	#ends after and within 50 kbs
          coords[name][scaff][1]=posend;coords[name][scaff][3]=hitend;coords[name][scaff][2].append([hitstart,hitend]);
        elif postart>coords[name][scaff][0] and postart<coords[name][scaff][1]:#starts inside current coordinates 
          if posend<coords[name][scaff][1]:
            coords[name][scaff][2].append([hitstart,hitend])  #if ending inside, just add alignemnt positions to deque
          elif posend>coords[name][scaff][1]: #if ending after current coordinates, extend
            coords[name][scaff][2][1]=posend;coords[name][scaff][2].append([hitstart,hitend]) 
  for i in coords:
    contest={};maxi=0
    for contig in coords[i]:
      sizer=disentangle(coords[i][contig][2])
      if contig not in contest:
        contest[contig]=0
      size=gargantua(sizer);
      contest[contig]=size
      if size>maxi:
        maxi=size
    for contig in contest:
      #if contest[contig]>0.7*maxi:
      out.write('%s\t%s\t%s\t%s\n' % (i,contig,max(0,coords[i][contig][0]-flank),coords[i][contig][1]+flank))
  out.close()
 
#Get coordinates, candidate regions and translate sequences (transcriptome analysis)
if mode=='transcriptome' or mode=='trans':
  print('*** Getting coordinates for candidate transcripts! ***')
  f=open('%s_tblastn' % args['abrev'])	#open input file
  dic={};transdic={}
  for i in f:				#get a dictionary of BUSCO matches vs candidate scaffolds
    if i.startswith('#'):
      pass
    else:
      line=i.strip().split()
      name=line[0];scaff=line[1];e_val=float(line[10]);leng=int(line[3])
      if name not in dic.keys() and e_val<=ev_cut:
        dic[name]=[scaff];maxi=leng
        transdic[scaff]=name
      elif e_val<=ev_cut and scaff not in dic[name] and len(dic[name])<3 and leng>=0.7*maxi:
        dic[name].append(scaff);transdic[scaff]=name

  scaff_list=[] #list of unique scaffolds with buscos matches
  for busco in dic:
    for scaff in dic[busco]:
       if scaff not in scaff_list:
         scaff_list.append(scaff)
  print('*** Extracting candidate transcripts! ***')
  f=open(args['genome']);check=0
  for i in f:
    if i.startswith('>'):
      i=i.strip().split();i=i[0][1:];
      if i in scaff_list:
        out=open('%s%s_.temp' % (i,args['abrev']),'w')
        out.write('>%s\n' % (i))
        check=1
      else:
        check=0
    elif check==1:
      out.write(i)
  out.close()
  if os.path.exists('%stranslated_proteins' % mainout)==False:
    os.system('mkdir %stranslated_proteins' % mainout)
  files=os.listdir('.');lista=[]
  for entry in files:
    if entry.endswith(args['abrev']+'_.temp'):
      lista.append(entry)
      
  print('Translating candidate transcripts !')    
  for entry in lista:
    #print(entry);a=input('press to continue')
    os.system('transeq -clean -frame 6 -trim -sequence %(scaffold)s -outseq %(translated_scaffold)s.fas' % {'scaffold':entry,'translated_scaffold':mainout+'translated_proteins/'+entry.split(args['abrev'])[0]+'_ts'})
  f2=open('%s/scores_cutoff' % (clade))	#open target scores file
  #Load dictionary of HMM expected scores and full list of groups
  score_dic={};
  for i in f2:
    i=i.strip().split()
    try:
      score_dic[i[0]]=float(i[1]); 	#float [1] = mean value; [2] = minimum value
    except:
      pass
  totalbuscos=len(list(score_dic.keys()))
  
#---------------------------AUGUSTUS steps START -------------------------------------------#

exitFlag=0

#Step-3 
#Extract candidate contigs/scaffolds from genome assembly (necessary because augustus doesn't handle multi-fasta files when running on a specific target region)
if mode=='genome' or mode=='augustus':
  #target_species=species_list[0]
  print('*** pre-Augustus scaffold extraction ***')
  coord=open('coordinates_%s' % args['abrev'])
  dic={};scaff_list=[]
  for i in coord:
    i=i.strip().split()
    if len(i)!=2:
      dic[i[0]]=[i[1],i[2],i[3]]
      if i[1] not in scaff_list:
        scaff_list.append(i[1])
  f=open(args['genome']);check=0
  for i in f:
    if i.startswith('>'):
      i=i.split();i=i[0][1:]
      if i in scaff_list:
        out=open('%s%s_.temp' % (i,args['abrev']),'w')
        out.write('>%s\n' % (i))
        check=1
      else:
        check=0
    elif check==1:
      out.write(i)
  out.close()
#################

#################


#Step-4 
#Augustus search on candidate regions using the pre-built Block profiles (msa2prfl.pl)
  print('*** Running Augustus prediction ***')
  if os.path.exists('%saugustus' % mainout)==False:
    os.system('mkdir %saugustus' % mainout)

  f=open('coordinates_%s' % args['abrev'])
  dic={};
  for i in f:
    i=i.strip().split('\t');name=i[0]
    if name not in dic:  
      dic[name]=[[i[1],i[2],i[3]]] #scaffold,start and end
    elif name in dic:
      dic[name].append([i[1],i[2],i[3]])
  strings=[]    
  for i in dic:
      if len(dic[i])>1:
        for z in range(0,len(dic[i])):
            subcommand = ('prfl/'+i+'.prfl', dic[i][z][1], dic[i][z][2], clade, target_species, dic[i][z][0]+args['abrev']+'_.temp', mainout+'augustus/'+i+'.out.'+str(z+1))
            #command='augustus --proteinprofile=%(clade)s/%(prot_profile)s --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' % {'prot_profile':'prfl/'+i+'.prfl','start_coord':dic[i][z][1],'end_coord':dic[i][z][2],'clade':clade,'species':target_species,'scaffold':dic[i][z][0]+args['abrev']+'_.temp','output':mainout+'augustus/'+i+'.out.'+str(z+1)}
            strings.append(subcommand)
      else:
        subcommand = ('prfl/'+i+'.prfl', dic[i][0][1], dic[i][0][2], clade, target_species, dic[i][0][0]+args['abrev']+'_.temp', mainout+'augustus/'+i+'.out')
        #command='augustus --proteinprofile=%(clade)s/%(prot_profile)s --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' % {'prot_profile':'prfl/'+i+'.prfl','start_coord':dic[i][0][1],'end_coord':dic[i][0][2],'clade':clade,'species':target_species,'scaffold':dic[i][0][0]+args['abrev']+'_.temp','output':mainout+'augustus/'+i+'.out'}
        strings.append(subcommand)
    
  print("Now splitting %i Augustus runs over %s cpus" % (len(strings), cpus))
  p = multiprocessing.Pool(int(cpus))
  rs = p.map_async(runAugustus, strings)
  p.close()
  while (True):
    if (rs.ready()): break
    remaining = rs._number_left
    print "Waiting for", remaining, "jobs to complete..."
    time.sleep(60)
  

#---------------------------AUGUSTUS steps END -------------------------------------------#
  
#---------------------------HMMER steps START -------------------------------------------#

if mode=='genome' or mode=='hmmer':  #should be augustus
  #STEP-1 EXTRACT AUGUSTUS PROTEINS
  print('*** Extracting predicted proteins ***')
  output_dir = mainout + 'augustus'
  files=os.listdir(output_dir)
  count=0;check=0
  #this sed command seems to be messed up, not sure why it even exists to delete lines 1-3, but I think I fixed it....Jon Palmer
  print('running this strange sed command....')
  for i in files:
    FileIn = mainout + 'augustus/' + i
    subprocess.call(['sed', '-i', '"1,3d"', FileIn])
    #os.system("sed -i '1,3d' %s/augustus/%s" % (output_dir,i)) 
  if os.path.exists(mainout+'augustus_proteins')==False:
    os.system('mkdir %saugustus_proteins' % mainout)  #there are too many of these to fix, shoudld just use os.makedirs.....
  
  for i in files:
    f=open(mainout+'augustus/'+i)
    if i.endswith('.out'):
      out=open('%saugustus_proteins/%s.fas' % (mainout,i[:-4]),'w')
    elif i.endswith(('.1','.2','.3')):
      out=open('%saugustus_proteins/%s.fas.%s' % (mainout,i[:-6],i[-1]),'w')    
    count=0;tr=0
    for line in f:
      if line.startswith('# start gene'):
        tr=1;
      elif tr==1:
        line=line.split();places=[line[0],line[3],line[4]];tr=0
      elif line.startswith('# protein'):
        line=line.strip().split('[')
        count+=1
        out.write('>g%s[%s:%s-%s]\n' % (count,places[0],places[1],places[2]))
        if line[1][-1]==']':
          line[1]=line[1][:-1]        
        out.write(line[1])
        check=1
      else:
        if line.startswith('# end'):
          check=0
          out.write('\n')
        elif check==1:
          line=line.split()[1]
          if line[-1]==']':
            line=line[:-1]
          out.write(line)
  out.close() 

#Run HMMer (genome mode)  
if mode=='genome' or mode=='hmmer':  
  print('*** Running HMMER to confirm orthology of predicted proteins ***')

  files=os.listdir(mainout+'augustus_proteins/')
  if os.path.exists(mainout+'hmmer_output')==False:
    os.system('mkdir %shmmer_output' % mainout)
  
  for i in files:
    if i.endswith('.fas'):
      f=open(mainout+'augustus_proteins/'+i)
      name=i[:-4]
      os.system('hmmsearch --domtblout %(output_file)s.out -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' % 
	    {'input_file':mainout+'/augustus_proteins/'+i,'db_size':Z,'cpu':cpus,'group_file':clade+'/hmms/'+name,'output_file':mainout+'hmmer_output/'+name})
    elif i.endswith(('.1','.2','.3')):
      f=open(mainout+'augustus_proteins/'+i)
      name=i[:-6]
      os.system('hmmsearch --domtblout %(output_file)s -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' % 
	    {'input_file':mainout+'/augustus_proteins/'+i,'db_size':Z,'cpu':cpus,'group_file':clade+'/hmms/'+name,'output_file':mainout+'hmmer_output/'+name+'.out.'+i[-1]}) 

#Run HMMer (transcriptome mode)
if mode=='trans' or mode=='transcriptome':  
  print('*** Running HMMER to confirm transcript orthology ***')
  files=os.listdir('%stranslated_proteins/' % mainout)
  if os.path.exists('%shmmer_output' % mainout)==False:
    os.system('mkdir %shmmer_output' % mainout)
  group='';grouplist=[]
  for i in files:
    if i.endswith('.fas'):
      f=open('%stranslated_proteins/%s' % (mainout,i))
      name=i[:-7];group=transdic[name]
      if group not in grouplist:
        grouplist.append(group)
        os.system('hmmsearch --domtblout %(output_file)s.out.1 -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' % 
	    {'input_file':mainout+'/translated_proteins/'+i,'db_size':Z,'cpu':cpus,'group_file':clade+'/hmms/'+group,'output_file':mainout+'hmmer_output/'+group})
      else:
        grouplist.append(group)
        os.system('hmmsearch --domtblout %(output_file)s.out.%(count)s -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' % 
	    {'input_file':mainout+'/translated_proteins/'+i,'db_size':Z,'cpu':cpus,'group_file':clade+'/hmms/'+group,'output_file':mainout+'hmmer_output/'+group,'count':str(grouplist.count(group))})

#OGS/Proteome module
if mode=='OGS':
  if os.path.exists(mainout+'hmmer_output')==False:
    os.system('mkdir %shmmer_output' % mainout)
  files=os.listdir(clade+'/hmms')
  f2=open('%s/scores_cutoff' % clade)	#open target scores file
  #Load dictionary of HMM expected scores and full list of groups
  score_dic={};
  for i in f2:
    i=i.strip().split()
    try:
      score_dic[i[0]]=float(i[1]); 	#[1] = mean value; [2] = minimum value
    except:
      pass
  totalbuscos=len(list(score_dic.keys()))
  for i in files:
    name=i[:-4]
    if name in score_dic:
      os.system('hmmsearch --domtblout %(output_file)s.out -o temp  -Z %(db_size)s --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' % 
	{'input_file':args['genome'],'db_size':Z,'cpu':cpus,'group_file':clade+'/hmms/'+name,'output_file':mainout+'hmmer_output/'+name})


###*******get list to be re-run
if mode=='genome' or mode=='hmmer':
  print('*** Parsing HMMER results ***')
  #Open the output file; if no name was specified the default name will be used    
  f2=open('%s/scores_cutoff' % clade)	#open target scores file
  #Load dictionary of HMM expected scores and full list of groups
  score_dic={};
  for i in f2:
    i=i.strip().split()
    try:
      score_dic[i[0]]=float(i[1]); 
    except:
      pass
  totalbuscos=len(list(score_dic.keys()))
  f=open('coordinates_%s' % args['abrev'])
  dic={};
  for i in f:
    i=i.strip().split('\t');name=i[0]
    if name not in dic:  
      dic[name]=[[i[1],i[2],i[3]]] #scaffold,start and end
    elif name in dic:
      dic[name].append([i[1],i[2],i[3]])
###*********


####Make summary

#Categorizing genes found in Complete; multi-copy and partial hits
leng_dic={};sd_dic={};
complete=[];frag=[];done=[]
cc=[];fcc=0;mcc=[];unique=[]
if mode=='genome' or mode!='OGS' or mode=='report' or mode=='hmmer':
  temp=os.listdir('%s/hmmer_output' % mainout)
  files=[]
  for i in temp:
    if i.endswith(('.out','.1','.2','.3')):
      files.append(i)
  f=open('%s/lengths_cutoff' % clade)
  for line in f:
    line=line.strip().split()
    leng_dic[line[0]]=float(line[3])
    sd_dic[line[0]]=float(line[2])
  for entry in files:
    f=open('%s/hmmer_output/%s' % (mainout,entry))
    hit_dic={}
    for line in f:
      if line.startswith('#'):
        pass
      else:
        line=line.strip().split()
        score=float(line[7]);group=line[3];prot=line[0];tlen=int(line[2]);qlen=int(line[5])
        if tlen>30*qlen:
           pass
        else:
           if prot not in hit_dic.keys() and score>=score_dic[group]:
             hit_dic[prot]=[[int(line[15]),int(line[16])]]
           elif score>=score_dic[group]:
             hit_dic[prot].append([int(line[15]),int(line[16])])
    length=measuring(hit_dic)  
    try:		#get maximum length of the putative gene in question
      if len(length)==1:
        length=length[0];
      else:
        length=max(length)+1
      sigma=abs(leng_dic[group]-length)/sd_dic[group]
      if sigma<=2:
        complete.append(entry);cc.append(group)
      elif sigma>2:
        frag.append(entry);
    except:
      pass
  #check the multi hits 
  for entry in complete:
    if entry.endswith('.out'):
      name=entry[:-4]
    else:
      name=entry[:-6]
    if name in done:
      if name not in mcc:
        mcc.append(name)
    done.append(name)
  for i in cc:
    if i not in mcc:
      unique.append(i)
  for entry in frag:
    if entry.endswith('.out'):
      name=entry[:-4]
    else:
      name=entry[:-6]
    if name not in done and entry not in complete:
      done.append(name);fcc+=1

if mode=='OGS':
  complete={};frag={};done=[];fcc=[]
  temp=os.listdir('%s/hmmer_output' % mainout)
  files=[]
  for i in temp:
    if i.endswith(('.out','.1','.2','.3')):
      files.append(i)
  f=open('%s/lengths_cutoff' % clade)
  for line in f:
    line=line.strip().split()
    leng_dic[line[0]]=float(line[3])
    sd_dic[line[0]]=float(line[2])
  for entry in files:
    f=open('%s/hmmer_output/%s' % (mainout,entry))
    hit_dic={}
    for line in f:
      if line.startswith('#'):
        pass
      else:
        line=line.strip().split()
        score=float(line[7]);group=line[3];prot=line[0];tlen=int(line[2]);qlen=int(line[5]);prediction=line[0]
        if group not in complete:
          complete[group]=[];frag[group]=[]
        if tlen>30*qlen:
           pass
        else:
           if prot not in hit_dic.keys() and score>=score_dic[group]:
             hit_dic[prot]=[[int(line[15]),int(line[16]),line[7]]]
           elif score>=score_dic[group]:
             hit_dic[prot].append([int(line[15]),int(line[16]),line[7]])
    lengths=measuring(hit_dic)  
    try:		#get maximum length of the putative gene in question
      if len(lengths)==1:
        length=lengths[0];
        sigma=abs(leng_dic[group]-length)/sd_dic[group]
        if sigma<=2:
          complete[group].append([list(hit_dic.keys())[lengths.index(length)],hit_dic[list(hit_dic.keys())[lengths.index(length)]][0][2],length]); 
        elif sigma>2:
          frag[group].append(list(hit_dic.keys())[lengths.index(length)]); 
      else:
        for length in lengths:
        #length=max(lengths)+1
          sigma=abs(leng_dic[group]-length)/sd_dic[group]
          if sigma<=2:
            complete[group].append([list(hit_dic.keys())[lengths.index(length)],hit_dic[list(hit_dic.keys())[lengths.index(length)]][0][2],length]); 
          elif sigma>2:
            frag[group].append(list(hit_dic.keys())[lengths.index(length)]); 
    except:
      pass
  #check the multi hits 
  for entry in complete:
    if len(complete[entry])==0:
      pass
    elif len(complete[entry])==1: #complete
      cc.append(entry)
    elif len(complete[entry])>1:
      mcc.append(entry)
  for entry in frag:
    if len(complete[entry])!=0:
      pass
    elif frag[entry]!=[]:
      fcc.append(entry)
      
#summarize results, print and write to output files
summary=open('short_summary_'+args['abrev'],'w');
if mode=='OGS':
  print('Total complete BUSCOs found in assembly (<2 sigma) :  %s\t(%s duplicated).' % (len(set(cc))+len(set(mcc)),len(mcc)))
  print('Total BUSCOs partially recovered (>2 sigma) :  %s' % (len(fcc)))
else:  
  print('Total complete BUSCOs found in assembly (<2 sigma) :  %s\t(%s duplicated).' % (len(set(unique)),len(mcc)))
  print('Total BUSCOs partially recovered (>2 sigma) :  %s' % (fcc))
print('Total groups searched: %s' % (totalbuscos))
try:
  if mode!='OGS':
    print('Total BUSCOs not found:  %s' % (totalbuscos-(len(set(cc))+fcc)))
  else: 
    print('Total BUSCOs not found:  %s' % (totalbuscos-(len(set(cc))+len(set(mcc))+len(fcc))))
except:
  print('Total BUSCOs not found:  %s' % (totalbuscos-(len(set(cc))+len(fcc))))


summary.write('#Summarized BUSCO benchmarking for file: %s\n#BUSCO was run in mode: %s\n\n' % (args['genome'],mode))
if mode!='OGS' and mode!='trans': 
  summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (shrink((len(set(cc))+len(set(mcc)))/totalbuscos),shrink(len(set(mcc))/totalbuscos),shrink(fcc/totalbuscos),shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),totalbuscos))
elif mode=='OGS':
  summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (shrink((len(set(cc))+len(set(mcc)))/totalbuscos),shrink(len(set(mcc))/totalbuscos),shrink(len(fcc)/totalbuscos),shrink((totalbuscos-(len(set(cc))+len(set(mcc))+len(fcc)))/totalbuscos),totalbuscos))
elif mode=='trans':  
  summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (shrink(len(set(cc))/totalbuscos),shrink(len(set(mcc))/totalbuscos),shrink(fcc/totalbuscos),shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),totalbuscos))

summary.write('Representing:\n')
if mode!='trans' and mode!='OGS':
  summary.write('\t%s\tComplete Single-copy BUSCOs\n' % (len(set(cc))))
  summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc))))
elif mode=='OGS':
  summary.write('\t%s\tComplete Single-copy BUSCOs\n' % (len(set(cc))+len(set(mcc))))
  summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc))))
elif mode=='trans':
  summary.write('\t%s\tComplete Single-copy BUSCOs\n' % (len(set(cc))-len(set(mcc))))
  summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc)))) 
if mode!='OGS':
  summary.write('\t%s\tFragmented BUSCOs\n' % (fcc))
  summary.write('\t%s\tMissing BUSCOs\n' % (totalbuscos-(len(set(cc))+fcc)))
elif mode=='OGS':
  summary.write('\t%s\tFragmented BUSCOs\n' % (len(fcc)))
  summary.write('\t%s\tMissing BUSCOs\n' % (totalbuscos-(len(set(cc))+len(set(mcc))+len(fcc))))


  
summary.write('\t%s\tTotal BUSCO groups searched\n' % (totalbuscos))
summary.close()
summary=open('full_table_%s' % args['abrev'],'w')
#write correct header
if mode=='genome' or mode=='report':
  summary.write('#BUSCO_group\tStatus\tScaffold\tStart\tEnd\tBitscore\tLength\n')
elif mode=='OGS':
  summary.write('#BUSCO_group\tStatus\tGene\tBitscore\tLength\n')
elif mode=='trans' or mode=='transcriptome':
  summary.write('#BUSCO_group\tStatus\tTranscript\tBitscore\tLength\n')
  
temp=os.listdir('%shmmer_output' % mainout);done=[]
files=[]
for i in temp:
  if i.endswith(('.out','.1','.2','.3')):
    files.append(i)

for i in files:
  if i.endswith('.out'):
    name=i[:-4];marker=0
  elif i.endswith(('.1','.2','.3')):
    name=i[:-6];marker=int(i[-1])-1
  f=open('%shmmer_output/%s'% (mainout,i));score=[]
  hit_dic={}
  for line in f:
    if line.startswith('#'):
      pass
    else:
      line=line.strip().split()
      score.append(float(line[7]));group=line[3];prot=line[0];tlen=int(line[2]);qlen=int(line[5]);prediction=line[0]
      if prot not in hit_dic.keys() and float(line[7])>=score_dic[group]:
         hit_dic[prot]=[[int(line[15]),int(line[16]),line[7]]]
      elif float(line[7])>=score_dic[group]:
         hit_dic[prot].append([int(line[15]),int(line[16]),line[7]])
  length=measuring(hit_dic);
  if mode=='genome' or mode=='report' or mode=='hmmer':
    if hit_dic=={}:
      pass
    elif i in complete and name not in mcc:
      summary.write('%s\tComplete\t%s\t%s\t%s\t%s\t%s\n' % (name,dic[group][marker][0],dic[group][marker][1],dic[group][marker][2],max(score),max(length)+1))    
    elif i in complete and name in mcc:
      summary.write('%s\tDuplicated\t%s\t%s\t%s\t%s\t%s\n' % (name,dic[group][marker][0],dic[group][marker][1],dic[group][marker][2],max(score),max(length)+1))    
    elif i in frag and name not in cc and name not in done:
      summary.write('%s\tFragmented\t%s\t%s\t%s\t%s\t%s\n' % (name,dic[group][marker][0],dic[group][marker][1],dic[group][marker][2],max(score),max(length)+1))    
  elif mode=='OGS':
    if hit_dic=={}:
      pass
    elif name in complete and name not in mcc and name in cc:
      summary.write('%s\tComplete\t%s\t%s\t%s\n' % (name,complete[name][0][0],max(score),max(length)+1))    
    elif name in mcc:
      for entry in complete[name]:
        summary.write('%s\tDuplicated\t%s\t%s\t%s\n' % (name,entry[0],entry[1],entry[2]+1))    
    elif name in fcc and name not in cc:
      summary.write('%s\tFragmented\t%s\t%s\t%s\n' % (name,frag[name][0],max(score),max(length)+1))    
  elif mode=='trans' or mode=='Transcriptome':
    if hit_dic=={}:
      pass
    elif i in complete and name not in mcc:
      summary.write('%s\tComplete\t%s\t%s\t%s\n' % (name,dic[group][marker],max(score),max(length)+1))    
    elif i in complete and name in mcc:
      summary.write('%s\tDuplicated\t%s\t%s\t%s\n' % (name,dic[group][marker],max(score),max(length)+1))    
    elif i in frag and name not in cc and name not in done:
      summary.write('%s\tFragmented\t%s\t%s\t%s\n' % (name,dic[group][marker],max(score),max(length)+1))   
summary.close();f.close()


f=open('full_table_%s' % args['abrev'],'r')
lista=[]
for i in f:
  i=i.strip().split()
  if i[0] not in lista:
    lista.append(i[0])
f.close()
out=open('missing_buscos_list_%s' % args['abrev'],'w')	#get final list of missing buscos
f=open('full_table_%s' % args['abrev'],'a')
for i in score_dic.keys():
  if i in lista:
    pass
  else:
    out.write(i+'\n')
    f.write('%s\tMissing\n' % (i))
out.close();f.close()        

#######retraining
exitFlag=0


if mode=='genome' or mode=='genome' or mode=='hmmer':
  if os.path.exists('%sselected' % mainout)==False:
    os.system('mkdir %sselected' % mainout)
  if os.path.exists('%sgffs'% mainout)==False:
    os.system('mkdir %sgffs' % mainout)
  if os.path.exists('%sgb' % mainout)==False:
    os.system('mkdir %sgb' % mainout)

  f=open('full_table_%s' % args['abrev']);
  lista=[];re_run=[]
  for line in f:
    status=line.split()[1]
    if status=='Complete':
      lista.append(line.split()[0])
    elif status=='Missing' or status=='Fragmented':
      re_run.append(line.split()[0]);

  files=os.listdir('%s/hmmer_output' % mainout);chosen=[]
  for i in files:
    if i.endswith('.out'):
      name=i[:-4]
      if name in lista:
        os.system('cp %s/hmmer_output/%s %s/selected/' % (mainout,i,mainout))
        chosen.append(i)
      
  for entry in chosen:
    f=open('%s/selected/%s' % (mainout,entry))
    out=open('%s/gffs/%s' % (mainout,entry),'w');choicy=''
    for line in f:
      if line.startswith('#'):
        pass
      elif choicy=='':
        choicy=line.split()[0].split('[')[0]
    f.close()
    f=open('%s/augustus/%s' % (mainout,entry));check=0
    for line in f:
      if line.startswith('# start gene'):
        name=line.strip().split()[-1]
        if name==choicy:
          check=1
      elif line.startswith('#'):
        check=0
      elif check==1:
        out.write(line)
  f.close()
  for entry in chosen:
    f=open('%s/gffs/%s' % (mainout,entry))      
    os.system('$AUGUSTUS_CONFIG_PATH/../scripts/gff2gbSmallDNA.pl %s/gffs/%s %s 1000 %s/gb/%s.raw.gb' % (mainout,entry,args['genome'],mainout,entry[:-4]))
  
  print('Training augustus gene predictor')  
  os.system('$AUGUSTUS_CONFIG_PATH/../scripts/new_species.pl --species=%s' % (args['abrev'])) #create new species config file from template
  os.system('cat %sgb/*.gb > training_set_%s' % (mainout,args['abrev']))    
  os.system('etraining --species=%s training_set_%s' % (args['abrev'],args['abrev'])) #train on new training set (complete single copy buscos)
  
  if args['long']==True:
     print('Optimizing augustus metaparameters, this may take around 20 hours')
     os.system('$AUGUSTUS_CONFIG_PATH/../scripts/optimize_augustus.pl --species=%s training_set_%s' % (args['abrev'],args['abrev'])) #train on new training set (complete single copy buscos)
     os.system('etraining --species=%s training_set_%s' % (args['abrev'],args['abrev'])) #train on new training set (complete single copy buscos)
  
  
  print('*** Re-running failed predictions with different constraints, total number %s ***' % len(re_run))
  done=[];target_species=args['abrev']
  strings=[];
  hammers=[];
  seds=[];
  ripped=[]
  for item in re_run:
    if item not in dic: #no coordinates found
      pass 
    elif len(dic[item])>1: #more than one target coordinate
      count=0
      for entry in dic[item]:
        count+=1;
        subcommand = ('prfl/'+item, entry[1], entry[2], clade, target_species, entry[0]+args['abrev']+'_.temp','output', mainout+'augustus/'+item+'.out.'+str(count))
        #command='augustus --proteinprofile=%(clade)s/%(prot_profile)s.prfl --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' % {'prot_profile':'prfl/'+item,'start_coord':entry[1],'end_coord':entry[2],'clade':clade,'species':target_species,'scaffold':entry[0]+args['abrev']+'_.temp','output':mainout+'augustus/'+item+'.out.'+str(count)}
        
        strings.append(subcommand)
        
        command = mainout+'augustus/'+item+'.out.'+str(count)
        
        seds.append(command)
        
        ripped.append(item+'.out.'+str(count));
        
        command = (mainout+'/augustus_proteins/'+item+'.fas.'+str(count), Z, clade+'/hmms/'+item, mainout+'hmmer_output/'+item+'.out.'+str(count))
        #command='hmmsearch --domtblout %(output_file)s -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' %   {'input_file':mainout+'/augustus_proteins/'+item+'.fas.'+str(count),'db_size':Z,'cpu':1,'group_file':clade+'/hmms/'+item,'output_file':mainout+'hmmer_output/'+item+'.out.'+str(count)}
        hammers.append(command)
      
    elif len(dic[item])==1:
      entry=dic[item][0]
      try:
        subcommand = ('prfl/'+item, entry[1], entry[2], clade, target_species, entry[0]+args['abrev']+'_.temp','output', mainout+'augustus/'+item+'.out')
        #command='augustus --proteinprofile=%(clade)s/%(prot_profile)s.prfl --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' %    {'prot_profile':'prfl/'+item,'start_coord':entry[1],'end_coord':entry[2],'clade':clade,'species':target_species,'scaffold':entry[0]+args['abrev']+'_.temp','output':mainout+'augustus/'+item+'.out'}
        strings.append(subcommand)
        
        command = mainout+'augustus/'+item+'.out'
        seds.append(command)
        
        ripped.append(item);name=item
        
        command = (mainout+'/augustus_proteins/'+name, Z, clade+'/hmms/'+name, mainout+'hmmer_output/'+name)
        #command='hmmsearch --domtblout %(output_file)s.out -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s.fas' % {'input_file':mainout+'/augustus_proteins/'+name,'db_size':Z,'cpu':1,'group_file':clade+'/hmms/'+name,'output_file':mainout+'hmmer_output/'+name}
        hammers.append(command)
      except:
        pass   
  #missing(mainout,args['abrev'],'missing_buscos_list_')  
  
  
  print("Now splitting %i Augustus runs over %s cpus" % (len(strings), cpus))
  p = multiprocessing.Pool(int(cpus))
  rs = p.map_async(runAugustus, strings)
  p.close()
  while (True):
    if (rs.ready()): break
    remaining = rs._number_left
    print "Waiting for", remaining, "jobs to complete..."
    time.sleep(60)
      
  print("Now splitting %i Sed calls over %s cpus" % (len(seds), cpus))
  p = multiprocessing.Pool(int(cpus))
  rs = p.map_async(runSEDS, seds)
  p.close()
  while (True):
    if (rs.ready()): break
    remaining = rs._number_left
    print "Waiting for", remaining, "jobs to complete..."
    time.sleep(10)
    
  
  print('Starting to run EXTRACT....')
  for entry in ripped:
      extract(mainout,entry)
  
  
  print("Now splitting %i HMMsearch runs over %s cpus" % (len(hammers), cpus))
  p = multiprocessing.Pool(int(cpus))
  rs = p.map_async(runHMMerSearch, hammers)
  p.close()
  while (True):
    if (rs.ready()): break
    remaining = rs._number_left
    print "Waiting for", remaining, "jobs to complete..."
    time.sleep(30)


###retraining and running over

#clean up temporary files
#no reason to use system here? why not use python?  os.system call to rm is dangerous?
if mode!='OGS':
  os.system('rm *%s_.temp' % args['abrev'])
  os.system('rm %s.nsq %s.nin %s.nhr'  % (args['abrev'],args['abrev'],args['abrev']))
  os.system('mv %s_tblastn run_%s' % (args['abrev'],args['abrev']))
  os.system('mv short_summary_%s run_%s' % (args['abrev'],args['abrev']))
  if mode!='trans':
    os.system('mv coordinates_%s run_%s' % (args['abrev'],args['abrev']))
  os.system('mv missing_buscos_list_%s run_%s' % (args['abrev'],args['abrev']))
  os.system('mv full_table_%s run_%s' % (args['abrev'],args['abrev']))
else:
  os.system('mv missing_buscos_list_%s run_%s' % (args['abrev'],args['abrev']))
  os.system('mv full_table_%s run_%s' % (args['abrev'],args['abrev']))
  os.system('mv short_summary_%s run_%s' % (args['abrev'],args['abrev']))
# Report run time per step
print('Total running time:  ',time.time() - start_time, "seconds")

#parse results and write final summary
#Categorizing genes found in Complete; multi-copy and partial hits
leng_dic={};sd_dic={};
complete=[];frag=[];done=[]
cc=[];fcc=0;mcc=[];unique=[]
if mode=='genome' or mode=='genome' or mode=='hmmer':
  temp=os.listdir('%s/hmmer_output' % mainout)
  files=[]
  for i in temp:
    if i.endswith(('.out','.1','.2','.3')):
      files.append(i)
  f=open('%s/lengths_cutoff' % clade)
  for line in f:
    line=line.strip().split()
    leng_dic[line[0]]=float(line[3])
    sd_dic[line[0]]=float(line[2])
  for entry in files:
    f=open('%s/hmmer_output/%s' % (mainout,entry))
    hit_dic={}
    for line in f:
      if line.startswith('#'):
        pass
      else:
        line=line.strip().split()
        score=float(line[7]);group=line[3];prot=line[0];tlen=int(line[2]);qlen=int(line[5])
        if tlen>30*qlen:
           pass
        else:
           if prot not in hit_dic.keys() and score>=score_dic[group]:
             hit_dic[prot]=[[int(line[15]),int(line[16])]]
           elif score>=score_dic[group]:
             hit_dic[prot].append([int(line[15]),int(line[16])])
    length=measuring(hit_dic)  
    try:		#get maximum length of the putative gene in question
      if len(length)==1:
        length=length[0];
      else:
        length=max(length)+1
      sigma=abs(leng_dic[group]-length)/sd_dic[group]
      if sigma<=2:
        complete.append(entry);cc.append(group)
      elif sigma>2:
        frag.append(entry);
    except:
      pass
  #check the multi hits 
  for entry in complete:
    if entry.endswith('.out'):
      name=entry[:-4]
    else:
      name=entry[:-6]
    if name in done:
      if name not in mcc:
        mcc.append(name)
    done.append(name)
  for i in cc:
    if i not in mcc:
      unique.append(i)
  for entry in frag:
    if entry.endswith('.out'):
      name=entry[:-4]
    else:
      name=entry[:-6]
    if name not in done and entry not in complete:
      done.append(name);fcc+=1

#summarize results, print and write to output files
  summary=open('short_summary_'+args['abrev'],'w');
  print('Total complete BUSCOs found in assembly (<2 sigma) :  %s\t(%s duplicated).' % (len(set(unique)),len(mcc)))
  print('Total BUSCOs partially recovered (>2 sigma) :  %s' % (fcc))
  print('Total groups searched: %s' % (totalbuscos))
  try:
    print('Total BUSCOs not found:  %s' % (totalbuscos-(len(set(cc))+fcc)))
  except:
    print('Total BUSCOs not found:  %s' % (totalbuscos-(len(set(cc))+len(fcc))))

  summary.write('#Summarized BUSCO benchmarking for file: %s\n#BUSCO was run in mode: %s\n\n' % (args['genome'],mode))
  summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (shrink(len(set(cc))/totalbuscos),shrink(len(set(mcc))/totalbuscos),shrink(fcc/totalbuscos),shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),totalbuscos))

  summary.write('Representing:\n')
  summary.write('\t%s\tComplete Single-Copy BUSCOs\n' % (len(set(cc))))
  summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc))))
  summary.write('\t%s\tFragmented BUSCOs\n' % (fcc))
  summary.write('\t%s\tMissing BUSCOs\n' % (totalbuscos-(len(set(cc))+fcc)))
  
  summary.write('\t%s\tTotal BUSCO groups searched\n' % (totalbuscos))
  summary.close()
  summary=open('full_table_'+args['abrev'],'w')
  summary.write('#BUSCO_group\tStatus\tScaffold\tStart\tEnd\tBitscore\tLength\n')

  
  temp=os.listdir('%shmmer_output' % mainout);done=[]
  files=[]
  for i in temp:
    if i.endswith(('.out','.1','.2','.3')):
      files.append(i)

  for i in files:
    if i.endswith('.out'):
      name=i[:-4];marker=0
    elif i.endswith(('.1','.2','.3')):
      name=i[:-6];marker=int(i[-1])-1
    f=open('%shmmer_output/%s'% (mainout,i));score=[]
    hit_dic={}
    for line in f:
      if line.startswith('#'):
        pass
      else:
        line=line.strip().split()
        score.append(float(line[7]));group=line[3];prot=line[0];tlen=int(line[2]);qlen=int(line[5]);prediction=line[0]
        if prot not in hit_dic.keys() and float(line[7])>=score_dic[group]:
           hit_dic[prot]=[[int(line[15]),int(line[16]),line[7]]]
        elif float(line[7])>=score_dic[group]:
           hit_dic[prot].append([int(line[15]),int(line[16]),line[7]])
    length=measuring(hit_dic);
    if hit_dic=={}:
      pass
    elif i in complete and name not in mcc:
      summary.write('%s\tComplete\t%s\t%s\t%s\t%s\t%s\n' % (name,dic[group][marker][0],dic[group][marker][1],dic[group][marker][2],max(score),max(length)+1))    
    elif i in complete and name in mcc:
      summary.write('%s\tDuplicated\t%s\t%s\t%s\t%s\t%s\n' % (name,dic[group][marker][0],dic[group][marker][1],dic[group][marker][2],max(score),max(length)+1))    
    elif i in frag and name not in cc and name not in done:
      summary.write('%s\tFragmented\t%s\t%s\t%s\t%s\t%s\n' % (name,dic[group][marker][0],dic[group][marker][1],dic[group][marker][2],max(score),max(length)+1))    
  summary.close()


  f=open('full_table_%s' % args['abrev'],'r');lista=[]
  for i in f:
    i=i.strip().split()
    if i[0] not in lista:
      lista.append(i[0])
  out=open('missing_buscos_list_%s' % args['abrev'],'w')	#get final list of missing buscos
  f=open('full_table_%s' % args['abrev'],'a')
  for i in score_dic.keys():
    if i in lista:
      pass
    else:
      out.write('%s\n' % i)
      f.write('%s\tMissing\n' % (i))
  out.close();f.close()        

  os.system('mv short_summary_%s run_%s' % (args['abrev'],args['abrev']))
  os.system('mv missing_buscos_list_%s run_%s' % (args['abrev'],args['abrev']))
  os.system('mv full_table_%s run_%s' % (args['abrev'],args['abrev']))
  os.system('mv training_set_%s run_%s' % (args['abrev'],args['abrev']))
