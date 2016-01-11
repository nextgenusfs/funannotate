#!/usr/bin/env python

import sys, os, multiprocessing, subprocess, csv, time, re, shutil, inspect
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

log_name = 'funannotate-p2g.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)


def tblastnFilter(input, query, cpus, output):
    global HitList
    HitList = []
    FNULL = open(os.devnull, 'w')
    if not os.path.exists(output):
        os.makedirs(output)
    #start by formatting blast db/dustmasker filtered format
    subprocess.call(['dustmasker', '-in', input, '-infmt', 'fasta', '-parse_seqids', '-outfmt', 'maskinfo_asn1_bin', '-out', 'genome_dust.asnb'], cwd = output, stdout = FNULL, stderr = FNULL)
    subprocess.call(['makeblastdb', '-in', input, '-dbtype', 'nucl', '-parse_seqids', '-mask_data', 'genome_dust.asnb', '-out', 'genome'], cwd = output, stdout = FNULL, stderr = FNULL)
    #okay, now run tblastn using uniprot proteins
    subprocess.call(['tblastn', '-num_threads', str(cpus), '-db', 'genome', '-query', query, '-max_target_seqs', '1', '-db_soft_mask', '11', '-threshold', '999', '-max_intron_length', MaxIntron, '-evalue', '1e-10', '-outfmt', '6', '-out', 'filter.tblastn.tab'], cwd = output, stdout = FNULL, stderr = FNULL)
    #now parse through results, generating a list for exonerate function
    with open(os.path.join(output, 'filter.tblastn.tab')) as input:
        reader = csv.reader(input, delimiter='\t')
        for cols in reader:
            hit = cols[0] + '::' + cols[1]
            if hit not in HitList:
                HitList.append(hit)

def runExonerate(input):
    FNULL = open(os.devnull, 'w')
    s = input.split('::')
    if s[0].startswith('sp|'):
        name = s[0].split("|")[1] + '_' + s[1]
    else:
        name = s[0].split()[0] + '_' + s[1]
    query = os.path.join(tmpdir, name+'.fa')
    with open(query, 'w') as output:
        rec = record_dict[s[0]]
        output.write(">%s\n%s\n" % (rec.id, rec.seq))
    scaffold = s[1] + '.fa'
    scaffold = os.path.join(tmpdir, scaffold)
    exonerate_out = 'exonerate_' + name + '.out'
    exonerate_out = os.path.join(tmpdir, exonerate_out)
    ryo = "AveragePercentIdentity: %pi\n"
    with open(exonerate_out, 'w') as output:
        subprocess.call(['exonerate', '--model', 'p2g', '--showvulgar', 'no', '--showalignment', 'no', '--showquerygff', 'no', '--showtargetgff', 'yes', '--maxintron', MaxIntron, '--percent', '80', '--ryo', ryo , query, scaffold], stdout = output, stderr = FNULL)

genome = os.path.abspath(sys.argv[2])
tmpdir = 'p2g_tmp'
proteins = os.path.abspath(sys.argv[1])
Output = sys.argv[3]
MaxIntron = sys.argv[4]
cpus = int(sys.argv[5])

lib.log.info("Running pre-filter tBlastn step")
tblastnFilter(genome, proteins, cpus, tmpdir)

#split genome fasta into individual scaffolds
if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)
with open(genome, 'rU') as input:
    for record in SeqIO.parse(input, "fasta"):
        SeqIO.write(record, os.path.join(tmpdir, record.id + ".fa"), "fasta")

#Now run exonerate on hits
lib.log.info("Polishing alignments with Exonerate")
record_dict = SeqIO.to_dict(SeqIO.parse(proteins, 'fasta'))
#print HitList
#print record_dict
p = multiprocessing.Pool(cpus)
rs = p.map_async(runExonerate, HitList)
p.close()
while (True):
    if (rs.ready()): break
    #remaining = rs._number_left
    #print "Waiting for", remaining, "exonerate jobs to complete..."
    #time.sleep(30)

#now collect all exonerate results into one
counter = 0
with open(Output, 'wb') as output:
    for root, dirs, files in os.walk(tmpdir):
        for file in files:
            if file.endswith('.out'):
                filename = os.path.join(root, file)
                with open(filename, 'rU') as readfile:
                    read_file = readfile.read()
                    line = read_file.splitlines()[3]
                    sub = read_file.splitlines()[3:]
                    if line.startswith("--"):
                        continue
                    counter += 1
                    for l in sub:
                        output.write(l+'\n')

#finally clean-up your mess
shutil.rmtree(tmpdir)

lib.log.info("Proteins mapped to genome finished: %i alignments generated" % counter)
