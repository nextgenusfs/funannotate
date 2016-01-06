#!/usr/bin/env python

import sys, multiprocessing, subprocess, os, time, glob, shutil
from itertools import izip_longest

perl = '/usr/bin/perl'
EVM = os.environ['EVM_HOME']
Partition = os.path.join(EVM, 'EvmUtils', 'partition_EVM_inputs.pl')
Commands = os.path.join(EVM, 'EvmUtils', 'write_EVM_commands.pl')
Execute = os.path.join(EVM, 'EvmUtils', 'execute_EVM_commands.pl')
Combine = os.path.join(EVM, 'EvmUtils', 'recombine_EVM_partial_outputs.pl')
Convert = os.path.join(EVM, 'EvmUtils', 'convert_EVM_outputs_to_GFF3.pl')

def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def worker(input):
    logfile = input + '.log'
    with open(logfile, 'w') as output:
        subprocess.call([perl, Execute, input], stdout = output, stderr = output)

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: worker(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))

FNULL = open(os.devnull, 'w')

#split partitions
print "Setting up EVM partitions"
subprocess.call([perl, Partition, '--genome', sys.argv[1], '--gene_predictions', sys.argv[2], '--protein_alignments', sys.argv[3], '--min_intron_length', '10', '--segmentSize', '100000', '--overlapSize', '10000', '--partition_listing', 'partitions_list.out'], stdout = FNULL, stderr = FNULL)

#generate commands
print "Generating EVM command list"
weights = os.path.abspath(sys.argv[4])
with open('commands.list', 'w') as output:
    subprocess.call([perl, Commands, '--genome', sys.argv[1], '--gene_predictions', sys.argv[2], '--protein_alignments', sys.argv[3], '--weights', weights, '--min_intron_length', '10', '--output_file_name', 'evm.out', '--partitions', 'partitions_list.out'], stdout = output, stderr = FNULL)


#count total lines
print "Now running commands in parallel"
num_lines = sum(1 for line in open('commands.list'))
print num_lines, "commands to run"
cpus = multiprocessing.cpu_count() - 14
print "Splitting over", cpus, "CPUs"
n = int(round(num_lines / cpus))

with open('commands.list', 'rU') as f:
    for i, g in enumerate(grouper(n, f, fillvalue=''), 1):
        with open('split_{0}'.format(i * n)+'.cmds', 'w') as fout:
            fout.writelines(g)

#now launch a process for each split file
files = ((os.path.join('.', f))
            for f in os.listdir('.') if f.endswith('cmds'))
p = multiprocessing.Pool(cpus)
rs = p.map_async(safe_run, files)
p.close()
while (True):
    if (rs.ready()): break
    remaining = rs._number_left
    print "Waiting for", remaining, "tasks to complete..."
    time.sleep(30)

#now combine the paritions
print "Combining EVM outputs that were partitioned"
subprocess.call([perl, Combine, '--partitions', 'partitions_list.out', '--output_file_name', 'evm.out'], stdout = FNULL, stderr = FNULL)

#now convert to GFF3
print "Converting EVM output to GFF3"
subprocess.call([perl, Convert, '--partitions', 'partitions_list.out', '--output', 'evm.out', '--genome', sys.argv[1]], stdout = FNULL, stderr = FNULL)

#now concatenate all GFF3 files together for a genome then
print "Now collecting all results"
with open('evm.all.gff3', 'w') as output:
    for root, dirs, files in os.walk('.'):
        for file in files:
            if file == 'evm.out.gff3':
                filename = os.path.join(root,file)
                with open(filename, 'rU') as readfile:
                    shutil.copyfileobj(readfile, output)

print "Cleaning up mess...." #remove all the folders in this directory, all
d='.'
dirs = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
for i in dirs:
    shutil.rmtree(i)
for file in os.listdir('.'):
    if file.endswith('cmds') or file.endswith('log'):
        os.remove(file)

print "holy shit, it worked???"