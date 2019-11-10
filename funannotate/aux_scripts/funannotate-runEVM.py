#!/usr/bin/env python

import sys
import multiprocessing
import subprocess
import os
import time
import shutil
from itertools import izip_longest
import funannotate.library as lib

# get EVM arguments, genome, protein, transcript, min_intron, weights all from command line
cpus = int(sys.argv[2])
tmpdir = sys.argv[3]
arguments = sys.argv[4:]  # logfile first, num cpus is second
Output = arguments[-1]
del arguments[-1]

log_name = sys.argv[1]
if os.path.isfile(log_name):
    os.remove(log_name)

# initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

# create output directory
if os.path.exists(tmpdir):
    shutil.rmtree(tmpdir)
os.makedirs(tmpdir)

perl = 'perl'
EVM = os.environ['EVM_HOME']
Partition = os.path.join(EVM, 'EvmUtils', 'partition_EVM_inputs.pl')
Commands = os.path.join(EVM, 'EvmUtils', 'write_EVM_commands.pl')
Execute = os.path.join(EVM, 'EvmUtils', 'execute_EVM_commands.pl')
Combine = os.path.join(EVM, 'EvmUtils', 'recombine_EVM_partial_outputs.pl')
Convert = os.path.join(EVM, 'EvmUtils', 'convert_EVM_outputs_to_GFF3.pl')

# need to pull out --genome genome.fasta from arguments list
genome_index = arguments.index('--genome')
genome_args = [arguments[genome_index], arguments[genome_index+1]]
# base commands
base_cmd1 = [perl, Partition, '--segmentSize', '100000',
             '--overlapSize', '10000', '--partition_listing', 'partitions_list.out']
base_cmd2 = [perl, Commands, '--output_file_name',
             'evm.out', '--partitions', 'partitions_list.out']
base_cmd5 = [perl, Convert, '--partitions',
             'partitions_list.out', '--output', 'evm.out']
# combined commands
cmd2 = base_cmd2 + arguments
cmd5 = base_cmd5 + genome_args
# need to remove weights for split partition command
index = arguments.index('--weights')
del arguments[index]
del arguments[index]
cmd1 = base_cmd1 + arguments
# remove intron from partitions command
del cmd1[-1]
del cmd1[-1]


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def worker(input):
    logfile = input + '.log'
    with open(logfile, 'w') as output:
        subprocess.call([perl, Execute, input], stdout=output, stderr=output)


def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        worker(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))


# split partitions
#lib.log.info("Setting up EVM partitions")
lib.runSubprocess(cmd1, tmpdir, lib.log)
#subprocess.call(cmd1, cwd = tmpdir, stdout = FNULL, stderr = FNULL)
# check output
lib.checkinputs(os.path.join(tmpdir, 'partitions_list.out'))

# generate commands
#lib.log.info("Generating EVM command list")
commands = os.path.join(tmpdir, 'commands.list')
with open(commands, 'w') as output:
    subprocess.call(cmd2, cwd=tmpdir, stdout=output, stderr=FNULL)

# count total lines
num_lines = sum(1 for line in open(commands))
# strange thing happens if you try to run with more cpus than commands
if num_lines < cpus:
    x = num_lines
else:
    x = cpus - 1
if x < 1:
    x = 1
lib.log.info("Running EVM commands with %i CPUs" % (x))
#print "Splitting over", cpus, "CPUs"
n = int(round(num_lines / x))

with open(commands, 'rU') as f:
    for i, g in enumerate(grouper(n, f, fillvalue=''), 1):
        with open(os.path.join(tmpdir, 'split_{0}'.format(i * n)+'.cmds'), 'w') as fout:
            fout.writelines(g)

# now launch a process for each split file
file_list = []
for file in os.listdir(tmpdir):
    if file.endswith('cmds'):
        f = os.path.join(tmpdir, file)
        file_list.append(f)

p = multiprocessing.Pool(cpus)
tasks = len(file_list)
results = []
for i in file_list:
    results.append(p.apply_async(safe_run, [i]))
while True:
    incomplete_count = sum(1 for x in results if not x.ready())
    if incomplete_count == 0:
        break
    sys.stdout.write("     Progress: %.2f%% \r" %
                     (float(tasks - incomplete_count) / tasks * 100))
    sys.stdout.flush()
    time.sleep(1)
p.close()
p.join()

# now combine the paritions
#lib.log.info("Combining partitioned EVM outputs")
partitioncmd = [perl, Combine, '--partitions',
                'partitions_list.out', '--output_file_name', 'evm.out']
lib.runSubprocess(partitioncmd, tmpdir, lib.log)

# now convert to GFF3
#lib.log.info("Converting EVM output to GFF3")
lib.runSubprocess(cmd5, tmpdir, lib.log)

# now concatenate all GFF3 files together for a genome then
lib.log.info("Converting to GFF3 and collecting all EVM results")
with open(Output, 'w') as out:
    for root, dirs, files in os.walk(tmpdir):
        for file in files:
            if file == 'evm.out.gff3':
                filename = os.path.join(root, file)
                with open(filename, 'rU') as readfile:
                    shutil.copyfileobj(readfile, out)
