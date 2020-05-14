#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import shutil
import argparse
from Bio import SeqIO
import funannotate.library as lib


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='funannotate-predict.py', usage="%(prog)s [options] -i genome.fasta",
                                     description='''Script that adds a proteome to the outgroups.''',
                                     epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help='Proteome in FASTA format')
    parser.add_argument('-s', '--species', required=True,
                        help='Species name "binomial in quotes"')
    parser.add_argument('-b', '--busco_db', default='dikarya', choices=['fungi', 'microsporidia', 'dikarya', 'ascomycota', 'pezizomycotina', 'eurotiomycetes', 'sordariomycetes', 'saccharomycetes', 'saccharomycetales', 'basidiomycota', 'eukaryota', 'protists',
                                                                        'alveolata_stramenophiles', 'metazoa', 'nematoda', 'arthropoda', 'insecta', 'endopterygota', 'hymenoptera', 'diptera', 'vertebrata', 'actinopterygii', 'tetrapoda', 'aves', 'mammalia', 'euarchontoglires', 'laurasiatheria', 'embryophyta'], help='BUSCO database to use')
    parser.add_argument('-c', '--cpus', default=2, type=int,
                        help='Number of CPUs to use')
    parser.add_argument('-d', '--database',
                        help='Path to funannotate database, $FUNANNOTATE_DB')
    args = parser.parse_args(args)

    if args.database:
        FUNDB = args.database
    else:
        try:
            FUNDB = os.environ["FUNANNOTATE_DB"]
        except KeyError:
            lib.log.error(
                'Funannotate database not properly configured, run funannotate setup.')
            sys.exit(1)

    parentdir = os.path.join(os.path.dirname(__file__))

    # get base name
    species = args.species.replace(' ', '_').lower()+'.'+args.busco_db
    OUTGROUPS = os.path.join(FUNDB, 'outgroups')

    # create log file
    log_name = species+'-add2outgroups.log'
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

    # check buscos, download if necessary
    if not os.path.isdir(os.path.join(FUNDB, args.busco_db)):
        lib.log.error("%s busco database is missing, install with funannotate setup -b %s" %
                      (args.busco_db, args.busco_db))
        sys.exit(1)

    ProtCount = lib.countfasta(args.input)
    lib.log.info('{0:,}'.format(ProtCount) + ' protein records loaded')

    # convert to proteins and screen with busco
    lib.log.info("Looking for BUSCO models with %s DB" % args.busco_db)
    BUSCODB = os.path.join(FUNDB, args.busco_db)
    BUSCO = os.path.join(parentdir, 'aux_scripts', 'funannotate-BUSCO2.py')
    cmd = [sys.executable, BUSCO, '-i', os.path.abspath(
        args.input), '-m', 'proteins', '--lineage', BUSCODB, '-o', species, '--cpu', str(args.cpus), '-f']
    lib.runSubprocess(cmd, '.', lib.log)

    # check that it ran correctly
    busco_results = os.path.join('run_'+species, 'full_table_'+species+'.tsv')
    if not lib.checkannotations(busco_results):
        lib.log.error("BUSCO failed, check logfile")
        sys.exit(1)
    nameChange = {}
    with open(busco_results, 'rU') as input:
        for line in input:
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            if cols[1] == 'Complete':
                if not cols[2] in nameChange:
                    nameChange[cols[2]] = cols[0]
                else:
                    lib.log.error(
                        "Duplicate ID found: %s %s. Removing from results" % (cols[2], cols[0]))
                    del nameChange[cols[2]]

    # output counts
    lib.log.info('{0:,}'.format(len(nameChange)) + ' BUSCO models found')

    # index the proteome for parsing
    SeqRecords = SeqIO.to_dict(SeqIO.parse(args.input, 'fasta'))

    # setup output proteome
    busco_out = os.path.join(OUTGROUPS, species+'_buscos.fa')
    with open(busco_out, 'w') as output:
        for k, v in list(nameChange.items()):
            rec = SeqRecords[k]
            output.write('>%s\n%s\n' % (v, rec.seq))
    lib.log.info("Results written to: %s" % busco_out)

    # clean up your mess
    shutil.rmtree('run_'+species)
    shutil.rmtree('tmp')


if __name__ == "__main__":
    main(sys.argv[1:])
