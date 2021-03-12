#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os
import subprocess
import shutil
import argparse
import io
from datetime import datetime
from goatools import obo_parser
from Bio import SeqIO
from natsort import natsorted
import pandas as pd
import funannotate.library as lib
import funannotate.resources as resources


def AnnotationFound(input):
    # test if the list of dictionaries has annotation
    test = False
    for x in input:
        if len(x) > 0:
            test = True
    return test


def main(args):
        # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='funannotate-compare.py', usage="%(prog)s [options] genome1.gbk genome2.gbk",
                                     description='''Funannotate comparative genomics.''',
                                     epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--input', nargs='+',
                        help='List of funannotate genome folders or GBK files')
    parser.add_argument('-o', '--out', default='funannotate_compare',
                        help='Name of output folder')
    parser.add_argument('--cpus', default=2, type=int,
                        help='Number of CPUs to utilize')
    parser.add_argument('--go_fdr', default=0.05, type=float,
                        help='P-value for FDR GO-enrichment')
    parser.add_argument('--heatmap_stdev', default=1.0, type=float,
                        help='Standard Deviation threshold for heatmap retention')
    parser.add_argument('--bootstrap', default=100, type=int,
                        help='Number of bootstraps to run with RAxML')
    parser.add_argument('--num_orthos', default=500, type=int,
                        help='Number of Single-copy orthologs to run with RAxML')
    parser.add_argument('--outgroup',
                        help='Name of species for RAxML outgroup')
    parser.add_argument('--eggnog_db', default='fuNOG', help='EggNog database')
    parser.add_argument('--run_dnds', choices=['estimate', 'full'],
                        help='Run dN/dS analysis with codeML for each ortholog (long runtime)')
    parser.add_argument('--proteinortho',
                        help='Pre-computed ProteinOrtho POFF')
    parser.add_argument('--ml_method', default='iqtree',
                        choices=['raxml', 'iqtree'], help='ML method')
    parser.add_argument('-d', '--database',
                        help='Path to funannotate database, $FUNANNOTATE_DB')
    parser.add_argument('--no-progress', dest='progress', action='store_false',
                        help='no progress on multiprocessing')
    args = parser.parse_args(args)

    parentdir = os.path.join(os.path.dirname(__file__))


    # setup funannotate DB path
    if args.database:
        FUNDB = args.database
    else:
        try:
            FUNDB = os.environ["FUNANNOTATE_DB"]
        except KeyError:
            print('ERROR: Funannotate database not properly configured, run funannotate setup.')
            sys.exit(1)

    # check database sources, so no problems later
    sources = [os.path.join(FUNDB, 'Pfam-A.clans.tsv'),
               os.path.join(FUNDB, 'interpro.tsv'),
               os.path.join(FUNDB, 'go.obo')]
    if not all([os.path.isfile(f) for f in sources]):
        print('ERROR: Database files not found in %s, run funannotate database and/or funannotate setup' % FUNDB)
        sys.exit(1)

    # remove slashes if they exist in output
    args.out = args.out.replace('/', '')

    # make output folder
    if not os.path.isdir(args.out):
        os.makedirs(args.out)
    go_folder = os.path.join(args.out, 'go_terms')
    protortho = os.path.join(args.out, 'protortho')
    phylogeny = os.path.join(args.out, 'phylogeny')
    ortho_folder = os.path.join(args.out, 'orthology')
    if os.path.isdir(go_folder):
        shutil.rmtree(go_folder)
        os.makedirs(go_folder)
    else:
        os.makedirs(go_folder)
    if not os.path.isdir(protortho):
        os.makedirs(protortho)
    if not os.path.isdir(phylogeny):
        os.makedirs(phylogeny)
    if not os.path.isdir(ortho_folder):
        os.makedirs(ortho_folder)

    # create log file
    log_name = os.path.join(args.out, 'funannotate-compare.log')
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

    if args.outgroup:
        if not os.path.isdir(os.path.join(FUNDB, 'outgroups')):
            lib.log.error("Outgroup folder is not properly configured")
            os._exit(1)
        files = [f for f in os.listdir(os.path.join(FUNDB, 'outgroups'))]
        files = [x.replace('_buscos.fa', '') for x in files]
        files = [x for x in files if not x.startswith('.')]
        if not args.outgroup in files:
            lib.log.error("%s is not found in outgroups" % args.outgroup)
            print(lib.list_columns(natsorted(files), cols=3))
        else:
            outgroup = True
            outgroup_species = os.path.join(
                FUNDB, 'outgroups', args.outgroup+'_buscos.fa')
            outgroup_name = args.outgroup
    else:
        outgroup = False
        outgroup_species = ''
        outgroup_name = ''

    # check dependencies and set path to proteinortho
    if args.run_dnds:
        programs = ['find_enrichment.py', 'mafft',
                    'trimal', 'proteinortho', 'ete3', 'phyml']
    else:
        programs = ['find_enrichment.py',
                    'mafft', 'trimal', 'proteinortho']
    if args.ml_method == 'raxml':
        programs = programs + ['raxmlHPC-PTHREADS']
    else:
        programs = programs + ['iqtree']
    lib.CheckDependencies(programs)

    # copy over html files
    if not os.path.isdir(os.path.join(args.out, 'css')):
        lib.copyDirectory(os.path.join(
            parentdir, 'html_template', 'css'), os.path.join(args.out, 'css'))
    if not os.path.isdir(os.path.join(args.out, 'js')):
        lib.copyDirectory(os.path.join(
            parentdir, 'html_template', 'js'), os.path.join(args.out, 'js'))

    # write input to logfile
    lib.log.debug("Input files/folders: %s\n" % args.input)

    # loop through each genome
    stats = []
    merops = []
    ipr = []
    cazy = []
    pfam = []
    eggnog = []
    busco = []
    gbkfilenames = []
    scinames = []
    signalp = []
    secmet = []
    sm_backbones = []
    transmembrane = []
    cogs = []
    names_seen = []
    num_input = len(args.input)
    if num_input == 0:
        lib.log.error("Error, you did not specify an input, -i")
        os._exit(1)
    lib.log.info("Now parsing %i genomes" % num_input)
    for i in range(0, num_input):
        # parse the input, I want user to give output folder for funannotate, put they might give a results folder, so do the best you can to check
        if not os.path.isdir(args.input[i]):
            if args.input[i].endswith('.gbk') or args.input[i].endswith('.gbff'):
                GBK = args.input[i]
            else:
                lib.log.error(
                    "Error, %s is not a funannotate folder and does not seem to be a GenBank file." % args.input[i])
                os._exit(1)
        else:  # split arguments into genomes and run a bunch of stats/comparisons
            # look for annotate_results folder
            GBK = ''
            # this means was not passed the whole folder
            if not os.path.isdir(os.path.join(args.input[i], 'annotate_results')):
                # set fun_dir up a directory to find other results if needed
                for file in os.listdir(args.input[i]):
                    if file.endswith('.gbk'):
                        GBK = os.path.join(args.input[i], file)
            else:  # whole folder is passed, now get the genbank file
                for file in os.listdir(os.path.join(args.input[i], 'annotate_results')):
                    if file.endswith('.gbk'):
                        GBK = os.path.join(
                            args.input[i], 'annotate_results', file)
        if not GBK:  # check this
            lib.log.error(
                "Error, was not able to find appropriate GenBank file in the annotate_results folder")
        gbkfilenames.append(GBK)
        # now run genome routines
        genomeStats = lib.genomeStats(GBK)
        if int(genomeStats[9].replace(',', '')) == 0:
            lib.log.error(
                "%s contains 0 gene models, exiting script" % genomeStats[0])
            sys.exit(1)
        stats.append(genomeStats)
        # this function will return list of dictionaries for each functional category
        functional = lib.getGBKannotation(GBK, FUNDB)
        # split those dictionaries and append to master list for each group of annotation
        pfam.append(functional[0])
        ipr.append(functional[1])
        cazy.append(functional[5])
        busco.append(functional[3])
        signalp.append(functional[7])
        transmembrane.append(functional[8])
        cogs.append(functional[6])
        secmet.append(functional[9])
        sm_backbones.append(functional[10])
        eggnog.append(functional[2])
        merops.append(functional[4])
        if stats[i][1]:
            name = stats[i][0].replace(' ', '_')+'_'+stats[i][1]
        else:
            name = stats[i][0].replace(' ', '_')
        if not name in names_seen:
            names_seen.append(name)
        else:
            if '-' in name:
                base, num = name.split('-')
                name = base+'-'+str(num+1)
            else:
                name = name+'-1'
        lib.parseGOterms(GBK, go_folder, name)
        lib.gb2proteinortho(GBK, protortho, name)
        scinames.append(name)

    # convert busco to dictionary
    busco = lib.busco_dictFlip(busco)

    # add species names to pandas table check for duplicate locus tag IDs
    names = []
    tags = []
    for i in stats:
        # locus tags
        tags.append(i[2])
        # naming
        sci_name = i[0]
        # isolate
        isolate_name = i[1]
        if '_' in sci_name:  # here I'm assuming that somebody used an abbreviated name and an underscore, this would be atypical I think
            names.append(sci_name)
        else:
            genus = sci_name.split(' ')[0]
            species = ' '.join(sci_name.split(' ')[1:])
            abbrev = genus[:1] + '.'
            if isolate_name:
                final_name = abbrev + ' ' + species + ' ' + isolate_name
            else:
                final_name = abbrev + ' ' + species
            names.append(final_name)

    if len(tags) != len(set(tags)):
        locus_tag_names = ['Genome\tlocus_tag']
        for x in range(0, len(stats)):
            value = scinames[x]+'\t'+stats[x][2]
            locus_tag_names.append(value)
        lib.log.error("Duplicate locus_tags found, each genome must have unique locus_tag ID\n%s" %
                      '\n'.join(locus_tag_names))
        lib.log.error(
            "You will either need to re-run prediction step with unique --names, or change one of the tags with find/replace (i.e. sed)")
        sys.exit(1)

    #Secondary metabolism#############################################
    # log raw data
    #lib.log.debug("Secondary metabolite raw data:\n%s" % secmet)
    if AnnotationFound(secmet):
        lib.log.info("Summarizing secondary metabolism gene clusters")
        if not os.path.isdir(os.path.join(args.out, 'secmet')):
            os.makedirs(os.path.join(args.out, 'secmet'))
        SM = {'NRPS': 'nonribosomal peptide synthase', 'PKS': 'polyketide synthase',
              'Hybrid': 'hybrid NRPS-PKS', 'Other': 'other backbone enzyme'}
        # first loop through results and add 'other' field to dictionary
        for i in range(0, len(secmet)):
            num_clusters = len(secmet[i])
            total = 0
            for k, v in sm_backbones[i].items():
                total += v
            others = num_clusters - total
            sm_backbones[i]['Other'] = others

        smdf = pd.DataFrame(sm_backbones)
        smdf['species'] = names
        smdf.set_index('species', inplace=True)

        # reorder the table NRPS, PKS, Hybrid, Other
        smorder = ['NRPS', 'PKS', 'Hybrid', 'Other']
        smdf = smdf[smorder]

        # draw plots for SM data
        # get totals for determining height of y-axis
        totals = smdf.sum(numeric_only=True, axis=1)
        max_num = max(totals)
        round_max = int(lib.roundup(max_num))
        diff = round_max - int(max_num)
        if diff < 100:
            ymax = round_max + 100
        else:
            ymax = round_max
        if round_max == 100 and diff > 30:
            ymax = max_num + 5

        # output to csv
        smdf.transpose().to_csv(os.path.join(args.out, 'secmet', 'SM.summary.results.csv'))

        # stackedbar graph
        if len(args.input) > 1:
            lib.drawStackedBar(smdf, 'Secondary Metabolism Clusters', SM, ymax, os.path.join(
                args.out, 'secmet', 'SM.graph.pdf'))

        # create html output
        with open(os.path.join(args.out, 'secmet.html'), 'w') as output:
            output.write(lib.HEADER)
            output.write(lib.SECMET)
            output.write(lib.FOOTER)
    else:
        lib.log.info("No secondary metabolite annotations found")
        # create html output
        with open(os.path.join(args.out, 'secmet.html'), 'w') as output:
            output.write(lib.HEADER)
            output.write(lib.MISSING)
            output.write(lib.FOOTER)

    #############################################
    #PFAM#############################################
    lib.log.info("Summarizing PFAM domain results")
    if not os.path.isdir(os.path.join(args.out, 'pfam')):
        os.makedirs(os.path.join(args.out, 'pfam'))

    # convert to counts
    pfamdf = lib.convert2counts(pfam)
    pfamdf.fillna(0, inplace=True)
    pfamdf['species'] = names
    pfamdf.set_index('species', inplace=True)

    # remove any "empty" genomes
    pfamdf = pfamdf[(pfamdf.T != 0).any()]

    # make an nmds
    if len(pfamdf.index) > 1:  # make sure number of species is at least two
        lib.distance2mds(pfamdf, 'braycurtis', 'PFAM',
                         os.path.join(args.out, 'pfam', 'PFAM.nmds.pdf'))

    # get the PFAM descriptions
    pfamdf2 = pfamdf.transpose().astype(int)
    PFAM = lib.pfam2dict(os.path.join(FUNDB, 'Pfam-A.clans.tsv'))
    pfam_desc = []
    for i in pfamdf2.index.values:
        pfam_desc.append(PFAM.get(i))
    pfamdf2['descriptions'] = pfam_desc
    # write to file
    pfamdf2.to_csv(os.path.join(args.out, 'pfam', 'pfam.results.csv'), encoding='utf-8')
    pfamdf2.reset_index(inplace=True)
    pfamdf2.rename(columns={'index': 'PFAM'}, inplace=True)
    pfamdf2['PFAM'] = '<a target="_blank" href="http://pfam.xfam.org/family/' + \
        pfamdf2['PFAM'].astype(str)+'">'+pfamdf2['PFAM']+'</a>'
    # create html output
    with open(os.path.join(args.out, 'pfam.html'), 'w') as output:
        try:
            pd.set_option('display.max_colwidth', None)
        except ValueError:
            pd.set_option('display.max_colwidth', 0)
        output.write(lib.HEADER)
        output.write(lib.PFAM)
        output.write(pfamdf2.to_html(
            index=False, escape=False, classes='table table-hover'))
        output.write(lib.FOOTER)

    ##################################################

    ####InterProScan##################################
    lib.log.info("Summarizing InterProScan results")
    if not os.path.isdir(os.path.join(args.out, 'interpro')):
        os.makedirs(os.path.join(args.out, 'interpro'))

    # convert to counts
    IPRdf = lib.convert2counts(ipr)
    IPRdf.fillna(0, inplace=True)  # fill in zeros for missing data
    IPRdf['species'] = names
    IPRdf.set_index('species', inplace=True)

    # some checking here of data, if genome is missing, i.e. counts are zero, drop it
    #print IPRdf
    #print len(IPRdf.columns)
    IPRdf = IPRdf[(IPRdf.T != 0).any()]
    #print len(IPRdf.index)

    # analysis of InterPro Domains
    # get IPR descriptions, we only need to get descriptions for terms in our study, limits memory footprint hopefully?
    uniqIPR = []
    for i in ipr:
        for x in i:
            if not x in uniqIPR:
                uniqIPR.append(x)
    uniqIPR = set(uniqIPR)
    lib.log.info("Loading InterPro descriptions")
    INTERPRO = lib.iprTSV2dict(os.path.join(FUNDB, 'interpro.tsv'), uniqIPR)
    # NMDS
    if len(IPRdf.index) > 1:  # count number of species
        if len(IPRdf.columns) > 1:  # count number of IPR domains
            lib.distance2mds(IPRdf, 'braycurtis', 'InterProScan', os.path.join(
                args.out, 'interpro', 'InterProScan.nmds.pdf'))

            # write to csv file
            ipr2 = IPRdf.transpose().astype(int)
            ipr_desc = []
            for i in ipr2.index.values:
                ipr_desc.append(INTERPRO.get(i))
            ipr2['descriptions'] = ipr_desc
            ipr2.to_csv(os.path.join(args.out, 'interpro',
                                     'interproscan.results.csv'))
            ipr2.reset_index(inplace=True)
            ipr2.rename(columns={'index': 'InterPro'}, inplace=True)
            ipr2['InterPro'] = '<a target="_blank" href="http://www.ebi.ac.uk/interpro/entry/' + \
                ipr2['InterPro'].astype(str)+'">'+ipr2['InterPro']+'</a>'

    # create html output
    with open(os.path.join(args.out, 'interpro.html'), 'w') as output:
        try:
            pd.set_option('display.max_colwidth', None)
        except ValueError:
            pd.set_option('display.max_colwidth', 0)
        output.write(lib.HEADER)
        output.write(lib.INTERPRO)
        if len(IPRdf.columns) > 1:
            if len(IPRdf.index) > 1:
                output.write(ipr2.to_html(
                    index=False, escape=False, classes='table table-hover'))
        output.write(lib.FOOTER)

    ##############################################

    ####MEROPS################################
    lib.log.info("Summarizing MEROPS protease results")
    if not os.path.isdir(os.path.join(args.out, 'merops')):
        os.makedirs(os.path.join(args.out, 'merops'))

    MEROPS = {'A': 'Aspartic Peptidase', 'C': 'Cysteine Peptidase', 'G': 'Glutamic Peptidase', 'M': 'Metallo Peptidase', 'N': 'Asparagine Peptide Lyase',
              'P': 'Mixed Peptidase', 'S': 'Serine Peptidase', 'T': 'Threonine Peptidase', 'U': 'Unknown Peptidase', 'I': 'Protease Inhibitors'}
    # convert to counts
    updatemerops = []
    for x in range(0, len(merops)):
        if None in merops[x]:
            lib.log.error(
                "Merops annotation was run with older database than currently installed, please re-run funannotate annotate for %s" % scinames[x])
            updatemerops.append({})
        else:
            updatemerops.append(merops[x])

    meropsdf = lib.convert2counts(updatemerops)
    meropsdf.fillna(0, inplace=True)
    meropsdf['species'] = names
    meropsdf.set_index('species', inplace=True)

    # make a simple table with just these numbers
    meropsA = meropsdf.filter(regex='^A').sum(numeric_only=True, axis=1)
    meropsC = meropsdf.filter(regex='^C').sum(numeric_only=True, axis=1)
    meropsG = meropsdf.filter(regex='^G').sum(numeric_only=True, axis=1)
    meropsM = meropsdf.filter(regex='^M').sum(numeric_only=True, axis=1)
    meropsN = meropsdf.filter(regex='^N').sum(numeric_only=True, axis=1)
    meropsP = meropsdf.filter(regex='^P').sum(numeric_only=True, axis=1)
    meropsS = meropsdf.filter(regex='^S').sum(numeric_only=True, axis=1)
    meropsT = meropsdf.filter(regex='^T').sum(numeric_only=True, axis=1)
    meropsU = meropsdf.filter(regex='^U').sum(numeric_only=True, axis=1)
    meropsI = meropsdf.filter(regex='^I').sum(numeric_only=True, axis=1)
    # get totals for determining height of y-axis
    totals = meropsdf.sum(numeric_only=True, axis=1)
    max_num = max(totals)
    round_max = int(lib.roundup(max_num))
    diff = round_max - int(max_num)
    if diff < 100:
        ymax = round_max + 100
    else:
        ymax = round_max
    if round_max == 100 and diff > 50:
        ymax = max_num + 10
    # recombine sums
    enzymes = ['A', 'C', 'G', 'M', 'N', 'P', 'S', 'T', 'U', 'I']
    meropsShort = pd.concat([meropsA, meropsC, meropsG, meropsM, meropsN,
                             meropsP, meropsS, meropsT, meropsU, meropsI], axis=1, keys=enzymes)
    meropsShort['species'] = names
    meropsShort.set_index('species', inplace=True)
    # remove any columns with no hits
    meropsShort = meropsShort.loc[:, (meropsShort != 0).any(axis=0)]
    meropsall = meropsdf.transpose()

    # write to file
    meropsdf.transpose().to_csv(os.path.join(
        args.out, 'merops', 'MEROPS.all.results.csv'))
    meropsShort.transpose().to_csv(os.path.join(
        args.out, 'merops', 'MEROPS.summary.results.csv'))

    if not meropsall.empty:
        # draw plots for merops data
        # stackedbar graph
        if len(args.input) > 1:
            lib.drawStackedBar(meropsShort, 'MEROPS families', MEROPS, ymax, os.path.join(
                args.out, 'merops', 'MEROPS.graph.pdf'))

        # drawheatmap of all merops families where there are any differences
        if len(args.input) > 1:
            stdev = meropsall.std(axis=1)
            meropsall['stdev'] = stdev
            if len(meropsall) > 25:
                df2 = meropsall[meropsall.stdev >= args.heatmap_stdev]
                lib.log.info("found %i/%i MEROPS familes with stdev >= %f" %
                             (len(df2), len(meropsall), args.heatmap_stdev))
            else:
                df2 = meropsall
                lib.log.info("found %i MEROPS familes" % (len(df2)))
            meropsplot = df2.drop('stdev', axis=1)
            if len(meropsplot) > 0:
                lib.drawHeatmap(meropsplot, 'BuPu', os.path.join(
                    args.out, 'merops', 'MEROPS.heatmap.pdf'), 6, False)

            meropsall = meropsall.astype(int)
            meropsall.reset_index(inplace=True)
            meropsall.rename(columns={'index': 'MEROPS'}, inplace=True)
            meropsall['MEROPS'] = '<a target="_blank" href="https://www.ebi.ac.uk/merops/cgi-bin/famsum?family=' + \
                meropsall['MEROPS'].astype(str)+'">'+meropsall['MEROPS']+'</a>'

    # create html output
    with open(os.path.join(args.out, 'merops.html'), 'w') as output:
        try:
            pd.set_option('display.max_colwidth', None)
        except ValueError:
            pd.set_option('display.max_colwidth', 0)
        output.write(lib.HEADER)
        output.write(lib.MEROPS)
        output.write(meropsall.to_html(
            escape=False, index=False, classes='table table-hover'))
        output.write(lib.FOOTER)

    #######################################################

    #####run CAZy routine#################################
    lib.log.info("Summarizing CAZyme results")
    if not os.path.isdir(os.path.join(args.out, 'cazy')):
        os.makedirs(os.path.join(args.out, 'cazy'))
    # convert to counts
    CAZydf = lib.convert2counts(cazy)

    # with CAZy there are 7 possible families
    CAZY = {'CBM': 'Carbohydrate-binding module', 'CE': 'Carbohydrate esterase', 'GH': 'Glycoside hydrolase',
            'GT': 'Glycosyltransferase', 'PL': 'Polysaccharide lyase', 'AA': 'Auxillary activities'}
    # make a simple table with just these numbers
    cazyAA = CAZydf.filter(regex='^AA').sum(numeric_only=True, axis=1)
    cazyGT = CAZydf.filter(regex='^GT').sum(numeric_only=True, axis=1)
    cazyPL = CAZydf.filter(regex='^PL').sum(numeric_only=True, axis=1)
    cazyCE = CAZydf.filter(regex='^CE').sum(numeric_only=True, axis=1)
    cazyCBM = CAZydf.filter(regex='^CBM').sum(numeric_only=True, axis=1)
    cazyGH = CAZydf.filter(regex='^GH').sum(numeric_only=True, axis=1)
    # get totals for determining height of y-axis
    totals = CAZydf.sum(numeric_only=True, axis=1)
    max_num = max(totals)
    round_max = int(lib.roundup(max_num))
    diff = round_max - int(max_num)
    if diff < 100:
        ymax = round_max + 100
    else:
        ymax = round_max
    if round_max == 100 and diff > 50:
        ymax = max_num + 10
    #print max_num, round_max, diff, ymax
    enzymes = ['AA', 'CBM', 'CE', 'GH', 'GT', 'PL']
    CAZyShort = pd.concat(
        [cazyAA, cazyCBM, cazyCE, cazyGH, cazyGT, cazyPL], axis=1, keys=enzymes)
    CAZydf['species'] = names
    CAZyShort['species'] = names
    CAZydf.set_index('species', inplace=True)
    CAZyShort.set_index('species', inplace=True)
    cazyall = CAZydf.transpose()

    # write to file
    CAZydf.transpose().to_csv(os.path.join(args.out, 'cazy', 'CAZyme.all.results.csv'))
    CAZyShort.transpose().to_csv(os.path.join(
        args.out, 'cazy', 'CAZyme.summary.results.csv'))

    # draw stacked bar graph for CAZY's
    if len(args.input) > 1:
        lib.drawStackedBar(CAZyShort, 'CAZyme families', CAZY,
                           ymax, os.path.join(args.out, 'cazy', 'CAZy.graph.pdf'))

    # if num of cazys greater than 25, drawheatmap of all CAZys that have standard deviation > X
    if len(args.input) > 1:
        stdev = cazyall.std(axis=1)
        cazyall['stdev'] = stdev
        if len(cazyall) > 25:
            df2 = cazyall[cazyall.stdev >= args.heatmap_stdev]
            lib.log.info("found %i/%i CAZy familes with stdev >= %f" %
                         (len(df2), len(cazyall), args.heatmap_stdev))
        else:
            df2 = cazyall
            lib.log.info("found %i CAZy familes" % (len(df2)))
        cazyplot = df2.drop('stdev', axis=1)
        if len(cazyplot) > 0:
            lib.drawHeatmap(cazyplot, 'YlOrRd', os.path.join(
                args.out, 'cazy', 'CAZy.heatmap.pdf'), 4, False)

        cazyall = cazyall.astype(int)
        cazyall.reset_index(inplace=True)
        cazyall.rename(columns={'index': 'CAZy'}, inplace=True)
        cazyall['CAZy'] = '<a target="_blank" href="http://www.cazy.org/' + \
            cazyall['CAZy'].astype(str)+'.html">'+cazyall['CAZy']+'</a>'

    # create html output
    with open(os.path.join(args.out, 'cazy.html'), 'w') as output:
        try:
            pd.set_option('display.max_colwidth', None)
        except ValueError:
            pd.set_option('display.max_colwidth', 0)
        output.write(lib.HEADER)
        output.write(lib.CAZY)
        output.write(cazyall.to_html(
            escape=False, index=False, classes='table table-hover'))
        output.write(lib.FOOTER)
    ########################################################

    ######COG families#####################
    if AnnotationFound(cogs):
        lib.log.info("Summarizing COG results")
        if not os.path.isdir(os.path.join(args.out, 'cogs')):
            os.makedirs(os.path.join(args.out, 'cogs'))

        COGSdf = lib.convert2counts(cogs)
        COGSdf['species'] = names
        COGSdf.set_index('species', inplace=True)
        COGSdf2 = COGSdf.div(COGSdf.sum(axis=1), axis=0).multiply(100)
        lib.donutplot(COGSdf2, resources.COGS, os.path.join(
            args.out, 'cogs', 'COGS.graph.pdf'))

        # make the csv file more informative
        LongNames = []
        for i in COGSdf.columns.tolist():
            if i in resources.COGS:
                LongNames.append(resources.COGS.get(i))
        COGSdf.columns = LongNames
        COGSdf = COGSdf.astype(int)
        COGSdf.transpose().to_csv(os.path.join(args.out, 'cogs', 'COGS.all.results.csv'))
        # create html output
        with open(os.path.join(args.out, 'cogs.html'), 'w') as output:
            try:
                pd.set_option('display.max_colwidth', None)
            except ValueError:
                pd.set_option('display.max_colwidth', 0)
            output.write(lib.HEADER)
            output.write(lib.COG)
            output.write(COGSdf.transpose().to_html(
                escape=False, index=True, classes='table table-hover'))
            output.write(lib.FOOTER)
    else:
        lib.log.info("No COG annotations found")
        # create html output
        with open(os.path.join(args.out, 'cogs.html'), 'w') as output:
            output.write(lib.HEADER)
            output.write(lib.MISSING)
            output.write(lib.FOOTER)

    ############################

    ####SignalP############################
    # flip the dict and just count number for each
    signalpDict = lib.busco_dictFlip(signalp)
    # log raw data
    lib.log.debug("SignalP raw data:\n%s" % signalpDict)
    if len(signalp[0]) > 1:
        lib.log.info("Summarizing secreted protein results")

        if not os.path.isdir(os.path.join(args.out, 'signalp')):
            os.makedirs(os.path.join(args.out, 'signalp'))
        sig = {}
        for i in range(0, len(names)):
            if names[i] not in sig:
                sig[names[i]] = len(signalpDict[i])
        sigdf = pd.DataFrame([sig])
        sigdf.transpose().to_csv(os.path.join(args.out, 'signalp', 'signalp.csv'))
        lib.drawbarplot(sigdf, os.path.join(
            args.out, 'signalp', 'signalp.pdf'))

        # create html output
        with open(os.path.join(args.out, 'signalp.html'), 'w') as output:
            output.write(lib.HEADER)
            output.write(lib.SIGNALP)
            output.write(lib.FOOTER)
    else:
        lib.log.info("No SignalP annotations found")
        lib.log.debug("%s\n" % signalp)
        # create html output
        with open(os.path.join(args.out, 'signalp.html'), 'w') as output:
            output.write(lib.HEADER)
            output.write(lib.MISSING)
            output.write(lib.FOOTER)

    ########################################################

    ####Transcription Factors############################
    lib.log.info("Summarizing fungal transcription factors")
    if not os.path.isdir(os.path.join(args.out, 'tfs')):
        os.makedirs(os.path.join(args.out, 'tfs'))
    # should be able to pull transcription factor counts from InterPro Domains, load into pandas df
    iprTF = os.path.join(parentdir, 'config', 'tf_interpro.txt')
    tf = pd.read_csv(iprTF, names=['InterPro', 'Description'])
    # convert to dictionary for all annotations later
    TFDict = tf.set_index('InterPro')['Description'].to_dict()
    iprall = IPRdf.transpose()
    iprall.reset_index(inplace=True)
    dfmerged = pd.merge(tf, iprall, left_on='InterPro',
                        right_on='index', how='left')
    dfmerged.drop('index', axis=1, inplace=True)
    dfmerged.fillna(0, inplace=True)
    dfmerged.to_csv(os.path.join(args.out, 'tfs',
                                 'transcription_factor_counts.csv'))
    dfmerged['TFs'] = dfmerged['InterPro'].map(
        str) + ': ' + dfmerged['Description']
    dfmerged.drop('InterPro', axis=1, inplace=True)
    dfmerged.drop('Description', axis=1, inplace=True)
    dfmerged.set_index('TFs', inplace=True)
    dfmerged.sort_index(axis=0, inplace=True)
    dfmerged = dfmerged.astype(int)
    #print dfmerged
    if len(dfmerged.columns) > 0:
        lib.drawHeatmap(dfmerged, 'Blues', os.path.join(
            args.out, 'tfs', 'TF.heatmap.pdf'), 6, True)
    else:
        lib.log.info("No transcription factor IPR domains found")
    # create html output
    with open(os.path.join(args.out, 'tf.html'), 'w') as output:
        output.write(lib.HEADER)
        output.write(lib.TF)
        output.write(lib.FOOTER)

    ########################################################

    ####GO Terms, GO enrichment############################
    if not os.path.isdir(os.path.join(args.out, 'go_enrichment')):
        os.makedirs(os.path.join(args.out, 'go_enrichment'))

    if len(args.input) > 1:
        # concatenate all genomes into a population file
        lib.log.info("Running GO enrichment for each genome")
        with open(os.path.join(go_folder, 'population.txt'), 'w') as pop:
            for file in os.listdir(go_folder):
                if not file.startswith('associations'):
                    file = os.path.join(go_folder, file)
                    with open(file) as input:
                        pop.write(input.read())

        # pass input and output folders to runner script
        cmd = [sys.executable, os.path.join(parentdir, 'aux_scripts', 'enrichment_parallel.py'), '-i', go_folder,
               '-o', os.path.join(args.out, 'go_enrichment'), '-d', FUNDB, '-c', str(args.cpus)]
        subprocess.call(cmd)

        # load into pandas and write to html
        with open(os.path.join(args.out, 'go.html'), 'w') as output:
            try:
                pd.set_option('display.max_colwidth', None)
            except ValueError:
                pd.set_option('display.max_colwidth', 0)
            pd.options.mode.chained_assignment = None  # turn off warning
            output.write(lib.HEADER)
            output.write(lib.GO)
            for f in os.listdir(os.path.join(args.out, 'go_enrichment')):
                if f.endswith('go.enrichment.txt'):
                    file = os.path.join(args.out, 'go_enrichment', f)
                    base = os.path.basename(file)
                    name = base.split('.go_enrichment.txt')[0]
                    #check goatools output, return is a tuple with True/False and header line #
                    goresult = lib.checkgoatools(file)
                    output.write(
                        '<h4 class="sub-header" align="left">GO Enrichment: '+name+'</h4>')
                    # goatools keeps changing output - which really sucks....trying now to parse the header, hopefully that doesnt change
                    # goatools changed output, empty files now have 9 lines instead of 3...
                    if goresult[0]:
                        # the get header row from tuple
                        df = pd.read_csv(file, sep='\t', skiprows=goresult[1])
                        df.columns = df.columns.str.replace(r'^# ', '')
                        df['enrichment'].replace('p', 'under', inplace=True)
                        df['enrichment'].replace('e', 'over', inplace=True)
                        df2 = df.loc[df['p_fdr'] < args.go_fdr]
                        df2.sort_values(by='enrichment', inplace=True)
                        if len(df2) > 0:
                            df2.to_csv(os.path.join(
                                args.out, 'go_enrichment', base+'.fdr_enriched.csv'), index=False)
                            # apparently goatools also changed the headers....arrggh...
                            df2['GO'] = '<a target="_blank" href="http://amigo.geneontology.org/amigo/search/ontology?q=' + \
                                df2['GO'].astype(str)+'">'+df2['GO']+'</a>'
                            # now output has all gene names in last column, drop this for generating HTML
                            df2.drop(
                                df2.columns[len(df2.columns)-1], axis=1, inplace=True)
                            output.write(df2.to_html(
                                escape=False, index=False, classes='table table-hover'))
                        else:
                            output.write(
                                '<table border="1" class="dataframe table table-hover">\n<th>No enrichment found</th></table>')
                    else:
                        output.write(
                            '<table border="1" class="dataframe table table-hover">\n<th>No enrichment found</th></table>')
            output.write(lib.FOOTER)

    ####################################################

    ##ProteinOrtho################################
    if not os.path.isdir(os.path.join(args.out, 'annotations')):
        os.makedirs(os.path.join(args.out, 'annotations'))
    scoCount = 0
    protOrthoTSV = os.path.join(args.out, 'protortho', 'funannotate.poff.tsv')
    if len(args.input) > 1:
        if not args.proteinortho:
            lib.log.info(
                "Running orthologous clustering tool, ProteinOrtho.  This may take awhile...")
            # setup protein ortho inputs, some are a bit strange in the sense that they use equals signs
            # generate list of files based on input order for consistency
            filelist = []
            for i in scinames:
                name = i+'.faa'
                filelist.append(name)
            # run diamond blastp for reciprocal hits, then follow with proteinortho for graph/clustering
            #lib.ReciprocalBlast(filelist, protortho, args.cpus)
            # setup command
            cmd = ['proteinortho', '-project=funannotate', '-synteny',
                   '-cpus='+str(args.cpus), '-singles', '-selfblast']
            cmd2 = cmd + filelist
            if not os.path.isfile(protOrthoTSV):
                lib.runSubprocess(cmd2, protortho, lib.log)
        else:
            shutil.copyfile(args.proteinortho, protOrthoTSV)

        # open poff in pandas to parse "easier" for stats, orthologs, etc
        df = pd.read_csv(protOrthoTSV, sep='\t', header=0)
        df.rename(columns=lambda x: x.replace('.faa', ''), inplace=True)
        # reorder table to it matches up with busco list of dicts
        newhead = [df.columns.values[0],
                   df.columns.values[1], df.columns.values[2]]
        newhead += scinames
        try:
            df = df[newhead]
        # means they were not found, likely need to then drop isolate name (I hope that catches them all)
        except KeyError:
            newhead = [i.rsplit('_', 1)[0] for i in newhead]
            for x in newhead:
                if not x in df.columns.values:
                    lib.log.error(
                        "Error: %s not found in ProteinOrtho results, exiting." % x)
                    sys.exit(1)
            df = df[newhead]
        scinames = newhead[3:]
        lib.log.debug(
            "There are %i entries in the proteinortho output" % len(df))
        #print(df)
        # now filter table to only single copy orthologs to use with phylogeny
        num_species = len(df.columns) - 3
        sco = df[(df['# Species'] == num_species) & (df['Genes'] == num_species)]
        sco_hits = sco.drop(sco.columns[0:3], axis=1)
        #print(sco_hits)
        # now cross reference with busco, as we want this for phylogeny
        keep = []
        sc_buscos = []
        #print(sco_hits)
        #print(busco)
        for index, row in sco_hits.iterrows():
            busco_check = []
            for i in range(0, num_species):
                if row[i] in busco[i]:
                    busco_check.append(busco[i].get(row[i]))
            busco_check = lib.flatten(busco_check)
            #print(row)
            #print(busco_check)
            # need to check if outgroup is passed and this model exists in that outgroup
            if len(set(busco_check)) == 1:
                if args.outgroup:
                    available_busco = []
                    with open(outgroup_species, 'r') as outfasta:
                        for line in outfasta:
                            if line.startswith('>'):
                                line = line.replace('\n', '')
                                name = line.replace('>', '')
                                available_busco.append(name)
                    if busco_check[0] in available_busco:
                        keep.append(index)
                        sc_buscos.append(busco_check[0])
                else:
                    keep.append(index)
        #print(keep)
        sco_final = sco_hits.loc[keep]
        #print(sco_final)
        lib.log.debug("There seem to be %i single copy orthologs" %
                      len(sco_final))
        # take dataframe and output the ortholog table.
        # trim down to just gene models
        dftrim = df.drop(df.columns[0:3], axis=1)
        # get rid of singletons in this dataset
        orthdf = df[(df['# Species'] > 1)]
        # trim to just gene models
        orth_hits = orthdf.drop(orthdf.columns[0:3], axis=1)
        lib.log.debug("There are a total of %i orthologs groups" %
                      len(orth_hits))
        if args.run_dnds:
            # load transcripts into index, first concatenation them all, remove if file found first
            AllTrans = os.path.join(ortho_folder, 'all.transcripts.fa')
            if os.path.isfile(AllTrans):
                os.remove(AllTrans)
            with open(AllTrans, 'w') as output2:
                for x in newhead[3:]:
                    with open(os.path.join(protortho, x+'.transcripts.fa')) as traninput2:
                        output2.write(traninput2.read())
            SeqTranscripts = SeqIO.index(AllTrans, 'fasta')

        # write orthologs output
        #lib.log.info("Writing ortholog summary to file")
        orthologstmp = os.path.join(
            args.out, 'orthology', 'orthology_groups.tmp')
        with open(orthologstmp, 'w') as output:
            if args.run_dnds:
                # calculate Ka/Ks for each ortholog while grabbing the buscos
                lib.log.info(
                    "Calculating dN/dS ratios for each ortholog group, %i orthologous groups" % len(orth_hits))
            # should be able to parse the pandas ortho dataframe now
            counter = 1
            for index, row in orth_hits.iterrows():
                ID = 'orth'+str(counter)
                counter += 1
                buscos = []
                eggs = []
                proteins = []
                for x in range(0, len(row)):
                    if row[x] != '*':
                        prots = row[x].split(',')
                        for y in prots:
                            proteins.append(y)
                            egghit = eggnog[x].get(y)
                            if not egghit in eggs:
                                eggs.append(egghit)
                            buscohit = busco[x].get(y)
                            if not buscohit in buscos:
                                buscos.append(buscohit)
                # clean up the None's that get added
                eggs = [x for x in eggs if x is not None]
                buscos = [x for x in buscos if x is not None]
                buscos = lib.flatten(buscos)

                # make sure no duplicated proteins
                proteins = set(proteins)

                # now grab proteins and transcripts for each ortholog group and create folder
                if args.run_dnds:
                    if not os.path.isdir(os.path.join(ortho_folder, ID)):
                        os.makedirs(os.path.join(ortho_folder, ID))
                    KaKstranscript = os.path.join(
                        ortho_folder, ID+'.transcripts.fa')
                    with open(KaKstranscript, 'w') as tranout:
                        for i in proteins:
                            SeqIO.write(SeqTranscripts[i], tranout, 'fasta')
                # write to output
                if len(eggs) > 0:
                    eggs = ', '.join(str(v) for v in eggs)
                else:
                    eggs = 'None'
                if len(buscos) > 0:
                    buscos = set(buscos)
                    buscos = ', '.join(str(v) for v in buscos)
                else:
                    buscos = 'None'
                # write now to file
                output.write("%s\t%s\t%s\t%s\n" %
                             (ID, eggs, buscos, ', '.join(proteins)))

    if args.run_dnds:
        # multiprocessing dN/dS on list of folders
        dNdSList = lib.get_subdirs(ortho_folder)
        if args.run_dnds == 'estimate':
            lib.log.debug("Running simple dN/dS estimate")
            lib.runMultiProgress(lib.rundNdSestimate, dNdSList, args.cpus,
                                 progress=args.progress)
        else:
            lib.log.debug(
                "Running exhasitve dN/dS ratio with Likelihood Ratio Tests")
            lib.runMultiProgress(lib.rundNdSexhaustive, dNdSList, args.cpus,
                                 progress=args.progress)

        # after all data is run, then parse result log files, return dictionary
        dNdSresults = lib.parsedNdS(ortho_folder)
    if len(args.input) > 1:
        orthologs = os.path.join(args.out, 'orthology', 'orthology_groups.txt')
        with open(orthologs, 'w') as output:
            with open(orthologstmp, 'r') as input:
                for line in input:
                    line = line.replace('\n', '')
                    cols = line.split('\t')
                    if args.run_dnds:
                        if cols[0] in dNdSresults:
                            dNdS = dNdSresults.get(cols[0])
                        else:
                            dNdS = ('NC', 'NC', 'NC')
                    else:
                        dNdS = ('NC', 'NC', 'NC')
                    if args.run_dnds == 'estimate':
                        output.write("%s\t%s (NC,NC)\t%s\t%s\t%s\n" % (
                            cols[0], dNdS[0], cols[1], cols[2], cols[3]))
                    else:
                        try:
                            output.write("%s\t%s (%f,%f)\t%s\t%s\t%s\n" % (cols[0], dNdS[0], round(
                                float(dNdS[1]), 4), round(float(dNdS[2]), 4), cols[1], cols[2], cols[3]))
                        except (ValueError, AttributeError):
                            output.write("%s\t%s (NA,NA)\t%s\t%s\t%s\n" % (
                                cols[0], dNdS[0], cols[1], cols[2], cols[3]))

        # cleanup
        os.remove(orthologstmp)

    if not os.path.isdir(os.path.join(args.out, 'stats')):
        os.makedirs(os.path.join(args.out, 'stats'))
    summary = []
    # get stats, this is all single copy orthologs
    if len(args.input) > 1:
        scoCount = len(sco_hits)
        for i in range(0, len(stats)):
            orthos = 0
            singletons = 0
            for index, row in orth_hits[scinames[i]].items():
                if row != '*':
                    add = row.count(',') + 1
                    orthos += add
            for index, row in dftrim.iterrows():
                if row[scinames[i]] != '*':
                    others = []
                    for y in range(0, len(row)):
                        others.append(row[y])
                    others = set(others)
                    if len(others) == 2:
                        singletons += 1
            stats[i].append("{0:,}".format(singletons))
            stats[i].append("{0:,}".format(orthos))
            stats[i].append("{0:,}".format(scoCount))

    else:
        scoCount = 0
        singletons = 0
        orthos = 0
        stats[i].append("{0:,}".format(singletons))
        stats[i].append("{0:,}".format(orthos))
        stats[i].append("{0:,}".format(scoCount))

    for i in range(0, len(stats)):
        summary.append(stats[i])

    # convert to dataframe for easy output
    header = ['species', 'isolate', 'locus_tag', 'Assembly Size', 'Largest Scaffold', 'Average Scaffold',
              'Num Scaffolds', 'Scaffold N50', 'Percent GC', 'Num Genes', 'Num Proteins', 'Num tRNA',
              'Unique Proteins', 'Prots atleast 1 ortholog', 'Single-copy orthologs']
    df = pd.DataFrame(summary, columns=header)
    df.set_index('species', inplace=True)
    df.transpose().to_csv(os.path.join(args.out, 'stats', 'genome.stats.summary.csv'))
    with open(os.path.join(args.out, 'stats.html'), 'w') as output:
        try:
            pd.set_option('display.max_colwidth', None)
        except ValueError:
            pd.set_option('display.max_colwidth', 0)
        output.write(lib.HEADER)
        output.write(lib.SUMMARY)
        output.write(df.to_html(classes='table table-condensed'))
        output.write(lib.FOOTER)
    ############################################

    # summarize all annotation for each gene in a table
    lib.log.info("Compiling all annotations for each genome")

    # get orthology into dictionary
    orthoDict = {}
    if len(args.input) > 1:
        with open(orthologs, 'r') as input:
            for line in input:
                line = line.replace('\n', '')
                col = line.split('\t')
                genes = col[-1].split(', ')
                for i in genes:
                    orthoDict[i] = col[0]

    # get GO associations into dictionary as well
    with lib.suppress_stdout_stderr():
        goLookup = obo_parser.GODag(os.path.join(FUNDB, 'go.obo'))
    goDict = {}
    go_errors = []
    with open(os.path.join(go_folder, 'associations.txt'), 'r') as input:
        for line in input:
            line = line.replace('\n', '')
            col = line.split('\t')
            gos = col[1].split(';')
            goList = []
            for i in gos:
                try:
                    description = i+' '+goLookup[i].name
                except KeyError:
                    go_errors.append(i)
                    #print '%s not found in go.obo, try to download updated go file' % i
                    description = i
                goList.append(description)
            goDict[col[0]] = goList
    lib.log.debug('GO terms not found in go.obo: {:}'.format(
        ','.join(set(go_errors))))
    iprDict = lib.dictFlipLookup(ipr, INTERPRO)
    pfamDict = lib.dictFlipLookup(pfam, PFAM)
    meropsDict = lib.dictFlip(updatemerops)
    cazyDict = lib.dictFlip(cazy)
    TMDict = lib.busco_dictFlip(transmembrane)

    # get Transcription factors in a dictionary
    TFLookup = {}
    for k, v in list(iprDict.items()):
        for x in v:
            IPRid = x.split(':')[0]
            if IPRid in TFDict:
                TFLookup[k] = TFDict.get(IPRid)

    header = ['GeneID', 'scaffold:start-end', 'strand', 'length', 'description', 'Ortho Group', 'EggNog', 'BUSCO', 'Secreted', 'TransMembrane',
              'Protease family', 'CAZyme family', 'Transcription factor', 'InterPro Domains', 'PFAM Domains', 'GO terms', 'SecMet Cluster', 'SMCOG']
    for y in range(0, num_input):
        outputname = os.path.join(
            args.out, 'annotations', scinames[y]+'.all.annotations.tsv')
        with open(outputname, 'w') as output:
            output.write("%s\n" % ('\t'.join(header)))
            with open(gbkfilenames[y], 'r') as input:
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
                            start = f.location.nofuzzy_start
                            end = f.location.nofuzzy_end
                            chr = record.id
                            location = str(chr)+':'+str(start)+'-'+str(end)
                            if f.location.strand == 1:
                                strand = '+'
                            else:
                                strand = '-'
                            if ID in iprDict:
                                IPRdomains = "; ".join(iprDict.get(ID))
                            else:
                                IPRdomains = ''
                            if ID in signalpDict[y]:
                                signalphit = signalpDict[y].get(ID)[0]
                            else:
                                signalphit = ''
                            if ID in TMDict[y]:
                                TMhit = TMDict[y].get(ID)[0]
                            else:
                                TMhit = ''
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
                            if ID in busco[y]:
                                buscogroup = busco[y].get(ID)[0]
                            else:
                                buscogroup = ''
                            if ID in goDict:
                                goTerms = "; ".join(goDict.get(ID))
                            else:
                                goTerms = ''
                            if ID in orthoDict:
                                orthogroup = orthoDict.get(ID)
                            else:
                                orthogroup = ''
                            if ID in TFLookup:
                                transfactor = TFLookup.get(ID)
                            else:
                                transfactor = ''
                            for k, v in list(f.qualifiers.items()):
                                if k == 'note':
                                    notes = v[0].split('; ')
                                    for i in notes:
                                        if i.startswith('EggNog:'):
                                            hit = i.replace('EggNog:', '')
                                            egg = hit
                                        if i.startswith('antiSMASH:'):
                                            cluster = i.replace(
                                                'antiSMASH:', '')
                                        if i.startswith('SMCOG:'):
                                            smcog = i

                            final_result = [ID, location, strand, str(length), description, orthogroup, egg, buscogroup, signalphit,
                                            TMhit, meropsdomains, cazydomains, transfactor, IPRdomains, pfamdomains, goTerms, cluster, smcog]
                            output.write("%s\n" % ('\t'.join(final_result)))
    ############################################

    # build phylogeny
    if not os.path.isfile(os.path.join(args.out, 'phylogeny', 'ML.phylogeny.pdf')):
        if outgroup:
            num_phylogeny = len(args.input) + 1
        else:
            num_phylogeny = len(args.input)
        if num_phylogeny > 3:
            lib.log.info("Inferring phylogeny using RAxML")
            folder = os.path.join(args.out, 'protortho')
            lib.ortho2phylogeny(folder, sco_final, args.num_orthos, busco, args.cpus, args.bootstrap,
                                phylogeny, outgroup, outgroup_species, outgroup_name, sc_buscos, args.ml_method)
            with open(os.path.join(args.out, 'phylogeny.html'), 'w') as output:
                output.write(lib.HEADER)
                output.write(lib.PHYLOGENY)
                output.write(lib.FOOTER)
        else:
            lib.log.info(
                "Skipping RAxML phylogeny as at least 4 taxa are required")
            with open(os.path.join(args.out, 'phylogeny.html'), 'w') as output:
                output.write(lib.HEADER)
                output.write(lib.NOPHYLOGENY)
                output.write(lib.FOOTER)

    ###########################################
    def addlink(x):
        if x.startswith('ENOG50'):
            baseurl = 'http://eggnog5.embl.de/'
        elif x.startswith('ENOG41'):
            baseurl = 'http://eggnog45.embl.de/'
        x = '<a target="_blank" href="{}#/app/results?target_nogs={}>{}</a>'.format(baseurl, x, x)
        return x

    def addlink2(x):
        x = '<a target="_blank" href="http://www.orthodb.org/?level=&species=&query='+x+'">'+x+'</a>'
        return x

    # building remaining HTML output
    if len(args.input) > 1:
        with open(os.path.join(args.out, 'orthologs.html'), 'w') as output:
            df = pd.read_csv(orthologs, sep='\t', header=None)
            orthtable = []
            for row in df.itertuples():
                t = row[3].split(', ')  # convert Eggnog to list
                if t[0] == 'None':
                    t = ['None']
                else:
                    t = [addlink(y) for y in t]
                try:
                    value = '; '.join(t)
                except TypeError:
                    value = 'None found'
                r = row[4].split(', ')  # convert BUSCO to list
                if r[0] == 'None':
                    r = ['None']
                else:
                    r = [addlink2(y) for y in r]
                try:
                    value2 = '; '.join(r)
                except TypeError:
                    value2 = 'None found'
                final = [row[0], row[1], row[2], value, value2, row[5]]
                orthtable.append(final)
            df2 = pd.DataFrame(orthtable)
            df2.columns = ['Index', 'Orthology Group',
                           'dN/dS ratio (LRTs M1/M2, M7/M8)', 'EggNog Ref', 'BUSCOs', 'Gene Names']
            df2.set_index('Index', inplace=True)
            try:
                pd.set_option('display.max_colwidth', None)
            except ValueError:
                pd.set_option('display.max_colwidth', 0)
            output.write(lib.HEADER)
            output.write(lib.ORTHOLOGS)
            output.write(df2.to_html(index=False, escape=False,
                                     classes='table table-hover'))
            output.write(lib.FOOTER)

    else:
        with open(os.path.join(args.out, 'orthologs.html'), 'w') as output:
            output.write(lib.HEADER)
            output.write(lib.ORTHOLOGS)
            output.write(lib.MISSING)
            output.write(lib.FOOTER)

    with open(os.path.join(args.out, 'citation.html'), 'w') as output:
        output.write(lib.HEADER)
        output.write(lib.CITATION)
        output.write(lib.FOOTER)

    # make the "homepage"
    date = datetime.now()
    d = list(date.timetuple())

    if d[3] > 12:
        hour = d[3] - 12
        m = 'pm'
    else:
        hour = d[3]
        m = 'am'

    d = [str(x) for x in d]

    with open(os.path.join(args.out, 'index.html'), 'w') as output:
        output.write(lib.HEADER)
        output.write(lib.INDEX)
        output.write('<p>Report generated on: ' +
                     d[1]+'/'+d[2]+'/'+d[0] + ' at '+str(hour)+':'+d[4] + ' '+m+'</p>')
        output.write(lib.FOOTER)

    lib.log.info("Compressing results to output file: %s.tar.gz" % args.out)
    lib.make_tarfile(args.out+'.tar.gz', args.out)
    lib.log.info("Funannotate compare completed successfully!")


if __name__ == "__main__":
    main(sys.argv[1:])
