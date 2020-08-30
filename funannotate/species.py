#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from natsort import natsorted
import json
import shutil
import funannotate.library as lib


def speciesAvailable(dir):
    # return dictionary of species name and path to info.json file
    Results = {}
    for f in os.listdir(dir):
        ff = os.path.join(dir, f)
        if os.path.isdir(ff) and lib.checkannotations(os.path.join(ff, 'info.json')):
            with open(os.path.join(ff, 'info.json')) as infile:
                data = json.load(infile)
            Results[f] = data
    return Results


def showAll(dir):
    Table = []
    TableHeader = ['Species', 'Augustus', 'GeneMark',
                   'Snap', 'GlimmerHMM', 'CodingQuarry', 'Date']
    for f in os.listdir(dir):
        ff = os.path.join(dir, f)
        if os.path.isdir(ff) and lib.checkannotations(os.path.join(ff, 'info.json')):
            with open(os.path.join(ff, 'info.json')) as infile:
                data = json.load(infile)
            sources = [f]
            for x in ['augustus', 'genemark', 'snap', 'glimmerhmm', 'codingquarry']:
                if x in data:
                    if len(data[x][0]) < 1:
                        sources.append('None')
                    else:
                        sourceFile = data[x][0]['source']
                        if ': ' in sourceFile:
                            sourceFile = sourceFile.split(':')[0]
                        sources.append(sourceFile)
            sources.append(data['augustus'][0]['date'])
        Table.append(sources)
    Table = natsorted(Table, key=lambda x: x[0])
    Table.insert(0, TableHeader)
    lib.print_table(Table, max_col_width=40)


def copyDir(src, dest):
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print(('Directory not copied. Error: %s' % e))
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print(('Directory not copied. Error: %s' % e))


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='species.py',
                                     description='''Script to show/update species training parameters.''',
                                     epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-s', '--species', help='Species name to show/update')
    parser.add_argument('-a', '--add', '--add-parameters',
                        dest='add', help='Parameter JSON file to add to database')
    parser.add_argument('-p', '--parameters', dest='parameters',
                        help='Parameter JSON file to add to database')
    parser.add_argument('-d', '--database',
                        help='Path to funannotate database, $FUNANNOTATE_DB')
    args = parser.parse_args(args)

    # setup funannotate DB path
    if args.database:
        FUNDB = args.database
    else:
        try:
            FUNDB = os.environ["FUNANNOTATE_DB"]
        except KeyError:
            print('Funannotate database not properly configured, run funannotate setup.')
            sys.exit(1)

    # process input here
    if args.parameters:  # just pretty-print JSON file
        with open(args.parameters) as input:
            table = json.load(input)
        print((json.dumps(table, indent=3)))
    elif args.species and args.add:  # have one to add to database
        SpFound = speciesAvailable(os.path.join(FUNDB, 'trained_species'))
        if not os.access(os.path.join(FUNDB, 'trained_species'), os.W_OK | os.X_OK):
            print(('ERROR: you do not have permissions to write to {:}'.format(
                os.path.join(FUNDB, 'trained_species'))))
            sys.exit(1)
        if args.species in SpFound:
            print(('ERROR: {:} is already in database, choose a different name or delete existing to continue'.format(
                args.species)))
            sys.exit(1)
        print(('Adding {:} to Database'.format(args.species)))
        newLoc = os.path.abspath(os.path.join(
            FUNDB, 'trained_species', args.species))
        if not os.path.isdir(newLoc):
            os.makedirs(newLoc)
        with open(args.add) as infile:
            data = json.load(infile)
        for x in data:
            if not 'path' in data[x][0]:
                continue
            newPath = os.path.join(
                newLoc, os.path.basename(data[x][0]['path']))
            if os.path.isdir(data[x][0]['path']):
                copyDir(data[x][0]['path'], newPath)
            elif os.path.isfile(data[x][0]['path']):
                shutil.copyfile(data[x][0]['path'], newPath)
            data[x][0]['path'] = os.path.abspath(newPath)
        # print new data to terminal
        print(('Following training data added for {:}'.format(args.species)))
        print((json.dumps(data, indent=3)))
        with open(os.path.join(newLoc, 'info.json'), 'w') as outfile:
            json.dump(data, outfile)

    elif args.species:  # look for in database and pretty-print JSON file
        SpFound = speciesAvailable(os.path.join(FUNDB, 'trained_species'))
        if args.species in SpFound:
            print((json.dumps(SpFound[args.species], indent=3)))
        else:
            print(('{:} not found in Funannotate trained species folder'.format(
                args.species)))
            print('Valid species are:')
            showAll(os.path.join(FUNDB, 'trained_species'))
    else:
        # just show all available species in the database and their training data
        showAll(os.path.join(FUNDB, 'trained_species'))
        # row_str = colour(row_str, header_format)
        print('\n')
        print((lib.colour('Options for this script:', 'bold')))
        print((lib.colour(' To print a parameter file to terminal:', 'none')))
        print((lib.colour('   funannotate species -p myparameters.json', 'dim')))
        print((lib.colour(
            ' To print the parameters details from a species in the database:', 'none')))
        print((lib.colour('   funannotate species -s aspergillus_fumigatus', 'dim')))
        print((lib.colour(' To add a new species to database:', 'none')))
        print((lib.colour(
            '   funannotate species -s new_species_name -a new_species_name.parameters.json\n', 'dim')))


if __name__ == "__main__":
    main(sys.argv[1:])
