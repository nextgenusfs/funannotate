#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import funannotate.library as lib
import funannotate.resources as resources


def main(args):
    # setup funannotate DB path
    try:
        FUNDB = os.environ["FUNANNOTATE_DB"]
    except KeyError:
        print('$FUNANNOTATE_DB not found, run funannotate setup and export ENV variable')
        sys.exit(1)
    if '--show-outgroups' in args:
        print('ERROR: `funannotate database --show-outgroups` was removed in 1.10. The BUSCO outgroups belonged to\n'
              '`funannotate compare`, which moved to funannotate-compare:\n'
              '  https://github.com/hyphaltip/funannotate-compare\n'
              'List them with `funannotate-compare outgroups`; download with `funannotate-compare setup -i outgroups`.')
        sys.exit(1)
    elif '--show-buscos' in args:
        print("-----------------------------")
        print("BUSCO DB tree: (# of models)")
        print("-----------------------------")
        print((resources.buscoTree))
    else:
        dbfile = os.path.join(FUNDB, 'funannotate-db-info.txt')
        db_list = [['Database', 'Type', 'Version',
                    'Date', 'Num_Records', 'Md5checksum']]
        if not os.path.isfile(dbfile):
            print('Database is not properly configured, re-run funannotate setup')
            sys.exit(1)
        with open(dbfile, 'r') as infile:
            for line in infile:
                line = line.rstrip()
                cols = line.split('\t')
                del cols[2]
                db_list.append(cols)
        msg = lib.bold_underline('Funannotate Databases currently installed:')
        print(('\n'+msg+'\n'))
        lib.print_table(db_list, alignments='LLLLRL', max_col_width=60)

        print((
            '\nTo update a database type:\n\tfunannotate setup -i DBNAME -d {:} --force\n'.format(FUNDB)))
        print('To see BUSCO tree type:\n\tfunannotate database --show-buscos\n')


if __name__ == "__main__":
    main(sys.argv[1:])
