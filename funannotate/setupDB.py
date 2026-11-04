#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import io
import re
import subprocess
import argparse
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
import datetime
import requests
import hashlib
import shutil
import tarfile
import json
import xml.etree.cElementTree as cElementTree
import funannotate.library as lib


def calcmd5(file):
    digest = hashlib.md5()
    with open(file, 'rb') as infile:
        for chunk in iter(lambda: infile.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def calcmd5remote(url, max_file_size=100*1024*1024):
    # Use requests instead of urlopen to handle 308 redirects (osf.io);
    # urllib only follows 308 from Python 3.11
    hash = hashlib.md5()
    total_read = 0
    with requests.get(url, stream=True, allow_redirects=True) as remote:
        remote.raise_for_status()
        for data in remote.iter_content(chunk_size=4096):
            total_read += 4096
            if not data or total_read > max_file_size:
                break
            hash.update(data)
    return hash.hexdigest()


def check4newDB(name, infoDB, info_key=None):
    # check remote md5 with stored in database
    # name is the DBURL key; info_key is the infoDB key when it differs
    if info_key:
        checkname = info_key
    elif '-' in name:
        checkname = name.split('-')[0]
    else:
        checkname = name
    if checkname not in infoDB:
        lib.log.error("%s not found in database" % name)
        return True
    else:
        oldmd5 = infoDB[checkname][5]
        newmd5 = calcmd5remote(DBURL.get(name))
        lib.log.debug("%s database, Old md5: %s; New md5: %s" %
                      (name, oldmd5, newmd5))
        if oldmd5 == newmd5:
            lib.log.info("%s database is current." % name)
            return False
        else:
            lib.log.info("%s database is out of date, updating." % name)
            return True


def download(url, name, wget=False):
    if wget:
        # download with wget
        cmd = ['wget', '-O', name, '--no-check-certificate', '-t', '2', '-c', url]
        subprocess.call(cmd)
    else:
        file_name = name
        try:
            # fix bug 1138
            # Use requests instead of urlopen to handle 308 redirects - was necessary for osf.io requests
            with requests.get(url, stream=True, allow_redirects=True) as r:
                r.raise_for_status()
                file_size = int(r.headers.get('content-length', 0))
                lib.log.info("Downloading:  {0} Bytes:  {1}".format(url, file_size))
                with open(file_name, 'wb') as f:
                    file_size_dl = 0
                    block_sz = 8192
                    for chunk in r.iter_content(chunk_size=block_sz):
                        if chunk:
                            file_size_dl += len(chunk)
                            f.write(chunk)
                            if file_size > 0:
                                p = float(file_size_dl) / file_size
                                status = r"{0} [{1:.2%}]".format(file_size_dl, p)
                                status = status + chr(8)*(len(status)+1)
                                sys.stdout.write(status)
                sys.stdout.flush()
        except requests.exceptions.RequestException as e:
            lib.log.error("Download failed: {}".format(e))
            raise


def wget(url, name):
    # download with wget
    cmd = ['wget', '-O', name, '--no-check-certificate', '--tries=3', '--retry-connrefused', '-c', url]
    subprocess.call(cmd)


def meropsDB(info, force=False, args={}):
    fasta = os.path.join(FUNDB, 'meropsscan.lib')
    filtered = os.path.join(FUNDB, 'merops.formatted.fa')
    database = os.path.join(FUNDB, 'merops.dmnd')
    if os.path.isfile(fasta) and args.update and not force:
        if check4newDB('merops', info):
            force = True
    if not os.path.isfile(fasta) or force or 'merops' not in info:
        lib.log.info('Downloading Merops database')
        for x in [fasta, filtered, database]:
            if os.path.isfile(x):
                os.remove(x)
        if args.wget:
            wget(DBURL.get('merops'), fasta)
        else:
            download(DBURL.get('merops'), fasta)
        md5 = calcmd5(fasta)
        # reformat fasta headers
        with open(filtered, 'w') as filtout:
            with io.open(fasta, encoding="utf8", errors='ignore') as infile:
                for line in infile:
                    if line.startswith('>'):
                        line = line.rstrip()
                        ID = line.split()[0]
                        family = line.split('#')[1]
                        filtout.write('{:} {:}\n'.format(ID, family))
                    else:
                        filtout.write(line)
        lib.log.info('Building diamond database')
        cmd = ['diamond', 'makedb', '--in',
               'merops.formatted.fa', '--db', 'merops']
        lib.runSubprocess(cmd, os.path.join(FUNDB), lib.log)
        num_records = lib.countfasta(filtered)
        info['merops'] = ('diamond', database, '12.5',
                          '2023-01-19', num_records, md5)
    type, name, version, date, records, checksum = info.get('merops')
    lib.log.info('MEROPS Database: version={:} date={:} records={:,}'.format(
        version, date, records))


def uniprotDB(info, force=False, args={}):
    '''
    download swissprot/uniprot database, format for diamond, and output date of database
    '''
    fasta = os.path.join(FUNDB, 'uniprot_sprot.fasta')
    database = os.path.join(FUNDB, 'uniprot.dmnd')
    versionfile = os.path.join(FUNDB, 'uniprot.release-date.txt')
    if os.path.isfile(fasta) and args.update and not force:
        if check4newDB('uniprot-release', info):
            force = True
    if not os.path.isfile(fasta) or force or 'uniprot' not in info:
        lib.log.info('Downloading UniProtKB/SwissProt database')
        for x in [fasta, fasta+'.gz', versionfile, database]:
            if os.path.isfile(x):
                os.remove(x)
        if args.wget:
            wget(DBURL.get('uniprot'), fasta+'.gz')
        else:
            download(DBURL.get('uniprot'), fasta+'.gz')
        subprocess.call(
            ['gunzip', '-f', 'uniprot_sprot.fasta.gz'], cwd=os.path.join(FUNDB))
        if args.wget:
            wget(DBURL.get('uniprot-release'), versionfile)
        else:
            download(DBURL.get('uniprot-release'), versionfile)
        md5 = calcmd5(versionfile)
        unidate = ''
        univers = ''
        with io.open(versionfile, encoding="utf8", errors='ignore') as infile:
            for line in infile:
                if line.startswith('UniProtKB/Swiss-Prot Release'):
                    rest, datepart = line.split(' of ')
                    unidate = datetime.datetime.strptime(
                        datepart.rstrip(), "%d-%b-%Y").strftime("%Y-%m-%d")
                    univers = rest.split(' ')[-1]
        lib.log.info('Building diamond database')
        cmd = ['diamond', 'makedb', '--in',
               'uniprot_sprot.fasta', '--db', 'uniprot']
        lib.runSubprocess(cmd, os.path.join(FUNDB), lib.log)
        num_records = lib.countfasta(
            os.path.join(FUNDB, 'uniprot_sprot.fasta'))
        info['uniprot'] = ('diamond', database, univers,
                           unidate, num_records, md5)
    type, name, version, date, records, checksum = info.get('uniprot')
    lib.log.info('UniProtKB Database: version={:} date={:} records={:,}'.format(
        version, date, records))


def dbcan_version_info(hmm_url, readme_head, download_date):
    """
    (version, date) recorded for the dbCAN HMM file that was downloaded.

    The readme.txt next to the V11 files now describes the newest dbCAN release
    ("release 14.0"), so its version would be recorded for V11 models. The version
    comes from the V<n> folder of the HMM link; the readme's version and date are
    used only when it agrees with that folder, else the date is the download date.
    """
    readme_vers = readme_head[0].split(' ')[-1].rstrip()
    m = re.search(r'/V(\d+)/', hmm_url or '')
    if not m:
        folder_vers = readme_vers
    else:
        folder_vers = '{}.0'.format(m.group(1))
    if folder_vers == readme_vers:
        try:
            date = datetime.datetime.strptime(
                readme_head[1].replace('# ', '').rstrip(), "%m/%d/%Y").strftime("%Y-%m-%d")
            return folder_vers, date
        except ValueError:
            pass
    return folder_vers, download_date


def dbCANDB(info, force=False, args={}):
    hmm = os.path.join(FUNDB, 'dbCAN.hmm')
    familyinfo = os.path.join(FUNDB, 'dbCAN-fam-HMMs.txt')
    versionfile = os.path.join(FUNDB, 'dbCAN.changelog.txt')
    if os.path.isfile(hmm) and args.update and not force:
        if check4newDB('dbCAN', info):
            force = True
    if not os.path.isfile(hmm) or force or 'dbCAN' not in info:
        lib.log.info('Downloading dbCAN database')
        for x in [os.path.join(FUNDB, 'dbCAN.tmp'), hmm, familyinfo, versionfile]:
            if os.path.isfile(x):
                os.remove(x)
        if args.wget:
            wget(DBURL.get('dbCAN'), os.path.join(FUNDB, 'dbCAN.tmp'))
            wget(DBURL.get('dbCAN-tsv'), familyinfo)
            wget(DBURL.get('dbCAN-log'), versionfile)
        else:
            download(DBURL.get('dbCAN'),
                     os.path.join(FUNDB, 'dbCAN.tmp'))
            download(DBURL.get('dbCAN-tsv'), familyinfo)
            download(DBURL.get('dbCAN-log'), versionfile)
        md5 = calcmd5(os.path.join(FUNDB, 'dbCAN.tmp'))
        num_records = 0
        dbdate = ''
        dbvers = ''
        with open(hmm, 'w') as out:
            with io.open(os.path.join(FUNDB, 'dbCAN.tmp'), encoding="utf8", errors='ignore') as input:
                for line in input:
                    if line.startswith('NAME'):
                        num_records += 1
                        line = line.replace('.hmm\n', '\n')
                    out.write(line)
        with io.open(versionfile, encoding="utf8", errors='ignore') as infile:
            head = [next(infile) for x in range(2)]
        dbvers, dbdate = dbcan_version_info(DBURL.get('dbCAN'), head, today)
        lib.log.info('Creating dbCAN HMM database')
        cmd = ['hmmpress', '-f', 'dbCAN.hmm']
        lib.runSubprocess(cmd, os.path.join(FUNDB), lib.log)
        info['dbCAN'] = ('hmmer3', hmm, dbvers, dbdate, num_records, md5)
        os.remove(os.path.join(FUNDB, 'dbCAN.tmp'))
    type, name, version, date, records, checksum = info.get('dbCAN')
    lib.log.info('dbCAN Database: version={:} date={:} records={:,}'.format(
        version, date, records))


def pfamDB(info, force=False, args={}):
    hmm = os.path.join(FUNDB, 'Pfam-A.hmm')
    familyinfo = os.path.join(FUNDB, 'Pfam-A.clans.tsv')
    versionfile = os.path.join(FUNDB, 'Pfam.version')
    if os.path.isfile(hmm) and args.update and not force:
        if check4newDB('pfam-log', info):
            force = True
    if not os.path.isfile(hmm) or force or 'pfam' not in info:
        for x in [hmm, hmm+'.gz', familyinfo, familyinfo+'.gz', versionfile, versionfile+'.gz']:
            if os.path.isfile(x):
                os.remove(x)
        lib.log.info('Downloading Pfam database')
        if args.wget:
            wget(DBURL.get('pfam'), hmm+'.gz')
            wget(DBURL.get('pfam-tsv'), familyinfo+'.gz')
            wget(DBURL.get('pfam-log'), versionfile+'.gz')
        else:
            download(DBURL.get('pfam'), hmm+'.gz')
            download(DBURL.get('pfam-tsv'), familyinfo+'.gz')
            download(DBURL.get('pfam-log'), versionfile+'.gz')
        subprocess.call(['gunzip', '-f', 'Pfam-A.hmm.gz'],
                        cwd=os.path.join(FUNDB))
        subprocess.call(['gunzip', '-f', 'Pfam-A.clans.tsv.gz'],
                        cwd=os.path.join(FUNDB))
        md5 = calcmd5(versionfile+'.gz')
        subprocess.call(['gunzip', '-f', 'Pfam.version.gz'],
                        cwd=os.path.join(FUNDB))
        num_records = 0
        pfamdate = ''
        pfamvers = ''
        with io.open(versionfile, encoding="utf8", errors='ignore') as input:
            for line in input:
                if line.startswith('Pfam release'):
                    pfamvers = line.split(': ')[-1].rstrip()
                if line.startswith('Pfam-A families'):
                    num_records = int(line.split(': ')[-1].rstrip())
                if line.startswith('Date'):
                    pfamdate = line.split(': ')[-1].rstrip()
        lib.log.info('Creating Pfam HMM database')
        cmd = ['hmmpress', '-f', 'Pfam-A.hmm']
        lib.runSubprocess(cmd, os.path.join(FUNDB), lib.log)
        info['pfam'] = ('hmmer3', hmm, pfamvers, pfamdate,  num_records, md5)
    type, name, version, date, records, checksum = info.get('pfam')
    lib.log.info('Pfam Database: version={:} date={:} records={:,}'.format(
        version, date, records))


def repeatDB(info, force=False, args={}):
    fasta = os.path.join(FUNDB, 'funannotate.repeat.proteins.fa')
    filtered = os.path.join(FUNDB, 'funannotate.repeats.reformat.fa')
    database = os.path.join(FUNDB, 'repeats.dmnd')
    if os.path.isfile(fasta) and args.update and not force:
        if check4newDB('repeats', info):
            force = True
    if not os.path.isfile(fasta) or force or 'repeats' not in info:
        lib.log.info('Downloading Repeat database')
        for x in [fasta, fasta+'.tar.gz', filtered, database]:
            if os.path.isfile(x):
                os.remove(x)
        if args.wget:
            wget(DBURL.get('repeats'), fasta+'.tar.gz')
        else:
            download(DBURL.get('repeats'), fasta+'.tar.gz')
        md5 = calcmd5(fasta+'.tar.gz')
        subprocess.call(
            ['tar', '-zxf', 'funannotate.repeat.proteins.fa.tar.gz'], cwd=os.path.join(FUNDB))
        with open(filtered, 'w') as out:
            with io.open(fasta, encoding="utf8", errors='ignore') as infile:
                for line in infile:
                    # this repeat fasta file has messed up headers....
                    if line.startswith('>'):
                        line = line.replace('#', '_')
                        line = line.replace('/', '-')
                        line = line.replace('&', '')
                    out.write(line)
        lib.log.info('Building diamond database')
        cmd = ['diamond', 'makedb', '--in', 'funannotate.repeats.reformat.fa',
               '--db', 'repeats', '-parse_seqids']
        lib.runSubprocess(cmd, os.path.join(FUNDB), lib.log)
        num_records = lib.countfasta(filtered)
        info['repeats'] = ('diamond', database, '1.0', today, num_records, md5)
    type, name, version, date, records, checksum = info.get('repeats')
    lib.log.info('Repeat Database: version={:} date={:} records={:,}'.format(
        version, date, records))


def goDB(info, force=False, args={}):
    goOBO = os.path.join(FUNDB, 'go.obo')
    if os.path.isfile(goOBO) and args.update and not force:
        if check4newDB('go-obo', info):
            force = True
    if not os.path.isfile(goOBO) or force or 'go' not in info:
        lib.log.info('Downloading GO Ontology database')
        for x in [goOBO]:
            if os.path.isfile(x):
                os.remove(x)
        if args.wget:
            wget(DBURL.get('go-obo'), goOBO)
        else:
            download(DBURL.get('go-obo'), goOBO)
        md5 = calcmd5(goOBO)
        num_records = 0
        version = ''
        with io.open(goOBO, encoding="utf8", errors='ignore') as infile:
            for line in infile:
                if line.startswith('data-version:'):
                    version = line.split(
                        ' ')[1].rstrip().replace('releases/', '')
                if line.startswith('[Term]'):
                    num_records += 1
        info['go'] = ('text', goOBO, version, version,  num_records, md5)
    type, name, version, date, records, checksum = info.get('go')
    lib.log.info('GO ontology version={:} date={:} records={:,}'.format(
        version, date, records))


def mibigDB(info, force=False, args={}):
    fasta = os.path.join(FUNDB, 'mibig.fa')
    database = os.path.join(FUNDB, 'mibig.dmnd')
    if os.path.isfile(fasta) and args.update and not force:
        if check4newDB('mibig', info):
            force = True
    if not os.path.isfile(fasta) or force or 'mibig' not in info:
        lib.log.info('Downloading MiBIG Secondary Metabolism database')
        for x in [fasta, database]:
            if os.path.isfile(x):
                os.remove(x)
        if args.wget:
            wget(DBURL.get('mibig'), fasta)
        else:
            download(DBURL.get('mibig'), fasta)
        md5 = calcmd5(fasta)
        version = os.path.basename(DBURL.get(
            'mibig')).split('_')[-1].replace('.fasta', '')
        lib.log.info('Building diamond database')
        cmd = ['diamond', 'makedb', '--in', 'mibig.fa', '--db', 'mibig']
        lib.runSubprocess(cmd, os.path.join(FUNDB), lib.log)
        num_records = lib.countfasta(fasta)
        info['mibig'] = ('diamond', database, version, today, num_records, md5)
    type, name, version, date, records, checksum = info.get('mibig')
    lib.log.info('MiBIG Database: version={:} date={:} records={:,}'.format(
        version, date, records))


def interproDB(info, force=False, args={}):
    iprXML = os.path.join(FUNDB, 'interpro.xml')
    iprTSV = os.path.join(FUNDB, 'interpro.tsv')
    if os.path.isfile(iprXML) and args.update and not force:
        if check4newDB('interpro', info):
            force = True
    if not os.path.isfile(iprXML) or force or 'interpro' not in info:
        lib.log.info('Downloading InterProScan Mapping file')
        for x in [iprXML, iprTSV, iprXML+'.gz']:
            if os.path.isfile(x):
                os.remove(x)
        if args.wget:
            wget(DBURL.get('interpro'), iprXML+'.gz')
            wget(DBURL.get('interpro-tsv'), iprTSV)
        else:
            download(DBURL.get('interpro'), iprXML+'.gz')
            download(DBURL.get('interpro-tsv'), iprTSV)
        md5 = calcmd5(iprXML+'.gz')
        subprocess.call(['gunzip', '-f', 'interpro.xml.gz'],
                        cwd=os.path.join(FUNDB))
        num_records = ''
        version = ''
        iprdate = ''
        for event, elem in cElementTree.iterparse(iprXML):
            if elem.tag == 'release':
                for x in list(elem):
                    if x.attrib['dbname'] == 'INTERPRO':
                        num_records = int(x.attrib['entry_count'])
                        version = x.attrib['version']
                        iprdate = x.attrib['file_date']
        try:
            iprdate = datetime.datetime.strptime(
                iprdate, "%d-%b-%y").strftime("%Y-%m-%d")
        except ValueError:
            iprdate = datetime.datetime.strptime(
                iprdate, "%d-%b-%Y").strftime("%Y-%m-%d")
        info['interpro'] = ('xml', iprXML, version, iprdate, num_records, md5)
    type, name, version, date, records, checksum = info.get('interpro')
    lib.log.info('InterProScan XML: version={:} date={:} records={:,}'.format(
        version, date, records))


def curatedDB(info, force=False, args={}):
    curatedFile = os.path.join(FUNDB, 'ncbi_cleaned_gene_products.txt')
    if os.path.isfile(curatedFile) and args.update and not force:
        if check4newDB('gene2product', info):
            force = True
    if not os.path.isfile(curatedFile) or force or 'gene2product' not in info:
        lib.log.info('Downloaded curated gene names and product descriptions')
        for x in [curatedFile]:
            if os.path.isfile(x):
                os.remove(x)
        if args.wget:
            wget(DBURL.get('gene2product'), curatedFile)
        else:
            download(DBURL.get('gene2product'), curatedFile)
        md5 = calcmd5(curatedFile)
        num_records = 0
        curdate = ''
        version = ''
        with io.open(curatedFile, encoding="utf8", errors='ignore') as infile:
            for line in infile:
                if line.startswith('#version'):
                    version = line.split(' ')[-1].rstrip()
                elif line.startswith('#Date'):
                    curdate = line.split(' ')[-1].rstrip()
                else:
                    num_records += 1
        curdate = datetime.datetime.strptime(
            curdate, "%m-%d-%Y").strftime("%Y-%m-%d")
        info['gene2product'] = ('text', curatedFile,
                                version, curdate, num_records, md5)
    type, name, version, date, records, checksum = info.get('gene2product')
    lib.log.info('Gene2Product: version={:} date={:} records={:,}'.format(
        version, date, records))


# BUSCO lineages funannotate can install (`setup -i busco -b NAME`, `predict --busco_db NAME`).
# NAME is the folder created in $FUNANNOTATE_DB. The dataset inside the download has its own
# name, set in downloads.json ("busco": {NAME: [url, dataset folder in the archive]}). They
# differ for two lineages (the odb9 datasets carry the taxonomic rank): sordariomycetes ->
# sordariomyceta_odb9, saccharomycetes -> saccharomyceta_odb9; and protists and
# alveolata_stramenophiles are `_ensembl` datasets, not odb9. The archive for sordariomycetes is
# named `sordariomyceta_odb9.tar-2.gz` on the db-mirror release. setup renames the dataset
# folder to NAME, so every command refers to a lineage by NAME; the dataset name is kept in
# NAME/dataset.cfg and shown by `lib.busco_lineage_label`.
BUSCO_LINEAGES = ['fungi', 'microsporidia', 'dikarya', 'ascomycota', 'pezizomycotina',
                  'eurotiomycetes', 'sordariomycetes', 'saccharomycetes', 'saccharomycetales',
                  'basidiomycota', 'eukaryota', 'protists', 'alveolata_stramenophiles', 'metazoa',
                  'nematoda', 'arthropoda', 'insecta', 'endopterygota', 'hymenoptera', 'diptera',
                  'vertebrata', 'actinopterygii', 'tetrapoda', 'aves', 'mammalia', 'euarchontoglires',
                  'laurasiatheria', 'embryophyta']


# OrthoDB 10 lineages (`<lineage>_odb10`) come from the BUSCO download server, not from the
# funannotate mirror. file_versions.tsv lists every dataset as
#   name <TAB> date <TAB> md5 <TAB> domain <TAB> type
# and the archive is lineages/<name>.<date>.tar.gz (it unpacks to a folder named <name>).
# The server redirects to S3, so downloads must follow redirects.
BUSCO_DATA_URL = 'https://busco-data.ezlab.org/v5/data/'
BUSCO_INDEX_URL = BUSCO_DATA_URL + 'file_versions.tsv'
ODB10_RE = re.compile(r'^[a-z0-9]+(?:_[a-z0-9]+)*_odb10$')
ODB12_RE = re.compile(r'_odb12(?:\.\d+)?$')


def busco_name_problem(name):
    """
    None when `name` can be installed (`all`, an odb9 lineage of downloads.json, or a
    `<lineage>_odb10` name), otherwise a message saying why not.
    """
    if name == 'all' or name in BUSCO_LINEAGES or ODB10_RE.match(name):
        return None
    if ODB12_RE.search(name):
        return ('{0}: OrthoDB 12 datasets cannot be used by funannotate\'s BUSCO (they have no '
                'lengths_cutoff file, and their HMM names do not match their score table). '
                'Use an _odb10 lineage such as fungi_odb10, or one of the odb9 lineages.'.format(name))
    return ('{0}: not a BUSCO lineage funannotate can install. Use one of: all, {1}, or a '
            '<lineage>_odb10 name such as fungi_odb10.'.format(name, ', '.join(BUSCO_LINEAGES)))


def parse_file_versions(text):
    """{dataset name: (date, md5, domain, type)} from the BUSCO file_versions.tsv."""
    index = {}
    for line in text.splitlines():
        cols = line.rstrip('\n').split('\t')
        if len(cols) >= 5 and cols[0]:
            index[cols[0]] = (cols[1], cols[2], cols[3], cols[4])
    return index


def odb10_lineage(name, index):
    """(url, md5) of the archive of the eukaryotic odb10 lineage `name`; ValueError if there is none."""
    if name not in index:
        raise ValueError('{0} is not in the BUSCO dataset list (file_versions.tsv); see '
                         '`busco --list-datasets` for the names'.format(name))
    date, md5, domain, kind = index[name]
    if kind != 'lineages' or domain != 'Eukaryota':
        raise ValueError('{0} is a {1} {2} dataset; funannotate annotates eukaryotes'.format(
            name, domain, kind))
    return '{0}lineages/{1}.{2}.tar.gz'.format(BUSCO_DATA_URL, name, date), md5


def fetch_text(url):
    response = requests.get(url, allow_redirects=True, timeout=120)
    response.raise_for_status()
    return response.text


def safe_extract(archive, dest):
    """Extract a tar.gz into `dest`, refusing members that would land outside it."""
    root = os.path.realpath(dest)
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar.getmembers():
            target = os.path.realpath(os.path.join(dest, member.name))
            if not (target == root or target.startswith(root + os.sep)):
                raise RuntimeError('unsafe path in archive {0}: {1}'.format(archive, member.name))
        tar.extractall(dest)


def busco_info_row(target, md5):
    """funannotate-db-info.txt row for an installed lineage: (type, path, dataset, created, BUSCOs, md5)."""
    cfg = lib.busco_lineage_info(target)
    return ('busco', target, cfg.get('name', os.path.basename(target)),
            cfg.get('creation_date', ''), cfg.get('number_of_BUSCOs', 0), md5)


def install_odb10_lineage(name, fundb, index, fetch, log=None):
    """
    Download, verify and install `<lineage>_odb10` into $FUNANNOTATE_DB/<name>.

    `fetch(url, path)` saves a URL to a file (following redirects). Returns the md5 of the
    archive; raises ValueError (unknown dataset) or RuntimeError (checksum, archive layout).
    """
    log = log or lib.log
    url, md5 = odb10_lineage(name, index)
    archive = os.path.join(fundb, name + '.tar.gz')
    if os.path.isfile(archive):
        os.remove(archive)
    fetch(url, archive)
    found = calcmd5(archive)
    if found != md5:
        os.remove(archive)
        raise RuntimeError('md5 of the {0} download is {1}, file_versions.tsv says {2}; '
                           'nothing was installed'.format(name, found, md5))
    safe_extract(archive, fundb)
    os.remove(archive)
    finish_busco_install(fundb, name, name, log)
    return md5


def finish_busco_install(fundb, name, dataset_folder, log=None):
    """
    Move the extracted dataset folder to $FUNANNOTATE_DB/NAME and check what was installed.

    Raises RuntimeError naming the folder that was expected when the archive did not
    contain it. Returns the dataset.cfg provenance (see lib.busco_lineage_info).
    """
    log = log or lib.log
    extracted = os.path.join(fundb, dataset_folder)
    target = os.path.join(fundb, name)
    if not os.path.isdir(extracted):
        raise RuntimeError(
            'BUSCO archive for {0} did not contain the folder "{1}" expected from '
            'downloads.json; nothing was installed'.format(name, dataset_folder))
    if extracted != target:
        os.rename(extracted, target)
    info = lib.busco_lineage_info(target)
    if not info:
        log.warning('BUSCO lineage {0}: no dataset.cfg in {1}'.format(name, target))
    elif info.get('name') != dataset_folder:
        log.warning('BUSCO lineage {0}: dataset.cfg names "{1}", downloads.json expects "{2}"'.format(
            name, info.get('name'), dataset_folder))
    log.info('BUSCO lineage {0}: dataset {1}'.format(name, lib.busco_lineage_label(target)))
    return info


def download_buscos(name, force=False, args={}, info=None):
    # name is a list of lineage names (see busco_name_problem)
    info = {} if info is None else info
    for i in name:
        problem = busco_name_problem(i)
        if problem:
            lib.log.error(problem)
            sys.exit(1)
    if 'all' in name:
        installList = list(BUSCO_LINEAGES)
        lib.log.info("Downloading all %i busco models" % len(installList))
    else:
        installList = name
        lib.log.info("Downloading busco models: %s" % ', '.join(installList))
    odb10_index = None
    for i in installList:
        target = os.path.join(FUNDB, i)
        if os.path.isdir(target) and not force:
            continue
        if os.path.isdir(target):
            shutil.rmtree(target)
        try:
            if ODB10_RE.match(i):
                if odb10_index is None:
                    odb10_index = parse_file_versions(fetch_text(BUSCO_INDEX_URL))
                md5 = install_odb10_lineage(
                    i, FUNDB, odb10_index,
                    lambda url, path: wget(url, path) if args.wget else download(url, path))
            elif i in busco_links:
                address = busco_links.get(i)[0]
                filename = os.path.join(FUNDB, i+'.tar.gz')
                foldername = os.path.join(FUNDB, busco_links.get(i)[1])
                if os.path.isfile(filename):
                    os.remove(filename)
                if os.path.isdir(foldername):
                    shutil.rmtree(foldername)
                if args.wget:
                    wget(address, filename)
                else:
                    download(address, filename)
                md5 = calcmd5(filename)
                safe_extract(filename, FUNDB)
                os.remove(filename)
                finish_busco_install(FUNDB, i, busco_links.get(i)[1])
            else:
                lib.log.error('No download link for BUSCO lineage "{0}" in downloads.json'.format(i))
                continue
        except (ValueError, RuntimeError, requests.exceptions.RequestException) as err:
            lib.log.error('BUSCO lineage {0}: {1}'.format(i, err))
            sys.exit(1)
        info['busco-' + i] = busco_info_row(target, md5)


def training_species(force):
    # copy over Augustus training data and generate JSON description file
    augustus_list = []
    for i in os.listdir(os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], 'species')):
        if not i.startswith('.'):
            augustus_list.append(os.path.join(
                os.environ["AUGUSTUS_CONFIG_PATH"], 'species', i))
    augustus_list = set(augustus_list)
    version = subprocess.Popen(['augustus', '--version'], stderr=subprocess.STDOUT,
                               stdout=subprocess.PIPE).communicate()[0].rstrip().decode('utf-8')
    version = version.split(' is ')[0]
    curdate = datetime.date.today().strftime("%Y-%m-%d")
    destDir = os.path.join(FUNDB, 'trained_species')
    if not os.path.isdir(destDir):
        os.makedirs(destDir)
    for sp_path in augustus_list:
        sp = os.path.basename(sp_path)
        spDir = os.path.join(destDir, sp)
        spAugDir = os.path.join(spDir, 'augustus')
        paramFile = os.path.join(spDir, 'info.json')
        if force:
            if os.path.isfile(paramFile):
                os.remove(paramFile)
            if os.path.isdir(spAugDir):
                shutil.rmtree(spAugDir)
        if not os.path.isdir(spAugDir):
            os.makedirs(spAugDir)
        if not os.path.isdir(spDir):
            os.makedirs(spDir)
        if os.path.isfile(paramFile):
            with open(paramFile) as json_file:
                data = json.load(json_file)
        else:
            data = {'augustus': [{'version': version,
                                  'source': 'augustus pre-trained',
                                  'date': curdate,
                                  'path': os.path.abspath(spAugDir)}],
                    'genemark': [{}],
                    'codingquarry': [{}],
                    'snap': [{}],
                    'glimmerhmm': [{}],
                    }
            with open(paramFile, 'w') as outfile:
                json.dump(data, outfile)
        # move augustus files
        src_files = os.listdir(sp_path)
        for f in src_files:
            ff = os.path.realpath(os.path.join(sp_path, f))
            if os.path.isfile(ff) and not os.path.isfile(os.path.join(spAugDir, f)):
                shutil.copyfile(ff, os.path.join(spAugDir, f))


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='funannotate-setup.py', usage="%(prog)s [options] -m all -d path/to/database",
                                     description='''Download/setup databases for funannotate''',
                                     epilog="""Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--install', nargs='+', default=['all'],
                        choices=['all', 'merops', 'uniprot', 'dbCAN', 'pfam', 'repeats', 'go', 'mibig',
                                 'interpro', 'busco_outgroups', 'gene2product', 'busco'],
                        help='Databases to download/install')
    parser.add_argument('-d', '--database', help='Path to database')
    parser.add_argument('-u', '--update', action='store_true',
                        help='Check if new DB is availabe and update')
    parser.add_argument('-l', '--local', action='store_true',
                        help='Use local json links')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite current database')
    parser.add_argument('-w', '--wget', action='store_true',
                        help='Use wget for download instead of python')
    parser.add_argument('-b', '--busco_db', default=['dikarya'], nargs='+',
                        metavar='LINEAGE',
                        help='BUSCO lineages to install: all, an odb9 lineage (%s), or a <lineage>_odb10 name such as fungi_odb10' % ', '.join(BUSCO_LINEAGES))
    args = parser.parse_args(args)
    for lineage in args.busco_db:
        problem = busco_name_problem(lineage)
        if problem:
            parser.error('argument -b/--busco_db: ' + problem)

    # create log file
    log_name = 'funannotate-setup.log'
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

    # look for environmental variable if -d not passed
    global FUNDB
    if args.database:
        FUNDB = args.database
    else:
        try:
            FUNDB = os.environ["FUNANNOTATE_DB"]
        except KeyError:
            lib.log.error(
                '$FUNANNOTATE_DB variable not found, specify DB location with -d,--database option')
            sys.exit(1)
    lib.log.info("Database location: %s" % FUNDB)

    # create directory if doesn't exist
    if not os.path.isdir(FUNDB):
        os.makedirs(FUNDB)

    global today
    today = datetime.datetime.today().strftime('%Y-%m-%d')

    installdbs = []
    if 'all' in args.install:
        installdbs = ['merops', 'uniprot', 'dbCAN', 'pfam', 'repeats', 'go',
                      'mibig', 'interpro', 'gene2product', 'busco']
    else:
        installdbs = args.install
    if 'busco_outgroups' in installdbs:
        lib.log.error('The BUSCO outgroups database is no longer installed by funannotate setup. '
                      'It belonged to `funannotate compare`, which moved to funannotate-compare '
                      '(https://github.com/hyphaltip/funannotate-compare): '
                      'run `funannotate-compare setup -i outgroups`.')
        sys.exit(1)

    # load download links from gitlab if possible, need to change this when merge into master
    if not args.local:
        try:
            lib.log.info('Retrieving download links from GitHub Repo')
            response = json.loads(requests.get("https://raw.githubusercontent.com/nextgenusfs/funannotate/master/funannotate/downloads.json").text)
        except:
            lib.log.error('Unable to download links from GitHub, using funannotate version specific links')
            with open(os.path.join(os.path.dirname(__file__), 'downloads.json')) as infile:
                response = json.load(infile)
    else:
        with open(os.path.join(os.path.dirname(__file__), 'downloads.json')) as infile:
            response = json.load(infile)

    global DBURL
    global busco_links
    DBURL = response['downloads']
    busco_links = response['busco']

    # if text file with DB info is in database folder, parse into Dictionary
    DatabaseFile = os.path.join(FUNDB, 'funannotate-db-info.txt')
    DatabaseInfo = {}
    if os.path.isfile(DatabaseFile):
        with open(DatabaseFile, 'r') as inDB:
            for line in inDB:
                line = line.rstrip()
                try:
                    db, type, name, version, date, records, md5checksum = line.split(
                        '\t')
                    DatabaseInfo[db] = (type, name, version,
                                        date, int(records), md5checksum)
                except ValueError:
                    pass

    if args.update and not args.force:
        lib.log.info("Checking for newer versions of database files")

    # install Augustus species into DB
    lib.log.info(
        'Parsing Augustus pre-trained species and porting to funannotate')
    training_species(args.force)

    for x in installdbs:
        if x == 'uniprot':
            uniprotDB(DatabaseInfo, args.force, args=args)
        elif x == 'merops':
            meropsDB(DatabaseInfo, args.force, args=args)
        elif x == 'dbCAN':
            dbCANDB(DatabaseInfo, args.force, args=args)
        elif x == 'pfam':
            pfamDB(DatabaseInfo, args.force, args=args)
        elif x == 'repeats':
            repeatDB(DatabaseInfo, args.force, args=args)
        elif x == 'go':
            goDB(DatabaseInfo, args.force, args=args)
        elif x == 'interpro':
            interproDB(DatabaseInfo, args.force, args=args)
        elif x == 'mibig':
            mibigDB(DatabaseInfo, args.force, args=args)
        elif x == 'gene2product':
            curatedDB(DatabaseInfo, args.force, args=args)
        elif x == 'busco':
            download_buscos(args.busco_db, args.force, args=args, info=DatabaseInfo)

    # output the database text file and print to terminal
    with open(DatabaseFile, 'w') as outDB:
        for k, v in list(DatabaseInfo.items()):
            data = '%s\t%s\t%s\t%s\t%s\t%i\t%s' % (
                k, v[0], v[1], v[2], v[3], v[4], v[5])
            outDB.write('{:}\n'.format(data))
    if args.database:
        lib.log.info('Funannoate setup complete. Add this to ~/.bash_profile or ~/.bash_aliases:\n\n\texport FUNANNOTATE_DB={:}\n'.format(
            os.path.abspath(args.database)))


if __name__ == "__main__":
    main(sys.argv[1:])
