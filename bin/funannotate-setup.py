#!/usr/bin/env python

import sys, os, subprocess, inspect, argparse, urllib2, datetime
import xml.etree.cElementTree as cElementTree
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='funannotate-setup.py', usage="%(prog)s [options] -m all -d path/to/database",
    description = '''Download/setup databases for funannotate''',
    epilog = """Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--install', nargs='+', default=['all'], choices=['all', 'merops', 'uniprot', 'dbCAN', 'pfam', 'repeats', 'go', 'mibig', 'interpro', 'busco_outgroups', 'gene2product'], help='Databases to download/install')
parser.add_argument('-d', '--database', required=True, help='Path to database')
parser.add_argument('-f', '--force', action='store_true', help='Overwrite current database')
args=parser.parse_args()

URL = { 'uniprot_sprot': 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',
        'uniprot-release': 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt',
        'merops': 'ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/meropsscan.lib',
        'dbCAN': 'http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt',
        'dbCAN-tsv': 'http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt',
        'dbCAN-log': 'http://csbl.bmb.uga.edu/dbCAN/download/readme.txt',
        'pfam': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam-A.hmm.gz',
        'pfam-tsv': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam-A.clans.tsv.gz',
        'pfam-log': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam.version.gz',
        'outgroups': 'https://uwmadison.box.com/shared/static/4pl3ngptpjjfs1cu4se6g27ei0wptsdt.gz',
        'repeats': 'https://uwmadison.box.com/shared/static/vcftxq6yuzc3u1nykiahxcqzk3jlvyzx.gz',
        'go-obo': 'http://purl.obolibrary.org/obo/go.obo', 
        'mibig': 'http://mibig.secondarymetabolites.org/MIBiG_prot_seqs_1.3.fasta',
        'interpro': 'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz',
        'gene2product': 'https://raw.githubusercontent.com/nextgenusfs/gene2product/master/ncbi_cleaned_gene_products.txt'}

def download(url, name):
    file_name = name
    try:
        u = urllib2.urlopen(url)
        f = open(file_name, 'wb')
        meta = u.info()
        file_size = int(meta.getheaders("Content-Length")[0])
        lib.log.info("Downloading: {0} Bytes: {1}".format(url, file_size))
        file_size_dl = 0
        block_sz = 8192
        while True:
            buffer = u.read(block_sz)
            if not buffer:
                break
            file_size_dl += len(buffer)
            f.write(buffer)
            p = float(file_size_dl) / file_size
            status = r"{0}  [{1:.2%}]".format(file_size_dl, p)
            status = status + chr(8)*(len(status)+1)
            sys.stdout.write(status)
        sys.stdout.flush()
        f.close()
    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            raise
        pass

def meropsDB(info, force=False):
    fasta = os.path.join(args.database, 'merops_scan.lib')
    filtered = os.path.join(args.database, 'merops.formatted.fa')
    database = os.path.join(args.database, 'merops.dmnd')
    if not os.path.isfile(fasta) or force:
        lib.log.info('Downloading Merops database')
        download(URL.get('merops'), fasta)
        #reformat fasta headers
        with open(filtered, 'w') as filtout:
            with open(fasta, 'rU') as infile:
                for line in infile:
                    if line.startswith('>'):
                        line = line.rstrip()
                        ID = line.split()[0]
                        family = line.split('#')[1]
                        filtout.write('{:} {:}\n'.format(ID, family))
                    else:
                        filtout.write(line)
        lib.log.info('Building diamond database')
        cmd = ['diamond', 'makedb', '--in', 'merops.formatted.fa', '--db', 'merops']
        lib.runSubprocess(cmd, os.path.join(args.database), lib.log)
        num_records = lib.countfasta(filtered)
        info['merops'] = ('diamond', database, '12.0', '2017-10-04', num_records)
    type, name, version, date, records = info.get('merops')
    lib.log.info('MEROPS Database: version={:} date={:} records={:,}'.format(version, date, records))

def uniprotDB(info, force=False):
    '''
    download swissprot/uniprot database, format for diamond, and output date of database
    '''
    fasta = os.path.join(args.database, 'uniprot_sprot.fasta')
    database = os.path.join(args.database, 'uniprot.dmnd')
    versionfile = os.path.join(args.database, 'uniprot.release-date.txt')
    if not os.path.isfile(fasta) or force:
        lib.log.info('Downloading UniProtKB/SwissProt database')
        download(URL.get('uniprot_sprot'), fasta+'.gz')
        subprocess.call(['gunzip', '-f', 'uniprot_sprot.fasta.gz'], cwd=os.path.join(args.database))
        download(URL.get('uniprot-release'), versionfile)
        unidate = None
        univers = None
        with open(versionfile, 'rU') as infile:
            for line in infile:
                if line.startswith('UniProtKB/Swiss-Prot Release'):
                    rest, datepart = line.split(' of ')
                    unidate = datetime.datetime.strptime(datepart.rstrip(), "%d-%b-%Y").strftime("%Y-%m-%d") 
                    univers = rest.split(' ')[-1]
        lib.log.info('Building diamond database')
        cmd = ['diamond', 'makedb', '--in', 'uniprot_sprot.fasta', '--db', 'uniprot']
        lib.runSubprocess(cmd, os.path.join(args.database), lib.log)
        num_records = lib.countfasta(os.path.join(args.database, 'uniprot_sprot.fasta'))
        info['uniprot'] = ('diamond', database, univers, unidate, num_records)
    type, name, version, date, records = info.get('uniprot')
    lib.log.info('UniProtKB Database: version={:} date={:} records={:,}'.format(version, date, records))
                
def dbCANDB(info, force=False):
    hmm = os.path.join(args.database, 'dbCAN.hmm')
    familyinfo = os.path.join(args.database, 'dbCAN-fam-HMMs.txt')
    versionfile = os.path.join(args.database, 'dbCAN.changelog.txt')
    if not os.path.isfile(hmm) or force:
        lib.log.info('Downloading dbCAN database')
        download(URL.get('dbCAN'), os.path.join(args.database,'dbCAN.tmp'))
        download(URL.get('dbCAN-tsv'), familyinfo)
        download(URL.get('dbCAN-log'), versionfile)
        num_records = 0
        dbdate = None
        dbvers = None
        with open(hmm, 'w') as out:
            with open(os.path.join(args.database,'dbCAN.tmp'), 'rU') as input:
                for line in input:
                    if line.startswith('NAME'):
                        num_records += 1
                        line = line.replace('.hmm\n', '\n')
                    out.write(line)
        with open(versionfile, 'rU') as infile:
            head = [next(infile) for x in xrange(2)]
        dbdate = head[1].replace('# ', '').rstrip()
        dbvers = head[0].split(' ')[-1].rstrip()
        dbdate = datetime.datetime.strptime(dbdate, "%m/%d/%Y").strftime("%Y-%m-%d") 
        lib.log.info('Creating dbCAN HMM database')
        cmd = ['hmmpress', 'dbCAN.hmm']
        lib.runSubprocess(cmd, os.path.join(args.database), lib.log)
        info['dbCAN'] = ('hmmer3', hmm, dbvers, dbdate, num_records)
        os.remove(os.path.join(args.database,'dbCAN.tmp'))
    type, name, version, date, records = info.get('dbCAN')
    lib.log.info('dbCAN Database: version={:} date={:} records={:,}'.format(version, date, records))

    
def pfamDB(info, force=False):
    hmm = os.path.join(args.database, 'Pfam-A.hmm')
    familyinfo = os.path.join(args.database, 'Pfam-A.clans.tsv')
    versionfile = os.path.join(args.database, 'Pfam.version')
    if not os.path.isfile(hmm) or force:
        lib.log.info('Downloading Pfam database')
        download(URL.get('pfam'), hmm+'.gz')
        subprocess.call(['gunzip', '-f', 'Pfam-A.hmm.gz'], cwd=os.path.join(args.database))
        download(URL.get('pfam-tsv'), familyinfo+'.gz')
        subprocess.call(['gunzip', '-f', 'Pfam-A.clans.tsv.gz'], cwd=os.path.join(args.database))
        download(URL.get('pfam-log'), versionfile+'.gz')
        subprocess.call(['gunzip', '-f', 'Pfam.version.gz'], cwd=os.path.join(args.database))
        num_records = 0
        pfamdate = None
        pfamvers = None
        with open(versionfile, 'rU') as input:
            for line in input:
                if line.startswith('Pfam release'):
                    pfamvers = line.split(': ')[-1].rstrip()
                if line.startswith('Pfam-A families'):
                    num_records = int(line.split(': ')[-1].rstrip())
                if line.startswith('Date'):
                    pfamdate = line.split(': ')[-1].rstrip()
        lib.log.info('Creating Pfam HMM database')
        cmd = ['hmmpress', 'Pfam-A.hmm']
        lib.runSubprocess(cmd, os.path.join(args.database), lib.log)
        info['pfam'] = ('hmmer3', hmm, pfamvers, pfamdate,  num_records)
    type, name, version, date, records = info.get('pfam')
    lib.log.info('Pfam Database: version={:} date={:} records={:,}'.format(version, date, records))

def repeatDB(info, force=False):
    fasta = os.path.join(args.database, 'funannotate.repeat.proteins.fa')
    filtered = os.path.join(args.database, 'funannotate.repeats.reformat.fa')
    database = os.path.join(args.database, 'repeats.dmnd')
    if not os.path.isfile(fasta) or force:
        lib.log.info('Downloading Repeat database')
        download(URL.get('repeats'), fasta+'.tar.gz')
        subprocess.call(['tar', '-zxf', 'funannotate.repeat.proteins.fa.tar.gz'], cwd=os.path.join(args.database))
        with open(filtered, 'w') as out:
            with open(fasta, 'rU') as infile:
                for line in infile:
                    #this repeat fasta file has messed up headers....
                    if line.startswith('>'):
                        line = line.replace('#', '_')
                        line = line.replace('/', '-')
                        line = line.replace('&', '')
                    out.write(line)
        lib.log.info('Building diamond database')
        cmd = ['diamond', 'makedb', '--in', 'funannotate.repeats.reformat.fa', '--db', 'repeats', '-parse_seqids']
        lib.runSubprocess(cmd, os.path.join(args.database), lib.log)
        num_records = lib.countfasta(filtered)
        info['repeats'] = ('diamond', database, '1.0', today, num_records)
    type, name, version, date, records = info.get('repeats')
    lib.log.info('Repeat Database: version={:} date={:} records={:,}'.format(version, date, records))
        
def outgroupsDB(info, force=False):
    OutGroups = os.path.join(args.database, 'outgroups')
    if not os.path.isdir(OutGroups):
        lib.log.info('Downloading pre-computed BUSCO outgroups')
        download(URL.get('outgroups'), os.path.join(args.database, 'busco_outgroups.tar.gz'))
        subprocess.call(['tar', '-zxf', 'busco_outgroups.tar.gz'], cwd=os.path.join(args.database))
        num_records = len([name for name in os.listdir(OutGroups) if os.path.isfile(os.path.join(OutGroups, name))])
        info['busco_outgroups'] = ('outgroups', OutGroups, '1.0', today,  num_records)
    type, name, version, date, records = info.get('merops')
    lib.log.info('BUSCO outgroups: version={:} date={:} records={:,}'.format(version, date, records))
        
def goDB(info, force=False):
    goOBO = os.path.join(args.database, 'go.obo')
    if not os.path.isfile(goOBO) or force:
        lib.log.info('Downloading GO Ontology database')
        download(URL.get('go-obo'), goOBO)
        num_records = 0
        version = None
        with open(goOBO, 'rU') as infile:
            for line in infile:
                if line.startswith('data-version:'):
                    version = line.split(' ')[1].rstrip().replace('releases/', '')
                if line.startswith('[Term]'):
                    num_records += 1
        info['go'] = ('text', goOBO, version, version,  num_records)
    type, name, version, date, records = info.get('go')
    lib.log.info('MEROPS Database: version={:} date={:} records={:,}'.format(version, date, records))
        
def mibigDB(info, force=False):
    fasta = os.path.join(args.database, 'mibig.fa')
    database = os.path.join(args.database, 'mibig.dmnd')
    if not os.path.isfile(fasta) or force:
        lib.log.info('Downloading MiBIG Secondary Metabolism database')
        download(URL.get('mibig'), fasta)
        version = os.path.basename(URL.get('mibig')).split('_')[-1].replace('.fasta', '')
        lib.log.info('Building diamond database')
        cmd = ['diamond', 'makedb', '--in', 'mibig.fa', '--db', 'mibig']
        lib.runSubprocess(cmd, os.path.join(args.database), lib.log)
        num_records = lib.countfasta(fasta)
        info['mibig'] = ('diamond', database, version, today, num_records)
    type, name, version, date, records = info.get('mibig')
    lib.log.info('MiBIG Database: version={:} date={:} records={:,}'.format(version, date, records))
        
def interproDB(info, force=False):
    iprXML = os.path.join(args.database, 'interpro.xml')
    if not os.path.isfile(iprXML) or force:
        lib.log.info('Downloading InterProScan Mapping file')
        download(URL.get('interpro'), iprXML+'.gz')
        subprocess.call(['gunzip', '-f', 'interpro.xml.gz'], cwd=os.path.join(args.database))
        num_records = None
        version = None
        iprdate = None
        for event, elem in cElementTree.iterparse(iprXML):
            if elem.tag == 'release':
                for x in elem.getchildren():
                    if x.attrib['dbname'] == 'INTERPRO':
                        num_records = int(x.attrib['entry_count'])
                        version = x.attrib['version']
                        iprdate = x.attrib['file_date']
        iprdate = datetime.datetime.strptime(iprdate, "%d-%b-%y").strftime("%Y-%m-%d")            
        info['interpro'] = ('xml', iprXML, version, iprdate, num_records)
    type, name, version, date, records = info.get('interpro')
    lib.log.info('InterProScan XML: version={:} date={:} records={:,}'.format(version, date, records))

def curatedDB(info, force=False):
    curatedFile = os.path.join(args.database, 'ncbi_cleaned_gene_products.txt')
    if not os.path.isfile(curatedFile) or force:
        lib.log.info('Downloaded curated gene names and product descriptions')
        download(URL.get('gene2product'), curatedFile)
        num_records = 0
        curdate = None
        version = None
        with open(curatedFile, 'rU') as infile:
            for line in infile:
                if line.startswith('#version'):
                    version = line.split(' ')[-1].rstrip()
                elif line.startswith('#Date'):
                    curdate = line.split(' ')[-1].rstrip()
                else:
                    num_records += 1
        curdate = datetime.datetime.strptime(curdate, "%m-%d-%Y").strftime("%Y-%m-%d")
        info['gene2product'] = ('text', curatedFile, version, curdate, num_records)
    type, name, version, date, records = info.get('gene2product')
    lib.log.info('Gene2Product: version={:} date={:} records={:,}'.format(version, date, records))

#create directory if doesn't exist
if not os.path.isdir(args.database):
    os.makedirs(args.database)

#create log file
log_name = 'funannotate-setup.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print("-------------------------------------------------------")
lib.SystemInfo()

#get version of funannotate
version = lib.get_version()
lib.log.info("Running %s" % version)

global today
today = datetime.datetime.today().strftime('%Y-%m-%d')

installdbs = []
if 'all' in args.install:
    installdbs = ['merops', 'uniprot', 'dbCAN', 'pfam', 'repeats', 'go', 'mibig', 'interpro', 'busco_outgroups', 'gene2product']
else:
    installdbs = args.install

#if text file with DB info is in database folder, parse into Dictionary
DatabaseFile = os.path.join(args.database, 'funannotate-db-info.txt')
DatabaseInfo = {}
if os.path.isfile(DatabaseFile):
    with open(DatabaseFile, 'rU') as inDB:
        for line in inDB:
            line = line.rstrip()
            db, type, name, version, date, records = line.split('\t')
            DatabaseInfo[db] = (type, name, version, date, int(records))

for x in installdbs:
    if x == 'uniprot':
        uniprotDB(DatabaseInfo, args.force)
    elif x == 'merops':
        meropsDB(DatabaseInfo, args.force)
    elif x == 'dbCAN':
        dbCANDB(DatabaseInfo, args.force)
    elif x == 'pfam':
        pfamDB(DatabaseInfo, args.force)
    elif x == 'repeats':
        repeatDB(DatabaseInfo, args.force)
    elif x == 'go':
        goDB(DatabaseInfo, args.force)
    elif x == 'interpro':
       interproDB(DatabaseInfo, args.force)
    elif x == 'mibig':
        mibigDB(DatabaseInfo, args.force)
    elif x == 'busco_outgroups':
        outgroupsDB(DatabaseInfo, args.force)
    elif x == 'gene2product':
        curatedDB(DatabaseInfo, args.force)
    
#output the database text file and print to terminal        
with open(DatabaseFile, 'w') as outDB:
    for k,v in DatabaseInfo.items():
        data = '%s\t%s\t%s\t%s\t%s\t%i' % (k, v[0], v[1], v[2], v[3], v[4])
        outDB.write('{:}\n'.format(data))

lib.log.info('Funannoate setup complete. Add this to ~/.bash_profile or ~/.bash_aliases:\n\n\texport FUNANNOTATE_DB={:}\n'.format(os.path.abspath(args.database)))
sys.exit(1)

