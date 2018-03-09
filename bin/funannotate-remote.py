#!/usr/bin/env python
from __future__ import division

import sys, os, subprocess, inspect, shutil, argparse, time, glob, requests, zipfile, urllib2
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

RUNIPRSCAN = os.path.join(parentdir, 'util', 'runIPRscan.py')
XMLCombine = os.path.join(parentdir, 'util', 'xmlcombine.py')

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-remote.py',
    description='''Script that adds functional annotation to a genome using remote searches.''',
    epilog="""Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--input', help='Folder from funannotate predict.')
parser.add_argument('-g', '--genbank', help='Annotated genome in GenBank format')
parser.add_argument('-m', '--methods', required=True, nargs='+', choices=['all', 'phobius', 'antismash', 'interproscan'], help='Method to run')
parser.add_argument('-o', '--out', help='Basename of output files')
parser.add_argument('-e', '--email', required=True, help='Email address for IPRSCAN server')
parser.add_argument('--force', action='store_true', help='Over-write output folder')
parser.add_argument('-a', '--antismash', default='fungi', choices=['fungi','plants'], help='antiSMASH server')
args=parser.parse_args()

def runIPRpython(Input):
    base = Input.split('/')[-1]
    base = base.split('.fa')[0]
    OUTPATH = os.path.join(IPROUT, base)
    subprocess.call([sys.executable, RUNIPRSCAN, '--goterms','--email', args.email, '--outfile', OUTPATH, '--input', Input], stderr=FNULL, stdout=FNULL)
    #now rename output files, just keep xml and tsv files?
    time.sleep(3) #make sure there is time for all files to show up
    os.rename(OUTPATH+'.xml.xml', OUTPATH+'.xml')
    os.rename(OUTPATH+'.tsv.txt', OUTPATH+'.tsv')
    os.rename(OUTPATH+'.svg.svg', OUTPATH+'.svg')
    os.rename(OUTPATH+'.sequence.txt', OUTPATH+'.fa')
    os.remove(OUTPATH+'.gff.txt')
    os.remove(OUTPATH+'.htmltarball.html.tar.gz')
    os.remove(OUTPATH+'.log.txt')
    os.remove(OUTPATH+'.out.txt')

def download(url, name):
    file_name = name
    u = urllib2.urlopen(url)
    f = open(file_name, 'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print("Downloading: {0} Bytes: {1}".format(url, file_size))
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
    f.close()

#create log file
log_name = 'funannotate-remote.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.SystemInfo()

#get version of funannotate
version = lib.get_version()
lib.log.info("Running %s" % version)

#need to do some checks here of the input
genbank = ''
Proteins = ''
if not args.input:
    #did not parse folder of funannotate results, so need either gb + gff or fasta + proteins, + gff and also need to have args.out for output folder
    if not args.out:
        lib.log.error("If you are not providing funannotate predict input folder, then you need to provide an output folder (--out)")
        sys.exit(1)
    else:
        outputdir = args.out
        #create outputdir and subdirs
        if not os.path.isdir(outputdir):
            os.makedirs(outputdir)
            os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            os.makedirs(os.path.join(outputdir, 'annotate_results'))
            os.makedirs(os.path.join(outputdir, 'logfiles')) 
    if not args.genbank:
        lib.log.error("You did not specifiy the apropriate input files, either: \n1) Funannotate input \n2) GenBank")
        sys.exit(1)
    else:
        #create output directories
        if not os.path.isdir(outputdir):
            os.makedirs(outputdir)
            os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            os.makedirs(os.path.join(outputdir, 'annotate_results'))
            os.makedirs(os.path.join(outputdir, 'logfiles'))
        else:
            lib.log.error("Output directory %s already exists" % (outputdir))
            if not os.path.isdir(os.path.join(outputdir, 'annotate_misc')):
                os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            if not os.path.isdir(os.path.join(outputdir, 'annotate_results')):
                os.makedirs(os.path.join(outputdir, 'annotate_results'))
            if not os.path.isdir(os.path.join(outputdir, 'logfiles')):
                os.makedirs(os.path.join(outputdir, 'logfiles'))                
        genbank = args.genbank
        Scaffolds = os.path.join(outputdir, 'annotate_misc', 'genome.scaffolds.fasta')
        Proteins = os.path.join(outputdir, 'annotate_misc', 'genome.proteins.fasta')
        Transcripts = os.path.join(outputdir, 'annotate_misc', 'genome.transcripts.fasta')
        GFF = os.path.join(outputdir, 'annotate_misc', 'genome.gff3')
        lib.log.info("Checking GenBank file for annotation")
        if not lib.checkGenBank(genbank):
            lib.log.error("Found no annotation in GenBank file, exiting")
            sys.exit(1)
        lib.gb2allout(genbank, GFF, Proteins, Transcripts, Scaffolds)
    
else:
    #should be a folder, with funannotate files, thus store results there, no need to create output folder
    if not os.path.isdir(args.input):
        lib.log.error("%s directory does not exist" % args.input)
        sys.exit(1)
    if os.path.isdir(os.path.join(args.input, 'update_results')): #funannotate results should be here
        inputdir = os.path.join(args.input, 'update_results')
        outputdir = args.input
    elif os.path.isdir(os.path.join(args.input, 'predict_results')):
    	inputdir = os.path.join(args.input, 'predict_results')
        outputdir = args.input
    else:
        inputdir = os.path.join(args.input) #here user specified the predict_results folder, or it is a custom folder

    #get files that you need
    for file in os.listdir(inputdir):
        if file.endswith('.gbk'):
            genbank = os.path.join(inputdir, file)
    #now create the files from genbank input file for consistency in gene naming, etc
    if not genbank:
        lib.log.error("Properly formatted 'funannotate predict' files do no exist in this directory")
        sys.exit(1)
    else:
        if 'predict_results' in inputdir or 'update_results' in inputdir: #if user gave predict_results folder, then set output to up one directory
            outputdir = lib.get_parent_dir(inputdir)
        else:
            if not args.out:
                outputdir = inputdir #output the results in the input directory
            else:
                outputdir = args.out
                if not os.path.isdir(outputdir):
                    os.makedirs(outputdir)
        #create output directories
        if not os.path.isdir(os.path.join(outputdir, 'annotate_misc')):
            os.makedirs(os.path.join(outputdir, 'annotate_misc'))
            os.makedirs(os.path.join(outputdir, 'annotate_results'))
        else:
            lib.log.error("Output directory %s already exists, will use any existing data.  If this is not what you want, exit, and provide a unique name for output folder" % (outputdir))
        lib.log.info("Parsing input files")
        Scaffolds = os.path.join(outputdir, 'annotate_misc', 'genome.scaffolds.fasta')
        Proteins = os.path.join(outputdir, 'annotate_misc','genome.proteins.fasta')
        Transcripts = os.path.join(outputdir, 'annotate_misc', 'genome.transcripts.fasta')
        GFF = os.path.join(outputdir, 'annotate_misc', 'genome.gff3')
        lib.log.debug("Generating files from %s" % genbank)
        lib.gb2allout(genbank, GFF, Proteins, Transcripts, Scaffolds)

#make sure logfiles directory is present, will need later
if not os.path.isdir(os.path.join(outputdir, 'logfiles')):
    os.makedirs(os.path.join(outputdir, 'logfiles'))

#get absolute path for all input so there are no problems later, not using Transcripts yet could be error? so take out here
Proteins = os.path.abspath(Proteins)
genbank = os.path.abspath(genbank)


if 'phobius' in args.methods or 'all' in args.methods:
    #run Phobius to predict secreted proteins and membrane, default is local if installed, otherwise remote
    phobius_out = os.path.join(outputdir, 'annotate_misc', 'phobius.results.txt')
    phobiusLog = os.path.join(outputdir, 'logfiles', 'phobius.log')
    lib.log.info("Predicting secreted and transmembrane proteins using Phobius")
    if not lib.checkannotations(phobius_out):
        if args.email:
            subprocess.call([os.path.join(parentdir, 'util', 'phobius-multiproc.py'), '-i', Proteins, '-o', phobius_out, '-e', str(args.email), '-l', phobiusLog])
        else:
            subprocess.call([os.path.join(parentdir, 'util', 'phobius-multiproc.py'), '-i', Proteins, '-o', phobius_out, '-l', phobiusLog])        

if 'interproscan' in args.methods or 'all' in args.methods:
    IPRCombined = os.path.join(outputdir, 'annotate_misc', 'iprscan.xml')
    #run interpro scan
    IPROUT = os.path.join(outputdir, 'annotate_misc', 'iprscan')
    PROTS = os.path.join(outputdir, 'annotate_misc', 'protein_tmp')
    for i in IPROUT,PROTS:
        if not os.path.exists(i):
            os.makedirs(i)
    #now run interproscan
    #split input into individual files
    lib.splitFASTA(Proteins, PROTS)
    #now iterate over list using pool and up to 25 submissions at a time
    proteins = []
    for file in os.listdir(PROTS):
        if file.endswith('.fa'):
            file = os.path.join(PROTS, file)
            proteins.append(file)
    
    num_files = len(glob.glob1(IPROUT,"*.xml"))
    num_prots = len(proteins)
    lib.log.info("Now running InterProScan search remotely using EBI servers on " + '{0:,}'.format(num_prots) + ' proteins')
    #build in a check before running (in case script gets stopped and needs to restart
    finished = []
    for file in os.listdir(IPROUT):
        if file.endswith('.xml'):
            base = file.split('.xml')[0]
            fasta_file = os.path.join(PROTS, base+'.fa')
            finished.append(fasta_file)

    finished = set(finished) #make sure no duplicates
    runlist = [x for x in proteins if x not in finished]
    if len(runlist) < num_prots:
        lib.log.info("Results found, querying remaining %i proteins" % len(runlist))
    #start up the list, max 25 at a time 
    lib.runMultiProgress(runIPRpython, runlist, 25)
    #clean up protein fasta files
    shutil.rmtree(PROTS)
    #now convert to single file and then clean up
    with open(IPRCombined, 'w') as output:
        subprocess.call([sys.executable, XMLCombine, IPROUT], stdout = output)
    if lib.checkannotations(IPRCombined):
        shutil.rmtree(IPROUT)

if 'antismash' in args.methods or 'all' in args.methods:
    if args.antismash == 'fungi':
        base_address = "https://fungismash.secondarymetabolites.org"
        job_parameters = {'email': args.email, 'smcogs': 'on', 'knownclusterblast': 'on', 'activesitefinder': 'on', 'subclusterblast': 'on'}
    elif args.antismash == 'plants':
        base_address = "https://plantismash.secondarymetabolites.org"
        job_parameters = {'email': args.email, 'knownclusterblast': 'on', 'subclusterblast': 'on'}
    version = requests.get(base_address+"/api/v1.0/version")
    as_vers = version.json()['antismash_generation']
    tax = version.json()['taxon']
    as_status = requests.get(base_address+"/api/v1.0/stats")
    queue = as_status.json()['queue_length']
    running = as_status.json()['running']
    lib.log.info("Connecting to antiSMASH %s v%s webserver" % (tax, as_vers))
    lib.log.info("Queue Length: %s; Jobs Running: %s" % (queue, running))
    lib.log.info("PLEASE to not abuse the webserver, be considerate!")
    if int(queue) > 10 and not args.force:
        lib.log.error("There are more than 10 antiSMASH jobs in queue, use --force to submit anyway")
        sys.exit(1)
    job_files = {'seq': open(genbank, 'rb')}

    lib.log.info("Uploading %s to webserver" % genbank)
    postjob = requests.post(base_address+"/api/v1.0/submit", files=job_files, data=job_parameters)
    jobid = postjob.json()['id']
    #now we can query the job every so often, not sure what is reasonable here, start with 2 minutes?
    lib.log.info("Waiting for results from job: %s" % jobid)
    while True:
        job_status = requests.get(base_address+"/api/v1.0/status/"+jobid)
        if job_status.json()['status'] == 'done':
            break
        time.sleep(120) #check every 2 minutes
    result_url = job_status.json()['result_url']
    base_url = result_url.replace('index.html', '')
    lib.log.info("antiSMASH v%s job finished" % (as_vers))
    lib.log.debug("%s" % job_status.json())
    #need to retrieve results, have to find link, seems like this might be first scaffold name?
    #after asking Kai Blin - there is no "easy" way to identify the output name, however, think I can grab the html file and parse it
    job_html = requests.get(base_address+result_url)
    for line in job_html.iter_lines():
        if 'Download GenBank summary file' in line:
             cols = line.split('a href="')
    for x in cols:
        if '.zip' in x:
            link = x.split('"')[0]
    baselink = link.replace('.zip', '')
    download_url = base_address+base_url+link
    download(download_url, 'antiSMASH.zip')
    #now unzip and move folder
    zipref = zipfile.ZipFile('antiSMASH.zip', 'r')
    zipref.extractall(outputdir)
    zipref.close()
    os.remove('antiSMASH.zip')
    lib.log.info("Results folder: %s/%s" % (outputdir, jobid))
    #now grab the GBK files from folder as you will need just that for annotation, place in annotate_misc folder for auto-detection
    anti_GBK = os.path.join(outputdir, jobid, baselink+'.final.gbk')
    final = os.path.join(outputdir, 'annotate_misc', 'antiSMASH.results.gbk')
    shutil.copyfile(anti_GBK, final)
    lib.log.info("Results GBK: %s" % final)
    
lib.log.info("Remote searches complete")
#move logfile
os.rename(log_name, os.path.join(outputdir, 'logfiles', log_name))