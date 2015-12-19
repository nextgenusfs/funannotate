import os, subprocess, logging, sys, argparse, inspect, csv

#get the working directory, so you can move back into DB folder to find the files you need
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
DB = os.path.join(parentdir, 'DB')


class colr:
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'
    
def which(name):
    try:
        with open(os.devnull) as devnull:
            subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True
    
def setupLogging(LOGNAME):
    global log
    if 'win32' in sys.platform:
        stdoutformat = logging.Formatter('%(asctime)s: %(message)s', datefmt='[%I:%M:%S %p]')
    else:
        stdoutformat = logging.Formatter(colr.GRN+'%(asctime)s'+colr.END+': %(message)s', datefmt='[%I:%M:%S %p]')
    fileformat = logging.Formatter('%(asctime)s: %(message)s')
    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)
    sth = logging.StreamHandler()
    sth.setLevel(logging.INFO)
    sth.setFormatter(stdoutformat)
    log.addHandler(sth)
    fhnd = logging.FileHandler(LOGNAME)
    fhnd.setLevel(logging.DEBUG)
    fhnd.setFormatter(fileformat)
    log.addHandler(fhnd)

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

'''
def runEggNog(file, output):
    
    #load in annotation dictionary
    EggNog = {}
    with open(os.path.join(DB,'fuNOG.annotations.tsv'), 'rU') as input:
        reader = csv.reader(input, delimiter='\t')
        for line in reader:
            EggNog[line[1]] = line[5]
    
    
    	#first run Hmmer, saving sig domains less than 1e-10
	hmmscan --cpu $4 --domtblout $2.out -E $3 $HMMDB/fuNOG.hmm $1
	#parse results with awk, just get values we care about, length 50% - 150% of target, and then get best match
	grep -v "^#" $2.out | gawk '{if($3/$6 >= 0.5 && $3/$6 <= 1.5) print $1"\t"$4"\t"$7"\t"$3"\t"$6"\t"$3/$6;}' | gsort -k2,2 -k3,3g | gsort -u -k2,2 > $2.filtered
	#now grab descriptions
	gawk 'BEGIN {FS=OFS="\t"} NR==FNR {v[$1]=$3;next} {print $2,$1": "v[$1],$3,$4,$5,$6}' $HMMDB/fuNOG.mapping.txt $2.filtered > $2.results.txt
	#clean up
	rm $2.out
	rm $2.filtered
	echo "Finished.  Ran HMMscan against fuNOG database with Evalue of $3, and then filtered results."
fi
'''