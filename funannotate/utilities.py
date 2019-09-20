import os
import subprocess


def which(name):
    try:
        with open(os.devnull) as devnull:
            subprocess.Popen([name, '--version'], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def CheckDependencies(input):
    missing = []
    for p in input:
        if which(p) == False:
            missing.append(p)
    if missing != []:
        error = ", ".join(missing)
        print("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
        sys.exit(1)
        
def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def softwrap(string, every=80):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)