#!/usr/bin/env python
from __future__ import division

#script to test out monitoring a folder and how I might do that with RunIprScan, as it keeps getting locked up if you just wait for it

#from stack overflow

import datetime, sys, glob, os, time, subprocess, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

#get executable for runiprscan
PATH2JAR = os.path.join(currentdir, 'util', 'RunIprScan-1.1.0', 'RunIprScan.jar')
RUNIPRSCAN_PATH = os.path.join(currentdir, 'util', 'RunIprScan-1.1.0')
IPROUT = 'iprscan'
if not os.path.isdir(IPROUT):
    os.makedirs(IRPOUT)

def update_progress(progress):
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

def runFunction(input, outputdir, num_complete):
    num_files = len(glob.glob1(outputdir,"*.xml"))
    while (num_files < num_complete):
        #launch process
        p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', input, '-m', 'palmer.jona@gmail.com', '-o', outputdir], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #time.sleep(120) #give it 2 minutes to generate something
        while p.poll() is None:
            #wait 1 minute, then check results again
            time.sleep(60)
            num_files = len(glob.glob1(outputdir,"*.xml"))
            pct = num_files / num_complete
            update_progress(pct)
            #monitor the output folder for recent changes in last 30 minutes
            now = datetime.datetime.now()
            ago = now - datetime.timedelta(minutes=30)  #if nothing happens in 30 minutes, relaunch service 
            file_list = []
            for path in glob.glob(IPROUT + "/*.xml"):
                (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(path)
                if datetime.datetime.fromtimestamp(mtime) > ago:
                    file_list.append(path)
            if not (file_list):
                try:
                    p.terminate()
                except OSError:
                    pass
                break
        num_files = len(glob.glob1(outputdir,"*.xml"))

  
runFunction(sys.argv[1], IPROUT, 10)

print("Finised, completely passed the while loops.....is this depending on previous function finishing?")