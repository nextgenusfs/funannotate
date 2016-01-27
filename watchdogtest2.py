#!/usr/bin/env python

#script to test out monitoring a folder and how I might do that with RunIprScan, as it keeps getting locked up if you just wait for it

#from stack overflow

import datetime, sys, glob, os, time

while True:
    now = datetime.datetime.now()
    ago = now - datetime.timedelta(minutes=60)

    file_list = []
    for path in glob.glob(sys.argv[1] + "/*.xml"):
        (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(path)

        if datetime.datetime.fromtimestamp(mtime) > ago:
            file_list.append(path)

    if not (file_list):
        ## process your list
        print("nothing has happened in last 10 secs")
    else:
        print("%i files changed in last hour" % len(file_list))

    time.sleep(10)