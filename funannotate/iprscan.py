#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import subprocess


def main(args):
    # little wrapper to run interproscan multiprocessing
    # just pass this onto iprscan-local.py
    parentdir = os.path.join(os.path.dirname(__file__))
    cmd = [sys.executable, os.path.join(
        parentdir, 'aux_scripts', 'iprscan-local.py')]
    cmd += args
    subprocess.call(cmd)


if __name__ == "__main__":
    main(sys.argv[1:])
