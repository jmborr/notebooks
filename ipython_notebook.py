#!/usr/bin/env python

import os
import argparse
import sys

parser = argparse.ArgumentParser(description='Instantiate a python notebook')
parser.add_argument('dir', help='directory containing the notebook.')
args=parser.parse_args()

if not os.path.exists(args.dir):
    raise IOError , 'project {0} not found'.format(args.dir)
    sys.exit(1)

os.system( ' ipython notebook --notebook-dir={0}'.format(args.dir) )
