#!/usr/env/python

from pdb import set_trace as tr
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Create PDB file tagging selected hydrogens.')
parser.add_argument('infile', help='input PDB file')
parser.add_argument('selection', type=str, help='One of: RNA, water')
parser.add_argument('outfile', help='output PDB file')
args=parser.parse_args()

# Check input first
if not os.path.exists(args.infile):
    raise IOError(args.infile)

sels=['RNA', 'water']
if args.selection not in sels:
    raise IOError('selection must be one of {0}'.format(','.join(sels)) )

if not os.path.dirname(args.outfile):
    raise IOError('directory {0} does not exists'.format(os.path.dirname(args.outfile)))

buf=''
for line in open(args.infile).readlines():
    bufline = line
    if line[0:5]=='ATOM ' and (line[13]=='H' or line[12]=='H'):
        resname = line[17:21]
        if (args.selection=='water' and resname=='TIP3') or \
           (args.selection=='RNA' and resname in ['ADE ', 'CYT ', 'GUA ', 'URA ']):
            bufline = line[:62] + '1.00' + line[66:] # Flag by setting B-factor as 1.00
    buf += bufline

open(args.outfile,'w').write(buf)

    

