import os
import argparse
from amber.amber14.cpptraj import hbond
from pdb import set_trace as tr

parser = argparse.ArgumentParser(description='Issue pymol command for the first lines of the input average hbond file')
parser.add_argument('avgin', type=str, help='input average hbond file')
parser.add_argument('--maxn', type=int, default=100, help='maximum number of bonds to look')
args = parser.parse_args()

inavg=hbond.AvgInfoFile(args.avgin) #create list of AvgInfo objects
for avginfo in inavg.records[0:args.maxn]:
    bondstring = str(avginfo.bond)
    print '{0}_onlyWAT.pdb'.format(bondstring)

# 25 lines for (water) acceptor
# 64 lines for (water) donor
