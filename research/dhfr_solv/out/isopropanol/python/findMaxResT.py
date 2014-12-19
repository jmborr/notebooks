import argparse
from pdb import set_trace as tr
from utilities.readingWritingFiles import read_column

parser = argparse.ArgumentParser(description='Extract, for each residue, the maximum residence time among the simulations with isopropanol ranging from 15% to 30%')
parser.add_argument('maxTfile',help='output file containing the maximum times')
args = parser.parse_args()

isos='0.15 0.20 0.25 0.30'.split() # isopropanol concentraitons
tdat='/projects2/research/dhfr_solv/out/isopropanol/_ISO_/solv/Prod/mdprod/isoResidenceTimes.dat' # template file
nres=161 # number of residues (includes the cofactor and ligand)
maxtimes=[0]*nres
for iso in isos:
    dat=tdat.replace('_ISO_',iso)
    times=read_column(dat, 2, isFloat=1)
    for i in range(nres):
        if maxtimes[i] < times[i]: maxtimes[i] = times[i]

buf='# residue maximum-residence-time\n'
for ires in range(nres):
    buf += '{0:3d} {1:6.0f}\n'.format(1+ires, maxtimes[ires])
open(args.maxTfile, 'w').write(buf)

