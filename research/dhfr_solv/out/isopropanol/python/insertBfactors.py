import argparse
from mypdb.smallThings import insertBfact
from utilities.readingWritingFiles import read_column

parser = argparse.ArgumentParser(description='insert Bfactor into a PDB file')
parser.add_argument('inpdb',help='input PDB file')
parser.add_argument('bfile',help='file containing the B-factors')
parser.add_argument('outpdb', help='output PDB file')
parser.add_argument('--icol', type=int, default=1, help='column in bfile for the B-factors. Default=1')
parser.add_argument('--byres', default='no', help='are the  b-factors in bfile per residue basis? Default=no')
parser.add_argument('--rescale', default='no', help='rescale to range [0,100]. Useful if bfile contains negative B-factors')
args = parser.parse_args()

# Are the B-factors only provided for the whole residue
if args.byres.lower()[0] == 'y':
    args.byres = True
else:
    args.byres = False

# Do we rescale?
if args.rescale.lower()[0] == 'y':
    args.rescale = True
else:
    args.rescale = False

blist = read_column(args.bfile, args.icol, isFloat=1)
buf = insertBfact(args.inpdb, blist, byres=args.byres, rescale=args.rescale)
open(args.outpdb,'w').write(buf)

