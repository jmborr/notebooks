import argparse
from amber.amber14.cpptraj import hbond
from copy import copy
from pdb import set_trace as tr

protocols=('highest', 'average')
L=54
Lmax=216 #maximun number for protein residue
np=Lmax/L #number of proteins
def analog_bond(bond1,bond2):
    '''Check if two bonds share analog protein residues
    Returns:
      true if bond1 and bond2 have analogous atoms, defined as atoms on different proteins
      but at the same residue
'''
    if bond1.acceptor.resn < Lmax:
        name1=bond1.acceptor.atomname
        name2=bond2.acceptor.atomname
        resn1=bond1.acceptor.resn 
        resn2=bond2.acceptor.resn
    else:
        name1=bond1.donorH.atomname
        name2=bond2.donorH.atomname
        resn1=bond1.donorH.resn
        resn2=bond2.donorH.resn
    if (name1==name2) and ((resn1>resn2 and resn1%L==resn2) or (resn2>resn1 and resn2%L==resn1)):
        return true
    return false

def analog_avginfo(avginfo1, avginfo2):
    '''Check if two average hbond lines share analog protein residue'''
    return analog_bond(avginfo1.bond, avginfo2.bond)

def generate_key(avginfo):
    '''Create the key identifying the protein atom of interest'''
    bond=avginfo.bond
    if bond.acceptor.resn < Lmax:
        analog_atom=copy(bond.acceptor)
    else:
        analog_atom=copy(bond.donorH)
    analog_atom.resn=analog_atom.resn%L
    return str(analog_atom)

parser = argparse.ArgumentParser(description='Prune average hbond file')
parser.add_argument('avgin', type=str, help='input average hbond file')
parser.add_argument('protocol', type=str, default=protocols[0],
                    help='selection protocol, one of {0}'.format(','.join(protocols)))
parser.add_argument('avgout', type=str, help='output average hbond file')
args = parser.parse_args()

inavg=hbond.AvgInfoFile(args.avgin) #create list of AvgInfo objects
pruned={}
#tr()
for avginfo in inavg.records:
    key=generate_key(avginfo)
    if key in pruned.keys():
        #update key according to protocol
        recorded=pruned[key]
        if args.protocol=='highest' and avginfo.frames>recorded.frames:
            pruned[key]=avginfo
        elif args.protocol=='average':
            recorded.frames += avginfo.frames
    else:
        pruned[key]=avginfo

#tr()
if args.protocol=='average':
    for recorded in pruned.values():
       recorded.frames = int(recorded.frames/np)
#tr()
#Output pruned records
outavg=hbond.AvgInfoFile()
outavg.header=inavg.header
outavg.records=pruned.values()
outavg.records.sort(key=lambda x: x.frames, reverse=True) #Sort by decreasing number of frames
open(args.avgout,'w').write(str(outavg))
