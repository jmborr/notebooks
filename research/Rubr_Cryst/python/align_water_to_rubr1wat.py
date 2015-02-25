import os
import argparse
import tempfile
from amber.amber14.cpptraj import hbond
from pdb import set_trace as tr

template='''parm rubr1wat.top
reference rubr1wat.crd
trajin _INCRD_
rms reference _REFMASK_
strip !:WAT
trajout _OUTPDB_ pdb
'''

template2='''parm rubr1wat.top
reference rubr1wat.crd
trajin _INCRD_
#bondinfo _ATOMMASK_ #Not neccessary since angleinfo also contains the bonded atoms
angleinfo _ATOMMASK_
'''

def genmask(atom,incrd):
    '''Create mask for atom, its bonded atoms, and its angle-atoms'''
    script=template2.replace('_INCRD_',incrd)
    atommask=':{0}@{1}'.format(atom.resn,atom.atomname)
    script=script.replace('_ATOMMASK_',atommask)
    handle,filename=tempfile.mkstemp(dir='/tmp',prefix='junk')
    handle,resultsfile=tempfile.mkstemp(dir='/tmp',prefix='junk')   
    open(filename,'w').write(script)
    os.system('cpptraj -i {0} &> {1}'.format(filename,resultsfile))
    neighbors=[]
    for line in open(resultsfile).readlines():
        if atommask in line and line.count('@')==3:
            triad=line.split()[-1]
            neighbors+=triad[1:-1].split(',')
    neighbors=list(set(neighbors)) #remove duplicates
    os.system('/bin/rm -f {0} {1}'.format(filename,resultsfile))
    return '@'+','.join(neighbors)

def align_water_to_rubr1wat(incrd,outpdb):
    nres=54 #number of residues in rubredoxin
    # Find hydrogen bond
    bondname = os.path.basename(args.incrd)[:-4] #assume name of input crd file is the bond string
    bond = hbond.Bond(bondname)

    if bond.acceptor.resname!='WAT':
        bond.acceptor.resn=bond.acceptor.resn%nres
    else:
        bond.donor.resn=bond.donor.resn%nres
        bond.donorH.resn=bond.donorH.resn%nres

    #Find protein atoms to serve as reference mask in the alignment
    donor=bond.donor
    if donor.resname!='WAT':
        if donor.atomname=='N' and donor.resn!=1:
            #donor is the nitrogen in the protein backbone and not the first residue
            refmask='(:{0}@H,N,CA)|(:{1}@C,O,CA)'.format(donor.resn,donor.resn-1)
        else:
            #Find first and second neighbors of the donor hydrogen
            refmask=genmask(donor,args.incrd)
    else:
        acceptor=bond.acceptor
        if acceptor.atomname=='O':
            #acceptor is the oxygen in the protein backbone, and not in the last residue
            refmask='(:{0}@C,O,CA)|(:{1}@H,N,CA)'.format(acceptor.resn,acceptor.resn+1)
        else:
            refmask=genmask(acceptor,args.incrd)

    script=template.replace('_INCRD_',args.incrd)
    script=script.replace('_REFMASK_',refmask)
    script=script.replace('_OUTPDB_',args.outpdb)
    handle,filename=tempfile.mkstemp(dir='/tmp',prefix='junk')
    open(filename,'w').write(script)
    os.system('cpptraj -i {0}'.format(filename))
    os.system('/bin/rm -f {0}'.format(filename))

parser = argparse.ArgumentParser(description='align the water of the input crd file to the reference rubredoxin, and save the file')
parser.add_argument('incrd', type=str, help='input average crd file')
parser.add_argument('outpdb', type=str, help='output pdb file containing only the aligned water molecule')
args = parser.parse_args()

align_water_to_rubr1wat(args.incrd,args.outpdb)
