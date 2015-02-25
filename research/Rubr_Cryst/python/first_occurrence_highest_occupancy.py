import os
import argparse
import tempfile
from amber.amber14.cpptraj import hbond
from pdb import set_trace as tr

template1='''parm top
trajin dcd
image
hbond H1 series donormask _DONORMASK_ donorhmask _DONORHMASK_ acceptormask _ACCEPTORMASK_
go
filter H1[solutehb] min 0.5 max 1.5
trajout junk.dcd dcd
'''

template2='''parm top
trajin junk.dcd
autoimage :_RESNW_
strip (:Na+)|(_EXTRA_WATER_)|(_EXTRA_PROTEIN_)
trajout _OUTPDB_ pdb start 1 stop 1 offset 1
trajout _OUTCRD_ start 1 stop 1 offset 1
'''

parser = argparse.ArgumentParser(description='Find and save snapshot of first occurence with highest occupancy')
parser.add_argument('avgin', type=str, help='input average hbond file')
parser.add_argument('outd', type=str, help='output directory where to store the PDB files')
parser.add_argument('--maxn', type=int, default=100, help='maximum number of bonds to look')
args = parser.parse_args()

inavg=hbond.AvgInfoFile(args.avgin) #create list of AvgInfo objects
for avginfo in inavg.records[0:args.maxn]:

    # First, create trajectory with only frames containing the hydrogen bond
    # find hydrogen bond mask
    bond = avginfo.bond
    acceptormask = ':{0}@{1}'.format(bond.acceptor.resn,bond.acceptor.atomname)
    script = template1.replace('_ACCEPTORMASK_', acceptormask)
    donormask = ':{0}@{1}'.format(bond.donor.resn,bond.donor.atomname)
    script = script.replace('_DONORMASK_',donormask)
    donorhmask = ':{0}@{1}'.format(bond.donorH.resn,bond.donorH.atomname)
    script = script.replace('_DONORHMASK_',donorhmask)
    # invoque cpptraj
    handle,filename=tempfile.mkstemp(prefix='junk', dir='/tmp')
    open(filename, 'w').write(script)
    os.system('cpptraj -i {0}'.format(filename))
    os.system('/bin/rm -f {0}'.format(filename))
    
    # Second, load the filtered out trajectory and select first frame
    # find extra water
    resnW=bond.acceptor.resn
    if bond.donor.resname == 'WAT': resnW=bond.donor.resn
    script = template2.replace('_RESNW_', str(resnW))
    extrawater =':WAT&!:{0}'.format(resnW)
    script = script.replace('_EXTRA_WATER_', extrawater)
    # find extra proteins
    nres=54 #54 residues per protein
    resnP=bond.acceptor.resn
    if bond.donor.resname != 'WAT': resnP=bond.donor.resn
    firstindex = 1 + nres * int(resnP/nres)
    lastindex = firstindex + nres - 1
    extraprotein =':1-216&!(:{0}-{1})'.format(firstindex,lastindex)
    script = script.replace('_EXTRA_PROTEIN_', extraprotein)
    # name
    name='{0}/{1}.pdb'.format(args.outd,str(bond))
    script = script.replace('_OUTPDB_', name)
    name='{0}/{1}.crd'.format(args.outd,str(bond))
    script = script.replace('_OUTCRD_', name)
    # invoque cpptraj
    handle,filename=tempfile.mkstemp(prefix='junk', dir='/tmp')
    open(filename, 'w').write(script)
    os.system('cpptraj -i {0}'.format(filename))
    os.system('/bin/rm -f {0}'.format(filename))
    
