import os
import argparse
from tempfile import mkstemp
from pdb import set_trace as tr

def execptraj(topfile, script):
    handle, scriptfile = mkstemp(dir='/tmp', prefix='junk', suffix='.ptraj')
    open(scriptfile,'w').write(script)
    os.system('ptraj {0} < {1}'.format(topfile, scriptfile) )
    os.system('/bin/rm {0}'.format(scriptfile))

def extract(topfile, indcdfile, protein_index, outdcdfile):
    '''Extract one of the unsolvated proteins into a dcd trajectory'''
    nresidue = 54  # 54 residues per rubredoxin
    firstindex = 1+nresidue*protein_index
    script = '''trajin {0}
strip !:{1}-{2}
trajout {3} charmm'''.format(indcdfile, firstindex, firstindex + nresidue-1, outdcdfile)
    execptraj(topfile, script)

def GenerateSingleTopfile(topfile, indcd, singletopfile):
    '''Create PDB of a single protein'''
    nresidue = 54
    script='''trajin {0} 1 1
strip !:1-{1}
trajout {2} pdb'''.format(indcd, nresidue, singletopfile)
    execptraj(topfile, script)
    os.system('/bin/mv {0}.1 {0}'.format(singletopfile))

def serialize_trajectories(topfile, singletrajs, serialtraj):
    '''Merge trajectories into one'''
    script=''
    for singletraj in singletrajs: script += 'trajin {0}\n'.format(singletraj)
    script += 'trajout {0} charmm'.format(serialtraj)
    execptraj(topfile, script)

def RMStoCentroid(topfile, serialtraj, mask, centroidPDB, averagePDB, outdcdfile, rmsfile):
    '''Find centroid, then RMS to it'''
    script = '''trajin {0} 1 40000 4
cluster out /tmp/cluster representative pdb average pdb averagelinkage clusters 1 rms !@H*'''.format(serialtraj)
    execptraj(topfile, script)
    os.system('/bin/mv /tmp/cluster.rep.c0 {0}'.format(centroidPDB))
    os.system('/bin/mv /tmp/cluster.avg.c0 {0}'.format(averagePBD))
    os.system('/bin/rm /tmp/cluster*')

    script = '''trajin {0}
reference {1}
rms reference out {2} {3}
trajout {4} charmm'''.format(serialtraj, averagePDB, rmsfile, mask, outdcdfile)
    execptraj(topfile, script)
    os.system('/bin/rm {0}'.format(serialtraj))

parser = argparse.ArgumentParser(description="""Extract and merge each of the protein trajectories, doing an RMS to centroid.
Example: python unsolv_serial.py topfile indcd outdcd""")
parser.add_argument('topfile', help='input PDB file of the solvated system')
parser.add_argument('indcd', help='input solvated trajectory of the four rubredoxins')
parser.add_argument('outdcd', help='output unsolvated trajectory, each frame one protein, RMS to centroid')
parser.add_argument('centroid', help='output PDB centroid conformation')
parser.add_argument('average', help='output PBD average conformation')
parser.add_argument('rms', help='output file storing the RMS of each conformation to average')
parser.add_argument('--mask', default='!@H*', help='mask to do the clustering search for the centroid. Default is heavy atoms only')
args=parser.parse_args()

# Extract each of the four proteins into a trajectory
nprotein = 4
singletrajs = []
tr()
for protein_index in range(nprotein):
    handle, outdcdfile = mkstemp(dir='/tmp',prefix='junk',suffix='.dcd') 
    extract(args.topfile, args.indcd, protein_index, outdcdfile)
    singletrajs.append(outdcdfile)

# Create temporary PDB containing one protein only
handle, singletopfile = mkstemp(dir='/tmp', prefix='junk', suffix='.pdb')
GenerateSingleTopfile(args.topfile, args.indcd, singletopfile)

# Merge the four trajectories
handle, serialtraj = mkstemp(dir='/tmp', prefix='junk', suffix='.dcd')
serialize_trajectories(singletopfile, singletrajs, serialtraj)
for singletraj in singletrajs: os.system('/bin/rm {0}'.format(singletraj))

# Find centroid and then do RMS of each snapshot to it
RMStoCentroid(singletopfile, serialtraj, args.mask, args.centroid, args.average, args.outdcd, args.rms)

# Final cleanup
os.system('/bin/rm {0}'.singletopfile)

