import os
import numpy
import argparse
from tempfile import mkstemp
from pdb import set_trace as tr

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Find centroid and average structure.
    Example: cluster_trajectory.py ./8styrene32.pdb ./equil_rms2first.dcd 1000 '@1-4133&@H*'""")
    parser.add_argument('topfile', help='input PDB file')
    parser.add_argument('trajfile', help='DCD trajectory')
    parser.add_argument('nsampling', type=int, help='number of frames to sample for the clustering')
    parser.add_argument('mask', help='ptraj mask defining the atoms for which to use for the clustering')
    parser.add_argument('centroid', help='output name of centroid PDB file')
    parser.add_argument('average', help='output name of average PDB file')
    args=parser.parse_args()


# Find number of frames in the trajectory
handle, tempfile = mkstemp(prefix='junk', dir='/tmp')
os.system('dumpdcd {0} | head -1 > {1}'.format(args.trajfile,tempfile))
tf = int(open(tempfile).readline().strip())
os.system('/bin/rm {0}'.format(tempfile))
dt = int(tf/args.nsampling)

script='''trajin {0}  1  {1} {2}
cluster out clustering representative pdb average pdb averagelinkage clusters 1 rms {3}
'''.format(args.trajfile, tf, dt, args.mask)

handle, scriptfile = mkstemp(prefix='junk', dir='/tmp')
open(scriptfile,'w').write(script)

os.system('module load amber/amber12; ptraj {0} < {1}'.format(args.topfile, scriptfile))
os.system('/bin/mv clustering.rep.c0 {0}'.format(args.centroid))
os.system('/bin/mv clustering.avg.c0 {0}'.format(args.average))
os.system('/bin/rm clustering.c0 clustering.txt {0}'.format(scriptfile))
