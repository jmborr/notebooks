import os
import numpy
import argparse
from tempfile import mkstemp
from utilities.readingWritingFiles import read_column
from utilities.path import which
from pdb import set_trace as tr

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Obtain (R)MSD averaged over different chunks of the trajectory.
    Example: atomic_fluct_t0average.py ./8styrene32.pdb ./equil_rms2first.dcd 30000 1000  10 '@1-4133&@H*'""")
    parser.add_argument('topfile', help='input PDB file')
    parser.add_argument('trajfile', help='DCD trajectory')
    parser.add_argument('span', type=int, help='number of frames in the trajectory')
    parser.add_argument('t', type=int, help='number of frames over which to calculate MSD')
    parser.add_argument('nt0', type=int, help='number of times to calculate MSD, for the average')
    parser.add_argument('mask', help='ptraj mask defining the atoms for which to calculate the MSD')
    parser.add_argument('outfile', help='filename to store average MSD')
    parser.add_argument('--seriesfile', type=str, help='filename storing MSD for each chunk')
    parser.add_argument('--reference', help='reference PDB file for extra RMSD step')
    parser.add_argument('--msd', type=str, default=None, help='Output MSD instead of RMSD')

    args=parser.parse_args()

# Check ptraj exists
if not which('ptraj'):
    raise Exception('ptraj not found')

template='''trajin _trajfile_  _t0_  _t0+t_
_REFERENCE_atomicfluct _mask_ out _outfn_
'''

template = template.replace('_mask_', args.mask)
template = template.replace('_trajfile_', args.trajfile)

if args.reference:
    template=template.replace('_REFERENCE_', 'reference {0}\nrms reference {1}\n'.format(args.reference,args.mask))
else:
    template=template.replace('_REFERENCE_', '')

dt0 = int( (args.span - args.t)/args.nt0 )
if dt0 < 1:
    raise ValueError("nt0 too large")


series='#it0 RMSD\n'
if args.msd:
    series='#it0 MSD\n'
av = 0

for it0 in range(args.nt0):
    print '{0} sampling REMAIN'.format(args.nt0-it0)
    # Define the ptraj script
    t0 = it0*dt0   # first frame
    tf = t0+args.t # last frame
    script = template.replace('_t0_', str(t0))
    script = script.replace('_t0+t_', str(tf))
    handle,outfn = mkstemp(prefix='junk', dir='/tmp')
    script = script.replace('_outfn_', outfn)

    # Run the ptraj job
    handle,scriptfile = mkstemp(prefix='junk', dir='/tmp')
    open(scriptfile,'w').write(script)
    os.system('ptraj {0} < {1}'.format(args.topfile, scriptfile))
    
    # Store results to memory
    rmsd = numpy.array( read_column(outfn, 2, isFloat=1) )
    if args.msd:
        msd = numpy.average(rmsd*rmsd)
        series += '{0:3d} {1:5.2f}\n'.format(it0, msd)
        av += msd
    else:
        rmsd = numpy.average(rmsd)
        series += '{0:3d} {1:5.2f}\n'.format(it0, rmsd)
        av += rmsd
       
    os.system('/bin/rm -f {0} {1}'.format(outfn,scriptfile))

# Write output files
open(args.outfile,'w').write('{0}\n'.format(av/args.nt0))
if args.seriesfile:
    open(args.seriesfile,'w').write(series)
if args.msd:
    print 'MSD = ', av/args.nt0
else:
    print 'RMSD = ', av/args.nt0
