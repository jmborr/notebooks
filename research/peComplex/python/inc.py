from pdb import set_trace as tr

import scatter.scatterMDA as scatterMDA
import MDAnalysis
import numpy
import argparse

#default values for optional arguments
dv={'asel':'*',
    'qvecs':'0.65 1.3 6.48',
    'nt':200,
    'dt':1,
    'ns':100
}
def pHelp(prefix, key):
    return prefix+', default='+str(dv[key])
parser = argparse.ArgumentParser(description='''Create I(Q,t)''')
parser.add_argument('inpsf', type=str, help='input PSF file')
parser.add_argument('indcd', type=str, help='input DCD file')
parser.add_argument('outIqt', type=str, help='output file containing I(Q,t)')
parser.add_argument('--asel', type=str, default=dv['asel'], help=pHelp('atom selection','asel'))
parser.add_argument('--qvecs', type=str, default=dv['qvecs'], help=pHelp('Q vectors','qvecs'))
parser.add_argument('--nt', type=int, default=dv['nt'], help=pHelp('number of time-points in I(Q,t)','nt'))
parser.add_argument('--dt', type=int, default=dv['dt'], help=pHelp('time unit in I(Q,t)','dt'))
parser.add_argument('--ns', type=int, default=dv['ns'], help=pHelp('number of t_0 sample-points','ns'))
args = parser.parse_args()

u = MDAnalysis.Universe(args.inpsf,args.indcd)
atom_selection = u.selectAtoms('name {0}'.format(args.asel))

print '{0} frames and {1} atoms'.format(len(u.trajectory),len(atom_selection))
#extract coordinates for selected atoms
frames=numpy.zeros([len(u.trajectory),len(atom_selection),3])
counter=0
for ts in u.trajectory:
    frames[counter]=atom_selection.coordinates()
    counter+=1

#calculate I(Q,t)
b=numpy.ones(len(atom_selection)) #scattering lengths
QQ=numpy.array([float(q) for q in args.qvecs.split()])

results = scatterMDA.II(frames, b, QQ, nsampling=args.ns, nt=args.nt, dt=args.dt)
open(args.outIqt,'w').write(results['buf'])

           
