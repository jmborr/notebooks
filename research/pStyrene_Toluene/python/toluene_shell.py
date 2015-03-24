import argparse
from amber.amber14.cpptraj.command import exec_cpptraj
from tempfile import mkstemp
from pdb import set_trace as tr

nsty_per_polymer=32 # number of residues in each styrene polymer
defaults={'spacing':0.2, 'minimum':3.0, 'maximum':10.0, 'npol':1}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Number of toluenes around a styrene monomer as function of distance.
    Example: radial.py 8styrene32.pdb equil_rms2first.dcd C3 C --npol 1 --spacing 0.2 --maximum 10.0""")
    parser.add_argument('top', type=str, help='topology or PDB file')
    parser.add_argument('dcd', type=str, help='DCD trajectory file')
    parser.add_argument('sty', type=str, help='atom name for styrene residue')
    parser.add_argument('tol', type=str, help='atom name in toluene molecule')
    parser.add_argument('outf', type=str, help='Output file with number of toluenes')
    parser.add_argument('--npol', type=int, default=defaults['npol'], 
                        help='number of styrene polymers in the simulation, default is {0}'.format(defaults['npol']))
    parser.add_argument('--spacing', type=float, default=defaults['spacing'],
                        help='binning, default is {0}'.format(defaults['spacing']))
    parser.add_argument('--minimum', type=float, default=defaults['minimum'],
                        help='minimum distance, default is {0}'.format(defaults['minimum']))
    parser.add_argument('--maximum', type=float, default=defaults['maximum'],
                        help='maximum distance, default is {0}'.format(defaults['maximum']))
    args=parser.parse_args()

handle, tmpfile = mkstemp(prefix='junk', dir='/tmp')
nsty = nsty_per_polymer * args.npol
solvmask = '(!:1-{0})&(@{1})'.format(nsty, args.tol)
solmask = '(:1-{0})&(@{1})'.format(nsty, args.sty)
lower = args.minimum
upper = lower + args.spacing
buf='# distance  number_of_tol/sty\n'

while upper < args.maximum:
    script='''trajin {0}
center {1}
watershell {2} {3} lower {4} upper {5} {6}'''.format(args.dcd, solmask, solmask, tmpfile, lower, upper, solvmask)

    exec_cpptraj(args.top, script)
    #tr()
    avgn1=0
    avgn2=0
    for line in open(tmpfile).readlines()[1:]:
        frame, n1, n2 = [float(x) for x in line.split()]
        avgn1 += n1
        avgn2 += n2
    avgn1 /= (frame*nsty)
    avgn2 /= (frame*nsty)
    buf += '{0:5.2f} {1:6.2f}\n{2:5.2f} {3:6.2f}\n'.format(lower, avgn1, upper, avgn2)
    lower = upper + args.spacing
    upper = lower + args.spacing

open(args.outf,'w').write(buf)

