import argparse
import numpy
from pdb import set_trace as tr

header = 'Time Temp Press Volume Density TotEng E_vdwl E_coul E_mol Lx'
items = header.split()
nitems =  len(items)
averages = numpy.array([0.0] * nitems)
footnote = 'Loop time'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Calculate thermodynamic averages from the log file.
    Example: thermo_averages.py equil.log""")
    parser.add_argument('log', type=str, help='log file')
    #parser.add_argument('--npol', type=int, default=defaults['npol'], 
    #                    help='number of styrene polymers in the simulation, default is {0}'.format(defaults['npol']))
    args=parser.parse_args()

#handle, tmpfile = mkstemp(prefix='junk', dir='/tmp')

#open(args.outf,'w').write(buf)
#tr()
handle = open(args.log)

# Find the header
line = handle.readline()
while header not in line:
    line = handle.readline()

# Read data
i = 0
line = handle.readline()
while footnote not in line and  line.strip():
    averages += numpy.array( [float(x) for x in line.split()] )    
    i += 1
    line = handle.readline()

averages /= i
print '#'+header
print ' '.join(['{0:9.3f}'.format(x) for x in averages])

