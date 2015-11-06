import numpy
import argparse
from pdb import set_trace as tr

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Obtain MSD for a given instrument
    Example: internalMSD.py msd msdcm msdint""")
    parser.add_argument('msd', help='(input) two column file containing time(ns) and total MSD(t)')
    parser.add_argument('msdcm', help='(input) two column file containing time(ns) and center of mass MSD(t)')
    parser.add_argument('msdint', help='(output) two column file containing time(ns) and internal MSD(t)')
    pargs=parser.parse_args()

    times = numpy.loadtxt(pargs.msd, usecols=(0,))
    othertimes = numpy.loadtxt(pargs.msdcm, usecols=(0,))
    if not (times==othertimes).all():
        raise 'times are different for the two input files'
    msd = numpy.loadtxt(pargs.msd, usecols=(1,))
    msdcm = numpy.loadtxt(pargs.msdcm, usecols=(1,))
    msdint = msd - msdcm
    buf='# time(ps)   MSDinteral\n'
    for i in range(len(times)):
        buf += '{0} {1}\n'.format(times[i], msdint[i])
    open(pargs.msdint, 'w').write(buf)

