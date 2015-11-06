import numpy
import argparse
from math import sqrt,exp,pi
from pdb import set_trace as tr

std = {'BASIS':460., 'HFBS':1940. } #standar deviations of the resolution function, in picoseconds

def resfunc(instrument):
    def func(x):
        ss = std[instrument]**2
        return 1/sqrt(2*pi*ss) * exp(-x*x/(2*ss))
    return func

def elastic(instrument, times):
    FWHM = std[instrument]*2.355
    times_reduced = times/FWHM
    minval = numpy.min( times_reduced[numpy.nonzero(times_reduced)] )
    times_reduced[ times == 0 ] = minval/1000  #substitue zeroes with smallish values
    return numpy.sin(times_reduced)/times_reduced #will yield error if times_reduces contains zeroes
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Obtain MSD for a given instrument
    Example: benedettoMSD.py msdt instrument""")
    parser.add_argument('msdt', help='two column file containing time(ns) and MSD(t)')
    parser.add_argument('instrument', help='Either "BASIS" or "HFBS"')
    parser.add_argument('--elastic', type=str, default='no', help='elastic method')
    pargs=parser.parse_args()

    times = numpy.loadtxt(pargs.msdt, usecols=(0,))
    msd_vals = numpy.loadtxt(pargs.msdt, usecols=(1,))
    func = resfunc(pargs.instrument)
    res_vals = [func(t) for t in times]

    if pargs.elastic.lower()[0]=='y':
        print 2*numpy.sum(elastic(pargs.instrument,times)*msd_vals*res_vals)
    else:
        print 2*numpy.sum(msd_vals*res_vals) #print MSD
