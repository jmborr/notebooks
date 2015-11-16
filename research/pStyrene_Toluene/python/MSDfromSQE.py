'''
Created on Nov 11, 2015

@author: jbq
'''

import os
import numpy
import argparse
import mantid.simpleapi as mtds
#import mantid as mtd
from scipy.stats import linregress
from pdb import set_trace as tr

parser = argparse.ArgumentParser(description="""Generate MSD from simulated S(Q,E)
Example: MSDfromSQE.py told3\n
Looks for file told3_sqw_HFBS.nxs in directory mantid_sqe/, calculates the elastic
intensity for HFBS resolution, and does the fit EIE(Q)~exp(-MSD Q^2/3)
Outputs file msd_from_sqe_HFBS.dat containing EIE(Q^2) and prints out the MSD value.""")
parser.add_argument('deutsche', type=str, help='deuteration scheme')
pargs=parser.parse_args()
deutscheme=pargs.deutsche

# Load S(Q,E) convolved with HFBS model resolution
sqe=mtds.LoadNexus('./mantid_sqe/{0}_sqw_HFBS.nxs'.format(deutscheme))

FWHM={'HFBS':0.8E-03, 'BASIS':3.4e-03} #in meV
iei=mtds.Integration(sqe, RangeLower=-FWHM['HFBS'], RangeUpper=FWHM['HFBS'], IncludePartialBins=1)
iei=mtds.Transpose(iei)

# Pick only small Q values when the gaussian e^(-MSD Q^2) has validity
#Qmin=0.36; Qmax = 1.51  #[0.36, 1.51] is the Q-range used in HFBS to obtain the MSD
Qmin=0.25; Qmax=0.95
iei=mtds.Rebin(iei,[Qmin,0.1,Qmax])

# Linear regression of IEI versus Q^2
q2=numpy.power(iei.dataX(0), 2) # Q^2
ieiY=iei.dataY(0)
#Prepend point for Q=0
mtds.LoadSassena('fqt_inc_{0}.h5'.format(deutscheme), OutputWorkspace='fq')
fq0=mtds.mtd['fq_fq0']
I_0=fq0.dataY(0)[0]
q2=numpy.insert(q2,0,[0.0,])
ieiY=numpy.insert(ieiY,0,[I_0,])
liei=numpy.log(ieiY)
slope, intercept, r_value, p_value, std_err = linregress(q2, liei)

#Print results and save to file
msd = -3*slope
open(os.getcwd()+'/{0}_msd_from_SQE.dat'.format(deutscheme), 'w').write(str(msd))
buf='# Q^2  IEI\n'
for i in range(len(q2)):
    buf += '{0} {1}\n'.format(q2[i], ieiY[i])
open(os.getcwd()+'/{0}_IEI_from_SQE.dat'.format(deutscheme), 'w').write(buf)
print 'Files generated:', '{0}_msd_from_SQE.dat'.format(deutscheme), '{0}_IEI_from_SQE.dat'.format(deutscheme)
