import numpy
from scipy.stats import linregress

datadir='/projects/development/ModeratorTzero/13087_ModeratorTzero_Extend/test'

#Load data and a mask that filters all detectors but a region ~ meters from the sample
Load(Filename=datadir+'/CORELLI_11326.nxs.h5', OutputWorkspace='CORELLI_11326')
LoadMask(Instrument='CORELLI', InputFile=datadir+'/roi.xml', OutputWorkspace='roi')
MaskDetectors(Workspace='CORELLI_11326', MaskedWorkspace='roi')

#Pile onto a single histogram and apply ModeratorTzero to the events workspace
SumSpectra(InputWorkspace='CORELLI_11326',OutputWorkspace='events')
ModeratorTzero(InputWorkspace='events', Emode='Elastic', OutputWorkspace='events_shifted')

#Apply ModeratorTzero to the histogram workspace
Rebin(InputWorkspace='events', Params=[1,4,17000], OutputWorkspace='histo')
ModeratorTzero(InputWorkspace='histo', Emode='Elastic', OutputWorkspace='histo_shifted')

#Visually compare with original spectrum
#Plot histo_smoothed and histo_shifted_smoothed together
SmoothData(InputWorkspace='histo', Npoints=20, OutputWorkspace='histo_smoothed')
SmoothData(InputWorkspace='histo_shifted', Npoints=20, OutputWorkspace='histo_shifted_smoothed')

#Visually compare application of ModeratorTzero to events and histogram workspaces:
#Plot histo_shifted_smoothed and histo_shifted_2_smoothed together. The plots
#should be "almost" identical.
Rebin(InputWorkspace='events_shifted', Params=[1,4,17000], OutputWorkspace='histo_shifted_2')
SmoothData(InputWorkspace='histo_shifted', Npoints=20, OutputWorkspace='histo_shifted_2_smoothed')

#Extract the t0 values by comparing the bin values before and after applying ModeratorTzero
h=mtd['histo']
tof_initial=h.dataX(0)
hs=mtd['histo_shifted']
tof_final=hs.dataX(0)
t0=tof_initial - tof_final

#Compare the extracted t0 values with the theoretical expectation.
#Given a time-of-flight tof, the incident energy Ei=a*(L1+L2/tof)^2 with
# a=5.22704e+06  (conversion factor)
# L1=20 (distance from source to sample, in meters)
# L2=2.56 (distance from sample to detectors, in meters)
Ei=5.22704e+06*(22.56/tof_final)**2
#The relation between t0 and incident energy stored in the CORELLI
#parameter file is t0 = 101.9*Ei^(-0.41) * exp(-Ei/282.0))
t0_theoretical=101.9*Ei**(-0.41) * numpy.exp(-Ei/282.0)
#A linear regression of t0 and t0_theoretical should give a slope~1
#and intercept ~0micro-seconds
slope,intercept=linregress(t0,t0_theoretical)[0:2]
print 'slope=',slope,'  intercept=',intercept
