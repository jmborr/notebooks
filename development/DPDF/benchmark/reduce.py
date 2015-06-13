#Script input parameters
wd='/projects/development/DPDF/benchmark'
dE=1.0 #1meV binning in energy
dtheta=0.15 #binning in theta angle. Advised is 0.2degrees
Qmin=0.0
Qmax=15.0
dQ=0.1 #0.1A^{-1} binning Q

#Load several event files into a sinle workspace. The nominal incidentenergy should be
#the same to avoid difference in energy resolution
event_files=['ARCS_56415_event.nxs',
    'ARCS_56416_event.nxs',
    'ARCS_56417_event.nxs',
    'ARCS_56418_event.nxs',
    'ARCS_56419_event.nxs',
    'ARCS_56420_event.nxs'    
    ]
event_files=[wd+'/'+x for x in event_files] #create absolute path directory
Load(Filename='+'.join(event_files),OutputWorkspace='events_sum')

#Load the vanadium file, assume to be preprocessed, meaning that for every detector
#all events whithin a particular wide wavelength range have been rebinned into a
#single histogram
Load(Filename='van56293.nxs', OutputWorkspace='van')

#Retrieve the mask from the vanadium workspace, and apply it to the data
MaskDetectors(Workspace='events_sum', MaskedWorkspace='van')

#Obtain incident energy as the mean of the nominal Ei values. There is one nominal value per events file
w=mtd['events_sum']
Ei=w.getRun()['EnergyRequest'].getStatistics().mean
Ei_min=-0.5*Ei #These are
Ei_max=0.95*Ei

#Convert to energy transger. The output workspace is S(detector-id,E)
Erange=','.join([str(x) for x in [Ei_min,dE,Ei_max]])
DgsReduction(SampleInputWorkspace='events_sum', EnergyTransferRange=Erange, OutputWorkspace='reduced')

#Convert to S(theta,E)
GenerateGroupingPowder(InputWorkspace='reduced', AngleStep=dtheta, GroupingFilename=wd+'/detectors_to_angle.xml')
GroupDetectors(InputWorkspace='reduced', MapFile=wd+'/detectors_to_angle.xml', OutputWorkspace='ste')

#Normalize by the vanadium intensity, but before that we need S(theta) for the vanadium. Recall
#every detector has all energies into a single bin, so we get S(theta) instead of S(theta,E)
GroupDetectors(InputWorkspace='van', MapFile=wd+'/detectors_to_angle.xml', OutputWorkspace='st_van')
Divide('ste','st_van',OutputWorkspace='ste_normalized')

#Mask the problematic angles
#MaskAngle()

#Convert S(\theta,E) to S(Q,E), then rebin to 
ConvertToMD(InputWorkspace='ste_normalized', QDimensions='|Q|', dEAnalysisMode='Direct', OutputWorkspace='sqe')
Qrange='|Q|,{0},{1}.{2}'.format(Qmin,Qmax,int( (Qmax-Qmin)/dQ ))
deltaErange='DeltaE,{0},{1}.{2}'.format(Ei_min, Ei_max, int( (Ei_max-Ei_min)/dE ))
BinMD(InputWorkspace='sqe', AxisAligned=1, AlignedDim0='|Q|,0,15,75',
    AlignedDim1=deltaErange, OutputWorkspace='sqe_binned')

#Slice the data by transforming to a Matrix2Dworkspace, with deltaE along the vertical axis
ConvertMDHistoToMatrixWorkspace(InputWorkspace='sqe_binned',
    Normalization='NumEventsNormalization', OutputWorkspace='sqe_binned_matrix')