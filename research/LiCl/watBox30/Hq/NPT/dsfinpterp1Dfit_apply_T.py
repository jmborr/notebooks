from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, DSFinterp, Fit

T=270

rootd='/projects/research/LiCl/EugeneData'
LoadNexus(Filename=rootd+'/elasticLine.nxs', OutputWorkspace='elasticLine')
LoadNexus(Filename=rootd+'/LiCl_{0}K.nxs'.format(T), OutputWorkspace='LiCl_{0}K'.format(T))

rootd='/projects/research/LiCl/watBox30/Hq/NPT'
parametervalues='0.390 0.400 0.402 0.404 0.406 0.408 0.410 0.412 0.414 0.416 0.420 0.422 0.424 0.426 0.428 0.430 0.432 0.440 0.450 0.460 0.470 0.480 0.490 0.500'
qdirs='Q390 Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450 Q460 Q470 Q480 Q490 Q500'

for qdir in qdirs.split():
	LoadSassena(Filename=rootd+'/{0}/T{1}/production/fqt_rms2first_inc.hd5'.format(qdir,T), TimeUnit=1.0, SortByQVectors=1, OutputWorkspace=qdir)
	#Set I(Q,t=0)=1 so that S(Q,E) is normalized to one
	SassenaFFT(InputWorkspace=qdir, FFTonlyRealPart=1, DetailedBalance=1, Temp=T)
	Scale(InputWorkspace='{0}_sqw'.format(qdir), Factor=1e-07, Operation='Multiply', OutputWorkspace='{0}_sqw'.format(qdir))
	Rebin(InputWorkspace='{0}_sqw'.format(qdir), Params=[-0.2, 0.0004, 0.2], OutputWorkspace='{0}_sqw'.format(qdir))
	
# Convolve with the elastic line
for qdir in qdirs.split():
	ConvolveWorkspaces(Workspace1='elasticLine', Workspace2='{0}_sqw'.format(qdir), OutputWorkspace='conv{0}'.format(qdir))

# Prepare the input
workspaces=''
for qdir in qdirs.split():
	workspaces += ' conv{0}'.format(qdir)
guess=0.43  # initial hydrogen charge
optimal_Hq=[]
print '# T   Q    Hq  Elastic QuasiE  Chi2'
# Fit each Q-value independently
for wi in range(4):
	#Fitting model: S(Q,E) = A*elastic_line + B*Simulated_structure_factor + Linear_background
	fit_string = 'name=TabulatedFunction,Workspace=elasticLine,WorkspaceIndex={0},Scaling=0.5,constraints=(0.0001<Scaling);'.format(wi)+\
		'name=DSFinterp1DFit,InputWorkspaces="{0}",ParameterValues="{1}",'.format(workspaces,parametervalues) +\
		'LoadErrors=0,LocalRegression=1,RegressionType=quadratic,RegressionWindow=6,' +\
		'WorkspaceIndex={0},Intensity=1.0,TargetParameter={1},'.format(wi,guess) +\
		'constraints=(0.0001<Intensity);' +\
		'name=LinearBackground,A0=0.0,A1=0.0'
	Fit(fit_string, InputWorkspace='LiCl_{0}K'.format(T), WorkspaceIndex=wi, StartX=-0.1, EndX=0.1, CreateOutput = 1 )
	# print the optimized K and the associated Chi-square of the fit
	ws=mtd['LiCl_{0}K_Parameters'.format(T)]
	optimal_Hq.append(ws.row(2)['Value'])
	print ' {0:3d} {1:3.1f} {2:6.4f} {3:6.4f} {4:6.4f} {5:6.3f}'.format(T, 0.3+0.2*wi, ws.row(2)['Value'], ws.row(0)['Value'], ws.row(1)['Value'], ws.row(5)['Value'])

# Now that we have obtained the optimal values of Hq for each Q-value,
# create the intermediate structure factor for each optimized Hq
# and save them into workspaces optimal_fqt_Q$Q
parametervalues='0.390 0.400 0.402 0.404 0.406 0.408 0.410 0.412 0.414 0.416 0.420 0.422 0.424 0.426 0.428 0.430 0.432 0.440 0.450 0.460 0.470 0.480 0.490 0.500'
qdirs='Q390 Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450 Q460 Q470 Q480 Q490 Q500'
parmlist=[float(x) for x in parametervalues.split()]
workspaces=['{0}_fqt'.format(x) for x in qdirs.split()]
print '#  Q   Hq'
for wi in range(4):
	Q=0.3+0.2*wi
	outparmlist=[ optimal_Hq[wi], ]
	outputworkspaces=[  'optimal_fqt_Q{0:3.1f}'.format(Q), ]
	for qdir in qdirs.split():
		ExtractSingleSpectrum(InputWorkspace='{0}_fqt.Re'.format(qdir), WorkspaceIndex=wi, OutputWorkspace='{0}_fqt'.format(qdir))
		Scale(InputWorkspace='{0}_fqt'.format(qdir), Factor=0.00000069, Operation='Multiply', OutputWorkspace='{0}_fqt'.format(qdir))
		Rebin(InputWorkspace='{0}_fqt'.format(qdir), Params=[5,1,9000], OutputWorkspace='{0}_fqt'.format(qdir))
	print '{0:3.1f} {1:5.3f}'.format(Q,optimal_Hq[wi])
	DSFinterp(Workspaces=workspaces, ParameterValues=parmlist, LoadErrors=0, LocalRegression=0,RegressionType='quadratic', RegressionWindow=6, OutputWorkspaces=outputworkspaces, TargetParameters=outparmlist)
	SaveAscii(InputWorkspace='optimal_fqt_Q{0:3.1f}'.format(Q), Filename='/tmp/optimal_fqt_Q{0:3.1f}'.format(Q))

	
