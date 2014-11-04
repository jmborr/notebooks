from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, DSFinterp, Fit

T=200

rootd='/projects/research/LiCl/EugeneData'
LoadNexus(Filename=rootd+'/elasticLine.nxs', OutputWorkspace='elasticLine')
LoadNexus(Filename=rootd+'/LiCl_{0}K.nxs'.format(T), OutputWorkspace='LiCl_{0}K'.format(T))

rootd='/projects/research/LiCl/watBox30/Hq'
qdirs='Q390 Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450 Q460 Q470 Q480 Q490 Q500'
for qdir in qdirs.split():
	LoadSassena(Filename=rootd+'/{0}/T{1}/production/fqt_inc.hd5'.format(qdir,T), TimeUnit=1.0, OutputWorkspace=qdir)
	SassenaFFT(InputWorkspace=qdir, FFTonlyRealPart=1, DetailedBalance=1, Temp=T)
	Scale(InputWorkspace='{0}_sqw'.format(qdir), Factor=1e-07, Operation='Multiply', OutputWorkspace='{0}_sqw'.format(qdir))
	Rebin(InputWorkspace='{0}_sqw'.format(qdir), Params=[-0.2, 0.0004, 0.2], OutputWorkspace='{0}_sqw'.format(qdir))
	
# Convolve with the elastic line
qdirs='Q390 Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450 Q460 Q470 Q480 Q490 Q500'
for qdir in qdirs.split():
	ConvolveWorkspaces(Workspace1='elasticLine', Workspace2='{0}_sqw'.format(qdir), OutputWorkspace='conv{0}'.format(qdir))

#####################
#### DSFinterp algorithm  ####
#####################
parametervalues='0.390 0.400 0.402 0.404 0.406 0.408 0.410 0.412 0.414 0.416 0.420 0.422 0.424 0.426 0.428 0.430 0.432 0.440 0.450 0.460 0.470 0.480 0.490 0.500'
parmlist=[float(x) for x in parametervalues.split()]
qdirs='Q390 Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450 Q460 Q470 Q480 Q490 Q500'
workspaces=[]
outputworkspaces=[]
for qdir in qdirs.split():
	workspaces.append('conv{0}'.format(qdir))
	outputworkspaces.append('int{0}'.format(qdir))
DSFinterp(Workspaces=workspaces, ParameterValues=parmlist, OutputWorkspaces=outputworkspaces, LoadErrors=0, LocalRegression=1,RegressionType='quadratic', RegressionWindow=6, TargetParameters=parmlist)

# Prepare the input
T=200
parametervalues='0.390 0.400 0.402 0.404 0.406 0.408 0.410 0.412 0.414 0.416 0.420 0.422 0.424 0.426 0.428 0.430 0.432 0.440 0.450 0.460 0.470 0.480 0.490 0.500'
#parametervalues='0.390 0.400 0.402 0.404 0.406 0.408 0.410 0.412 0.414 0.416 0.420 0.422 0.424 0.426 0.428 0.430 0.432 0.440 0.450'
#parametervalues='0.390 0.400 0.410 0.420 0.430 0.440 0.450 0.460 0.470 0.480 0.490 0.500'
qdirs='Q390 Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450 Q460 Q470 Q480 Q490 Q500'
#qdirs='Q390 Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450'
#qdirs='Q390 Q400 Q410 Q420 Q430 Q440 Q450 Q460 Q470 Q480 Q490 Q500'
workspaces=''
for qdir in qdirs.split():
	workspaces += ' conv{0}'.format(qdir)
guess=0.425  # initial charge
for wi in range(4):
	#Fitting model: S(Q,E) = A*elastic_line + B*Simulated_structure_factor + Linear_background
	fit_string = 'name=TabulatedFunction,Workspace=elasticLine,WorkspaceIndex={0},Scaling=0.5,constraints=(0.0001<Scaling);'.format(wi)+\
		'name=DSFinterp1DFit,InputWorkspaces="{0}",ParameterValues="{1}",'.format(workspaces,parametervalues) +\
		'LoadErrors=0,LocalRegression=1,RegressionType=quadratic,RegressionWindow=6,' +\
		'WorkspaceIndex={0},Intensity=1.0,TargetParameter={1},'.format(wi,guess) +\
		'constraints=(0.0001<Intensity);' +\
		'name=LinearBackground,A0=0.0,A1=0.0'
	# Fit for WorkspaceIndex=8 corresponding to highest measured Q  (Q=19nm^(-1) or L~3Anstromgs)
	Fit(fit_string, InputWorkspace='LiCl_{0}K'.format(T), WorkspaceIndex=wi, StartX=-0.1, EndX=0.1, CreateOutput = 1 )
	# print the optimized K and the associated Chi-square of the fit
	ws=mtd['LiCl_{0}K_Parameters'.format(T)]
	print 'Q=',0.3+0.2*wi, 'Hq=', ws.row(2)['Value'], 'Elastic=', ws.row(0)['Value'], 'Quasielastic=', ws.row(1)['Value'],'Chi2=', ws.row(5)['Value']