from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, DSFinterp, Fit

rootd='/projects/research/LiCl/EugeneData'
LoadNexus(Filename=rootd+'/elasticLine.nxs', OutputWorkspace='elasticLine')
LoadNexus(Filename=rootd+'/LiCl_290K.nxs', OutputWorkspace='LiCl_290K')


rootd='/projects/research/LiCl/watBox30/Hq'
qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q421 Q422 Q423 Q424 Q425 Q426 Q427 Q428 Q429 Q430 Q431 Q432 Q433 Q434 Q435 Q436 Q437 Q438 Q439 Q44 Q46 Q48 Q50'
#qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q44 Q46 Q48 Q50'
for qdir in qdirs.split():
	LoadSassena(Filename=rootd+'/{0}/T290/production/fqt_inc.hd5'.format(qdir), TimeUnit=1.0, OutputWorkspace=qdir)
	# Rebin in momentum transfer and normalize to one
	Transpose(InputWorkspace='{0}_fqt.Re'.format(qdir), OutputWorkspace='{0}_fqt.Re'.format(qdir))
	Rebin(InputWorkspace='{0}_fqt.Re'.format(qdir), Params=[0.2,0.2,1.0], OutputWorkspace='{0}_fqt.Re'.format(qdir))
	Transpose(InputWorkspace='{0}_fqt.Re'.format(qdir), OutputWorkspace='{0}_fqt.Re'.format(qdir))
	factor = 1./22 # renormalize after rebin in momentum transfer
	if qdir in 'Q32 Q34 Q36 Q38 Q40 Q42 Q44 Q46 Q48 Q50': factor = 1./2 # These simulations had different number of momentum 
	Scale(InputWorkspace='{0}_fqt.Re'.format(qdir), Factor=factor, Operation='Multiply', OutputWorkspace='{0}_fqt.Re'.format(qdir))
	SassenaFFT(InputWorkspace=qdir, FFTonlyRealPart=1, DetailedBalance=1, Temp=290.0)
	Scale(InputWorkspace='{0}_sqw'.format(qdir), Factor=1e-07, Operation='Multiply', OutputWorkspace='{0}_sqw'.format(qdir))
	Rebin(InputWorkspace='{0}_sqw'.format(qdir), Params=[-0.2, 0.0004, 0.2], OutputWorkspace='{0}_sqw'.format(qdir))
	#Let's remove the points close to the elastic line
	
# Convolve with the elastic line
qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q421 Q422 Q423 Q424 Q425 Q426 Q427 Q428 Q429 Q430 Q431 Q432 Q433 Q434 Q435 Q436 Q437 Q438 Q439 Q44 Q46 Q48 Q50'
#qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q44 Q46 Q48 Q50'
for qdir in qdirs.split():
	ConvolveWorkspaces(Workspace1='elasticLine', Workspace2='{0}_sqw'.format(qdir), OutputWorkspace='conv{0}'.format(qdir))

#####################
#### DSFinterp algorithm  ####
#####################
#parametervalues='0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50'
parametervalues='0.32 0.34 0.36 0.38 0.40 0.42 0.421 0.422 0.423 0.424 0.425 0.426 0.427 0.428 0.429 0.430 0.431 0.432 0.433 0.434 0.435 0.436 0.437 0.438 0.439 0.44 0.46 0.48 0.50'
parmlist=[float(x) for x in parametervalues.split()]
#qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q44 Q46 Q48 Q50'
qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q421 Q422 Q423 Q424 Q425 Q426 Q427 Q428 Q429 Q430 Q431 Q432 Q433 Q434 Q435 Q436 Q437 Q438 Q439 Q44 Q46 Q48 Q50'
workspaces=[]
outputworkspaces=[]
for qdir in qdirs.split():
	workspaces.append('conv{0}'.format(qdir))
	outputworkspaces.append('int{0}'.format(qdir))
DSFinterp(Workspaces=workspaces, ParameterValues=parmlist, OutputWorkspaces=outputworkspaces, LoadErrors=0, LocalRegression=1,RegressionType='quadratic', RegressionWindow=6, TargetParameters=parmlist)

# Prepare the input
parametervalues='0.32 0.34 0.36 0.38 0.40 0.42 0.421 0.422 0.423 0.424 0.425 0.426 0.427 0.428 0.429 0.430 0.431 0.432 0.433 0.434 0.435 0.436 0.437 0.438 0.439 0.44 0.46 0.48 0.50'
qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q421 Q422 Q423 Q424 Q425 Q426 Q427 Q428 Q429 Q430 Q431 Q432 Q433 Q434 Q435 Q436 Q437 Q438 Q439 Q44 Q46 Q48 Q50'
#parametervalues='0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50'
#qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q44 Q46 Q48 Q50'
workspaces=''
for qdir in qdirs.split():
	workspaces += ' conv{0}'.format(qdir)
guess=0.40  # initial charge
for wi in range(4):
	#Fitting model: S(Q,E) = A*elastic_line + B*Simulated_structure_factor + Linear_background
	fit_string = 'name=TabulatedFunction,Workspace=elasticLine,WorkspaceIndex={0},Scaling=1,constraints=(0.0001<Scaling);'.format(wi)+\
		'name=DSFinterp1DFit,InputWorkspaces="{0}",ParameterValues="{1}",'.format(workspaces,parametervalues) +\
		'LoadErrors=0,LocalRegression=1,RegressionType=quadratic,RegressionWindow=6,' +\
		'WorkspaceIndex={0},Intensity=1.0,TargetParameter={1},'.format(wi,guess) +\
		'constraints=(0.0001<Intensity);' +\
		'name=LinearBackground,A0=0.0,A1=0.0'
	# Fit for WorkspaceIndex=8 corresponding to highest measured Q  (Q=19nm^(-1) or L~3Anstromgs)
	Fit(fit_string, InputWorkspace='LiCl_290K', WorkspaceIndex=wi, StartX=-0.1, EndX=0.1, CreateOutput = 1 )
	# print the optimized K and the associated Chi-square of the fit
	ws=mtd['LiCl_290K_Parameters']
	print 'Q=',0.3+0.2*wi, 'Hq=', ws.row(2)['Value'], 'Elastic=', ws.row(0)['Value'], 'Quasielastic=', ws.row(1)['Value'],'Chi2=', ws.row(5)['Value']