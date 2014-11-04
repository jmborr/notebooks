from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, DSFinterp, Fit

T=290

rootd='/projects/research/LiCl/EugeneData'
LoadNexus(Filename=rootd+'/elasticLine.nxs', OutputWorkspace='elasticLine')
LoadNexus(Filename=rootd+'/LiCl_290K.nxs', OutputWorkspace='LiCl_290K')

rootd='/projects/research/LiCl/watBox30/Hq'
qdirs='Q32 Q34 Q36 Q38 Q40 Q42 Q421 Q422 Q423 Q424 Q425 Q426 Q427 Q428 Q429 Q430 Q431 Q432 Q433 Q434 Q435 Q436 Q437 Q438 Q439 Q44 Q46 Q48 Q50'
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

# Convolve with the elastic line
for qdir in qdirs.split():
	ConvolveWorkspaces(Workspace1='elasticLine', Workspace2='{0}_sqw'.format(qdir), OutputWorkspace='conv{0}'.format(qdir))

#####################
#### DSFinterp algorithm  ####
#####################
# Input
parametervalues='0.32 0.34 0.36 0.38 0.40 0.42 0.421 0.422 0.423 0.424 0.425 0.426 0.427 0.428 0.429 0.430 0.431 0.432 0.433 0.434 0.435 0.436 0.437 0.438 0.439 0.44 0.46 0.48 0.50'
parmlist=[float(x) for x in parametervalues.split()]
workspaces=[]
for qdir in qdirs.split():
	workspaces.append('conv{0}'.format(qdir))
# Output
outparmlist=[]
outputworkspaces=[]
x=390
while x < 461:
	outparmlist.append(x*0.001)
	outputworkspaces.append('intQ{0}'.format(x))
	x+=2
N=len(outparmlist)
#Run the algorithm
DSFinterp(Workspaces=workspaces, ParameterValues=parmlist, OutputWorkspaces=outputworkspaces, LoadErrors=0, LocalRegression=1,RegressionType='quadratic', RegressionWindow=6, TargetParameters=outparmlist)

for wi in range(4):
	for i in range(N):
		Hq=outparmlist[i]
		workspace=outputworkspaces[i]
		fit_string='name=TabulatedFunction,Workspace=elasticLine,WorkspaceIndex={0},Scaling=0.8;'.format(wi) +\
		'name=TabulatedFunction,Workspace={0},WorkspaceIndex={1},Scaling=0.2;'.format(workspace,wi) +\
		'name=LinearBackground,A0=0.0,A1=0.0'
		Fit(fit_string, InputWorkspace='LiCl_{0}K'.format(T), WorkspaceIndex=wi, StartX=-0.1, EndX=0.1, CreateOutput=1)
		# print the optimized K and the associated Chi-square of the fit
		ws=mtd['LiCl_{0}K_Parameters'.format(T)]
		print 'Q=',0.3+0.2*wi, 'Hq=',  Hq, 'Chi2=', ws.row(4)['Value']
