from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, DSFinterp, Fit

T=230

rootd='/projects/research/LiCl/EugeneData'
LoadNexus(Filename=rootd+'/elasticLine.nxs', OutputWorkspace='elasticLine')
LoadNexus(Filename=rootd+'/LiCl_{0}K.nxs'.format(T), OutputWorkspace='LiCl_{0}K'.format(T))

rootd='/projects/research/LiCl/watBox30/Vickie'
qdirs='Q0.417 Q0.421023 Q0.42117 Q0.425276 Q0.429528 Q0.42971 Q0.43229 Q0.433444 Q0.4340506 Q0.436612 Q0.4378223 Q0.438391 Q0.4422006 Q0.4455 Q0.45 Q0.4545'
for qdir in qdirs.split():
	LoadSassena(Filename=rootd+'/T{0}/{1}/fqt_inc.hd5'.format(T,qdir), TimeUnit=1.0, OutputWorkspace=qdir)
	SortByQVectors(InputWorkspace=qdir)
	SassenaFFT(InputWorkspace=qdir, FFTonlyRealPart=1, DetailedBalance=1, Temp=T)
	Scale(InputWorkspace='{0}_sqw'.format(qdir), Factor=1e-07, Operation='Multiply', OutputWorkspace='{0}_sqw'.format(qdir))
	Rebin(InputWorkspace='{0}_sqw'.format(qdir), Params=[-0.2, 0.0004, 0.2], OutputWorkspace='{0}_sqw'.format(qdir))

# Convolve with the elastic line
for qdir in qdirs.split():
	ConvolveWorkspaces(Workspace1='elasticLine', Workspace2='{0}_sqw'.format(qdir), OutputWorkspace='conv{0}'.format(qdir))

#####################
#### DSFinterp algorithm  ####
#####################
# Input
parametervalues='0.417 0.421023 0.42117 0.425276 0.429528 0.42971 0.43229 0.433444 0.4340506 0.436612 0.4378223 0.438391 0.4422006 0.4455 0.45 0.4545'
parmlist=[float(x) for x in parametervalues.split()]
workspaces=[]
for qdir in qdirs.split():
	workspaces.append('conv{0}'.format(qdir))
# Output
outparmlist=[]
outputworkspaces=[]
x=418
while x < 454:
	outparmlist.append(x*0.001)
	outputworkspaces.append('intQ{0}'.format(x))
	x+=1
N=len(outparmlist)
#Run the algorithm
DSFinterp(Workspaces=workspaces, ParameterValues=parmlist, OutputWorkspaces=outputworkspaces, LoadErrors=0, LocalRegression=1,RegressionType='quadratic', RegressionWindow=6, TargetParameters=outparmlist)


for wi in range(4):
	f0_Scaling=0.131464
	f1_Scaling=1.66077
	for i in range(N):
		Hq=outparmlist[i]
		workspace=outputworkspaces[i]
		fit_string='name=TabulatedFunction,Workspace=elasticLine,WorkspaceIndex={0},Scaling={1};'.format(wi, f0_Scaling) +\
		'name=TabulatedFunction,Workspace={0},WorkspaceIndex={1},Scaling={2};'.format(workspace,wi, f1_Scaling) +\
		'name=LinearBackground,A0=0.0,A1=0.0'
		Fit(fit_string, InputWorkspace='LiCl_{0}K'.format(T), WorkspaceIndex=wi, StartX=-0.1, EndX=0.1, CreateOutput=1)
		# print the optimized K and the associated Chi-square of the fit
		ws=mtd['LiCl_{0}K_Parameters'.format(T)]
		print 'Q=',0.3+0.2*wi, 'Hq=',  Hq, 'Chi2=', ws.row(4)['Value']
