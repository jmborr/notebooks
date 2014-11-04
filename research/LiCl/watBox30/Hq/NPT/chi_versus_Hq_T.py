from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, DSFinterp, Fit

T=290

rootd='/projects/research/LiCl/EugeneData'
LoadNexus(Filename=rootd+'/elasticLine.nxs', OutputWorkspace='elasticLine')
LoadNexus(Filename=rootd+'/LiCl_{0}K.nxs'.format(T), OutputWorkspace='LiCl_{0}K'.format(T))

rootd='/projects/research/LiCl/watBox30/Hq/NPT'

parametervalues='0.390 0.400 0.402 0.404 0.406 0.408 0.410 0.412 0.414 0.416 0.420 0.422 0.424 0.426 0.428 0.430 0.432 0.440 0.450 0.460 0.470 0.480 0.490 0.500'
qdirs='Q390 Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450 Q460 Q470 Q480 Q490 Q500'

for qdir in qdirs.split():
	LoadSassena(Filename=rootd+'/{0}/T{1}/production/fqt_rms2first_inc.hd5'.format(qdir,T), TimeUnit=1.0, SortByQVectors=1, OutputWorkspace=qdir)
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
parmlist=[float(x) for x in parametervalues.split()]
workspaces=[]
for qdir in qdirs.split():
	workspaces.append('conv{0}'.format(qdir))
# Output
outparmlist=[]
outputworkspaces=[]
x=390
while x < 500:
	outparmlist.append(x*0.001)
	outputworkspaces.append('intQ{0}'.format(x))
	x+=2
N=len(outparmlist)
#Run the algorithm
DSFinterp(Workspaces=workspaces, ParameterValues=parmlist, OutputWorkspaces=outputworkspaces, LoadErrors=0, LocalRegression=1,RegressionType='quadratic', RegressionWindow=6, TargetParameters=outparmlist)

nQ=4
# Calculate fits for the original structure factors
buf=['{0:5.3f}'.format(x) for x in parmlist]
average_Chi2=[0]*len(parmlist)
for wi in range(nQ):
	for i in range(len(workspaces)):
		Hq=parmlist[i]
		workspace=workspaces[i]
		fit_string='name=TabulatedFunction,Workspace=elasticLine,WorkspaceIndex={0},Scaling=0.8;'.format(wi) +\
		'name=TabulatedFunction,Workspace={0},WorkspaceIndex={1},Scaling=0.2;'.format(workspace,wi) +\
		'name=LinearBackground,A0=0.0,A1=0.0'
		Fit(fit_string, InputWorkspace='LiCl_{0}K'.format(T), WorkspaceIndex=wi, StartX=-0.1, EndX=0.1, CreateOutput=1)
		# print the optimized K and the associated Chi-square of the fit
		ws=mtd['LiCl_{0}K_Parameters'.format(T)]
		chi2=ws.row(4)['Value']
		print 'Q=',0.3+0.2*wi, 'Hq=',  Hq, 'Chi2=', chi2
		buf[i]+=' {0:7.3f}'.format(chi2)
		average_Chi2[i]+=chi2
for i in range(len(workspaces)): buf[i]+=' {0:7.3f}'.format(average_Chi2[i]/nQ)
buf='#Hq Chi2(Q=0.3) Chi2(Q=0.5) Chi2(Q=0.7) Chi2(Q=0.9) average_Chi2\n'+'\n'.join(buf)+'\n'
open('{0}/chi_versus_Hq_T{1}_conv.dat'.format(rootd,T),'w').write(buf)

# Calculate fits for the interpolated structure factors
buf=['{0:5.3f}'.format(x) for x in outparmlist]
average_Chi2=[0]*len(outparmlist)
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
		chi2=ws.row(4)['Value']
		print 'Q=',0.3+0.2*wi, 'Hq=',  Hq, 'Chi2=', chi2
		buf[i]+=' {0:7.3f}'.format(chi2)
		average_Chi2[i]+=chi2
for i in range(len(outputworkspaces)): buf[i]+=' {0:7.3f}'.format(average_Chi2[i]/nQ)
buf='#Hq Chi2(Q=0.3) Chi2(Q=0.5) Chi2(Q=0.7) Chi2(Q=0.9)\n'+'\n'.join(buf)+'\n'
open('{0}/chi_versus_Hq_T{1}_int.dat'.format(rootd,T),'w').write(buf)
