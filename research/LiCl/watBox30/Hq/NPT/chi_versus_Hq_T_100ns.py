from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, DSFinterp, Fit
import os

T=210
T2T={200:'201', 210:'210'} #There isn't an equivalence between experimental and simulated temperatures for T=201K
rootd='/projects/research/LiCl/EugeneData/NSEdata'
LoadAscii(Filename=rootd+'/NSE_{0}K.dat'.format(T2T[T]), OutputWorkspace='NSE_{0}K'.format(T2T[T]),Separator='Tab',CommentIndicator='#',Unit='Time')
ScaleX(InputWorkspace='NSE_{0}K'.format(T2T[T]), OutputWorkspace='NSE_{0}K'.format(T2T[T]), Factor=1000, Operation='Multiply' ) #from ns to ps
 #These are linear interpolations of Eugene's data so that we have one point every 10ps, as in the simulations
datafile={200:'NSE_201K_linearInterp.dat', 210:'NSE_210K_linearInterp.dat'}
maxtimes={200:26200, 210:8100}  #maximum time (ps) for each experimental dataset
LoadAscii(Filename=rootd+'/'+datafile[T], OutputWorkspace='LiCl_{0}K'.format(T),Separator='Space',CommentIndicator='#',Unit='Time')

rootd='/projects/research/LiCl/watBox30/Hq/NPT'

# Not all 100ns simulations were successfully carried out, thus not all parametervalues have their corresponding structure factor. Hence the prefix 'nominal'
nominal_parametervalues='0.400 0.402 0.404 0.406 0.408 0.410 0.412 0.414 0.416 0.420 0.422 0.424 0.426 0.428 0.430 0.432 0.440 0.450 0.460 0.470 0.480 0.490 0.500'.split()
nominal_qdirs='Q400 Q402 Q404 Q406 Q408 Q410 Q412 Q414 Q416 Q420 Q422 Q424 Q426 Q428 Q430 Q432 Q440 Q450 Q460 Q470 Q480 Q490 Q500'.split()
#nominal_parametervalues='0.400 0.410 0.420 0.430 0.440 0.450 0.460 0.470 0.480 0.490 0.500'.split()
#nominal_qdirs='Q400 Q410 Q420 Q430 Q440 Q450 Q460 Q470 Q480 Q490 Q500'.split()
parametervalues=''
qdirs=''
for i in range( len(nominal_qdirs) ):
	qdir=nominal_qdirs[i]
	sassenafile=rootd+'/{0}/T{1}/production100ns/fqt_rms2first_inc.hd5'.format(qdir,T)
	if os.path.exists(sassenafile):
		LoadSassena(Filename=sassenafile, TimeUnit=10.0, SortByQVectors=1, OutputWorkspace=qdir)
		if mtd['{0}_fqt.Re'.format(qdir)].dataX(0)[-1] > maxtimes[T]:
			Rebin(InputWorkspace='{0}_fqt.Re'.format(qdir), Params=[-5,10.0,maxtimes[T]+5], OutputWorkspace='{0}_iqt'.format(qdir))
			factor=1.0/mtd['{0}_fq0'.format(qdir)].dataY(0)[0] # Normalize I(Q,t=0)=1
			Scale(InputWorkspace='{0}_iqt'.format(qdir), OutputWorkspace='{0}_iqt'.format(qdir), Factor=factor, Operation='Multiply' )
			parametervalues += ' '+nominal_parametervalues[i]
			qdirs += ' '+qdir

#####################
#### DSFinterp algorithm  ####
#####################
# Input
parmlist=[float(x) for x in parametervalues.split()]
workspaces=[]
outputworkspaces=[]
for qdir in qdirs.split():
	workspaces.append('{0}_iqt'.format(qdir))
	outputworkspaces.append('int{0}'.format(qdir))

#####################
#### DSFinterp1DFit Fit     ####
#####################
parmlist_str=' '.join([ str(x) for x in parmlist])
fit_string   ='name=DSFinterp1DFit,InputWorkspaces="{0}",ParameterValues="{1}",'.format(' '.join(workspaces), parmlist_str)
#fit_string+='LoadErrors=0,LocalRegression=0,RegressionType=linear,RegressionWindow=4,'  #  FITTING AT 200K
fit_string+='LoadErrors=0,LocalRegression=1,RegressionType=linear,RegressionWindow=4,'
fit_string+='WorkspaceIndex=0,Intensity=1.0,TargetParameter=0.424;'
fit_string+='name=FlatBackground,A0=0.1,constraints=(0<A0)'
Fit(fit_string, InputWorkspace='NSE_{0}K'.format(T2T[T]), WorkspaceIndex=0, StartX=0, EndX=maxtimes[T],CreateOutput=1) 


#####################
#### DSFinterp algorithm  ####
#####################
# Input
outparmlist=[]
outputworkspaces=[]
x=410
while x < 440:
	outparmlist.append(x*0.001)
	outputworkspaces.append('intQ{0}'.format(x))
	x+=1
N=len(outparmlist)
#Run the algorithm
DSFinterp(Workspaces=workspaces, ParameterValues=parmlist, OutputWorkspaces=outputworkspaces, LoadErrors=0, LocalRegression=1,RegressionType='linear', RegressionWindow=4, TargetParameters=outparmlist)

nQ=1
# Calculate fits for the interpolated structure factors
buf=['{0:5.3f}'.format(x) for x in outparmlist]
average_Chi2=[0]*len(outparmlist)
for wi in range(nQ):
	for i in range(N):
		Hq=outparmlist[i]
		workspace=outputworkspaces[i]
		fit_string    ='name=TabulatedFunction,Workspace={0},WorkspaceIndex={1},Scaling=1.0;'.format(workspace,wi)
		fit_string +='name=FlatBackground,A0=0.1,constraints=(0<A0)'
		#Fit(fit_string, InputWorkspace='NSE_{0}K'.format(T2T[T]), WorkspaceIndex=wi, StartX=0, EndX=maxtimes[T], CreateOutput=1)
		Fit(fit_string, InputWorkspace='NSE_{0}K'.format(T2T[T]), WorkspaceIndex=wi, StartX=21, EndX=8100, CreateOutput=1)
		# print the optimized K and the associated Chi-square of the fit
		ws=mtd['NSE_{0}K_Parameters'.format(T2T[T])]
		chi2=ws.row(2)['Value']
		print 'Q=',0.45, 'Hq=',  Hq, 'Chi2=', chi2, 'I=',ws.row(0)['Value'], 'A0=',ws.row(1)['Value']
		buf[i]+=' {0:7.3f} {1:5.3} {2:5.3}'.format(chi2,ws.row(0)['Value'], ws.row(1)['Value'])
		average_Chi2[i]+=chi2
buf='#Hq     Chi2    I    A0\n'+'\n'.join(buf)+'\n'
open('{0}/chi_versus_Hq_T{1}_100ns_int.dat'.format(rootd,T),'w').write(buf)


