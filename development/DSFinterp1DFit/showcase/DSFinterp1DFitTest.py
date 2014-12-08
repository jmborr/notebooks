''' Comparison between structure factors derived from quasielastic data and molecular dynamics (MD) simulations.
- System: a octa-methyl silsesquioxane molecule ( http://goo.gl/LaCjAk )
- Fitting parameter of interest: barrier governing the rotation of the methyl groups around their three-fold axis.
- From the quasielastic data, the barrier has been estimated at K=0.068Kcal/mol
- Each MD simulation is carried out with a different value of the barrier, from K=0.02Kcal/mol up to K=0.15Kcal/mol
The DSFinterp1DFit function will load the estructure factors S(Q,E,K_i) derived from the simulations, and construct an S(Q,E,K) to fit against the
experimental data, with K one more fitting parameter.
The fit should return a value of K~0.064Kcal/mol and Chi-square~6.9
'''

testdir = '/tmp/DSFinterp1DFitTest'  # PLEASE UPDATE THIS VARIABLE IF NECCESSARY
testdir = '/projects/development/DSFinterp1DFit/showcase'
from mantid.simpleapi import LoadNexus, LoadSassena, Transpose, Rebin, Scale, SassenaFFT, ConvolveWorkspaces, Fit, mtd

LoadNexus(Filename=testdir+'/elastic.nxs', OutputWorkspace='elastic')  # resolution file, which is also the elastic line
LoadNexus(Filename=testdir+'/exp200K.nxs', OutputWorkspace='exp200K')  # Quasielastic data (taken at T=200K)

# Load S(Q,E,K_i) dervied from simulations (at  T=200K) for different values of the dihedral barrier K
workspaces=''
parametervalues='0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15'
for K in parametervalues.split():
	LoadNexus(Filename=testdir+'/simSMK{0}.nxs'.format(K), OutputWorkspace='simSMK{0}'.format(K))
	workspaces += ' simSMK{0}'.format(K)
	
guess=0.11  # initial value of the barrier
#Fittin model: S(Q,E) = A*elastic_line + B*Simulated_structure_factor + Linear_background
fit_string = 'name=TabulatedFunction,Workspace=elastic,WorkspaceIndex=8,Scaling=1,constraints=(0.0001<Scaling);'+\
	'name=DSFinterp1DFit,InputWorkspaces="{0}",ParameterValues="{1}",'.format(workspaces,parametervalues) +\
	'LoadErrors=0,LocalRegression=1,RegressionType=quadratic,RegressionWindow=6,' +\
	'WorkspaceIndex=8,Intensity=1.0,TargetParameter={0},'.format(guess) +\
	'constraints=(0.0001<Intensity);' +\
	'name=LinearBackground,A0=0.0,A1=0.0'

# Fit for WorkspaceIndex=8 corresponding to highest measured Q  (Q=19nm^(-1) or L~3Anstromgs)
Fit(fit_string, InputWorkspace='exp200K', WorkspaceIndex=8, StartX=-0.13, EndX=0.1, CreateOutput = 1 )
# print the optimized K and the associated Chi-square of the fit
print 'Kopt=', mtd['exp200K_Parameters'].row(2)['Value'], 'Chi2=', mtd['exp200K_Parameters'].row(5)['Value']