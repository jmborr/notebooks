# Fit 
expdir='/SNSlocal/projects/jbq/BaSO4/expdata'
simdir='/SNSlocal/projects/jbq/BaSO4/andrew_simulation'
e_range=(-0.1,0.0004,0.5)
LoadNexus(Filename='%s/elastic.nxs'%expdir, OutputWorkspace='elastic')
LoadNexus(Filename='%s/elastic_extended.nxs'%expdir, OutputWorkspace='elastic_extended')
LoadDaveGrp(Filename='%s/q300.dat'%(expdir), OutputWorkspace='exp300', XAxisUnits='DeltaE', YAxisUnits='Empty', IsMicroEV=1)  # experiment

# Load different components (notice the time step is 0.5ps here)
simdir="/SNSlocal/projects/jbq/BaSO4/boundWater/"
for component in ("", "YY", "NY", "YN", "NN"):
    LoadSassena(Filename='%s/fqt%s_inc.h5'%(simdir,component),OutputWorkspace='sim%s'%component,TimeUnit=0.5) # load simulated data
    SumSpectra(InputWorkspace='sim%s_fqt.Re'%component, OutputWorkspace='sim%s_fqt.Re'%component) # sum all Q-components
    Rebin(InputWorkspace='sim%s_fqt.Re'%component, OutputWorkspace='sim%s_fqt.Re'%component, Params=[-5000,0.5,5000])  # Eliminate wiggly tails
    SassenaFFT(InputWorkspace='sim%s'%component, FFTonlyRealPart=1, DetailedBalance=1, Temp=300)
    Scale(InputWorkspace='sim%s_sqw'%component,Factor=1.0e-08, Operation='Multiply',OutputWorkspace='sim%s_sqw'%component) # More manageable for the fitting procedure

# Fit sim_sqw to data
fitstr = 'name=TabulatedFunction,Workspace=elastic,WorkspaceIndex=0,Scaling=0.278159,Shift=6.35095e-05,ties=(XScaling=1);' +\
    '(composite=Convolution,FixResolution=true,NumDeriv=true;' +\
    'name=TabulatedFunction,Workspace=elastic_extended,WorkspaceIndex=0,Scaling=1;' +\
    'name=TabulatedFunction,Workspace=sim_sqw,WorkspaceIndex=0,Scaling=29.6129,Shift=-5.92427e-05,ties=(XScaling=1));' +\
    'name=LinearBackground,A0=0.00204131,A1=-0.00296182'
Fit(fitstr, InputWorkspace='exp300', WorkspaceIndex=0,
    CreateOutput=1, Output="sim_fit", OutputCompositeMembers=1, startX=e_range[0], endX=e_range[-1])

# For each component, evaluate the function (Fit with zero iterations) with the optimized parameters of the previous fit
fitstr_template = 'name=TabulatedFunction,Workspace=elastic,WorkspaceIndex=0,Scaling=0.296008,Shift=6.05556e-05,ties=(XScaling=1);' +\
    '(composite=Convolution,FixResolution=true,NumDeriv=true;' +\
    'name=TabulatedFunction,Workspace=elastic_extended,WorkspaceIndex=0,Scaling=1;' +\
    'name=TabulatedFunction,Workspace=_COMPONENT_,WorkspaceIndex=0,Scaling=3.27381,Shift=-4.18335e-05,ties=(XScaling=1));' +\
    'name=LinearBackground,A0=0.000179823,A1=0.000738576'
for component in ("YY", "NY", "YN", "NN"):
    sqwName = 'sim%s_sqw'%component
    fitstr = fitstr_template.replace("_COMPONENT_", sqwName)
    Fit(fitstr, InputWorkspace='exp300', WorkspaceIndex=0,
        Output='sim%s_fit'%component, CreateOutput=1, MaxIterations=0, startX=e_range[0], endX=e_range[-1])







