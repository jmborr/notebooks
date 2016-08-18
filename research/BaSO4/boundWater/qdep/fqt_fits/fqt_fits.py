from os.path import join as pjoin

saveDir="/SNSlocal/projects/jbq/BaSO4/boundWater/qdep/fqt_fits"
# Load experimental data containing Q-values 0.3, 0.5,..,1.9 (nine values)
expdir = "/SNSlocal/projects/jbq/BaSO4/expdata/qdep"
LoadNexus(Filename='%s/elastic.nxs'%expdir, OutputWorkspace='resolution')
LoadNexus(Filename='%s/exp300.nxs'%(expdir), OutputWorkspace='data')
e_range=(-0.1,0.0004,0.5)

# Load simulated data containing Q-values 0.3, 0.5,..,1.1  (five values)
# Load total and components (notice the time step is 0.5ps here)
simdir="/SNSlocal/projects/jbq/BaSO4/boundWater/qdep"
for component in ("", "YY", "NY", "YN", "NN"):
    LoadSassena(Filename='%s/fqt%s_inc.h5'%(simdir,component),OutputWorkspace='sim%s'%component,TimeUnit=0.5) # load simulated data
    Rebin(InputWorkspace='sim%s_fqt.Re'%component, OutputWorkspace='sim%s_fqt.Re'%component, Params=[-5000,0.5,5000])  # Eliminate wiggly tails
    SassenaFFT(InputWorkspace='sim%s'%component, FFTonlyRealPart=1, DetailedBalance=1, Temp=300)
    # More manageable for the fitting procedure
    Scale(InputWorkspace='sim%s_sqw'%component,Factor=1.0e-09, Operation='Multiply',OutputWorkspace='sim%s_sqw'%component)
    Rebin(InputWorkspace='sim%s_sqw'%component, Params=[-1,0.0004,1], OutputWorkspace='sim%s_sqw'%component)

for wi in range(5):
    ######################
    ## Fit sim_sqw to data
    ######################    
    fitstr=\
    """(composite=Convolution,FixResolution=false,NumDeriv=true;
        name=TabulatedFunction,Workspace=resolution,WorkspaceIndex={0},
        Scaling=1,Shift=0,XScaling=1,ties=(Scaling=1,XScaling=1);
        (name=DeltaFunction,Height=1,Centre=0,ties=(Centre=0);
         name=TabulatedFunction,Workspace=sim_sqw,WorkspaceIndex={0}
        ,Scaling=1,Shift=0,XScaling=1,ties=(Shift=0,XScaling=1)
        )
       );
       name=LinearBackground,A0=0,A1=0""".format(wi)
    Fit(fitstr, InputWorkspace='data', WorkspaceIndex=wi,
        CreateOutput=1, Output="fit{0}".format(wi),
        OutputCompositeMembers=1, ConvolveMembers=1,
        startX=e_range[0], endX=e_range[-1])
    ######################
    ## Adjust components to the fit
    ######################    
    ## Retrieve parameters from the fit
    ws = mtd["fit{0}_Parameters".format(wi)]
    shift = ws.row(1)['Value']
    height = ws.row(3)['Value']
    scaling = ws.row(5)['Value']
    a0 = ws.row(8)['Value']
    a1 = ws.row(9)['Value']
    # For each component, evaluate the function (Fit with zero iterations)
    # with the optimized parameters of the previous fit
    single_elastic = ExtractSingleSpectrum(InputWorkspace="resolution", WorkspaceIndex=wi)
    for component in ("YY", "NY", "YN", "NN"):
        fitstr=\
        """(composite=Convolution,FixResolution=false,NumDeriv=true;
            name=TabulatedFunction,Workspace=resolution,WorkspaceIndex={0},
            Scaling=1,Shift={1},XScaling=1,ties=(Scaling=1,XScaling=1);
            (name=DeltaFunction,Height={2},Centre=0,ties=(Centre=0);
             name=TabulatedFunction,Workspace={3},WorkspaceIndex={0}
             ,Scaling={4},Shift=0,XScaling=1,ties=(Shift=0,XScaling=1)
            )
           );
           name=LinearBackground,A0={5},A1={6}""".format(wi,shift,height,"sim%s_sqw"%component,scaling,a0,a1)
        Fit(fitstr, InputWorkspace='DATA', WorkspaceIndex=0,
            Output="fit{0}{1}".format(wi,component),
            CreateOutput=1, OutputCompositeMembers=1, ConvolveMembers=1,
            startX=e_range[0], endX=e_range[-1], MaxIterations=0)

##########
## Prepare workspaces containing data for Andrew
## e.g. fitQ0.3 contains 8 spectra, in this order: (1) Data, (2) Model, (3) Resolution, (4) Background,
##  (5) BB component, (6) UU component, (7) BU component (8) UB component
##########
for index in range(5):
    fitAllWs = "fit{0}_Workspace".format(index)
    Q = 0.3 + 0.2*index
    fitAndrewWs = "fitQ{0:3.1f}".format(Q)
    ExtractSpectra(InputWorkspace=fitAllWs, OutputWorkspace=fitAndrewWs, WorkspaceIndexList='0,1,2,3,5')
    for component in ('YY', 'NN', 'YN', 'NY'):
        fitWs = "fit{0}{1}_Workspace".format(index,component)
        ExtractSingleSpectrum(InputWorkspace=fitWs, WorkspaceIndex=4, OutputWorkspace="junk")
        AppendSpectra(InputWorkspace1=fitAndrewWs, InputWorkspace2="junk", OutputWorkspace=fitAndrewWs)
    # Save to Ascii file
    SaveAscii(InputWorkspace=fitAndrewWs,
        Filename=pjoin(saveDir,fitAndrewWs+".dat"),
        CommentIndicator='#', Separator="Space")


