root='/SNSlocal/projects/jbq/nanodiamond/nodiamond'

# Load
LoadNexus(Filename='{0}/expdata/SQE_RNA_plus_Water.nxs'.format(root), OutputWorkspace='exp')
LoadNexus(Filename='{0}/expdata/Resolution_RNA_plus_Water.nxs'.format(root), OutputWorkspace='res')
ScaleX(InputWorkspace='exp', OutputWorkspace='exp', Factor=0.001, Operation='Multiply')
ScaleX(InputWorkspace='res', OutputWorkspace='res', Factor=0.001, Operation='Multiply')
Scale(InputWorkspace='exp', OutputWorkspace='exp', Factor=1e11, Operation='Multiply')
Scale(InputWorkspace='res', OutputWorkspace='res', Factor=1e11, Operation='Multiply')

LoadSassena(Filename='{0}/equil/fqt_inc_RNA.h5'.format(root), TimeUnit=1.0, SortByQVectors=1, OutputWorkspace='fqt')
Rebin(InputWorkspace='fqt_fqt.Re', OutputWorkspace='fqt_fqt.Re', Params=[-5000,1,5000])
SassenaFFT(InputWorkspace='fqt',Temp=300, DetailedBalance=1, FFTonlyRealpart=1)
ws=mtd['fqt_sqw']
for wi in range(ws.getNumberHistograms()):
    avgY=(ws.dataY(wi)[4998] + ws.dataY(wi)[5002]) / 2
    for ix in (4999, 5000, 5001):
        ws.dataY(wi)[ix] = avgY
Scale(InputWorkspace='fqt_sqw', OutputWorkspace='fqt_sqw', Factor=1/avgY, Operation='Multiply')

LoadSassena(Filename='{0}/equil/fqt_inc_HeavyAtoms_RNA.h5'.format(root), TimeUnit=1.0, SortByQVectors=1, OutputWorkspace='fqtHA')
Rebin(InputWorkspace='fqtHA_fqt.Re', OutputWorkspace='fqtHA_fqt.Re', Params=[-5000,1,5000])
SassenaFFT(InputWorkspace='fqtHA',Temp=300, DetailedBalance=1, FFTonlyRealpart=1)
ws=mtd['fqtHA_sqw']
for wi in range(ws.getNumberHistograms()):
    avgY=(ws.dataY(wi)[4998] + ws.dataY(wi)[5002]) / 2
    for ix in (4999, 5000, 5001):
        ws.dataY(wi)[ix] = avgY
Scale(InputWorkspace='fqtHA_sqw', OutputWorkspace='fqtHA_sqw', Factor=1/avgY, Operation='Multiply')

# Iterative fitting
ws=mtd['exp']
f0_f0_Shift=-2.709844e-06
f0_f1_f0_Height=0.481817
f0_f1_f1_Scaling=75.954263
for wi in range(ws.getNumberHistograms()-1, -1, -1):
    fitString='(composite=Convolution,FixResolution=false,NumDeriv=true;'+\
        'name=TabulatedFunction,Workspace=res,WorkspaceIndex={0},Scaling=1,Shift={1},ties=(Scaling=1);'.format(wi, f0_f0_Shift)+\
        '(name=DeltaFunction,Height={0};'.format(f0_f1_f0_Height)+\
        'name=TabulatedFunction,Workspace=fqt_sqw,WorkspaceIndex={0},Scaling={1},Shift=0.0,ties=(Shift=0.0)));'.format(wi+1, f0_f1_f1_Scaling)+\
        'name=LinearBackground,A0=0.0,A1=0.0'
    outputname='fit_{0}'.format(wi)
    Fit(fitString, InputWorkspace='exp', WorkspaceIndex=wi, startX=-0.2, endX=0.2, CreateOuptut=1, Output=outputname)
    wfit = mtd[outputname+'_Parameters']
    f0_f0_Shift = wfit.row(1)['Value']
    f0_f1_f0_Height = wfit.row(2)['Value']
    f0_f1_f1_Scaling = wfit.row(3)['Value']
    chi2 = wfit.row(7)['Value']
    Q= 0.5 + wi * 0.2
    print '{0} {1} {2}'.format(wi, Q, chi2)
