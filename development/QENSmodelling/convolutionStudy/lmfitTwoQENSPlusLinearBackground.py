'''Use of lmfit + scipy + numpy
Convolution is "helped"
'''
import mantid.simpleapi as ms
import numpy as np
import lmfit as lf
import sys
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pdb import set_trace as tr

data_dir="/projects/development/QENSmodelling/convolutionStudy"  #substitute with your own directory

def createTableWorkspace(result, workspaceName):
    '''Create Mantid table workspace with results from the fit'''            
    table=CreateEmptyTableWorkspace(OutputWorkspace=workspaceName)
    table.addColumn('str', 'Name')
    table.addColumn('float','Value')
    table.addColumn('float','Error')
    for parname in result.params:
        p = result.params[parname]
        table.addRow((p.name, p.value, p.stderr))
    table.addRow(('chisqr', result.redchi, 0.0))
    return table

def quasiElasticSignal(x, A, B, f1, f2):
    '''Two Lorentzians
         /     f1           f2   \
(1.0/pi)*|A*-------- + B*---------|
         |    2    2       2    2 |
         \  f1  + x      f2  + x  /
    '''
    y = (1./np.pi)*(A*f1/(x*x+f1*f1)+B*f2/(x*x+f2*f2))
    return y

def linearBackground(x, a, b):
    y = a*x+b
    return y
   
def convolve(x, resolution, qeCallback, parameters):
    '''Convolve usin numpy.convolution
    Extend the domain to avoid boundary effects with the convolution
    '''
    A= parameters['F1'].value
    B = parameters['F2'].value
    f1 = parameters['f1'].value
    f2 = parameters['f2'].value
    left = x-x[-1]  # negative energies domain
    right = x-x[0]  # positive energies domain
    xx = np.concatenate([left, right[1:]]) #Avoid common value between left and right
    qe = qeCallback(xx, A, B, f1, f2) # evaluate on the extended domain
    y = np.convolve(resolution, qe, mode='valid')
    #qe = qeCallback(x, A, B, f1, f2)
    #y = np.convolve(resolution, qe, mode='same')
    #np.savetxt("quasielastic.dat", zip(xx,qe))
    #np.savetxt("convolution.dat", zip(x,y/2184.6677557517305))
    return y

def residuals(parameters, dataX, dataY, dataE, resY):
    residuals.counter += 1
    newResolution = resY(parameters['e'].value+dataX)
    model =  parameters['E'].value*newResolution  # Elastic signal
    model += convolve( dataX, newResolution, quasiElasticSignal, parameters ) # QuasiElastic signal
    model += linearBackground(dataX, parameters['a'].value, parameters['b'].value) # Background
    #print 'residuals called',residuals.counter
    #parameters.pretty_print()
    #np.savetxt("model.dat", zip(dataX,model))
    return (dataY-model)/dataE
residuals.counter = 0

params = lf.Parameters()
params.add('E', value=1, min=0) #elastic
params.add('F1', value=1, min=0) #first lorentzian
params.add('F2', value=1, min=0) #second lorentzian
params.add('e', value=0) # resolution shift
params.add('f1', value=0.1, min=0) # FWHM first lorentzian
params.add('f2', value=0.2, min=0) # FWHM second lorentzian
params.add('a', value=0) # slope of the background
params.add('b', value=0) # intercept of the background

skewToFlag = {'sym':'', 'asym':'Asym-'} #differentiation used in the file names
skewToseeds = {'sym' : ('0.0-7220', '0.0-7240', '0.1-7220', '0.1-7240', '0.5-7220', '0.5-7240', '1.0-7220', '1.0-7240'),
               'asym': ('0.0-7220', '0.1-7220', '0.5-7220', '1.0-7220')}
for skew in ('sym', 'asym'):
    for seed in skewToseeds[skew]:
        ## Load resolution and signal files (reduced files)
        resolution_file = os.path.join(data_dir, skew, 'QENS_Synthetic_Resolution-{0}{1}.dat'.format(skewToFlag[skew],seed))
        resolution = ms.LoadAscii(Filename=resolution_file, Unit='DeltaE', Separator='Space')
        signal_file = os.path.join(data_dir, skew, 'QENS_Synthetic_Signal-{0}{1}.dat'.format(skewToFlag[skew],seed))
        signal = ms.LoadAscii(Filename=signal_file, Unit='DeltaE', Separator='Space')
        dataX = signal.dataX(0)
        dataY = signal.dataY(0)
        dataE = signal.dataE(0)
        ## reset parameters
        for (name,value) in {'E':1, 'F1':1, 'F2':1, 'e':0, 'f1':0.1, 'f2':0.2, 'a':0, 'b':0}.items():
            params[name].value=value
        '''Generate the interpolator for the resolution function to treat shifts, resY'''
        rX = resolution.dataX(0)
        rY = resolution.dataY(0)
        npad=9 #padding
        dx0 = (rX[-1]-rX[0])/len(rX)
        leftX = rX[0] - dx0*np.arange(npad,0,-1)
        rightX = rX[-1] + dx0*np.arange(1,1+npad)
        xx = np.concatenate( [leftX, rX, rightX])
        yy = np.concatenate( [rY[0]*np.ones(npad), rY, rY[-1]*np.ones(npad)] )
        resY = interp1d(xx, yy, kind='linear', bounds_error=False,fill_value=rY.min())
        pargs = (dataX, dataY, dataE, resY ) #positional parameters for residuals function
        minimizer = lf.Minimizer(residuals, params, fcn_args=pargs)
        result = minimizer.minimize()
        modelY = dataY + result.residual
        ms.CreateWorkspace(OutputWorkspace="fit-{0}{1}".format(skewToFlag[skew],seed),
                        Nspec=3,
                        DataX = np.tile(dataX,3),
                        DataY = np.concatenate([dataY,dataY+result.residual,result.residual]),
                        dataE = np.concatenate([dataE, np.zeros(len(dataE)),dataE]),
                        UnitX='DeltaE' )
        table = createTableWorkspace(result, "params-{0}{1}".format(skewToFlag[skew],seed))
        lf.report_fit(result.params)
        #save some ascii files
        np.savetxt("data.dat", zip(dataX,dataY))
        np.savetxt("model.dat", zip(dataX,modelY))
