import os
import numpy
import argparse
import mantid.simpleapi as mtds
#import mantid as mtd
from pdb import set_trace as tr

# Create resolution function
def resolfunc(instrument,npsec,template):
    FWHM={'HFBS':0.8E-03, 'BASIS':3.4e-03} #in meV
    s=FWHM[instrument]/2.355 # standard deviation
    X=template.dataX(0)
    Y= 1/(s*numpy.sqrt(2*numpy.pi)) * numpy.exp(-X*X/(2*s*s))
    X=numpy.tile(X, nspec)
    Y=numpy.tile(Y, nspec)
    ws_res_hfbs='resolution_'+instrument
    ws=mtds.CreateWorkspace(DataX=X,
        DataY=Y,
        Nspec=nspec,
        UnitX='EnergyTransfer',
        OutputWorkspace=ws_res_hfbs)
    return ws

def lose_tails(ws_sas, keep):
    ws_fqt=ws_sas+'_fqt.Re'
    tmax = dt/2 + dt*int(keep * max(mtds.mtd[ws_fqt].dataX(0))/dt)
    mtds.Rebin(InputWorkspace=ws_fqt,
               Params=[-tmax, dt, tmax],
               OutputWorkspace=ws_fqt
               )
    return ws_fqt

parser = argparse.ArgumentParser(description="""Generate S(Q,E) from simulated I(Q,t)
and convolve with model resolution function for HFBS instrument.\n
Example: genSQE.py 300 told3\n
Creates directory mantid_sqe/ and stores nexus files for HFBS or BASIS
resolution function, simulated S(Q,E), and S(Q,E) convolved with the resolution.""")
parser.add_argument('temperature', type=str, help='temperature')
parser.add_argument('deutsche', type=str, help='deuteration scheme')
parser.add_argument('--instrument', type=str, default='HFBS', help='instrument. Either HFBS or BASIS. Default is BASIS')
parser.add_argument('--outSQE', type=str, default='yes', help='Do we save S(Q,E) before convolution? Default=yes')
parser.add_argument('--outSQEconv', type=str, default='yes', help='Do we save S(Q,E) after convolution? Default=yes')
parser.add_argument('--outRes', type=str, default='yes', help='Do we save model resolution function for the instrument? Default=yes')
pargs=parser.parse_args()
deutscheme=pargs.deutsche
T=pargs.temperature
instr=pargs.instrument

# Load I(Q,t)
dt=1.0 # 1 picosecond
ws_sas=deutscheme
mtds.LoadSassena(Filename='fqt_inc_{0}.h5'.format(deutscheme),
                 OutputWorkspace=ws_sas,
                 TimeUnit=dt,
                 SortByQVectors=1
                 )
    

# Rebin losing tails with poor statistics
ws_fqt=lose_tails(ws_sas, 0.9)

# Fourier transform to get S(Q,E)
mtds.SassenaFFT(InputWorkspace=ws_sas,
                FFTonlyRealPart=1,
                DetailedBalance=1,
                Temp=float(T)
                )
ws_sqe=ws_sas+'_sqw'
wsqe=mtds.mtd[ws_sqe]

#Rescale so that sum over energies is I(Q,t=0)
I0=mtds.mtd[ws_sas+'_fq0'].dataY(0)[0]
sums=mtds.Integration(wsqe)
wsqe=mtds.Scale(wsqe, Factor=I0/sums.dataY(0)[0], Operation='Multiply')

# Convolve with resolution function
nspec=wsqe.getNumberHistograms()
wreshfbs=resolfunc(instr,nspec,wsqe) # Create HFBS resolution function
ws_hfbs=ws_sqe+'_'+instr
whfbs=mtds.ConvolveWorkspaces(wsqe, wreshfbs, OutputWorkspace=ws_hfbs)

#Save to files
savedir=os.getcwd()+'/mantid_sqe'
os.system('mkdir -p {0}'.format(savedir))

files_saved=list()
if pargs.outSQE.lower()[0]=='y':
    fname_sqe='{0}.nxs'.format(ws_sqe)
    mtds.SaveNexus(wsqe, savedir+'/'+fname_sqe)  #not convolved
    files_saved.append(fname_sqe)

if pargs.outSQEconv.lower()[0]=='y':
    fname_convolved='{0}.nxs'.format(ws_hfbs)
    mtds.SaveNexus(whfbs,savedir+'/'+fname_convolved)  #convolved
    files_saved.append(fname_convolved)

if pargs.outRes.lower()[0]=='y':
    fname_resolution='{0}.nxs'.format(wreshfbs.name())
    mtds.SaveNexus(wreshfbs, savedir+'/'+fname_resolution) #resolution
    files_saved.append(fname_resolution)

print 'Files saved:', ', '.join(files_saved)



