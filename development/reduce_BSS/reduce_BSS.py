import os
import sys
import shutil
import numpy

mantid_root = "/opt/Mantid"
mantid_bin = sys.path.append(os.path.join(mantid_root, "bin"))

sys.path.insert(0,"/mnt/software/lib/python2.6/site-packages/matplotlib-1.2.0-py2.6-linux-x86_64.egg/")
from matplotlib import *                                                                                      
use("agg")                                                                                                    
import matplotlib.pyplot as plt                                                                           
numpy.seterr(all='ignore')                                                                                    
import warnings
warnings.filterwarnings('ignore',module='numpy') 

from mantid.simpleapi import *

nexus_file=sys.argv[1]
output_directory=sys.argv[2]

filename = os.path.split(nexus_file)[-1]
run_number = filename.split('_')[1]

autows = "__auto_ws"
autows_monitor = autows + "_monitor"

dave_grp_filename = os.path.join(output_directory, "BASIS_" + run_number + "_1run.dat")
processed_filename = os.path.join(output_directory, "BSS_" + run_number + "_silicon111sqw.nxs")

Load(Filename=nexus_file, OutputWorkspace=autows)
data=mtd[autows].extractY()[0:2520*4]
LoadMask(Instrument='BASIS', OutputWorkspace='BASIS_MASK', InputFile='/SNS/BSS/shared/autoreduce/BASIS_Mask.xml')
MaskDetectors(Workspace=autows, MaskedWorkspace='BASIS_MASK')
ModeratorTzeroLinear(InputWorkspace=autows,OutputWorkspace=autows)
LoadParameterFile(Workspace=autows, Filename=os.path.join(mantid_root, 'instrument', 'BASIS_silicon_111_Parameters.xml'))
LoadNexusMonitors(Filename=nexus_file, OutputWorkspace=autows_monitor)
ModeratorTzeroLinear(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor)
Rebin(InputWorkspace=autows_monitor,OutputWorkspace=autows_monitor,Params='10')
ConvertUnits(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, Target='Wavelength')
OneMinusExponentialCor(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, C='0.20749999999999999', C1='0.001276')
Scale(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, Factor='9.9999999999999995e-07')
ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='Wavelength', EMode='Indirect')
RebinToWorkspace(WorkspaceToRebin=autows, WorkspaceToMatch=autows_monitor, OutputWorkspace=autows)
Divide(LHSWorkspace=autows, RHSWorkspace=autows_monitor,  OutputWorkspace=autows)
ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='DeltaE', EMode='Indirect')
CorrectKiKf(InputWorkspace=autows, OutputWorkspace=autows,EMode='Indirect')

#RenameWorkspace(InputWorkspace=autows,OutputWorkspace='bss20200_silicon111_red')
Rebin(InputWorkspace=autows, OutputWorkspace=autows, Params='-0.12,0.0004,0.12')
#GroupDetectors(InputWorkspace=autows, OutputWorkspace=autows, MapFile='/SNS/BSS/shared/autoreduce/BASIS_Grouping.xml', Behaviour='Sum')
QAxisBinning='0.2,0.2,2.0'
autows_sqw = SofQW3(InputWorkspace=autows, OutputWorkspace=autows+'_sqw', QAxisBinning=QAxisBinning, EMode='Indirect', EFixed='2.082')
SaveDaveGrp(Filename=dave_grp_filename, InputWorkspace=autows+'_sqw', ToMicroEV=True)
SaveNexus(Filename=processed_filename, InputWorkspace=autows+'_sqw')

# Save experiment log file
logname = os.path.join(output_directory, 'experiment_log.csv')
filemode = 'new'
if os.path.exists(logname):
    filemode = 'fastappend'
comment = mtd[autows].getComment()
AddSampleLog(autows, LogName='Comment', LogText=comment, LogType='String')
ExportExperimentLog(InputWorkspace=autows, OutputFilename=logname, FileMode=filemode,
                    SampleLogTitles = 'Run number,Title,Comment,StartTime,EndTime,Duration,ProtonCharge,Mean Sensor A, Min Sensor A, Max Sensor A,  Mean Sensor B, Min Sensor B, Max Sensor B, Wavelength, Chopper 1, Chopper 2, Chopper 3, Slit S1t, Slit S1b, Slit S1l, Slit S1r',
                    SampleLogNames = 'run_number,run_title,Comment,start_time,end_time,duration,gd_prtn_chrg,SensorA,SensorA,SensorA,SensorB,SensorB,SensorB,LambdaRequest,Speed1,Speed2,Speed3,s1t,s1b,s1l,s1r',
                    SampleLogOperation = '0,0,0,0,0,0,0,average,min,max,average,min,max,0,0,0,0,0,0,0,0',
                    FileFormat = 'comma (csv)',
)

# Make Figures
fig = plt.gcf()
fig.set_size_inches(8.0,16.0)

# Instrument Figure
ax1 = plt.subplot2grid((7,2), (0,0), colspan=2, rowspan=2)
bss=numpy.zeros((91,113))
try:
    b1=data[0:2520].reshape(56,45).transpose()[::-1,::-1]
    bss[0:45,0:56]=b1
    b2=data[2520:5040].reshape(56,45).transpose()[::-1,::-1]
    bss[46:91,0:56]=b2
    b3=data[5040:7560].reshape(56,45).transpose()[::-1,::-1]
    bss[0:45,57:113]=b3
    b4=data[7560:10080].reshape(56,45).transpose()[::-1,::-1]
    bss[46:91,57:113]=b4
except:
    pass
plt.imshow(numpy.log(bss))
plt.axis('off')

# Monitor Figure
ax1 = plt.subplot2grid((7,2), (2,0))
ax1.ticklabel_format(style = 'sci', axis='x',scilimits=(0,0))
mon = mtd[autows_monitor]
x = mon.readX(0)
y = mon.readY(0)                                                                                          
plt.plot(x[1:],y)
plt.xlabel('TOF ($\mu$s)')
plt.ylabel('Intensity')
plt.title('Monitor')
plt.axis('on')

# Spectra Figures
Qm,dQ,QM = [float(x) for x in QAxisBinning.split(',')]
nQ = int( (QM-Qm)/dQ )
for i in range(nQ):
    irow=(i+1)/2+2
    icol=(i+1)%2
    ax1 = plt.subplot2grid((7,2), (irow,icol) )
    ax1.ticklabel_format(style = 'sci', axis='x',scilimits=(0,0)) 
    x = autows_sqw.readX(i)
    y = autows_sqw.readY(i)                                                                                          
    plt.plot(x[1:],y)
    plt.xlabel('Energy ($\mu$eV)')
    plt.ylabel('Intensity')
    plt.yscale('log')
    plt.title('Q={0} '.format(Qm + dQ/2 + i*dQ)+"$\AA^{-1}$")

plt.tight_layout(1.08)
plt.show()
plt.savefig(processed_filename+'.png', bbox_inches='tight')
plt.close()
