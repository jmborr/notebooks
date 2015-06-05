import os
import sys

rootd='/projects/development/BASISmask/t11811/test'

# Input files
nexus_file='BSS_47884_event.nxs'
mask_file='BASIS_Mask.xml'
parameters_file='BASIS_silicon_111_Parameters.xml'

# Output files
run_number = nexus_file.split('_')[1]
dave_grp_filename = "BASIS_" + run_number + "_1run.dat"
processed_filename = "BSS_" + run_number + "_silicon111sqw.nxs"


autows = "__auto_ws"
autows_monitor = autows + "_monitor"
Load(Filename=os.path.join(rootd,nexus_file), OutputWorkspace=autows)
data=mtd[autows].extractY()[0:2520*4]
LoadMask(Instrument='BASIS', OutputWorkspace='BASIS_MASK', InputFile=os.path.join(rootd,mask_file))
MaskDetectors(Workspace=autows, MaskedWorkspace='BASIS_MASK')
ModeratorTzeroLinear(InputWorkspace=autows,OutputWorkspace=autows)
LoadParameterFile(Workspace=autows, Filename=os.path.join(parameters_file))
LoadNexusMonitors(Filename=os.path.join(rootd,nexus_file), OutputWorkspace=autows_monitor)
MonTemp=CloneWorkspace(autows_monitor)
ModeratorTzeroLinear(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor)
Rebin(InputWorkspace=autows_monitor,OutputWorkspace=autows_monitor,Params='10')
ConvertUnits(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, Target='Wavelength')
OneMinusExponentialCor(InputWorkspace=autows_monitor, 
                       OutputWorkspace=autows_monitor, C='0.20749999999999999', C1='0.001276')
Scale(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, Factor='9.9999999999999995e-07')
ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='Wavelength', EMode='Indirect')
RebinToWorkspace(WorkspaceToRebin=autows, WorkspaceToMatch=autows_monitor, OutputWorkspace=autows)
Divide(LHSWorkspace=autows, RHSWorkspace=autows_monitor,  OutputWorkspace=autows)
ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='DeltaE', EMode='Indirect')
CorrectKiKf(InputWorkspace=autows, OutputWorkspace=autows,EMode='Indirect')

Rebin(InputWorkspace=autows, OutputWorkspace=autows, Params='-0.12,0.0004,0.12')
QAxisBinning='0.2,0.2,2.0'
SofQW3(InputWorkspace=autows, OutputWorkspace=autows+'_sqw', QAxisBinning=QAxisBinning,
       EMode='Indirect', EFixed='2.082')
ClearMaskFlag(Workspace=autows+'_sqw')
SaveDaveGrp(Filename=os.path.join(rootd,dave_grp_filename), InputWorkspace=autows+'_sqw', ToMicroEV=True)
SaveNexus(Filename=os.path.join(rootd,processed_filename), InputWorkspace=autows+'_sqw')

# compare the old and new workspaces
wold = LoadNexus(Filename=os.path.join(rootd,'old','BSS_47884_silicon111sqw.nxs'),OutputWorkspace='wold')
wnew = LoadNexus(Filename=os.path.join(rootd,'BSS_47884_silicon111sqw.nxs'),OutputWorkspace='wnew')
plotSpectrum(wold, range(wold.getNumberHistograms()))
plotSpectrum(wnew, range(wold.getNumberHistograms()))