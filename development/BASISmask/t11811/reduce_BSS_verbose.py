import os
import sys
import os

rootd='/projects/development/BASISmask/t11811/test'

# Input files
nexus_file='BSS_47884_event.nxs'
mask_file='BASIS_Mask.xml'
parameters_file='BASIS_silicon_111_Parameters.xml'

# Output files
run_number = nexus_file.split('_')[1]
dave_grp_filename = "BASIS_" + run_number + "_1run.dat"
processed_filename = "BSS_" + run_number + "_silicon111sqw.nxs"

def reduce(filename, suffix='', idf=None):
    autows = "auto_ws"+suffix      
    autows_monitor = autows + "_monitor"
    if idf:
        os.system('/bin/cp /tmp/BASIS_Definition_20140101-.xml /home/jbq/repositories/mantidproject/mantid/Code/Mantid/instrument/')
    Load(Filename=os.path.join(rootd,filename), OutputWorkspace=autows)
    #if idf:
    #    os.system('/bin/rm /home/jbq/repositories/mantidproject/mantid/Code/Mantid/instrument/BASIS_Definition_20140101-.xml')
    #return None
    
    data=mtd[autows].extractY()[0:2520*4]
    LoadMask(Instrument='BASIS', OutputWorkspace='BASIS_MASK', InputFile=os.path.join(rootd,mask_file))
    MaskDetectors(Workspace=autows, MaskedWorkspace='BASIS_MASK')
    ModeratorTzeroLinear(InputWorkspace=autows,OutputWorkspace=autows)
    LoadParameterFile(Workspace=autows, Filename=os.path.join(parameters_file))
    LoadNexusMonitors(Filename=os.path.join(rootd,filename), OutputWorkspace=autows_monitor)
    MonTemp=CloneWorkspace(autows_monitor)
    ModeratorTzeroLinear(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor)
    ConvertUnits(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, Target='Wavelength')
    OneMinusExponentialCor(InputWorkspace=autows_monitor, 
                        OutputWorkspace=autows_monitor, C='0.20749999999999999', C1='0.001276')
    Scale(InputWorkspace=autows_monitor, OutputWorkspace=autows_monitor, Factor='9.9999999999999995e-07')
    CloneWorkspace(InputWorkspace=autows, OutputWorkspace=autows+'.bak')
    ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='Wavelength', EMode='Indirect')
    
    #if idf:
    #    os.system('/bin/rm /home/jbq/repositories/mantidproject/mantid/Code/Mantid/instrument/BASIS_Definition_20140101-.xml')
    #return None

    
    RebinToWorkspace(WorkspaceToRebin=autows, WorkspaceToMatch=autows_monitor, OutputWorkspace=autows)
    Divide(LHSWorkspace=autows, RHSWorkspace=autows_monitor,  OutputWorkspace=autows)
    ConvertUnits(InputWorkspace=autows, OutputWorkspace=autows, Target='DeltaE', EMode='Indirect')
    CorrectKiKf(InputWorkspace=autows, OutputWorkspace=autows,EMode='Indirect')
    Rebin(InputWorkspace=autows, OutputWorkspace=autows, Params='-0.12,0.0004,0.12')
    QAxisBinning='0.2,0.2,2.0'
    SofQW3(InputWorkspace=autows, OutputWorkspace=autows+'_sqw', QAxisBinning=QAxisBinning,
        EMode='Indirect', EFixed='2.082')
    ClearMaskFlag(Workspace=autows+'_sqw')
    #SaveDaveGrp(Filename=os.path.join(rootd,dave_grp_filename), InputWorkspace=autows+'_sqw', ToMicroEV=True)
    #SaveNexus(Filename=os.path.join(rootd,processed_filename), InputWorkspace=autows+'_sqw')

reduce(nexus_file)
reduce(nexus_file, suffix='new', idf='/tmp/BASIS_Definition_20140101-.xml')
#CheckWorkspacesMatch(Workspace1='auto_ws', Workspace2='auto_wsnew',
#    Tolerance=0.001, CheckType=1, CheckAxes=1, CheckSpectraMap=0, CheckMasking=1, CheckAllData=0)