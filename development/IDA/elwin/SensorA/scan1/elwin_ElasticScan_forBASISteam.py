###################################################################
# This script to be run in the Script Window of Mantid will update
# reduced files bss41137_silicon111_sqw to bss41277_silicon111_sqw.
# These files represent an elastic scan from 20K to 300K.
#
# Once updated, go to directory /tmp to load the files in the ELWIN
# feature of the Indirect Data Analsys interface. Go to tutorial at
# http://www.mantidproject.org/IDA:Elwin#Tutorial to see how to do this.
####################################################################
target_dir='/projects/development/IDA/elwin/SensorA/scan1'
FirstRunNumber=41137
LastRunNumber=41277
for runNumber in range(FirstRunNumber, 1+LastRunNumber):
	name='bss{0:05d}_silicon111_sqw'.format(runNumber)
	LoadNexus(Filename='/SNS/BSS/IPTS-11889/shared/autoreduce/{0}.nxs'.format(name), OutputWorkspace=name)
	LoadInstrument(Workspace=name, InstrumentName='BASIS')
	ClearMaskFlag(Workspace=name)
	SaveNexus(InputWorkspace=name, FileName='/{0}/{1}.nxs'.format(target_dir,name))
	DeleteWorkspace(Workspace=name)
