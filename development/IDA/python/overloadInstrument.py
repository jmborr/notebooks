#
#   Meant to be used with MantidNightly in bac.sns.gov
#

# Overwrite instrument definition
rootd='/projects/development/IDA/elwin/scan1'
for i in range(41137, 41277):
	wsname='bss{0}_silicon111_sqw'.format(i)
	Load(Filename=rootd+'/{0}.nxs'.format(wsname), OutputWorkspace=wsname)
	LoadInstrument(wsname, InstrumentName='BASIS', RewriteSpectraMap=1)
	outfile=rootd+'/bss{0}_sqw.nxs'.format(i)
	SaveNexus(wsname, outfile)
	DeleteWorkspace(wsname)
# Overwrite instrument definition
rootd='/projects/development/IDA/elwin/scan2'
for i in range(40996, 41136):
	wsname='bss{0}_silicon111_sqw'.format(i)
	Load(Filename=rootd+'/{0}.nxs'.format(wsname), OutputWorkspace=wsname)
	LoadInstrument(wsname, InstrumentName='BASIS', RewriteSpectraMap=1)
	outfile=rootd+'/bss{0}_sqw.nxs'.format(i)
	SaveNexus(wsname, outfile)
	DeleteWorkspace(wsname)

# Add temperature series to the sample logs
rootd='/projects/development/IDA/elwin/scan1'
T=20.0
for i in range(41137, 41277):
	wsname='bss{0}_sqw'.format(i)
	Load(Filename=rootd+'/{0}.nxs'.format(wsname), OutputWorkspace=wsname)
	AddSampleLog(wsname,LogName='Temperature', LogText=str(T), LogType='Number Series') #only one value in the series
	SaveNexus(wsname, rootd+'/{0}.nxs'.format(wsname))
	DeleteWorkspace(wsname)
	T += 2.0
# Add temperature series value to the sample logs
rootd='/projects/development/IDA/elwin/scan2'
T=300.0
for i in range(40996, 41136):
	wsname='bss{0}_sqw'.format(i)
	Load(Filename=rootd+'/{0}.nxs'.format(wsname), OutputWorkspace=wsname)
	AddSampleLog(wsname,LogName='Temperature', LogText=str(T), LogType='Number Series') #only one value in the series
	SaveNexus(wsname, rootd+'/{0}.nxs'.format(wsname))
	DeleteWorkspace(wsname)
	T -= 2.0
	