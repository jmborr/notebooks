#
#   Meant to be used with MantidNightly in bac.sns.gov
#
Load(Filename='/SNS/BSS/IPTS-10679/shared/autoreduce/bss35181_silicon111_sqw.nxs', OutputWorkspace='bss35181_silicon111_sqw')
LoadInstrument('bss35181_silicon111_sqw',InstrumentName='BASIS',RewriteSpectraMap=1)
SaveNexus('bss35181_silicon111_sqw', '/tmp/bss35181_silicon111_sqw.nxs')
