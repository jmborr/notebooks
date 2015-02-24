#Convert ASCII experimental data to Nexus
rootd='/SNSlocal/projects/jbq/nanodiamond/expdata'
for temperature in '220 250 260 280 290 300 310'.split():
	filename='{0}/ND+_RNA_+D2O_{1}_KASCII.txt'.format(rootd,temperature)
	wname='exp{0}'.format(temperature)
	LoadAscii(Filename=filename, Separator='Space', Unit='Energy', Version=1, Outputworkspace=wname)
	ScaleX(InputWorkspace=wname, Factor=0.001, Operation='Multiply', Outputworkspace=wname) #Ascii Units are microeV
	outfilename='{0}/{1}.nxs'.format(rootd,wname)
	SaveNexus(InputWorkspace=wname, Filename=outfilename)

#Special for the resolution function
filename='{0}/Resolution_+ND+_RNA_+D2OASCII.txt'.format(rootd)
wname='resolution'
LoadAscii(Filename=filename, Separator='Space', Unit='Energy', Version=1,Outputworkspace=wname)
ScaleX(InputWorkspace=wname, Factor=0.001, Operation='Multiply', Outputworkspace=wname) #Ascii Units are microeV
outfilename='{0}/resolution.nxs'.format(rootd,wname)
SaveNexus(InputWorkspace=wname, Filename=outfilename)

