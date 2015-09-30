rootd='/SNSlocal/projects/jbq/crowDif/experiment/data'
run2rT={26518:'r000_T300', 26587:'r000_T015', 26588:'r200_T300', 26653:'r200_T015', 26654:'r010_T300', 26723:'r010_T015', 26728:'empty_can', 26733:'r100_T300', 26802:'r100_T015', 26804:'r025_T300', 26873:'r025_T015' }

for (run_number, rT) in run2rT.items():
    LoadDaveGrp(Filename=rootd+'/BASIS_{0}_1run_divided.dat'.format(run_number), IsMicroEV=1, OutputWorkspace=rT)
    SaveNexus(InputWorkspace=rT, Filename=rootd+'/{0}.nxs'.format(rT))

#print the renaming convention
for (run_number, rT) in run2rT.items():
    print 'BASIS_{0}_1run_divided.dat --> {1}.nxs'.format(run_number,rT)
