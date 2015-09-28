rootd='/SNSlocal/projects/jbq/pStyrene_Toluene/experiment/SNS/BASIS/data'
run2T={50911:10, 50912:97, 50913:137, 50914:177, 50915:217, 50916:257, 50917:297, 50918:337}

for run, T in run2T.items():
    wsname = 'BASIS_{0}'.format(T)
    LoadDaveGrp(Filename=rootd+'/BASIS_{0}_1run_divided.dat'.format(run), isMicroEV=1, OutputWorkspace=wsname)
    SaveNexus(InputWorkspace=wsname, Filename=rootd+'/BASIS_{0}.nxs'.format(T))
    