import os
 
#Reduced files from /SNS/BSS/IPTS-11889
 
data_dir = "/SNS/users/ele/Desktop/ElasticScanScript/ReductionRegular/"
van_file = data_dir+"/BASIS_32264_sqw.nxs"
 
save_dir = "/SNS/users/ele/Desktop/ElasticScanScript/ReductionRegular/divided/"
 
run_min = 41137
run_max = 41277
 
for i in range(run_min, run_max+1):
#for i in range(run_min, run_min+1):
    file_name = "BASIS_" + str(i) + "_sqw.nxs"
    data_file = data_dir + "/" + file_name
    print data_file
    if (os.path.isfile(data_file)):
        data_ws = Load(data_file)
        van_ws = Load(van_file)
 
        div_ws = data_ws / van_ws
         
        ####################################################
        # No need for the two lines below because the files have been recently reduced #
        ####################################################
        #LoadInstrument(Workspace=div_ws.name(), InstrumentName='BASIS')
        #ClearMaskFlag(Workspace=div_ws.name())
 
        out_name = "BASIS_" + str(i) + "_divided.dat"
        SaveDaveGrp(div_ws, save_dir+"/"+out_name, True)
        out_name = "BASIS_" + str(i) + "_divided_sqw.nxs"
        SaveNexus(InputWorkspace=div_ws, Filename= save_dir+"/"+out_name)
    else:
        print "Skipped", file_name
