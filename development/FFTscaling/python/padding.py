
wd='/projects/development/FFTscaling'


LoadSassena(Filename=wd+'/T450.100ps/fqt_inc_styH.h5',TimeUnit=0.01,SortByQVectors=1,OutputWorkspace='10fs')
LoadSassena(Filename=wd+'/T450.1ns/fqt_inc_styH.h5',TimeUnit=0.1,SortByQVectors=1,OutputWorkspace='100fs')
LoadSassena(Filename=wd+'/T450.10ns/fqt_inc_styH.h5',TimeUnit=1.0,SortByQVectors=1,OutputWorkspace='1000fs')

SassenaFFT(InputWorkspace='10fs',FFTonlyRealpart=1,DetailedBalance=1,Temp=450)
SassenaFFT(InputWorkspace='100fs',FFTonlyRealpart=1,DetailedBalance=1,Temp=450)
SassenaFFT(InputWorkspace='1000fs',FFTonlyRealpart=1,DetailedBalance=1,Temp=450)


#Study the effects of padding the data with one zero
ExtractSingleSpectrum(InputWorkspace='1000fs_fqt.Re',WorkspaceIndex=9,OutputWorkspace='IQt_notpadded')
FFT(InputWorkspace='IQt_notpadded',OutputWorkspace='SQE_notpadded')
ExtractSingleSpectrum(InputWorkspace='SQE_notpadded',WorkspaceIndex=3,OutputWorkspace='SQE_notpadded')
Rebin(InputWorkspace='IQt_notpadded',Params=[-10000.5,1,10000.5],OutputWorkspace='IQt_padded')
FFT(InputWorkspace='IQt_padded',OutputWorkspace='SQE_padded')
ExtractSingleSpectrum(InputWorkspace='SQE_padded',WorkspaceIndex=3,OutputWorkspace='SQE_padded')
