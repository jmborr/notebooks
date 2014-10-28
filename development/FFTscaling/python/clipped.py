wd='/projects/development/FFTscaling'


LoadSassena(Filename=wd+'/T450.100ps/fqt_inc_styH.h5',TimeUnit=0.01,SortByQVectors=1,OutputWorkspace='10fs')
LoadSassena(Filename=wd+'/T450.1ns/fqt_inc_styH.h5',TimeUnit=0.1,SortByQVectors=1,OutputWorkspace='100fs')
LoadSassena(Filename=wd+'/T450.10ns/fqt_inc_styH.h5',TimeUnit=1.0,SortByQVectors=1,OutputWorkspace='1000fs')

SassenaFFT(InputWorkspace='10fs',FFTonlyRealpart=1,DetailedBalance=1,Temp=450)
SassenaFFT(InputWorkspace='100fs',FFTonlyRealpart=1,DetailedBalance=1,Temp=450)
SassenaFFT(InputWorkspace='1000fs',FFTonlyRealpart=1,DetailedBalance=1,Temp=450)

#Study the effects of cutting the long-time tails of I(Q,t)
L=10000; dL = 10; iQ=9
while L > 9900:
	iqt='IQt_L{0}'.format(L);    sqe='SQE_L{0}'.format(L)
	ExtractSingleSpectrum(InputWorkspace='1000fs_fqt.Re',WorkspaceIndex=iQ,OutputWorkspace=iqt)
	Rebin(InputWorkspace=iqt,Params=[-L+0.5,1,L+0.5],OutputWorkspace=iqt)
	FFT(InputWorkspace=iqt,OutputWorkspace=sqe)
	L -= dL

#Study the effects of cutting the long-time tails of I(Q,t)
L=10000; dL = 100; iQ=9
while L > 9000:
	iqt='IQt_L{0}'.format(L);    sqe='SQE_L{0}'.format(L)
	ExtractSingleSpectrum(InputWorkspace='1000fs_fqt.Re',WorkspaceIndex=iQ,OutputWorkspace=iqt)
	Rebin(InputWorkspace=iqt,Params=[-L+0.5,1,L+0.5],OutputWorkspace=iqt)
	FFT(InputWorkspace=iqt,OutputWorkspace=sqe)
	L -= dL

#Study the effects of cutting the long-time tails of I(Q,t)
L=10000; dL = 1000; iQ=9
while L > 4000:
	iqt='IQt_L{0}'.format(L);    sqe='SQE_L{0}'.format(L)
	ExtractSingleSpectrum(InputWorkspace='1000fs_fqt.Re',WorkspaceIndex=iQ,OutputWorkspace=iqt)
	Rebin(InputWorkspace=iqt,Params=[-L+0.5,1,L+0.5],OutputWorkspace=iqt)
	FFT(InputWorkspace=iqt,OutputWorkspace=sqe)
	L -= dL
