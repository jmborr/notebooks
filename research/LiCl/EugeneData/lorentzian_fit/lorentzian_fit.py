rootdir='/projects/research/LiCl/EugeneData'
LoadNexus(Filename=rootdir+'/elasticLine.nxs', OutputWorkspace='elasticLine')

nQ=4 # four values of momentum transfer modulus
Iel=0.137311
Ilo=0.598152
fwhm=0.0533234
chi2=0

for T in [290, 280, 270, 260, 250, 240, 230, 220, 210, 200, 190]:
	print T
	LoadNexus(Filename=rootdir+'/LiCl_{0}K.nxs'.format(T), OutputWorkspace='LiCl_{0}K'.format(T))
	iw=3
	while iw >= 0:
		Q=0.3+0.3*iw # Q-value
		fit_string    = 'name=TabulatedFunction,Workspace=elasticLine,WorkspaceIndex={0},Scaling={1};'.format(iw,Iel)
		fit_string += '(composite=Convolution,FixResolution=true,NumDeriv=true;name=TabulatedFunction,'
		fit_string += 'Workspace=elasticLine,WorkspaceIndex={0},Scaling=1;'.format(iw)
		fit_string += 'name=Lorentzian,Amplitude={0},PeakCentre=0.0,FWHM={1});'.format(Ilo, fwhm) 
		fit_string += 'name=LinearBackground,A0=0.0,A1=0.0'
		Fit(fit_string, InputWorkspace='LiCl_{0}K'.format(T), WorkspaceIndex=iw, StartX=-0.1, EndX=0.1, CreateOutput = 1 )
		ws=mtd['LiCl_{0}K_Parameters'.format(T)]
		Iel=ws.row(0)['Value']
		Ilo=ws.row(2)['Value']
		fwhm=ws.row(4)['Value']
		chi2=ws.row(7)['Value']
		print 'T=',T,'Q=',Q,'fwhm=',fwhm,'chi2=',chi2
		RenameWorkspace(InputWorkspace='LiCl_{0}K_Parameters'.format(T) ,OutputWorkspace='LiCl_{0}K_index{1}_Parameters'.format(T,iw))
		RenameWorkspace(InputWorkspace='LiCl_{0}K_Workspace'.format(T) ,OutputWorkspace='LiCl_{0}K_index{1}_Workspace'.format(T,iw))
		RenameWorkspace(InputWorkspace='LiCl_{0}K_NormalisedCovarianceMatrix'.format(T) ,OutputWorkspace='LiCl_{0}K_index{1}_NormalisedCovarianceMatrix'.format(T,iw))
		iw -= 1