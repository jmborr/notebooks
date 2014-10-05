import numpy

def normalizeHistogram(inputworkspace, outputworkspace=None, RangeLower=None, RangeUpper=None):
	'''Set to unity the integral of the histogram. NOTE: the bin width is used as dx in the integral'''
	kwargs = {}
	if RangeLower: kwargs['RangeLower']=RangeLower
	if RangeUpper: kwargs['RangeUpper']=RangeUpper
	if not outputworkspace: outputworkspace = inputworkspace
	ConvertToHistogram(InputWorkspace=inputworkspace,OutputWorkspace='junk_normalizeHistogram')
	Integration(InputWorkspace='junk_normalizeHistogram',OutputWorkspace='junk2_normalizeHistogram', **kwargs)
	# Recall that Integration algorithm "fails" to multiply by the bin width
	ws=mtd['junk_normalizeHistogram']
	dX = ws.dataX(0)[1] - ws.dataX(0)[0]
	Scale(InputWorkspace='junk2_normalizeHistogram', Factor=dX,Operation='Multiply',OutputWorkspace='junk2_normalizeHistogram')
	Divide(LHSWorkspace=inputworkspace,RHSWorkspace='junk2_normalizeHistogram',OutputWorkspace=outputworkspace)
	DeleteWorkspace(Workspace='junk_normalizeHistogram')
	DeleteWorkspace(Workspace='junk2_normalizeHistogram')

def setToUnity(iwsg):
	'''Set I(Q,t=0)=1'''
	ExtractSingleSpectrum(InputWorkspace='{0}_fq0'.format(iwsg),WorkspaceIndex=0,OutputWorkspace='junk')
	Transpose(InputWorkspace='junk', OutputWorkspace='junk')
	Divide(LHSWorkspace='{0}_fqt.Re'.format(iwsg), RHSWorkspace='junk', OutputWorkspace='{0}_fqt.Re'.format(iwsg))
	DeleteWorkspace(Workspace='junk')

def removeFlatBackground(sqw, outputworkspace=None):
	'''S(Q,E) goes to a constant as E goes to infinity, instead of zero'''
	import numpy
	if not outputworkspace:
		outputworkspace = sqw
	else:
		CloneWorkspace(InputWorkspace=sqw, OutputWorkspace=outputworkspace)
	percent=0.05
	ws=mtd[outputworkspace]
	for ix in range( ws.getNumberHistograms() ):
		print "len=",len(ws.dataY(ix))
		N=len(ws.dataY(ix))*percent
		print 'N=',N
		average = numpy.sum(ws.dataY(ix)[-N:])/N
		print 'average=',average
		ws.dataY(ix)[:] -= average

def rebinInQ(wsname, params):
	'''Rebin along the vertical axes'''
	Transpose(InputWorkspace=wsname, OutputWorkspace=wsname)
	Rebin(InputWorkspace=wsname, Params=params, OutputWorkspace=wsname)
	Transpose(InputWorkspace=wsname, OutputWorkspace=wsname)

def clipIQ0(wsname, outname=None, mode='quadratic', **kwargs):
	'''Remove the peak of IQt. We don't know the shape of I(Q,t) in the range [0,dt], 
	where dt is the time between simulation frames. This produces a high frequency background in S(Q,E).
	The background can be removed if the extrapolate I(Q,t=0) from a fit to the nearby points I(Q,t>0).
	Since I(Q,t=0) must be one, one can always renormalize the intensities so that I(Q,0)=1, although
	this is not done in this routine
	'''
	if outname:
		CloneWorkspace(InputWorkspace=wsname, OutputWorkspace=outname)
	else:
		outname=wsname
	ws = mtd[outname] # print 'wsname=', wsname
	nk = ws.getNumberHistograms() #print 'nk=',nk
	if mode=='quadratic':
		'''fit I(Q,0) with quadratic function'''
		nt=4
		if 'nt' in kwargs.keys(): nt=kwargs['nt']
		for k in range(nk):
			Y = ws.dataY(k)
			hY = numpy.argmax(Y) # index with maximum of I(Q,t)
			X = ws.dataX(k)
			x = X[hY+1:hY+1+nt]
			if len(X) != len(Y):
				#Histogram data, the maximum may be duplicated if I(Q,t) is even. Find the highest index then
				hY = numpy.where( Y==Y[hY] )[-1]
				x = ( X[hY+1:hY+1+nt] + X[hY+2:hY+2+nt] ) / 2.0
			y = Y[hY+1:hY+1+nt]
			p = numpy.polyfit( x, y, 2)
			Y[ numpy.where( Y==Y[hY] ) ] = p[-1] # substitute the maximum intensity with the intercept of the quadratic fit


def clipSQ0(wsname, outname=None, mode='truncate', **kwargs):
	'''Remove the peak of SQE
	mode
		truncate: clipp the peak, leaving a flat central area
		lorentzian: fit a lorentzian and substite the peak for the lorentzian values
	'''
	import numpy
	if outname:
		CloneWorkspace(InputWorkspace=wsname, OutputWorkspace=outname)
	else:
		outname=wsname
	ws = mtd[outname] # print 'wsname=', wsname
	nk = ws.getNumberHistograms() #print 'nk=',nk
	if mode=='truncate':
		ratio = 0.1
		if 'ratio' in kwargs.keys(): ratio = kwargs['ratio']
		for k in range(nk):
			Y = ws.dataY(k)
			nx = len(Y) # print 'nx=',nx
			x0 = numpy.argmax(Y)  # index of histogram maximum print 'x0=',x0
			x0l = [x0,] # list of index where to clip SQE
			y0 = Y[x0] # value of histogram maximum print 'y0=', y0
			# Search adjacent points x < x0
			search = 1
			y0left = y0
			x = x0 - 1
			while search and x >= 0:
				if Y[x] > y0*ratio:
					x0l.append(x)
					y0left = Y[x]
					x = x - 1
				else:
					yleft = Y[x]
					search = 0
			search = 1
			y0right = y0
			x = x0 + 1
			while search and x <= nx:
				if Y[x] > y0*ratio:
					x0l.append(x)
					y0right = Y[x]
					x = x + 1
				else:
					yright = Y[x]
					search = 0
			# clip values
			yfinal = (yleft+yright) / 2 # print 'yfinal=', yfinal  print 'x0l=',str(x0l)
			for x in x0l:
				Y[x] = yfinal
	elif mode == 'lorentzian':
		def lorentz(x, intensity, FWHM):
			''' Lorentzian function '''
			return intensity / (FWHM**2 + x*x)
		from scipy.optimize import curve_fit
		emin = 0.0034 # assumed Basis resolution
		emax = 5*emin
		if 'emin' in kwargs.keys(): emin = kwargs['emin']
		if 'emax' in kwargs.keys(): emax = kwargs['emax']
		for k in range(nk):
			Y = ws.dataY(k)
			nY = len(Y)
			X = ws.dataX(k)
			nX = len(X)
			indexes = numpy.intersect1d( numpy.where(X>emin)[0], numpy.where(X<emax)[0] )
			y = Y[indexes]
			x = X[indexes]
			if nX != nY:
				x = (X[indexes] + X[indexes+1])/2  # deal with histogram data
			popt, pcov = curve_fit(lorentz, x, y, p0=[y[0],emin]) # Fit (x,y) with a Lorentzian1D
			# Evaluate the lorentzian in [-emin, emin] and substitute the data with it
			x0 = X[indexes[0]+1]  #extend the region one bin beyond the emin threshold
			indexes = numpy.intersect1d( numpy.where(X>-x0)[0], numpy.where(X<x0)[0] )
			x = X[indexes]
			if nX != nY:
				x = (X[indexes] + X[indexes+1])/2  # deal with histogram data
			y = lorentz(x, *popt)
			Y[indexes] = y	

def cancel_negative_intensities(inputworkspace, outputworkspace=None):
	'''Navigate the structure factor and cancel negative intensities by taking positive intensities from other energy values.
	This transformation preserves the integrated intensity'''
	if outputworkspace: 
		CloneWorkspace(InputWorkspace=inputworkspace, OutputWorkspace=outputworkspace)
	else:
		outputworkspace=inputworkspace
	ws = mtd[outputworkspace] # print 'wsname=', wsname
	nk = ws.getNumberHistograms() #print 'nk=',nk
	for k in range(nk):
		print 'k=',k
		Y = ws.dataY(k)
		en = len(Y) # number of energy points
		eh = numpy.argmax(Y) # index where maximum of SQE is located
		#Smooth the negative energies
		area=0.0 #intensity "to the left" of ie
		for ie in range(0,eh):
                        #print 'k=',k,'area=',area,'ie=',ie
			# Negative intensity !
			if Y[ie] < 0.0:
                                #print 'ie',ie,'Y[ie]=',Y[ie]
				# smooth Y[ie] grabbing intensity from Y[je<ie]	
				if area > 0.0:
					for je in range(0,ie):
						if Y[je] > 0.0:
							# intensity at Y[je] not enough to set Y[ie] to zero
							if Y[je] < abs(Y[ie]):
								Y[ie] += Y[je]
								area -= Y[je]
								Y[je] = 0.0
							else:
								Y[je] += Y[ie]
								area  += Y[ie]
								Y[ie] = 0.0
						if Y[ie] == 0.0:
							break
			# If there was not enough area in Y[je<ie], grab intensity from Y[ie<je]
			if Y[ie] < 0.0:
				for je in range(ie+1,eh):
					if Y[je] > 0.0:
						# intensity at Y[je] not enough to set Y[ie] to zero
						if Y[je] < abs(Y[ie]):
							Y[ie] += Y[je]
							Y[je] = 0.0
						else:
							Y[je] += Y[ie]
							Y[ie] = 0.0
					if Y[ie] == 0.0:
						break
			else:
				area += Y[ie] # Y[ie] could have been positive from the start
		#Smooth the positive energies
		area=0.0 #intensity "to the right" of ie
		for ie in range(en-1, eh, -1):
                        #print 'k=',k,'area=',area,'ie=',ie
			# Negative intensity !
			if Y[ie] < 0.0:
				# smooth Y[ie] grabbing intensity from Y[ie<je]	
				if area > 0.0:
					for je in range(en-1,ie,-1):
						if Y[je] > 0.0:
							# intensity at Y[je] not enough to set Y[ie] to zero
							if Y[je] < abs(Y[ie]):
								Y[ie] += Y[je]
								area -= Y[je]
								Y[je] = 0.0
							else:
								Y[je] += Y[ie]
								area  += Y[ie]
								Y[ie] = 0.0
						if Y[ie] == 0.0:
							break
			# If there was not enough area in Y[je<ie], grab intensity from Y[je<ie]
			if Y[ie] < 0.0:
				for je in range(ie-1,eh,-1):
					if Y[je] > 0.0:
						# intensity at Y[je] not enough to set Y[ie] to zero
						if Y[je] < abs(Y[ie]):
							Y[ie] += Y[je]
							Y[je] = 0.0
						else:
							Y[je] += Y[ie]
							Y[ie] = 0.0
					if Y[ie] == 0.0:
						break
			else:
				area += Y[ie] # Y[ie] could have been positive from the start

rootd='/SNSlocal/projects/jbq/POSS/ibu8t8'
erange=[-0.13, 0.0004, 0.07] # energy range without contamination from Bragg peaks
LoadDaveGrp(Filename='{0}/experiment/13K.grp'.format(rootd),isMicroEV=1,OutputWorkspace='resolution')
Rebin(InputWorkspace='resolution', Params=erange, OutputWorkspace='resolution')
normalizeHistogram('resolution')

for T in '190 230 250 270 300 320 330 340 350 360 370'.split():
	LoadDaveGrp(Filename='{0}/experiment/{1}K.grp'.format(rootd,T),isMicroEV=1,OutputWorkspace='e{0}K'.format(T))
	Rebin(InputWorkspace='e{0}K'.format(T), Params=erange, OutputWorkspace='e{0}K'.format(T))
	normalizeHistogram('e{0}K'.format(T))

rootd='/SNSlocal/projects/jbq/POSS/ibu8t8/ibu8t8c'	
for T in [180, 200, 220, 240, 260, 300, 320, 340, 360, 380]:
	LoadSassena(Filename='{0}/{1}K/fqt_inc.h5'.format(rootd,T),TimeUnit=1.0, SortByQvectors=1, OutputWorkspace='{0}Kc'.format(T))
	setToUnity('{0}Kc'.format(T))  #Set I(Q,t=0)=1
	clipIQ0('{0}Kc_fqt.Re'.format(T)) #Fit I(Q,t=0) from the nearby points. This removes the high frequency background for SQE
	Rebin(InputWorkspace='{0}Kc_fqt.Re'.format(T),Params=[-2000.5,1,2000.5],OutputWorkspace='{0}Kc_fqt.Re'.format(T))  #clip tails with poor statistics
	SassenaFFT(InputWorkspace='{0}Kc'.format(T), FFTonlyRealpart=1, DetailedBalance=1, Temp=float(T))
	rebinInQ('{0}Kc_sqw'.format(T), [0.2, 0.2, 2.0])  #same Q binning as in experiments
	clipSQ0( '{0}Kc_sqw'.format(T), outname='s{0}K'.format(T), mode='lorentzian') #Try to remove the elastic line
	rebin(InputWorkspace='s{0}K'.format(T), Params=[-2.0002, 0.0004,2.0002], OutputWorkspace='s{0}K'.format(T))
	normalizeHistogram( 's{0}K'.format(T) ) # Integral of S(Q,E) over E set to one. OK since the energy range is [-2, 2]meV, much wider than erange

#Output elastic scan ( here just S(Q,E=0) )
for T in  [180, 200, 220, 240, 260, 300, 320, 340, 360, 380]:
	buf='# Q  S(Q,E=0)\n'
	ws=mtd['{0}Kc_sqw'.format(T)]
	n=len(ws.dataX(0))/2  # center of histogram
	print 'n=',n
	for ix in range( ws.getNumberHistograms() ):
		Q=0.3+0.2*ix
		buf+='{0} {1}\n'.format(Q,ws.dataY(ix)[n])
	open('{0}/{1}K/es.dat'.format(rootd,T),'w').write(buf)
	
# Output elastic scan for a given energy range
dE=0.0034  # 3.4ueV
rootd='/SNSlocal/projects/jbq/POSS/ibu8t8/ibu8t8c'
for T in  [180, 200, 220, 240, 260, 300, 320, 340, 360, 380]:
	buf='# Q  S(Q,E=0)\n'
	ConvertToHistogram(InputWorkspace='{0}Kc_sqw'.format(T),OutputWorkspace='junk')
	Integration(InputWorkspace='junk', OutputWorkspace='junk2',RangeLower=-dE, RangeUpper=dE)
	ws=mtd['junk2']
	for ix in range( ws.getNumberHistograms() ):
		Q=0.3+0.2*ix
		buf+='{0} {1}\n'.format(Q,ws.dataY(ix)[0])
	open('{0}/{1}K/es_{2}.dat'.format(rootd,T,dE*1000),'w').write(buf)


# DSFinterp, then normalize
# The result is a quasi-elastic normalized signal interpolated from the simulation
inparams = [180, 200, 220, 240, 260, 300, 320, 340, 360, 380]
inputw = [ 's{0}K'.format(T) for T in inparams ]
targetsimparams=range(180,382,2); print 'number of interpolated temperatures=', len(targetsimparams)
outw = [ 'd{0}K'.format(T) for T in targetsimparams]
DSFinterp(Workspaces=inputw,LoadErrors=0,ParameterValues=inparams,LocalRegression=1,RegressionWindow=5,RegressionType='quadratic',TargetParameters=targetsimparams,OutputWorkspaces=outw)
for ws in outw: normalizeHistogram(ws) 

# Convolve with resolution and limit to experimental energy range.
Rebin(InputWorkspace='resolution', Params=[-0.07, 0.0004, 0.07], OutputWorkspace='junk')  # need same negative and positive domain
normalizeHistogram('junk')
for T in [180, 200, 220, 240, 260, 300, 320, 340, 360, 380]:
	ConvolveWorkspaces(Workspace1='s{0}K'.format(T), Workspace2='junk', OutputWorkspace='s{0}Kconv'.format(T))
	Rebin(InputWorkspace='s{0}Kconv'.format(T),Params=erange,OutputWorkspace='s{0}Kconv'.format(T))
for ws in outw:
	ConvolveWorkspaces(Workspace1=ws, Workspace2='junk', OutputWorkspace=ws+'conv')
	Rebin(InputWorkspace=ws+'conv',Params=erange,OutputWorkspace=ws+'conv')
DeleteWorkspace(Workspace='junk')


# Calculate EISF. Assumed that simulated clipped SQE are purely quasielastic
#  Model is Convolution( Resolution,  a*DeltaDirac + b Ssim ) + LinearBackground
elastic_I=0.5; quasie_I=500.0;  buf = '# Q  T   EISF   Chi2\n'
for index in range(8,-1,-1):
	Q = 0.3+0.2*index
	for T in [190, 230, 250, 270, 300, 320, 330, 340, 350, 360, 370]:
		fit_string    ='name=TabulatedFunction,Workspace=resolution,WorkspaceIndex={0},Scaling={1};'.format(index, elastic_I)
		fit_string +='name=TabulatedFunction,Workspace=d{0}Kconv,WorkspaceIndex={1},Scaling={2};'.format(T,index, quasie_I)
		fit_string +='name=LinearBackground,A0=0.0,A1=0.0'
		Fit(fit_string, InputWorkspace='e{0}K'.format(T), WorkspaceIndex=index, CreateOutput = 1, Output='e{0}K_EISF'.format(T) )
		ws=mtd['e{0}K_EISF_Parameters'.format(T0)]  # Workspace resulting from the fitting containing the optimal value of the fitting parameters and the Chi2
		elastic_I = ws.row(0)['Value'];  quasie_I = ws.row(1)['Value']; Chi2=ws.row(4)['Value']
		EISF = elastic_I / (elastic_I + quasie_I)
		dbuf = '{0:3.1f} {1:3f} {2:4.2f} {3:5.2f}\n'.format(Q, T, EISF, Chi2)
		buf = dbuf + buf; #print dbuf,
	buf = '&\n'+ buf
open('/tmp/EISF.dat','w').write(buf)


# DSFinterp1DFit
inparams = [180, 200, 220, 240, 260, 300, 320, 340, 360, 380]
inputw = ' '.join( [ 's{0}Kc_clipped'.format(T) for T in inparams ] )
inparamstr = ' '.join( str(T) for T in inparams )
T0=190; elastic_I=0.7; quasie_I=89.0; buf=''
for index in range(8,-1,-1):
	myfunc   ='name=TabulatedFunction,Workspace=resolution,WorkspaceIndex={0},Scaling={1};'.format(index, elastic_I)
	myfunc+='name=LinearBackground,A0=0,A1=0;'
	myfunc+='name=DSFinterp1DFit,LoadErrors=0,LocalRegression=1,RegressionType=quadratic,RegressionWindow=5,'
	myfunc+='InputWorkspaces={0},ParameterValues={1},'.format(inputw,inparamstr)
	myfunc+='WorkspaceIndex={0},Intensity={1},TargetParameter={2}'.format(index, quasie_I, T0)
	Fit(Function=myfunc, InputWorkspace='e{0}K'.format(T0), WorkspaceIndex=index, CreateOutput=1, Output='e{0}K_DSFinterp1DFit'.format(T0))
	ws = mtd['e{0}K_DSFinterp1DFit_Parameters'.format(T0)]
	elastic_I = ws.row(0)['Value'];  quasie_I = ws.row(3)['Value'];  Topt = ws.row(4)['Value'];  Chi2 = ws.row(5)['Value']
	Q = 0.3+0.2*index; dbuf = '{0:3.1f} {1:3f} {2:5.2f}\n'.format(Q, Topt, Chi2); print dbuf; buf = dbuf + buf
buf ='#Q Topt\n' + buf
open('/tmp/{0}_Topt.dat'.format(T0), 'w').write(buf)
	
# Chi-square plots with DSFinterp

T0=250; elastic_I=0.5; quasie_I=500.0; buf=''; chi2values=[];  tlist=targetparams[:]; tlist.reverse()
for index in range(8,-1,-1):
	Q = 0.3+0.2*index
	for T in tlist:
		fit_string    ='name=TabulatedFunction,Workspace=resolution,WorkspaceIndex={0},Scaling={1};'.format(index, elastic_I)
		fit_string +='name=TabulatedFunction,Workspace=d{0}Kc_clipped,WorkspaceIndex={1},Scaling={2};'.format(T,index, quasie_I)
		fit_string +='name=LinearBackground,A0=0.0,A1=0.0'
		Fit(fit_string, InputWorkspace='e{0}K'.format(T0), WorkspaceIndex=index, CreateOutput = 1, Output='e{0}K_DSFinterp'.format(T0) )
		ws=mtd['e{0}K_DSFinterp_Parameters'.format(T0)]  # Workspace resulting from the fitting containing the optimal value of the fitting parameters and the Chi2
		elastic_I = ws.row(0)['Value'];  quasie_I = ws.row(1)['Value']; Chi2=ws.row(4)['Value']
		dbuf = '{0:3f} {1:5.2f}\n'.format(T, Chi2) #dbuf = '{0:3.1f} {1:3f} {2:5.2f}\n'.format(Q, T, Chi2);
		buf = dbuf + buf; #print dbuf,
	buf = '&\n'+ buf
buf ='#Q T Chi2\n' + buf
open('/tmp/{0}K_DSFinterp.dat'.format(T0), 'w').write(buf)

targetparams=range(180,382,2); print 'number of interpolated temperatures=', len(targetparams)
outw = [ 'd{0}Kc_clipped'.format(T) for T in targetparams]
	

# Texp vs Tsim Chi2

# Fit interpolated experimental structure factors versus interpolated simulated structure factors
inparams = [190, 230, 250, 270, 300, 320, 330, 340, 350, 360, 370]
inputw = [ 'e{0}K'.format(T) for T in inparams ]
targetparams=range(190,372,2); print 'number of interpolated temperatures=', len(targetparams)
outw = [ 'de{0}K'.format(T) for T in targetparams]
DSFinterp(Workspaces=inputw,LoadErrors=1,ParameterValues=inparams,LocalRegression=0,RegressionWindow=5,RegressionType='quadratic',TargetParameters=targetparams,OutputWorkspaces=outw)
for ws in outw: normalizeHistogram(ws)

elastic_I=0.5; quasie_I=0.5; buf=''; 
Txs=targetparams[:]; Tss=targetsimparams[:]
for index in range(9):
	buf = '#Texp Tsim Chi2\n'
	Q = 0.3+0.2*index
	for Tx in Txs:
		for Ts in Tss:
			print "Q=",Q
			fit_string    ='name=TabulatedFunction,Workspace=resolution,WorkspaceIndex={0},Scaling={1};'.format(index, elastic_I)
			fit_string +='name=TabulatedFunction,Workspace=d{0}Kconv,WorkspaceIndex={1},Scaling={2};'.format(Ts,index, quasie_I)
			fit_string +='name=LinearBackground,A0=0.0,A1=0.0'
			Fit(fit_string, InputWorkspace='de{0}K'.format(Tx), WorkspaceIndex=index, CreateOutput = 1, Output='junk')
			ws=mtd['junk_Parameters']  # Workspace resulting from the fitting containing the optimal value of the fitting parameters and the Chi2
			Chi2=ws.row(4)['Value']
			dbuf = '{0:3f} {1:3f} {2:5.2f}\n'.format(Tx, Ts, Chi2); print dbuf,
			buf += dbuf
	open(rootd+'/analysis/Texp_Tsim_Chi2_Q{0:3.1f}.dat'.format(Q),'w').write(buf)


