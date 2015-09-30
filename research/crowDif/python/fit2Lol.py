#Fit to single lorentzian
rs=('r000', 'r010', 'r025', 'r100', 'r200')
#Initial guesses for Q=0.3
initGuess={'r000':{'Height':0.07, 'Amplitude':0.8, 'FWHM':0.032},
    'r010':{'Height':0.07, 'Amplitude':0.8, 'FWHM':0.033},
    'r025':{'Height':0.08, 'Amplitude':0.8, 'FWHM':0.033},
    'r100':{'Height':0.11, 'Amplitude':0.7, 'FWHM':0.030},
    'r200':{'Height':0.15, 'Amplitude':0.7, 'FWHM':0.027}
    }

#for r in ('r000', 'r010', 'r025', 'r100', 'r200'):
for r in ('r200',):
    resolution='{0}_T015'.format(r)
    ig=initGuess[r]
    nW=9 #number of workspaces
    for iw in range(nW):
        Q=0.3+0.2*iw
        outfit='{0}_Q{1}'.format(r,Q)
        h=ig['Height']
        a=ig['Amplitude']
        f=ig['FWHM']
        function  ='name=LinearBackground,A0=0.0,A1=0.0;'
        function +='(composite=Convolution,FixResolution=false,NumDeriv=true;'
        function +='name=TabulatedFunction,Workspace={0},WorkspaceIndex={1},Scaling=1,Shift=0.001,ties=(Scaling=1);'.format(resolution,iw)
        function +='(name=DeltaFunction,Height={0};'.format(h)
        function +='name=Lorentzian,Amplitude={0},PeakCentre=0,FWHM={1},ties=(PeakCentre=0)))'.format(a,f)
        Fit(Function=function, InputWorkspace='{0}_T300'.format(r), WorkspaceIndex=iw,
            StartX=-0.12, EndX=0.12, Output=outfit, CreateOutput=1)
        ws=mtd['{0}_Parameters'.format(outfit)] #workspace with fitting parameters
        #Initial guess for next workspace
        h=ws.row(4)['Value']
        a=ws.row(5)['Value']
        f=ws.row(7)['Value']
        chi=ws.row(8)['Value']
        #print 'iw={0}, Q={1}, h={2}, a={3}, f={4}, chi={5}'.format(iw,Q,h,a,f,chi)
        print '{0} {1}'.format(Q*Q,f)