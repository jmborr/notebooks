from mantid.simpleapi import LoadNexus, DSFinterp, mtd

# Input structure factors S(Q,E,K) and associated values for the dihedral barrier K.
# These structure factors were derived from molecular dynamics simulations of a small crystal
# of octa-methyl silsesquioxanes molecules (http://en.wikipedia.org/wiki/Silsesquioxane).
# The dihedral barrier referes to the energy barrier to methyl rotations in the molecules
workspaces=[]
parametervalues=[]
for K in '0.040  0.050  0.060  0.070  0.080  0.090  0.100  0.101  0.102  0.103  0.104  0.105  0.106  0.107  0.108  0.109  0.110  0.111  0.112  0.113  0.114  0.115  0.116  0.117  0.118  0.120'.split():
	LoadNexus(Filename='simK{0}.nxs'.format(K), OutputWorkspace='simK{0}'.format(K))
	workspaces.append('simK{0}'.format(K))
	parametervalues.append(float(K))

# Output structure factors for the dihedral values of interest
outputworkspaces = [ 'intK{0}'.format(K) for K in '0.040 0.040 0.050 0.055 0.060 0.065 0.070 0.075 0.080 0.085 0.090 0.095 0.100 0.105 0.110 0.115 0.120'.split()]
targetparameters = [ float(x) for x in '0.040 0.045 0.050 0.055 0.060 0.065 0.070 0.075 0.080 0.085 0.090 0.095 0.100 0.105 0.110 0.115 0.120'.split()]

# Interpolation
DSFinterp(Workspaces=workspaces, LoadErrors=False, ParameterValues=parametervalues,LocalRegression=True,RegressionWindow=6, RegressionType='quadratic',TargetParameters=targetparameters,OutputWorkspaces=outputworkspaces)

# Plot the interpolated structure factors versus energy E and dihedral barrier K.
# Q is kept constant in this plot at 0.9A^(-1)
try:
  from mpl_toolkits.mplot3d import axes3d
  import matplotlib.pyplot as plt
  from matplotlib import cm
except Exception as e:
  print str(e)
  print 'You need to install certain matplotlib libraries to see the plots!'
import numpy as np

X = mtd['intK0.040'].dataX(7)  # energy values
dX=0.004 # energy bin, in meV
nX = len(X)
nY = 17
dY = 0.005
Y = 0.040+0.005*np.arange(nY)  # dihedral barrier values
X, Y = np.meshgrid(X, Y)
Z = []
for K in '0.040 0.040 0.050 0.055 0.060 0.065 0.070 0.075 0.080 0.085 0.090 0.095 0.100 0.105 0.110 0.115 0.120'.split():
  Z +=	list( mtd['intK{0}'.format(K)].dataY(7) ) # Histogram with index=7 corresponds to Q=0.9A^(-1)
Z = np.array(Z).reshape(nY,nX)
Z = np.log(Z) # logarithmic scale

print "\nPlotting the Structure factors...\n"
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)

cset = ax.contourf(X, Y, Z, zdir='z', offset=-9, cmap=cm.coolwarm)
#cset = ax.contourf(X, Y, Z, zdir='x', offset=-0.18, cmap=cm.coolwarm)
cset = ax.contourf(X, Y, Z, zdir='y', offset=0.18, cmap=cm.coolwarm)

ax.set_xlabel('E(meV)')
ax.set_xlim(-0.18, 0.1)
ax.set_ylabel('K(Kcal/mol)')
ax.set_ylim(0.040, 0.18)
ax.set_zlabel('S(E,K|Q)')
ax.set_zlim(-9, 0)


plt.show()
