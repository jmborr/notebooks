# imports for compatibility with python 3.5
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import numpy as np
import MDAnalysis as mda
import scatter.scatterMDA as sc
from pdb import set_trace as tr

Q = sc.Qlist(0.45, 1.0, 0.05) # same Q-values as in experiment

# Find scattering lenghts of the hydrogens
topology = "q300.pdb"
trajectory = "q300.dcd"
all =  mda.Universe(topology, trajectory)
hydr = all.select_atoms("name H*")
b = sc.scatL(hydr)['incoherent'] # Scattering lengths

# Extract coordinates of hydrogens every 2.0ps (skip=20
# since frames in original trajectory distant by 0.1ps)
skip=20
frames = sc.getFrames(all, hydr, skip=skip, format='fac')
# There are 5000 frames, each separated by 2ps

# Load flag that tells at each frame whether each hydrogen
# is bound to the surface
bFlag = np.loadtxt("300.boundFlag.Hydr.dat", dtype=np.int32)
bFlag = bFlag[1::skip] # select every 2ps

""" Calculate IIv2(Q,t).
  Arguments:
   nsampling: number of sampling t0 points for each time lapse "t"
   nt: Calculate II for this number of "t"
   dt: separation between consecutive "t"'s, in number of frames
   begt: value of first "t", in number of frames
   c1: extra correlation, of shape (#frames, #atoms). Will multiply
       each self-interference term by c1(t0}(C1(t0+t)
"""
factors = sc.IIv2(frames, b, Q, bFlag, nsampling=5000, nt=5000, dt=1, begt=0)
# Save structure factors to file
sc.saveIISassenaFormat(factors['sf'], Q, "fqt_incBound.h5")

factors = sc.IIv2(frames, b, Q, 1-bFlag, nsampling=5000, nt=5000, dt=1, begt=0)
# Save structure factors to file
sc.saveIISassenaFormat(factors['sf'], Q, "fqt_incUnbound.h5")



