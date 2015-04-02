#!/usr/bin/env python

from pdb import set_trace as tr
import numpy


if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='calculate first and second order correlations')
  parser.add_argument('--infile',help='input file with dihedral angles')
  parser.add_argument('--outfile',help='output file')
  args = parser.parse_args()

  phis=(numpy.pi/180.)*numpy.loadtxt(args.infile, comments='#', usecols=(1,))

  # First order Legendre Polynomial autocorrelation <cos(phi(t)-phi(t+t0)>_t0
  # We use cos(a-b)=cos(a)cos(b)+sin(a)sin(b)
  cosphis=numpy.cos(phis)
  sinphis=numpy.sin(phis)
  c1=numpy.correlate(cosphis,cosphis,mode='same')+numpy.correlate(sinphis,sinphis,mode='same')
  c1/=max(c1)

  # Second order Legendre Polynomial autocorrelation  <0.5*(3*cos^2(phi(t)-phi(t+t0))-1)>_t0
  # We use the relation cos^2(x)=(1+cos(2x))/2 and cos(a-b)=cos(a)cos(b)+sin(a)sin(b)
  cos_2phis=numpy.cos(2*phis)  # cos(2*phi)
  sin_2phis=numpy.cos(2*phis)  # sin(2*phi)
  c2=numpy.correlate(cos_2phis,cos_2phis,mode='same')+numpy.correlate(sin_2phis,sin_2phis,mode='same') #<cos(2phi(t)-2phi(t+t0)>_t0
  c2=(1+c2)/2  # <cos^2(phi(t)-phi(t+t0)>_t0
  c2=(3*c2-1)/2 
  c2/=max(c2)

  #ouput
  n=len(phis)
  nh=n/2
  dt=1.0  #picoseconds
  buf='# time   c1    c2\n'
  for i in range(nh,n):
    buf+='{0:4.0f} {1:7.4f} {2:7.4f}\n'.format( dt*(i-nh), c1[i], c2[i] )
    
  outfile='/home/jbq/Downloads/junk.dat'
  open(args.outfile,'w').write(buf)
