import sys
sys.path.append("/usr/local/pizza/pizza-2Jul14/src")  # for pizza.py modules
import data
sys.path.append("/SNSlocal/projects/jbq/peComplex/python") # for xyz module
import xyzReader
import argparse
import numpy
from pdb import set_trace as tr

parser = argparse.ArgumentParser(description='''Create output dump file from XYZ components
Assumed the following files exists in working directory:
  head.xyz  Ncions.xyz  Pcions.xyz  poly.xyz  tail.xyz
''')
parser.add_argument('indata', type=str, help='input LAMMPS DATA file')
parser.add_argument('outdump', type=str, help='output LAMMPS DUMP file')
parser.add_argument('--nframes', type=int, default=-1, help='number of frames. All frames if not specified')

args = parser.parse_args()

# Template for timestep item
timestep_item_tpl='''ITEM: TIMESTEP
_TIMESTEP_
'''

# Number of atoms item
# Load xyz files
poly_filehandle=open("poly.xyz")
poly=xyzReader.read_xyz(poly_filehandle)
head_filehandle=open("head.xyz")
head=xyzReader.read_xyz(head_filehandle)
tail_filehandle=open("tail.xyz")
tail=xyzReader.read_xyz(tail_filehandle)
ncions_filehandle=open("Ncions.xyz")
ncions=xyzReader.read_xyz(ncions_filehandle)
pcions_filehandle=open("Pcions.xyz")
pcions=xyzReader.read_xyz(pcions_filehandle)
natoms=len(poly.atomtypes)+len(head.atomtypes)+len(tail.atomtypes)+len(ncions.atomtypes)+len(pcions.atomtypes)
natoms_items='''ITEM: NUMBER OF ATOMS
{0}
'''.format(natoms)

# Create box item
d=data.data(args.indata)
x=d.headers['xlo xhi']
y=d.headers['ylo yhi']
z=d.headers['zlo zhi']
LxH=(x[1]-x[0])/2
LyH=(y[1]-y[0])/2
LzH=(z[1]-z[0])/2
LxyzH=numpy.array([LxH,LyH,LzH])
box_item='''ITEM: BOX BOUNDS pp pp pp
{0:7.4f} {1:7.4f}
{2:7.4f} {3:7.4f}
{4:7.4f} {5:7.4f}
'''.format(x[0],x[1],y[0],y[1],z[0],z[1])

# Template for Atoms item
atoms_item_tpl='''ITEM: BOX BOUNDS pp pp pp
_ATOMS_
'''

def print_atominfo(iframe,iatom,mol,species,index,prevxyz,previjk):
    '''Helper funcion'''
    xyz=numpy.array(species.coords[index])
    r=' '.join(['{0:10.6f}'.format(x) for x in species.coords[index]])
    s=' 0 0 0' # original cell coordinates
    if iframe:
        # calculate new cell coordinates
        xyz_old=prevxyz[iatom-1]
        ijk_old=previjk[iatom-1]
        ijk=ijk_old-numpy.modf((xyz-xyz_old)/LxyzH)[1]
        s=' '.join(['{0}'.format(int(u)) for u in ijk])
        previjk[iatom-1]=ijk #update previous cell coordinates
    prevxyz[iatom-1]=xyz #update previous atomic wrapped coordinates
    return '{0} {1} {2} {3} {4}\n'.format(iatom,mol,species.atomtypes[index],r,s)

# Iterate over frames.
out_filehandle=open(args.outdump,'w')
iframe=0
prevxyz=numpy.zeros([natoms,3]) #previous frame wrapped coordinates
previjk=numpy.zeros([natoms,3],dtype=int) #previous frame cell coordinates
while ncions:
    print iframe
    buf =timestep_item_tpl.replace('_TIMESTEP_',str(iframe))  # add timestep item
    buf+=natoms_items  # add number of atoms items
    buf+=box_item # add box item
    buf+='ITEM: ATOMS id mol type x y z ix iy iz\n'
    iatom=1
    mol=1 # molecule number
    #Add polyelectrolite
    for i in range(len(poly.atomtypes)):
        buf+=print_atominfo(iframe,iatom,mol,poly,i,prevxyz,previjk)
        iatom+=1
    mol+=1
    #Add surfactant
    tail_length=11
    for i in range(len(head.atomtypes)):
        buf+=print_atominfo(iframe,iatom,mol,head,i,prevxyz,previjk)
        iatom+=1
        for j in range(tail_length):
            buf+=print_atominfo(iframe,iatom,mol,tail,i*tail_length+j,prevxyz,previjk)
            iatom+=1
        mol+=1 #next surfactant molecule
    #Add negative counterions
    for i in range(len(ncions.atomtypes)):
        buf+=print_atominfo(iframe,iatom,mol,ncions,i,prevxyz,previjk)
        iatom+=1
        mol+=1
    #Add positive counterions
    for i in range(len(pcions.atomtypes)):
        buf+=print_atominfo(iframe,iatom,mol,pcions,i,prevxyz,previjk)
        iatom+=1
        mol+=1
    #Write to file
    out_filehandle.write(buf)
    # Read next frame
    poly=xyzReader.read_xyz(poly_filehandle)
    head=xyzReader.read_xyz(head_filehandle)
    tail=xyzReader.read_xyz(tail_filehandle)
    ncions=xyzReader.read_xyz(ncions_filehandle)
    pcions=xyzReader.read_xyz(pcions_filehandle)
    iframe+=1
    if args.nframes>0 and iframe>=args.nframes:
        sys.exit(0)
