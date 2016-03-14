import sys
import os
sys.path.append("/SNSlocal/sw/pizza/pizza-2Jul14/src")  # for pizza.py modules
import data as pizzadata
sys.path.append("/SNSlocal/projects/jbq/PEblock/python") # for xyzReader.py module
import xyzReader
import argparse
import numpy
from pdb import set_trace as tr

parser = argparse.ArgumentParser(description='''Create output dump file from XYZ components
Assumed the following files exists in source directory:
  poly1.dat poly2.dat polyCharge.dat head.dat tail.dat  Ncions.xyz  Pcions.xyz
''')
parser.add_argument('sourcedir', type=str, help='directory where the .dat files reside')
parser.add_argument('indata', type=str, help='input LAMMPS DATA file')
parser.add_argument('nchains', type=str, help='number of polymer chains')
parser.add_argument('outdump', type=str, help='output LAMMPS DUMP file')
parser.add_argument('--nframes', type=int, default=-1, help='number of frames. All frames if not specified')

args = parser.parse_args()

# Template for timestep item
timestep_item_tpl='''ITEM: TIMESTEP
_TIMESTEP_
'''

# Number of atoms item
# Load xyz files
poly1_filehandle=open(os.path.join(args.sourcedir,"poly1.dat"))
poly1=xyzReader.read_xyz(poly1_filehandle)
poly2_filehandle=open(os.path.join(args.sourcedir,"poly2.dat"))
poly2=xyzReader.read_xyz(poly2_filehandle)
polyCharge_filehandle=open(os.path.join(args.sourcedir,"polyCharge.dat"))
polyCharge=xyzReader.read_xyz(polyCharge_filehandle)
head_filehandle=open(os.path.join(args.sourcedir,"head.dat"))
head=xyzReader.read_xyz(head_filehandle)
tail_filehandle=open(os.path.join(args.sourcedir,"tail.dat"))
tail=xyzReader.read_xyz(tail_filehandle)
ncions_filehandle=open(os.path.join(args.sourcedir,"Ncions.dat"))
ncions=xyzReader.read_xyz(ncions_filehandle)
pcions_filehandle=open(os.path.join(args.sourcedir,"Pcions.dat"))
pcions=xyzReader.read_xyz(pcions_filehandle)
natoms=len(poly1.atomtypes)+len(polyCharge.atomtypes)+len(poly2.atomtypes)+\
        len(head.atomtypes)+len(tail.atomtypes)+len(ncions.atomtypes)+len(pcions.atomtypes)
natoms_items='''ITEM: NUMBER OF ATOMS
{0}
'''.format(natoms)

# Create box item
d=pizzadata.data(os.path.join(args.sourcedir,args.indata))
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
    '''Helper funcion
    species: xyzReader for each particle type
    index: current index of the species reader (restarts at 0 every frame)
    iatom: absolute atom number (starts at 1, does not restart every frame)
    mol: molecule number (restarts at 1 every frame)
    '''
    xyz=numpy.array(species.coords[index]) #cartesian (wrapped) coordinates
    r=' '.join(['{0:10.6f}'.format(x) for x in species.coords[index]])
    s=' 0 0 0' # original cell coordinates, because input trajectory is wrapped
    if iframe:
        # calculate new cell coordinates
        xyz_old=prevxyz[iatom-1]  #previous cartesian (wrapped) coordinates
        ijk_old=previjk[iatom-1]  #previous cell coordinates
        ijk=ijk_old-numpy.modf((xyz-xyz_old)/LxyzH)[1]
        s=' '.join(['{0}'.format(int(u)) for u in ijk]) #new cell coordinates
        previjk[iatom-1]=ijk #update previous cell coordinates
    prevxyz[iatom-1]=xyz #update previous atomic wrapped coordinates
    return '{0} {1} {2} {3} {4}\n'.format(iatom,mol,species.atomtypes[index],r,s)

# Iterate over frames.
out_filehandle=open(args.outdump,'w')
iframe=0
prevxyz=numpy.zeros([natoms,3]) #wrapped coordinates of the previous frame
previjk=numpy.zeros([natoms,3],dtype=int) #cell coordinates of the previous frame
while ncions:
    print iframe
    buf =timestep_item_tpl.replace('_TIMESTEP_',str(iframe))  # add timestep item
    buf+=natoms_items  # add number of atoms items
    buf+=box_item # add box item
    buf+='ITEM: ATOMS id mol type x y z ix iy iz\n'
    iatom=1
    mol=1 # molecule number
    ####################
    #Add polyelectrolite chains (poly1+poly2+poly3)
    #One chain is made of two blocks. The first block is made of 6 repeats, and
    #Each repeat contains 4 poly1 particles and 1 polyCharge partile. The
    # second block contains 30 poly2 particles
    ####################
    chain_index=1
    nrepeats=6         # 6 repeats in the first block
    npoly_in_repeat=4  # 4 poly1 particles in each repeat
    poly1_index=0      # must be cleared every frame
    npoly2_in_chain=30 # 30 poly2 particles on each chain
    poly2_index=0      # must be cleared every frame
    for chain_index in range(args.nchains):
        for repeat_index in range(nrepeats):
            for j in range(npoly_in_repeat):
                buf+=print_atominfo(iframe,iatom,mol,poly1,poly1_index,prevxyz,previjk)
                poly1_index+=1
                iatom+=1
            buf+=print_atominfo(iframe,iatom,mol,polyCharge,polyCharge_index,prevxyz,previjk)
            iatom+=1
        for j in range(npoly2_in_chain):
            buf+=print_atominfo(iframe,iatom,mol,poly2,poly2_index,prevxyz,previjk)
            poly2_index+=1
            iatom+=1
        mol+=1 # next polymer chain
    ###############    
    #Add surfactant (head+tail)
    #One surfactant contains one head and 11 tail particles.
    ###############    
    tail_length=11
    tail_index=0
    for i in range(len(head.atomtypes)):
        buf+=print_atominfo(iframe,iatom,mol,head,i,prevxyz,previjk)
        iatom+=1
        for j in range(tail_length):
            buf+=print_atominfo(iframe,iatom,mol,tail,tail_index,prevxyz,previjk)
            tail_index+=1
            iatom+=1
        mol+=1 #next surfactant molecule
    #########################        
    #Add negative counterions
    #########################
    for i in range(len(ncions.atomtypes)):
        buf+=print_atominfo(iframe,iatom,mol,ncions,i,prevxyz,previjk)
        iatom+=1
        mol+=1
    #########################        
    #Add positive counterions
    #########################    
    for i in range(len(pcions.atomtypes)):
        buf+=print_atominfo(iframe,iatom,mol,pcions,i,prevxyz,previjk)
        iatom+=1
        mol+=1
    ##################################
    #Write to file and Read next frame
    ##################################
    out_filehandle.write(buf)
    poly=xyzReader.read_xyz(poly_filehandle)  #advance one frame
    polyCharge=xyzReader.read_xyz(poly_filehandle)  #advance one frame
    poly2=xyzReader.read_xyz(poly_filehandle)  #advance one frame
    head=xyzReader.read_xyz(head_filehandle)
    tail=xyzReader.read_xyz(tail_filehandle)
    ncions=xyzReader.read_xyz(ncions_filehandle)
    pcions=xyzReader.read_xyz(pcions_filehandle)
    iframe+=1
    if args.nframes>0 and iframe>=args.nframes:
        sys.exit(0)
