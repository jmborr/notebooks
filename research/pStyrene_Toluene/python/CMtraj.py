import os
import argparse
import numpy
import pytraj
import parmed
from utilities.path import which
from tempfile import mkstemp, mkdtemp
from utilities.readingWritingFiles import write_from_numpy
from pdb import set_trace as tr
from copy import deepcopy

if __name__ == '__main__':
    argsdef={ 'aggregrate':('all', 'byres'),
    }
    parser = argparse.ArgumentParser(description="""create PBD and trajectory files of center of mass(es).
    Example: CMtrajectory.py inpdb indcd mask  outpdb outdcd""")
    parser.add_argument('topfile', help='input topology as a PDB file')
    parser.add_argument('trajfile', help='input trajectory')
    parser.add_argument('mask', type=str, help='amber mask for residues to be considered')
    parser.add_argument('outtop', type=str, help='PDB file that will serve as "topology" of the output trajectory')
    parser.add_argument('outtraj', type=str, help='trajectory containing the center of mass of the residues')
    parser.add_argument('--aggregate', type=str, default=argsdef['aggregrate'][0], help='CoM options:'+','.join(argsdef['aggregrate'])+'. Default=all')
    p_args=parser.parse_args()

    traj = pytraj.iterload(p_args.trajfile, p_args.topfile)  # faster than pytraj.load
    #traj = pytraj.iterload(p_args.trajfile, p_args.topfile, frame_slice=(0,256,3)) #from first frame to frame 255, skipping two frames
    #traj = pytraj.load(p_args.trajfile, p_args.topfile, mask=p_args.mask, frame_indices=numpy.arange(1,9)) #much slower than iterload
    ptop = parmed.load_file(p_args.topfile) # parmed topology object
    masktop = traj.top[p_args.mask] #topology of the mask atoms
    CMselection = list() #Used only for the output PDB file
    if p_args.aggregate == 'byres':
        nres = masktop.n_residues
        CoM = numpy.empty((nres,traj.n_frames,3))
        ires = 0
        for residue in masktop.residues:
            try:
                resnum = residue.resnum
            except:
                resnum = int(residue.name)
            presidue = ptop.residues[resnum-1]
            first_atom = presidue[0] # selected atom to represent the whole residue in the output PDB file
            CMselection.append('@{0}'.format(first_atom.number))
            #Next, calculate CM for given residue in all frames
            v = pytraj.center_of_mass(traj,mask='({0})&(:{1})'.format(p_args.mask,resnum))
            CoM[ires] = v
            ires += 1
        CoM = CoM.transpose((1,0,2))  #shape=(n_frames, n_residues, 3)
        CMselection = ','.join(CMselection) # mask containing the first atoms of each considered residue
    else:
        CoM = pytraj.center_of_mass(traj)
        CoM = CoM.reshape((traj.n_frames, 1, 3))
        first_residue = ptop.residues[0]
        first_atom = first_residue[0]
        CMselection = '@{0}'.format(first_atom.number)
    outtop = traj.top[CMselection]
    for i in range(outtop.n_atoms):
        atom = outtop[i]
        atom.set_mol(0)
    stripped_traj = traj[CMselection]
    pytraj.write_traj(filename=p_args.outtop, traj=stripped_traj, frame_indices=[0,], top=outtop, overwrite=True)
    pytraj.write_traj(p_args.outtraj, CoM, top=outtop, overwrite=True)
