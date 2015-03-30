import argparse
from pdb import set_trace as tr

def read_residue(handle, firstline):
    '''Read one residue'''
    residue=[firstline,]
    resid_curr = firstline[23:26]
    line = handle.readline()
    resid = line[23:26]
    while resid == resid_curr:
        residue.append(line)
        line = handle.readline()
        if not line:
            return residue, line
        resid = line[23:26]
    return residue, line

def check_residues_same_length(inX):
    '''Check that residues of the same type have the same number of atoms'''
    all = {}
    handle =  open(inX)
    line = handle.readline()
    # scan mutated.1 file and load the residues into memory
    while line:
        if line[0:4] != 'ATOM':
            line = handle.readline()
            next
        resname = line[17:20]
        if resname not in all.keys():
            all[resname]=[]
        residue, line = read_residue(handle, line)
        all[resname].append(residue)
    # check the lenghts of each residue type
    for resname, residues in all.items():
        nres = len(residues)
        if nres > 1:
            l = len(residues[0]) # lenght of first residue
            for ires in range(1,nres):
                if len(residues[ires]) != l:
                    print 'Residues of different lenght!'
                    print '#############################'
                    print residues[0]
                    print 'has ',len(residues[0]),'residues'
                    print '#############################'
                    print residues[ires]
                    print 'has ',len(residues[ires]),'residues'
                    raise
    print 'Residue lenghts are harmonious'
    nres=0
    for resname, residues in all.items():
        print resname, [len(residue) for residue in residues]
        nres += len(residues)
    print 'Number of residues in file:', nres

# Translation from mutated.1 to final_snapshot_rubr1 for atoms with different names
#rosseta={
#    'CYM':{ 'D':'HN', 'DA':'HA', 'DB2':'HB2', 'DB3':'HB3'},
#    'ILE':{ 'D':'H', 'DA':'HA', 'DB':'HB',}
#    'ALA':{ 'D':'H', 'DA':'HA', 'DB1':'HB1', 'DB2':'HB2', 'DB3':'HB3'},
#}

def translate_deuterium_to_hydrogen(inX):
    '''Change "D" to "H" for all deuterium atoms. Deuterium atoms begin with "D" in file mutated.1.pdb'''
    buf=''
    for line in open(inX).readlines():
        if line[0:4] == 'ATOM':
            atom_name = line[12:16]
            if atom_name.strip()[0] == 'D':
                atom_name = atom_name.replace('D', 'H', 1)
                line = line[:12] + atom_name + line[16:]
        buf += line
    return buf

parser = argparse.ArgumentParser(description="""Translate deuterium to hydrogen.
Example: python mutated.1_to_mutated2.py mutated.1.pdb mutated.2.pdb""")
parser.add_argument('inX', type=str, help='input mutated.1 file')
parser.add_argument('outX', type=str, help='output mutated.2 file')
args=parser.parse_args()

check_residues_same_length(args.inX)
buf = translate_deuterium_to_hydrogen(args.inX)
open(args.outX, 'w').write(buf)
