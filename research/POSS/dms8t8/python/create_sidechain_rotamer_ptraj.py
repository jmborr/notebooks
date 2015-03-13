import os
from pdb import set_trace as tr

indexes='''11 10 19 15
14 11 15 12
6  3  9  7
3  1  7  5
4  9  16 13
5  2  8  6
8  14 20 16
13 4  10 8'''

template='''trajin dms8t8_rms2first.dcd
dihedral d1 @O_IXO_ @Si_IXSi_ @O_IXC1_ @Si_IXC2_ out _DIHEDRALFILE_ time 1.0 range360
'''                    

#tr()
os.system('/bin/mkdir -p /SNSlocal/projects/jbq/POSS/dms8t8/cpptraj/dihedrals/O-Si-O-Si')
for row in indexes.split('\n'):
    ixs=row.split()
    dihedral='O{0}-Si{1}-O{2}-Si{3}'.format(*ixs)
    dihedralfile='dihedrals/O-Si-O-Si/{0}.dat'.format(dihedral)
    print dihedral
    buf=template
    buf=buf.replace('_DIHEDRALFILE_',dihedralfile)
    buf=buf.replace('_IXO_',ixs[0])
    buf=buf.replace('_IXSi_',ixs[1])
    buf=buf.replace('_IXC1_',ixs[2])
    buf=buf.replace('_IXC2_',ixs[3])
    open('/SNSlocal/projects/jbq/POSS/dms8t8/cpptraj/dihedrals/O-Si-O-Si/{0}.ptraj'.format(dihedral),'w').write(buf)
    
