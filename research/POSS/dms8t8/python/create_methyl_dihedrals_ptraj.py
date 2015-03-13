import os
from pdb import set_trace as tr

indexes='''16 13 11 14
16 13 12 17
9  7  5  1
9  7  6  6
7  5  2  28
7  5  1  25
10 8  7  7
10 8  8  12
20 16 16 22
20 16 15 21
8  6  4  36
8  6  3  31
15 12 10 47
15 12 9  44
19 15 13 38
19 15 14 41'''

template='''trajin dms8t8_rms2first.dcd
dihedral d1 @O_IXH1_ @Si_IXC1_ @C_IXC2_ @H_IXH2_ out _DIHEDRALFILE_ time 1.0 range360
'''                    

#tr()
os.system('/bin/mkdir -p /SNSlocal/projects/jbq/POSS/dms8t8/cpptraj/dihedrals/O-Si-C-H')
for row in indexes.split('\n'):
    ixs=row.split()
    dihedral='O{0}-Si{1}-C{2}-H{3}'.format(*ixs)
    print dihedral
    dihedralfile='dihedrals/O-Si-C-H/{0}.dat'.format(dihedral)
    buf=template
    buf=buf.replace('_DIHEDRALFILE_',dihedralfile)
    buf=buf.replace('_IXH1_',ixs[0])
    buf=buf.replace('_IXC1_',ixs[1])
    buf=buf.replace('_IXC2_',ixs[2])
    buf=buf.replace('_IXH2_',ixs[3])
    open('/SNSlocal/projects/jbq/POSS/dms8t8/cpptraj/dihedrals/O-Si-C-H/{0}.ptraj'.format(dihedral),'w').write(buf)
    
