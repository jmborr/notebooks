import os
from pdb import set_trace as tr

indexes='''11 15 12 9
10 19 15 13
14 20 16 16
4  10 8  8
3  9  7  6
2  8  6  4
1  7  5  2
9  16 13 11'''

template='''trajin dms8t8_rms2first.dcd
dihedral d1 @Si_IXSi_ @O_IXC1_ @Si_IXC2_ @C_IXC3_ out _DIHEDRALFILE_ time 1.0 range360
'''                    

#tr()
os.system('/bin/mkdir -p /SNSlocal/projects/jbq/POSS/dms8t8/cpptraj/dihedrals/Si-O-Si-C')
for row in indexes.split('\n'):
    ixs=row.split()
    dihedral='Si{0}-O{1}-Si{2}-C{3}'.format(*ixs)
    print dihedral
    dihedralfile='dihedrals/Si-O-Si-C/{0}.dat'.format(dihedral)
    buf=template
    buf=buf.replace('_DIHEDRALFILE_',dihedralfile)
    buf=buf.replace('_IXSi_',ixs[0])
    buf=buf.replace('_IXC1_',ixs[1])
    buf=buf.replace('_IXC2_',ixs[2])
    buf=buf.replace('_IXC3_',ixs[3])
    open('/SNSlocal/projects/jbq/POSS/dms8t8/cpptraj/dihedrals/Si-O-Si-C/{0}.ptraj'.format(dihedral),'w').write(buf)
    
