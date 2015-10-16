import sys
sys.path.append('/home/jbq/repositories/ContactMapAnalysis')
import ContactMapAnalysis as CMA
import argparse
from pdb import set_trace as tr
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='save all residence times as a single column file')
    parser.add_argument('hdf5',help='HDF5 file containing the residence times')
    parser.add_argument('--outfile',help='save to file. Otherwise print to standard output"')
    args = parser.parse_args()

    rtl = CMA.loadResidenceTimesListFromFile(args.hdf5)
    buf=''
    for ires in range(1,rtl.n):
        if rtl.data[ires]:
            L=rtl.data[ires]  #list of isopropanol residence times to this particular residue
            buf+='\n'.join([str(x) for x in L])+'\n'
    if args.outfile:
        open(args.outfile,'w').write(buf)
    else:
        print buf
