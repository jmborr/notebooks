import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Remove negative times from I(Q,t).
    Example: nonnegative.py fqt_head.dat fqtNN_head.dat""")
    parser.add_argument('infile', help='input fqt file, in ASCII format')
    parser.add_argument('outfile', help='output file')
    args=parser.parse_args()

    buf=''
    for line in open(args.infile,'r').readlines():
        if line[0]!='-':
            buf+=line
    open(args.outfile,'w').write(buf)
