names=('poly1', 'poly2', 'polyCharge', 'head', 'tail', 'Ncions', 'Pcions')
ljdir1s=("Low_surf_conc", "Philic_BCP2", "Phobic_BCP2", "Strong_Phobic_BCP2")
ljdir2s=("lj10_prod", "lj12_prod", "lj13_prod", "lj11_prod")
rootd='/SNSlocal/projects/jbq/PEblock'
N=len(ljdir1s)
for i in range(N):
    for name in names:
        inname='/'.join([ rootd, ljdir1s[i], ljdir2s[i], 'fqt_inc_{0}.h5'.format(name) ])
        LoadSassena(Filename=inname, TimeUnit=50000, SortByQvectors=1, OutputWorkspace='iqt')
        alpha= 1.0/mtd['iqt_fq0'].dataY(0)[0]
        Scale(InputWorkspace='iqt_fqt.Re', Factor=alpha, Operation='Multiply', OutputWorkspace='iqt_fqt.Re')
        outname= '/'.join([ rootd, ljdir1s[i], ljdir2s[i], 'fqt_inc_{0}.dat'.format(name) ])
        SaveAscii(InputWorkspace='iqt_fqt.Re', Filename=outname, Separator='Space', CommentIndicator='#', Version=1)
        DeleteWorkspace('iqt')
        print outname
