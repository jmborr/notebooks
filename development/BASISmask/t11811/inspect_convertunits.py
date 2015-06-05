def distance_to_sample(workspace, detectorID):
    i=workspace.getInstrument()
    s=i.getSample()
    d=w.getDetector(detectorID)
    return d.getDistance(s)

w=mtd['auto_ws']
wn=mtd['auto_wsnew']
N=w.getNumberHistograms()

def getTwoTheta(workspace, detectorID):
    i=workspace.getInstrument()
    s=i.getSample()
    d=w.getDetector(detectorID)
    return d.getTwoTheta( s.getPos(), V3D(0,0,1) )
    
delta_dID=23
dID=0
while dID < N:
    #print distance_to_sample(w, dID) - distance_to_sample(wn, dID),
    #print getTwoTheta(w, dID) - getTwoTheta(wn, dID), 
    #print dID, w.getDetector(dID).getPos() - wn.getDetector(dID).getPos() #physical positions, not neutronic
    print dID, w.getDetector(dID).getNumberParameter('Efixed'), wn.getDetector(dID).getNumberParameter('Efixed') 
    dID += delta_dID

print N, wn.getNumberHistograms()