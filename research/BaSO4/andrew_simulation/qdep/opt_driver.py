from pdb import set_trace as tr

from camm.simulation.src.beamline.assemblemodel import modelB_freeE_C
from camm.simulation.src.molmec.ffupdate.ff_update import getParams
expdir='/projects/research/BaSO4/expdata/qdep'
simdir='/projects/research/BaSO4/andrew_simulation/qdep'

tr()
parms,model=getParams('params.in')
model='b0=1.3211; b1=0.00; e0.0=0.99; e0.1=0.99; e0.2=0.99; e0.3=0.99; e0.4=0.99; c0=2.3; eshift=0.0' #from params.in
resolution='%s/resolution.nxs'%expdir
convolved='%s/res_x_sqt_inc_T230_water.nxs'%simdir
assembled='%s/assembled_inc_T230_water.nxs'%simdir
modelB_freeE_C(model, resolution, convolved, assembled, expdata=None, costfile='results.out', derivdata=None, derivexclude=[], doshift=True)
