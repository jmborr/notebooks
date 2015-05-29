import os

rootd='Here write the path where you downloaded the test files'
#rootd='/projects/development/IDA/elwin/t11810/testd'

wsnames=('BASIS_43180_sqw',
    'BASIS_43181_sqw',
    'BASIS_43182_sqw',
    'BASIS_43183_sqw',
    'BASIS_43184_sqw',
    'BASIS_43185_sqw',
    'BASIS_43186_sqw',
    'BASIS_43187_sqw',
    'BASIS_43188_sqw',
    'BASIS_43189_sqw'
    )

for wsname in wsnames:
    Load(Filename=os.path.join(rootd,wsname)+'.nxs', OutputWorkspace=wsname)

GroupWorkspaces(InputWorkspaces=','.join(wsnames), OutputWorkspace='IDA_Elwin_Input', Version=1)

ElasticWindowMultiple(InputWorkspaces='IDA_Elwin_Input',
    Range1Start=-0.0035,
    Range1End=0.0035,
    Range2Start='-0.035',
    Range2End='-0.0315',
    SampleEnvironmentLogName='SetpointLP1',
    OutputInQ='BASIS_43180_elwin_eq',
    OutputInQSquared='BASIS_43180_elwin_eq2',
    OutputELF='BASIS_43180_elwin_elf',
    OutputELT='BASIS_43180_elwin_elt'
    )
MSDFit(InputWorkspace='BASIS_43180_elwin_eq2',
    XStart=0.01,
    Xend = 2.89,
    SpecMin=0,
    SpecMax=9,
    OutputWorkspace='BASIS_43180_elwin_msd'
    )

plotSpectrum('BASIS_43180_elwin_msd_A1', 0, error_bars=True)