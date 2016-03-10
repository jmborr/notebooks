import os
#NOTE: replace this directory with the directory where you
#      saved the fit_results_example.txt
saveDir="/home/jbq/Downloads"

Load(Filename=os.path.join(saveDir,'fit_results_example.txt'), OutputWorkspace='poinData')
ConvertToHistogram(InputWorkspace='poinData', OutputWorkspace='histogramData')
""" 1. Open Interfaces -> DynamicPDF -> tdst the displayCurveFit widget
    2. Load any of the pointData and histogramData workspaces.
       The first point should be shown at X=1.05, corresponding to the first
       center bin of histogramData and the first point in pointData.
"""
