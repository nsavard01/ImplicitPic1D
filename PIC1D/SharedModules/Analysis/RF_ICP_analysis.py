# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *
print('Loaded modules')
    
#%%

print('Loading data')



Exp_ICP_test_lowDensity = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_lowDensity/')
Exp_ICP_test_lowDensity_implicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_lowDensity_implicit/')
Exp_ICP_test_lowDensity_fullyImplicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_lowDensity_fullyImplicit/')

Exp_ICP_test_highDensity = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_implicit/')
Exp_ICP_test_highDensity_regRF = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_regRF/')
Exp_ICP_test_highDensity_regRF_6400delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_regRF_6400delT/')

Exp_ICP_test_highDensity_6400delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT/')
Exp_ICP_test_highDensity_6400delT_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_doublePart/')
Exp_ICP_test_highDensity_6400delT_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_quadPart/')
Exp_ICP_test_highDensity_6400delT_x8Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_x8Part/')
Exp_ICP_test_highDensity_6400delT_x16Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_x16Part/')
Exp_ICP_test_highDensity_6400delT_x32Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_x32Part/')

Exp_ICP_test_highDensity_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_doublePart/')
Exp_ICP_test_highDensity_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_quadPart/')

Exp_ICP_test_highDensity_lowPressure_6400delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT/')
Exp_ICP_test_highDensity_lowPressure_12800delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_12800delT/')
Exp_ICP_test_highDensity_lowPressure_6400delT_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_doublePart/')
Exp_ICP_test_highDensity_lowPressure_6400delT_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_quadPart/')
Exp_ICP_test_highDensity_lowPressure_6400delT_x8Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_x8Part/')

Exp_ICP_test_highDensity_lowPressure_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_quadPart/')
Exp_ICP_test_highDensity_lowPressure_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_doublePart/')
Exp_ICP_test_highDensity_lowPressure = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure/')


NGP_ICP_test_lowDensity_fullyImplicit = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_test_lowDensity_fullyImplicit/')