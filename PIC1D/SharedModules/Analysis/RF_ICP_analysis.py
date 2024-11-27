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

Exp_ICP_test_highDensity_implicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_implicit/')
Exp_ICP_test_highDensity_fullyImplicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_fullyImplicit/')