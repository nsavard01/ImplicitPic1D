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

Exp_test_fortran = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_test_fortran/')

test_cpp = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/test_Exp_cpp/')