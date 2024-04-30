# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *
import scipy.integrate as integ
print('Loaded modules')
    
#%%

print('Loading data')




def Z_func(x):
    return complex(0,1) * np.sqrt(np.pi) * np.exp(-x**2) * (1 + scipy.special.erf(complex(0,1) * x))

def innerFunc(x, A):
    return np.exp(-x**2) / (x - A)

def Z_small(psi):
    longTerm = 2 * psi * (1 - (2/3)*(psi**2) + (4/15)*(psi**4) - (8/105) * (psi**6))
    return complex(0,1) * np.sqrt(np.pi) * np.exp(-psi**2) - longTerm

def Z_high(psi):
    sigma = 0
    if (psi.imag > 1/abs(psi.real)):
        sigma = 0
    elif (abs(psi.imag) > 1/abs(psi.real)):
        sigma = 1
    elif (psi.imag < -1/abs(psi.real)):
        sigma = 2
    longTerm = 1 + 0.5/(psi**2) + 0.75 / (psi**4) + (15/8)/(psi**6)
    return complex(0,1) * np.sqrt(np.pi) * sigma * np.exp(-psi**2) - (1/psi) * longTerm

def root(omega, k):
    a = 1 + (1/k**2) * (1 + Z_func(omega/k/np.sqrt(2)) * omega/k/np.sqrt(2))
    return a


Exp_IAW_Chacon = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IAW_Chacon/')
Exp_IASW_Chacon = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon/')
Exp_IASW_Chacon2013 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013/')
Exp_IASW_Chacon2013_10thPart = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_10thPart/')
Exp_IASW_Chacon2013_HalfPart = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_HalfPart/')
Exp_IASW_Chacon2013_1024Cells_halfTime = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_1024Cells/')

NGP_IASW_Chacon2013_Curv = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_curv/')
NGP_IASW_Chacon2013_Curv_nonSmooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_curv_nonSmooth/')
NGP_IASW_Chacon2013_Curv_Smooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_curv_Smooth/')
NGP_IASW_Chacon2013_uniform_nonSmooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_uniform_nonSmooth/')
NGP_IASW_Chacon2013_uniform_Smooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_uniform_Smooth/')


CIC_IASW_Chacon = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon/')
CIC_IASW_Chacon2013_Curv = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv/')
CIC_IASW_Chacon2013_Curv_2delT_10thPart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_2delT_10thPart/')
CIC_IASW_Chacon2013_Curv_2delT_JFNK_1e3Res = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_2delT_JFNK_1e3Res/')
CIC_IASW_Chacon2013_Curv_2delT_JFNK_1e1Res = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_2delT_JFNK_1e1Res/')
CIC_IASW_Chacon2013_Curv_5delT_JFNK_1e1Res = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_5delT_JFNK_1e1Res/')


CIC_IASW_Chacon2013_Curv_1delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_1delT/')
CIC_IASW_Chacon2013_Curv_1delT_neg1Error = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_1delT_neg1Error/')
CIC_IASW_Chacon2013_Curv_1delT_changeMove = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_1delT_changeMove/')
CIC_IASW_Chacon2013_Curv_1delT_neg1Error_changeMove = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_1delT_neg1Error_changeMove/')
CIC_IASW_Chacon2013_Curv_1delT_pos1Error_changeMove = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_1delT_pos1Error_changeMove/')
CIC_IASW_Chacon2013_Curv_1delT_pos1Error = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_1delT_pos1Error/')
CIC_IASW_Chacon2013_Curv_test = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_test/')


CIC_IASW_Chacon_2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon_2delT/')