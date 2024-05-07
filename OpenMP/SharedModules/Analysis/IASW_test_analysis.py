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


NGP_IASW_Chacon2013_Curv_Smooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_curv_Smooth/')
NGP_IASW_Chacon2013_Curv_noSmooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_curv_noSmooth/')


CIC_IASW_Chacon2013_Curv_noSmooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_curv_noSmooth/')
CIC_IASW_Chacon2013_Curv_Smooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_curv_Smooth/')
CIC_IASW_Chacon2013_curv_Smooth_2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_curv_Smooth_2delT/')
CIC_IASW_Chacon2013_curv_Smooth_10delT_anal_lowRes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_curv_Smooth_10delT_anal_lowRes/')