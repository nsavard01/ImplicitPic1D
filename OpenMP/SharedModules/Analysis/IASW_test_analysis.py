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


Exp_IASW_Chacon2013_32Threads = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads/')
Exp_IASW_Chacon2013_8Threads = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads/')
Exp_IASW_Chacon2013_8Threads_try2 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_try2/')
Exp_IASW_Chacon2013_8Threads_200PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_200PPC/')
Exp_IASW_Chacon2013_8Threads_100PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_100PPC/')
Exp_IASW_Chacon2013_8Threads_100PPC_0p2delT = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_100PPC_0p2delT/')
Exp_IASW_Chacon2013_8Threads_50PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_50PPC/')


NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT/')
NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes/')
NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_1000PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_1000PPC/')
NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_100PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_100PPC/')
NGP_IASW_Chacon2013_Curv_nonSmooth_8Threads_2delT_lowerRes_100PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_nonSmooth_8Threads_2delT_lowerRes_100PPC/')
NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_50PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_50PPC/')
NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_Res = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_Res/')
NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_JFNK = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_JFNK/')
NGP_IASW_Chacon2013_Curv_Smooth_8Threads_5delT_JFNK = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_5delT_JFNK/')
# NGP_IASW_Chacon2013_Curv_noSmooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_curv_noSmooth/')

CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes/')
CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowRes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowRes/')
CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT/')
CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_8000PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_8000PPC/')
CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_100PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_100PPC/')
CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_50PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_50PPC/')
CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_50PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_50PPC/')
CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_100PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_100PPC/')
# CIC_IASW_Chacon2013_Curv_noSmooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_curv_noSmooth/')
# CIC_IASW_Chacon2013_Curv_Smooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_curv_Smooth/')
# CIC_IASW_Chacon2013_curv_Smooth_2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_curv_Smooth_2delT/')
# CIC_IASW_Chacon2013_curv_Smooth_10delT_anal_lowRes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_curv_Smooth_10delT_anal_lowRes/')

plt.interactive(False)
plt.plot(Exp_IASW_Chacon2013_8Threads.grid, Exp_IASW_Chacon2013_8Threads.getDensity('ion',20), '-', label = 'Exp.')
plt.plot(NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.grid, NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.getDensity('ion',20), linestyle = '--', marker = '.', label = r'NGP')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.grid, CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.getDensity('ion',20), linestyle = '--', marker = '.', label = r'CIC')
plt.xlim(Exp_IASW_Chacon2013_8Threads.grid[0], Exp_IASW_Chacon2013_8Threads.grid[-1])
plt.xlabel('Distance (m)')
plt.ylabel(r'Ion Density (m$^{-3}$)')
plt.legend(loc = 'best')
plt.savefig('IASW/IASW_Exp_Imp.pdf')
plt.close()

plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.grid, CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$varepsilon_a = 10^{-6}$, 37.7 s')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowRes.grid, CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowRes.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$varepsilon_a = 10^{-4}$, 27.7 s')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes.grid, CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$varepsilon_a = 10^{-2}$, 18.6 s')
plt.xlim(Exp_IASW_Chacon2013_8Threads.grid[0], Exp_IASW_Chacon2013_8Threads.grid[-1])
plt.xlabel('Distance (m)')
plt.ylabel(r'Ion Density (m$^{-3}$)')
plt.legend(loc = 'best')
plt.savefig('IASW/IASW_CIC_ResTime.pdf')
plt.close()

plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.globDiag['time(s)'].values, CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.globDiag['gaussError'], linestyle = '--', marker = 'o', label = r'$varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.globDiag['time(s)'].values, CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT.globDiag['gaussError'], linestyle = '--', marker = 'o', label = r'$varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes.grid, CIC_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes.getDensity('ion',20), linestyle = '--', marker = 'o', label = r'$varepsilon_a = 10^{-2}$, 18.6 s')
plt.xlim(Exp_IASW_Chacon2013_8Threads.grid[0], Exp_IASW_Chacon2013_8Threads.grid[-1])
plt.xlabel('Distance (m)')
plt.ylabel(r'Ion Density (m$^{-3}$)')
plt.legend(loc = 'best')
plt.savefig('IASW/IASW_CIC_ResTime.pdf')
plt.close()