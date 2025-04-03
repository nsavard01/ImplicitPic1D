# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""
import matplotlib.pyplot as plt
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


Exp_IASW_Chacon2013_4Threads = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_4Threads/')
Exp_IASW_Chacon2013_4Threads_200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_4Threads_200PPC/')
Exp_IASW_Chacon2013_4Threads_100PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_4Threads_100PPC/')
# Exp_IASW_Chacon2013_8Threads_100PPC_0p2delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_100PPC_0p2delT/')
# Exp_IASW_Chacon2013_8Threads_50PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_50PPC/')
# Exp_IASW_Chacon2013_8Threads_50PPC_400Cells_0p2delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_50PPC_400Cells_0p2delT/')
Exp_IASW_Chacon2013_4Threads_40PPC_400Cells_0p2delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_4Threads_40PPC_400Cells_0p2delT/')

Exp_IAW_Chen2011 = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_IAW_Chen2011/')

NGP_IASW_Chacon2013_Curv_NoSmooth_32Threads_2delT_try = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_NoSmooth_32Threads_2delT_try/')
NGP_density_uniform_sin = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_NoSmooth_32Threads_2delT_try', 'ion', 20)

NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT/')
NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1/')
NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_1000PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_1000PPC/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_100PPC/')
# NGP_IASW_Chacon2013_Curv_nonSmooth_8Threads_2delT_lowerRes_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_nonSmooth_8Threads_2delT_lowerRes_100PPC/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_50PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_50PPC/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_Res = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_Res/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_JFNK = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_JFNK/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_5delT_JFNK = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_5delT_JFNK/')
# # NGP_IASW_Chacon2013_Curv_noSmooth = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_curv_noSmooth/')
#

NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT/')
NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_0p1delT_epsneg8 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_0p1delT_epsneg8/')

CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_4delT_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_4delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_3delT_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_3delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_1delT_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_1delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p5delT_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p5delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p25delT_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p25delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_250PPC_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_250PPC_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3/')


CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT/')

Exp_Norm_density = getAveDensityFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads', 'ion', 20)
NGP_Norm_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_32Threads_2delT', 'ion', 20)
CIC_Norm_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT', 'ion', 20)
CIC_Norm_noSmooth_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_NoSmooth_32Threads_2delT', 'ion', 20)
Exp_100PPC_density = getAveDensityFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_100PPC', 'ion', 20)
Exp_200PPC_density = getAveDensityFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_200PPC', 'ion', 20)
Exp_lowRes_density = getAveDensityFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_30PPC_400Cells_0p2delT', 'ion', 20)
CIC_250PPC_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_250PPC_epsNeg3', 'ion', 20)
CIC_100PPC_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_100PPC_epsNeg3', 'ion', 20)
CIC_30PPC_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_30PPC_epsNeg3', 'ion', 20)
CIC_40PPC_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_40PPC_epsNeg3', 'ion', 20)
CIC_40PPC_noSmooth_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_NoSmooth_32Threads_2delT_40PPC_epsNeg3', 'ion', 20)
Exp_lowResOther_density = getAveDensityFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_40PPC_400Cells_0p2delT', 'ion', 20)
NGP_HighRes_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT', 'ion', 20)
NGP_HigherRes_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_0p1delT_epsneg8', 'ion', 20)
CIC_HighRes_density = getAveDensityFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT', 'ion', 20)

# Exp_Norm_phi = getAvePhiFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads', 20)
# NGP_Norm_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_32Threads_2delT', 20)
# CIC_Norm_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT', 20)
# CIC_Norm_noSmooth_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_NoSmooth_32Threads_2delT', 20)
# Exp_100PPC_phi = getAvePhiFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_100PPC', 20)
# Exp_200PPC_phi = getAvePhiFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_200PPC', 20)
# Exp_lowRes_phi = getAvePhiFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_30PPC_400Cells_0p2delT',  20)
# CIC_250PPC_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_250PPC_epsNeg3', 20)
# CIC_100PPC_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_100PPC_epsNeg3', 20)
# CIC_30PPC_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_30PPC_epsNeg3',  20)
# CIC_40PPC_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_40PPC_epsNeg3',  20)
# CIC_40PPC_noSmooth_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_NoSmooth_32Threads_2delT_40PPC_epsNeg3', 20)
# Exp_lowResOther_phi = getAvePhiFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_40PPC_400Cells_0p2delT', 20)
# NGP_HighRes_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT', 20)
# NGP_HigherRes_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_0p1delT_epsneg8', 20)
# CIC_HighRes_phi = getAvePhiFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT', 20)

Exp_Norm_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads', 20)
NGP_HighRes_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT', 20)
CIC_HighRes_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT', 20)
CIC_Norm_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT', 20)
Exp_100PPC_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_100PPC', 20)
Exp_200PPC_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_200PPC', 20)
Exp_lowResOther_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_32Threads_40PPC_400Cells_0p2delT', 20)
CIC_40PPC_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_40PPC_epsNeg3',  20)
CIC_250PPC_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_250PPC_epsNeg3', 20)
CIC_100PPC_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_32Threads_2delT_100PPC_epsNeg3', 20)
CIC_40PPC_noSmooth_EField = getAveEFieldFiles('Y:/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_NoSmooth_32Threads_2delT_40PPC_epsNeg3', 20)


plt.plot(Exp_IASW_Chacon2013_4Threads.grid*1e2, Exp_Norm_density, '-', label = 'Expl.')
plt.plot(NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid*1e2, NGP_HighRes_density, linestyle = '--', label = r'INGP')
plt.plot(CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid*1e2, CIC_HighRes_density, linestyle = '--', label = r'ICIC')
# plt.plot(NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid * 1e2, NGP_HighRes_density, linestyle = '--', marker = '.', label = 'NGP High Res')
# plt.plot(CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid * 1e2, CIC_HighRes_density, linestyle = '--', marker = '.', label = 'CIC High Res')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Ion Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best',fontsize=12)
plt.tight_layout()
plt.savefig('IASW/IASW_Exp_Imp.pdf')
plt.savefig('IASW/IASW_Exp_Imp.png')
plt.close()



CIC_times = [CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p25delT_epsNeg3.totPotTime, \
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p5delT_epsNeg3.totPotTime,\
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_1delT_epsNeg3.totPotTime,\
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.totPotTime, \
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_3delT_epsNeg3.totPotTime, \
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_4delT_epsNeg3.totPotTime]
CIC_timeSteps = [0.25, 0.5, 1, 2, 3, 4]

plt.plot(CIC_timeSteps, CIC_times, linestyle = '--', marker = 'o')
plt.xlabel('Time Step $\omega_{pe}^{-1}$', fontsize = 14)
plt.ylabel(r'Total Potential Solver Time (s)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('IASW/CIC_time_vs_timeStep.pdf')
plt.savefig('IASW/CIC_time_vs_timeStep.png')
plt.close()

plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.grid*1e2, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$\varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.grid*1e2, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$\varepsilon_a = 10^{-3}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.grid*1e2, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$\varepsilon_a = 10^{-1}$')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Ion Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('IASW/IASW_CIC_diffEps.pdf')
plt.savefig('IASW/IASW_CIC_diffEps.png')
plt.close()

plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['time(s)'].values, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['gaussError'].values, linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['time(s)'].values, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['gaussError'].values, linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-3}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['time(s)'].values, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['gaussError'].values, linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-1}$')
plt.xlabel('Time Stamp (s)', fontsize = 14)
plt.ylabel(r'Error Gauss Law (a.u.)', fontsize = 14)
plt.yscale('log')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('IASW/IASW_CIC_GaussError_Res.pdf')
plt.savefig('IASW/IASW_CIC_GaussError_Res.png')
plt.close()

plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['time(s)'].values, np.abs(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['TotalEnergy(J/m^2)'].values - CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['TotalEnergy(J/m^2)'].values[0])/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['TotalEnergy(J/m^2)'].values[0], linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['time(s)'].values, np.abs(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['TotalEnergy(J/m^2)'].values - CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['TotalEnergy(J/m^2)'].values[0])/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['TotalEnergy(J/m^2)'].values[0], linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-3}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['time(s)'].values, np.abs(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['TotalEnergy(J/m^2)'].values-CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['TotalEnergy(J/m^2)'].values[0])/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['TotalEnergy(J/m^2)'].values[0], linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-1}$')
plt.xlabel('Time Stamp (s)', fontsize = 14)
plt.ylabel(r'Relative Energy Error', fontsize = 14)
plt.yscale('log')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('IASW/IASW_CIC_EnergyError_Res.pdf')
plt.savefig('IASW/IASW_CIC_EnergyError_Res.png')
plt.close()


plt.plot(Exp_IASW_Chacon2013_4Threads_100PPC.grid*1e2, Exp_100PPC_density, '--', label = r'100 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_IASW_Chacon2013_4Threads_200PPC.grid*1e2, Exp_200PPC_density, '--', label = r'200 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_IASW_Chacon2013_4Threads.grid*1e2, Exp_Norm_density, '--', label = r'2000 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_IASW_Chacon2013_4Threads_40PPC_400Cells_0p2delT.grid*1e2, Exp_lowResOther_density, '--', label = r'40 PPC, 400 Cells, $\frac{0.2}{\omega_{pe}}$')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Ion Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'center right',fontsize=11)
plt.tight_layout()
plt.savefig('IASW/IASW_Exp_ChangePPC_delX.pdf')
plt.savefig('IASW/IASW_Exp_ChangePPC_delX.png')
plt.close()



plt.plot(CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid*1e2, CIC_HighRes_density, linestyle = '-', label = r'2000 PPC, 512 Cells')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.grid*1e2, CIC_Norm_density, linestyle = '--', marker = '.', label = r'2000 PPC')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_250PPC_epsNeg3.grid*1e2, CIC_250PPC_density, linestyle = '--', marker = '.', label = r'250 PPC')
# plt.plot(NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3.grid, NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'NGP,188 PPC')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.grid*1e2, CIC_100PPC_density, linestyle = '--', marker = '.', label = r'100 PPC')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.grid*1e2, CIC_40PPC_density, linestyle = '--', marker = '.', label = r'40 PPC')
# plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.grid*1e2, CIC_40PPC_noSmooth_density, linestyle = '--', marker = '.', label = r'40 PPC no smooth')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Ion Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best',fontsize=12)
plt.tight_layout()
plt.savefig('IASW/IASW_CIC_ChangePPC.pdf')
plt.savefig('IASW/IASW_CIC_ChangePPC.png')
plt.close()



plt.plot(CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid*1e2, CIC_HighRes_density, linestyle = '--', marker = '.', label = r'2000 PPC, 512 Cells')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.grid*1e2, CIC_40PPC_noSmooth_density, linestyle = '--', marker = '.', label = r'40 PPC')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Ion Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best',fontsize=12)
plt.tight_layout()
plt.savefig('IASW/IASW_CIC_40PPC_smooth.pdf')
plt.savefig('IASW/IASW_CIC_40PPC_smooth.png')
plt.close()

plt.plot(CIC_HighRes_EField[0]*1e2, CIC_HighRes_EField[1], linestyle = '--', marker = '.', label = r'2000 PPC, 512 Cells')
plt.plot(CIC_40PPC_noSmooth_EField[0]*1e2, CIC_40PPC_noSmooth_EField[1], linestyle = '--', marker = '.', label = r'40 PPC')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Field (V/m)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best',fontsize=12)
plt.tight_layout()
plt.savefig('IASW/IASW_CIC_40PPC_smooth_EField.pdf')
plt.savefig('IASW/IASW_CIC_40PPC_smooth_EField.png')
plt.close()

plt.figure()
xi = np.arange(64) + 1
plt.plot(xi, np.diff(NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.grid) * 1e2, 'o')
plt.xlim(1, 64)
plt.xlabel(r'Cell Number', fontsize = 14)
plt.ylabel(r'Cell Size (cm)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('IASW/IASW_cell_map.pdf')
plt.savefig('IASW/IASW_cell_map.png')
plt.close()


y = CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['TotalMomentum(kg/m/s)'].values[0]
z = CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.globDiag['TotalMomentum(kg/m/s)'].values[0]
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['time(s)'].values, (CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['TotalMomentum(kg/m/s)'].values - y)/y, linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.globDiag['time(s)'].values, (CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.globDiag['TotalMomentum(kg/m/s)'].values-z)/z, linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-3}$')


dev_sqr_exp = 0.0
for i in range(Exp_IASW_Chacon2013_4Threads_40PPC_400Cells_0p2delT.Nx-1):
    dev_sqr_exp += (Exp_lowResOther_density[i] - np.interp(Exp_IASW_Chacon2013_4Threads_40PPC_400Cells_0p2delT.grid[i], Exp_IASW_Chacon2013_4Threads.grid, Exp_Norm_density))**2
print('exp 40 PPC', 1e-13 * np.sqrt(dev_sqr_exp/(Exp_IASW_Chacon2013_4Threads_40PPC_400Cells_0p2delT.Nx-1)))

dev_sqr_exp = 0.0
for i in range(Exp_IASW_Chacon2013_4Threads_100PPC.Nx-1):
    dev_sqr_exp += (Exp_100PPC_density[i] - Exp_Norm_density[i])**2
print('exp 100 PPC', 1e-13 * np.sqrt(dev_sqr_exp/(Exp_IASW_Chacon2013_4Threads_100PPC.Nx-1)))

dev_sqr_exp = 0.0
for i in range(Exp_IASW_Chacon2013_4Threads_200PPC.Nx-1):
    dev_sqr_exp += (Exp_200PPC_density[i] - Exp_Norm_density[i])**2
print('exp 200 PPC', 1e-13 * np.sqrt(dev_sqr_exp/(Exp_IASW_Chacon2013_4Threads_200PPC.Nx-1)))



dev_sqr_impl = 0.0
for i in range(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.Nx):
    dev_sqr_impl += (CIC_40PPC_density[i] - np.interp(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.grid[i], CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid, CIC_HighRes_density))**2
print('CIC 40 PPC', 1e-13 * np.sqrt(dev_sqr_impl/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.Nx))

dev_sqr_impl = 0.0
for i in range(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.Nx):
    dev_sqr_impl += (CIC_100PPC_density[i] - np.interp(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.grid[i], CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid, CIC_HighRes_density))**2
print('CIC 100 PPC', 1e-13 * np.sqrt(dev_sqr_impl/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.Nx))

dev_sqr_impl = 0.0
for i in range(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.Nx):
    dev_sqr_impl += (CIC_250PPC_density[i] - np.interp(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.grid[i], CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid, CIC_HighRes_density))**2
print('CIC 250 PPC', 1e-13 * np.sqrt(dev_sqr_impl/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.Nx))

dev_sqr_impl = 0.0
for i in range(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.Nx):
    dev_sqr_impl += (CIC_Norm_density[i] - np.interp(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.grid[i], CIC_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid, CIC_HighRes_density))**2
print('CIC 2000 PPC', 1e-13 * np.sqrt(dev_sqr_impl/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.Nx))


dev_sqr_impl = 0.0
for i in range(NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.Nx-1):
    dev_sqr_impl += (NGP_Norm_density[i] - np.interp(NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.grid[i], NGP_IASW_Chacon2013_512cells_NoSmooth_32Threads_2delT.grid, NGP_HighRes_density))**2
print('NGP 2000 PPC', 1e-13 * np.sqrt(dev_sqr_impl/(NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.Nx-1)))

plt.figure()
plt.plot(Exp_Norm_EField[0]*1e2, Exp_Norm_EField[1], linestyle = '--', marker = '.', label = r'Expl.')
plt.plot(NGP_HighRes_EField[0]*1e2, NGP_HighRes_EField[1], linestyle = '--', marker = '.', label = r'INGP')
plt.plot(CIC_HighRes_EField[0]*1e2, CIC_HighRes_EField[1], linestyle = '--', marker = '.', label = r'ICIC')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Field (V/m)', fontsize = 14)
plt.xlim(3.72 ,11.15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best',fontsize=12)
plt.tight_layout()
plt.savefig('IASW/IASW_Exp_Imp_EField.pdf')
plt.savefig('IASW/IASW_Exp_Imp_EField.png')
plt.close()


plt.figure()
plt.plot(Exp_100PPC_EField[0]*1e2, Exp_100PPC_EField[1], '--', marker = '.', label = r'100 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_200PPC_EField[0]*1e2, Exp_200PPC_EField[1], '--', marker = '.', label = r'200 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_Norm_EField[0]*1e2, Exp_Norm_EField[1], linestyle = '-', label = r'2000 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_lowResOther_EField[0]*1e2, Exp_lowResOther_EField [1], linestyle = '--', marker = '.', label = r'40 PPC, 400 Cells, $\frac{0.2}{\omega_{pe}}$')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Field (V/m)', fontsize = 14)
plt.xlim(3.72, 11.15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best',fontsize=11)
plt.tight_layout()
plt.savefig('IASW/IASW_Exp_ChangePPC_delX_EField.pdf')
plt.savefig('IASW/IASW_Exp_ChangePPC_delX_EField.png')
plt.close()

plt.figure()
plt.plot(CIC_HighRes_EField[0]*1e2, CIC_HighRes_EField[1], linestyle = '-', label = r'2000 PPC, 512 Cells')
plt.plot(CIC_Norm_EField[0]*1e2, CIC_Norm_EField[1], linestyle = '--', marker = '.', label = r'2000 PPC')
plt.plot(CIC_250PPC_EField[0]*1e2, CIC_250PPC_EField[1], linestyle = '--', marker = '.', label = r'250 PPC')
# plt.plot(NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3.grid, NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'NGP,188 PPC')
plt.plot(CIC_100PPC_EField[0]*1e2, CIC_100PPC_EField[1], linestyle = '--', marker = '.', label = r'100 PPC')
plt.plot(CIC_40PPC_EField[0]*1e2, CIC_40PPC_EField[1], linestyle = '--', marker = '.', label = r'40 PPC')
# plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_40PPC_epsNeg3.grid*1e2, CIC_40PPC_noSmooth_density, linestyle = '--', marker = '.', label = r'40 PPC no smooth')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Field (V/m)', fontsize = 14)
plt.xlim(3.72, 11.15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best',fontsize=12)
plt.tight_layout()
plt.savefig('IASW/IASW_CIC_ChangePPC_EField.pdf')
plt.savefig('IASW/IASW_CIC_ChangePPC_EField.png')
plt.close()
