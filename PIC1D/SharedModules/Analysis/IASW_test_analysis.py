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


Exp_IASW_Chacon2013_4Threads = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_4Threads/')
Exp_IASW_Chacon2013_4Threads_200PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_4Threads_200PPC/')
Exp_IASW_Chacon2013_4Threads_100PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_4Threads_100PPC/')
# Exp_IASW_Chacon2013_8Threads_100PPC_0p2delT = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_100PPC_0p2delT/')
# Exp_IASW_Chacon2013_8Threads_50PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_50PPC/')
# Exp_IASW_Chacon2013_8Threads_50PPC_400Cells_0p2delT = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_8Threads_50PPC_400Cells_0p2delT/')
Exp_IASW_Chacon2013_4Threads_30PPC_400Cells_0p2delT = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IASW_Chacon2013_4Threads_30PPC_400Cells_0p2delT/')

Exp_IAW_Chen2011 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_IAW_Chen2011/')

NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT/')
NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1/')
NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_1000PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_1000PPC/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_100PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_100PPC/')
# NGP_IASW_Chacon2013_Curv_nonSmooth_8Threads_2delT_lowerRes_100PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_nonSmooth_8Threads_2delT_lowerRes_100PPC/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_50PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_lowerRes_50PPC/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_Res = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_Res/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_JFNK = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_2delT_JFNK/')
# NGP_IASW_Chacon2013_Curv_Smooth_8Threads_5delT_JFNK = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chacon2013_Curv_Smooth_8Threads_5delT_JFNK/')
# # NGP_IASW_Chacon2013_Curv_noSmooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_IASW_Chen2013_curv_noSmooth/')
#

CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_1e18/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_4delT_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_4delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_3delT_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_3delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_1delT_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_1delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p5delT_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p5delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p25delT_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p25delT_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3/')
CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_30PPC_epsNeg3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_30PPC_epsNeg3/')

plt.plot(Exp_IASW_Chacon2013_4Threads.grid, Exp_IASW_Chacon2013_4Threads.getDensity('ion',20), '-', label = 'Exp.')
plt.plot(NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.grid, NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.getDensity('ion',20), linestyle = '--', marker = '.', label = r'NGP')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.grid, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.getDensity('ion',20), linestyle = '--', marker = '.', label = r'CIC')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1])
plt.xlabel('Distance (m)')
plt.ylabel(r'Ion Density (m$^{-3}$)')
plt.legend(loc = 'best')
plt.savefig('IASW/IASW_Exp_Imp.pdf')
plt.close()

CIC_times = [CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p25delT_epsNeg3.totPotTime, \
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_0p5delT_epsNeg3.totPotTime,\
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_1delT_epsNeg3.totPotTime,\
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.totPotTime, \
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_3delT_epsNeg3.totPotTime, \
             CIC_IASW_Chacon2013_Curv_Smooth_4Threads_4delT_epsNeg3.totPotTime]
CIC_timeSteps = [0.25, 0.5, 1, 2, 3, 4]

plt.plot(CIC_timeSteps, CIC_times, linestyle = '--', marker = 'o')
plt.xlabel('Time Step $\omega_{pe}^{-1}$')
plt.ylabel(r'Total Potential Solver Time (s)')
plt.savefig('IASW/CIC_time_vs_timeStep.pdf')
plt.close()

plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.grid, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$\varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.grid, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$\varepsilon_a = 10^{-3}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.grid, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.getDensity('ion',20), linestyle = '--', marker = '.', label = r'$\varepsilon_a = 10^{-1}$')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1])
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'Ion Density (m$^{-3}$)', fontsize = 14)
plt.legend(loc = 'best', fontsize = 12)
plt.savefig('IASW/IASW_CIC_diffEps.pdf')
plt.close()

plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['time(s)'].values, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['gaussError'].values, linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['time(s)'].values, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['gaussError'].values, linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-3}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['time(s)'].values, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['gaussError'].values, linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-1}$')
plt.xlabel('Time Stamp (s)', fontsize = 14)
plt.ylabel(r'Error Gauss Law (a.u.)', fontsize = 14)
plt.yscale('log')
plt.legend(loc = 'best', fontsize = 12)
plt.savefig('IASW/IASW_CIC_GaussError_Res.pdf')
plt.close()

plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['time(s)'].values, np.abs(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['TotalEnergy(J/m^2)'].values - CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['TotalEnergy(J/m^2)'].values[0])/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT.globDiag['TotalEnergy(J/m^2)'].values[0], linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-6}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['time(s)'].values, np.abs(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['TotalEnergy(J/m^2)'].values - CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['TotalEnergy(J/m^2)'].values[0])/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.globDiag['TotalEnergy(J/m^2)'].values[0], linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-3}$')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['time(s)'].values, np.abs(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['TotalEnergy(J/m^2)'].values-CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['TotalEnergy(J/m^2)'].values[0])/CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg1.globDiag['TotalEnergy(J/m^2)'].values[0], linestyle = '--', marker = 'o', label = r'$\varepsilon_a = 10^{-1}$')
plt.xlabel('Time Stamp (s)', fontsize = 14)
plt.ylabel(r'Relative Energy Error', fontsize = 14)
plt.yscale('log')
plt.legend(loc = 'best', fontsize = 12)
plt.savefig('IASW/IASW_CIC_EnergyError_Res.pdf')
plt.close()


plt.plot(Exp_IASW_Chacon2013_4Threads_100PPC.grid, Exp_IASW_Chacon2013_4Threads_100PPC.getDensity('ion',20), '--', label = r'100 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_IASW_Chacon2013_4Threads_200PPC.grid, Exp_IASW_Chacon2013_4Threads_200PPC.getDensity('ion',20), '--', label = r'200 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_IASW_Chacon2013_4Threads.grid, Exp_IASW_Chacon2013_4Threads.getDensity('ion',20), '--', label = r'2000 PPC, 512 Cells, $\frac{0.1}{\omega_{pe}}$')
plt.plot(Exp_IASW_Chacon2013_4Threads_30PPC_400Cells_0p2delT.grid, Exp_IASW_Chacon2013_4Threads_30PPC_400Cells_0p2delT.getDensity('ion',20), '--', label = r'30 PPC, 400 Cells, $\frac{0.2}{\omega_{pe}}$')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1])
plt.xlabel('Distance (m)')
plt.ylabel(r'Ion Density (m$^{-3}$)')
plt.legend(loc = 'best')
plt.savefig('IASW/IASW_Exp_ChangePPC_delX.pdf')
plt.close()


plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.grid, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'2000 PPC')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3.grid, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'187.5 PPC')
# plt.plot(NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3.grid, NGP_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_188PPC_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'NGP,188 PPC')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.grid, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_100PPC_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'100 PPC')
plt.plot(CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_30PPC_epsNeg3.grid, CIC_IASW_Chacon2013_Curv_Smooth_4Threads_2delT_30PPC_epsNeg3.getDensity('ion',20), linestyle = '--', marker = '.', label = r'30 PPC')
plt.xlim(Exp_IASW_Chacon2013_4Threads.grid[0], Exp_IASW_Chacon2013_4Threads.grid[-1])
plt.xlabel('Distance (m)')
plt.ylabel(r'Ion Density (m$^{-3}$)')
plt.legend(loc = 'best')
plt.savefig('IASW/IASW_CIC_ChangePPC.pdf')
plt.close()