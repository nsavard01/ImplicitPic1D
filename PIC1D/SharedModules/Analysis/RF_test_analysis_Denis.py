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
exit()

def readBenchmarkResults():
    firstIndx = 0
    data = np.loadtxt('TurnerBenchmark/turner_benchmark_results.dat')
    case1 = data[0:129,:]
    firstIndx = firstIndx + 129
    case2 = data[firstIndx:firstIndx+257,:]
    firstIndx = firstIndx + 257
    case3 = data[firstIndx:firstIndx+513, :]
    firstIndx = firstIndx+513
    case4 = data[firstIndx:firstIndx+513]
    return case1,case2,case3,case4

def readBenchmarkResultsRefined():
    data = np.loadtxt('TurnerBenchmark/turner_benchmark_refined_results.dat')
    indx = np.where(data[:,0] == 0)[0]
    case1 = data[indx[0]:indx[1], :]
    case2 = data[indx[1]:indx[2], :]
    case3 = data[indx[2]:indx[3], :]
    case4 = data[indx[3]:,:]
    return case1, case2, case3, case4

case1, case2, case3, case4 = readBenchmarkResults()
case1_refined, case2_refined, case3_refined, case4_refined = readBenchmarkResultsRefined()

Exp_RFBenchmark_Denis = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis/')
Exp_RFBenchmark_Denis_doubleTime = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_doubleTime/')
Exp_RFBenchmark_Denis_doubleTime_halfPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_doubleTime_halfPart/')
Exp_RFBenchmark_Denis_doubleTime_halfPart_PCG = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_doubleTime_halfPart_PCG/')
Exp_RFBenchmark_Denis_doubleTime_quartPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_doubleTime_quartPart/')
Exp_RFBenchmark_Denis_tripleTime = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_tripleTime/')
Exp_RFBenchmark_Denis_moreCells = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells/')
Exp_RFBenchmark_Denis_moreCells_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_doublePart/')
Exp_RFBenchmark_Denis_moreCells_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_quadPart/')
Exp_RFBenchmark_Denis_moreCells_40xPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_40xPart/')


Exp_RFBenchmark_Denis_lowerDen_doubleTime_1024Cells_128PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_lowerDen_doubleTime_1024Cells_128PPC/')
Exp_RFBenchmark_Denis_lowerDen_1024Cells_128PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_lowerDen_1024Cells_128PPC/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_halfPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_halfPart/')
Exp_RFBenchmark_Denis_moreCells_lowerDen = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_doublePart/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_quadPart/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_10xPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_10xPart/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_10xPart_halfTime = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_10xPart_halfTime/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_20xPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_20xPart/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_20xPart_3000cells_thirdsTime = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_20xPart_3000cells_thirdsTime/')

# ----------------------- NGP ---------------------------------

NGP_RFBenchmark_Denis_262PPC_128cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_262PPC_128cells_noSmoothing/')
NGP_RFBenchmark_Denis_524PPC_128cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_524PPC_128cells_noSmoothing/')


NGP_RFBenchmark_Denis_131PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_131PPC_128cells_noSmoothing_sinCenter/')
NGP_RFBenchmark_Denis_131PPC_128cells_noSmoothing_sinCenter_Curv = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_131PPC_128cells_noSmoothing_sinCenter_Curv/')
NGP_RFBenchmark_Denis_262PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_262PPC_128cells_noSmoothing_sinCenter/')
NGP_RFBenchmark_Denis_524PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_524PPC_128cells_noSmoothing_sinCenter/')
NGP_RFBenchmark_Denis_1048PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_1048PPC_128cells_noSmoothing_sinCenter/')

NGP_RFBenchmark_Denis_131PPC_128cells_Smoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_131PPC_128cells_Smoothing_sinCenter/')
NGP_RFBenchmark_Denis_262PPC_128cells_Smoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_262PPC_128cells_Smoothing_sinCenter/')
NGP_RFBenchmark_Denis_524PPC_128cells_Smoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_524PPC_128cells_Smoothing_sinCenter/')
NGP_RFBenchmark_Denis_1048PPC_128cells_Smoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_1048PPC_128cells_Smoothing_sinCenter/')

NGP_RFBenchmark_Denis_131PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_131PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_262PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_262PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_524PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_524PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_1048PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_1048PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_1048PPC_256cells_noSmoothing_PCG = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_1048PPC_256cells_noSmoothing_PCG/')
NGP_RFBenchmark_Denis_2096PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_2096PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_2096PPC_256cells_noSmoothing_OGTime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_2096PPC_256cells_noSmoothing_OGTime/')

NGP_RFBenchmark_Denis_1048PPC_1024cells_noSmoothing_OGTime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_1048PPC_1024cells_noSmoothing_OGTime/')

NGP_RFBenchmark_Denis_262PPC_128cells_noSmoothing_normTime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_262PPC_128cells_noSmoothing_normTime/')

NGP_RFBenchmark_Denis_lowerDen_131PPC_128cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_131PPC_128cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_262PPC_128cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_262PPC_128cells_noSmoothing/')


NGP_RFBenchmark_Denis_lowerDen_524PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_524PPC_128cells_noSmoothing_sinCenter/')
NGP_RFBenchmark_Denis_lowerDen_1048PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_128cells_noSmoothing_sinCenter/')
NGP_RFBenchmark_Denis_lowerDen_2096PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_2096PPC_128cells_noSmoothing_sinCenter/')

NGP_RFBenchmark_Denis_lowerDen_1048PPC_128cells_noSmoothing_sinCenter_50 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_128cells_noSmoothing_sinCenter_50/')

NGP_RFBenchmark_Denis_lowerDen_131PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_131PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_262PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_262PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_524PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_524PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_1048PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_256cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_1048PPC_256cells_noSmoothing_PCG = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_256cells_noSmoothing_PCG/')
NGP_RFBenchmark_Denis_lowerDen_2096PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_2096PPC_256cells_noSmoothing/')

NGP_RFBenchmark_Denis_lowerDen_131PPC_512cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_131PPC_512cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_262PPC_512cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_262PPC_512cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_524PPC_512cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_524PPC_512cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_1048PPC_512cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_512cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_1048PPC_512cells_noSmoothing_halfTime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_512cells_noSmoothing_halfTime/')
NGP_RFBenchmark_Denis_lowerDen_1048PPC_512cells_noSmoothing_OGTime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_512cells_noSmoothing_OGTime/')

NGP_RFBenchmark_Denis_lowerDen_524PPC_1024cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_524PPC_1024cells_noSmoothing/')
NGP_RFBenchmark_Denis_lowerDen_1048PPC_1024cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_1024cells_noSmoothing/')

NGP_RFBenchmark_Denis_lowerDen_524PPC_2000cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_1048PPC_1024cells_noSmoothing/')

NGP_RFBenchmark_Denis_lowerDen_524PPC_2000cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_Denis_lowerDen_524PPC_2000cells_noSmoothing/')

#-------------------------- CIC -----------------

CIC_RFBenchmark_Denis_131PPC_128cells = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_131PPC_128cells/')
CIC_RFBenchmark_Denis_262PPC_128cells = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_262PPC_128cells/')
CIC_RFBenchmark_Denis_262PPC_128cells_NoSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_262PPC_128cells_NoSmoothing/')
CIC_RFBenchmark_Denis_262PPC_128cells_normTime = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_262PPC_128cells_normTime/')


CIC_RFBenchmark_Denis_131PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_131PPC_256cells_noSmoothing/')
CIC_RFBenchmark_Denis_262PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_262PPC_256cells_noSmoothing/')
CIC_RFBenchmark_Denis_524PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_524PPC_256cells_noSmoothing/')
CIC_RFBenchmark_Denis_1048PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_1048PPC_256cells_noSmoothing/')
CIC_RFBenchmark_Denis_2096PPC_256cells_noSmoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_2096PPC_256cells_noSmoothing/')

CIC_RFBenchmark_Denis_131PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_131PPC_128cells_noSmoothing_sinCenter/')
CIC_RFBenchmark_Denis_262PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_262PPC_128cells_noSmoothing_sinCenter/')
CIC_RFBenchmark_Denis_524PPC_128cells_noSmoothing_sinCenter = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_Denis_524PPC_128cells_noSmoothing_sinCenter/')

#------------- plots --------------------

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells, 'He+', '131 PPC')
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_doublePart, 'He+', '262 PPC')
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_quadPart, 'He+', '524 PPC')
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_40xPart, 'He+', '5240 PPC')
plotAveDensity(Exp_RFBenchmark_Denis_doubleTime_halfPart, 'He+', r'128 PPC, x2 $\Delta t$, 1024 cells')
plt.xlim(Exp_RFBenchmark_Denis_moreCells.grid[0], Exp_RFBenchmark_Denis_moreCells.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/Exp_diffCases.pdf')
plt.savefig('Denis_RF/Exp_diffCases.png')
plt.close()

plt.figure()
plotAveEPF(Exp_RFBenchmark_Denis_moreCells, 'e', label = '131 PPC')
plotAveEPF(Exp_RFBenchmark_Denis_moreCells_40xPart, 'e', label = '5240 PPC')
plotAveEPF(Exp_RFBenchmark_Denis_doubleTime_halfPart, 'e', label = r'128 PPC, x2 $\Delta t$, 1024 cells')
plt.ylim(0,None)
plt.xlim(0.1,150)
plt.xlabel('Energy (eV)', fontsize = 14)
plt.ylabel(r'EPDF (eV$^{-3/2}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/EPF_Exp_comp.pdf')
plt.savefig('Denis_RF/EPF_Exp_comp.png')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_40xPart, 'He+', 'Ref')
plotAveDensity(NGP_RFBenchmark_Denis_262PPC_128cells_noSmoothing, 'He+', '262 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_524PPC_128cells_noSmoothing, 'He+', '524 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells.grid[0], Exp_RFBenchmark_Denis_moreCells.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/NGP_128cells_even.pdf')
plt.savefig('Denis_RF/NGP_128cells_even.png')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_40xPart, 'He+', 'Ref')
plotAveDensity(NGP_RFBenchmark_Denis_131PPC_256cells_noSmoothing, 'He+', '131 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_262PPC_256cells_noSmoothing, 'He+', '262 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_524PPC_256cells_noSmoothing, 'He+', '524 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_1048PPC_256cells_noSmoothing, 'He+', '1048 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells.grid[0], Exp_RFBenchmark_Denis_moreCells.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/NGP_256cells_even.pdf')
plt.savefig('Denis_RF/NGP_256cells_even.png')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_40xPart, 'He+', 'Ref')
plotAveDensity(NGP_RFBenchmark_Denis_131PPC_128cells_noSmoothing_sinCenter, 'He+', '131 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_262PPC_128cells_noSmoothing_sinCenter, 'He+', '262 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_524PPC_128cells_noSmoothing_sinCenter, 'He+', '524 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_1048PPC_128cells_noSmoothing_sinCenter, 'He+', '1048 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells.grid[0], Exp_RFBenchmark_Denis_moreCells.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/NGP_128cells_noneven.pdf')
plt.savefig('Denis_RF/NGP_128cells_noneven.png')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_40xPart, 'He+', 'Ref')
plotAveDensity(NGP_RFBenchmark_Denis_131PPC_128cells_Smoothing_sinCenter, 'He+', '131 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_262PPC_128cells_Smoothing_sinCenter, 'He+', '262 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_524PPC_128cells_Smoothing_sinCenter, 'He+', '524 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_1048PPC_128cells_Smoothing_sinCenter, 'He+', '1048 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells.grid[0], Exp_RFBenchmark_Denis_moreCells.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/NGP_128cells_noneven_smooth.pdf')
plt.savefig('Denis_RF/NGP_128cells_noneven_smooth.png')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_40xPart, 'He+', 'Ref')
plotAveDensity(CIC_RFBenchmark_Denis_131PPC_128cells_noSmoothing_sinCenter, 'He+', '131 PPC')
plotAveDensity(CIC_RFBenchmark_Denis_262PPC_128cells_noSmoothing_sinCenter, 'He+', '262 PPC')
plotAveDensity(CIC_RFBenchmark_Denis_524PPC_128cells_noSmoothing_sinCenter, 'He+', '524 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells.grid[0], Exp_RFBenchmark_Denis_moreCells.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/CIC_128cells_noneven.pdf')
plt.savefig('Denis_RF/CIC_128cells_noneven.png')
plt.close()


# ------------------- lower pressure -------------------------------

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_lowerDen, 'He+', '131 PPC')
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_lowerDen_doublePart, 'He+', '262 PPC')
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_lowerDen_quadPart, 'He+', '524 PPC')
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_lowerDen_10xPart, 'He+', '1310 PPC')
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_lowerDen_20xPart, 'He+', r'2620 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells_lowerDen.grid[0], Exp_RFBenchmark_Denis_moreCells_lowerDen.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/Exp_diffCases_lowDen.pdf')
plt.savefig('Denis_RF/Exp_diffCases_lowDen.png')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_lowerDen_20xPart, 'He+', 'Ref')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_131PPC_256cells_noSmoothing, 'He+', '131 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_262PPC_256cells_noSmoothing, 'He+', '262 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_524PPC_256cells_noSmoothing, 'He+', '524 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_1048PPC_256cells_noSmoothing, 'He+', '1048 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_2096PPC_256cells_noSmoothing, 'He+', '2096 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells_lowerDen.grid[0], Exp_RFBenchmark_Denis_moreCells_lowerDen.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/NGP_lowPress_256cells_even.pdf')
plt.savefig('Denis_RF/NGP_lowPress_256cells_even.png')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_lowerDen_20xPart, 'He+', 'Ref')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_131PPC_512cells_noSmoothing, 'He+', '131 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_262PPC_512cells_noSmoothing, 'He+', '262 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_524PPC_512cells_noSmoothing, 'He+', '524 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_1048PPC_512cells_noSmoothing, 'He+', '1048 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells_lowerDen.grid[0], Exp_RFBenchmark_Denis_moreCells_lowerDen.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/NGP_lowPress_512cells_even.pdf')
plt.savefig('Denis_RF/NGP_lowPress_512cells_even.png')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_Denis_moreCells_lowerDen_20xPart, 'He+', 'Ref')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_524PPC_128cells_noSmoothing_sinCenter, 'He+', '524 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_1048PPC_128cells_noSmoothing_sinCenter, 'He+', '1048 PPC')
plotAveDensity(NGP_RFBenchmark_Denis_lowerDen_2096PPC_128cells_noSmoothing_sinCenter, 'He+', '2096 PPC')
plt.xlim(Exp_RFBenchmark_Denis_moreCells_lowerDen.grid[0], Exp_RFBenchmark_Denis_moreCells_lowerDen.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('Denis_RF/NGP_lowPress_128cells_noneven.pdf')
plt.savefig('Denis_RF/NGP_lowPress_128cells_noneven.png')
plt.close()