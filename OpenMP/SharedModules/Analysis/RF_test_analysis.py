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


#case 1
Exp_RFBenchmark_case1 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1/')
# Exp_RFBenchmark_case1_0p05delT_10xParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_0p05delT_10xParticles/')
# Exp_RFBenchmark_case1_10xParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_10xParticles/')
# Exp_RFBenchmark_case1_20xParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_20xParticles/')
# Exp_RFBenchmark_case1_halfTime_doubleCells_doubleParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_halfTime_doubleCells_doubleParticles/')
# Exp_RFBenchmark_case1_halfCells_doubleTime_halfParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_halfCells_doubleTime_halfParticles/')
# Exp_RFBenchmark_case1_34Cells_512PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_34Cells_512PPC/')
# Exp_RFBenchmark_case1_34Cells_1024PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_34Cells_1024PPC/')
# Exp_RFBenchmark_case1_34Cells_0p5delT_1024PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_34Cells_0p5delT_1024PPC/')
# Exp_RFBenchmark_case1_34Cells_0p5delT_512PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_34Cells_0p5delT_512PPC/')
#
# Exp_RFBenchmark_case2 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2/')
# Exp_RFBenchmark_case2_10xParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_10xParticles/')
#
# Exp_RFBenchmark_case2_40amu = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_40amu/')
# Exp_RFBenchmark_case2_40amu_10xParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_40amu_10xParticles/')
#
# Exp_VassBenchmark_Ar = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar/')
# Exp_VassBenchmark_Ar_noRam = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_noRam/')
# Exp_VassBenchmark_Ar_noRam_10xParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_noRam_10xParticles/')
# Exp_VassBenchmark_Ar_512Cells_200PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_512Cells_200PPC/')
# Exp_VassBenchmark_Ar_1024Cells_200PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_1024Cells_200PPC/')
# Exp_VassBenchmark_Ar_1024Cells_800PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_1024Cells_800PPC/')
# Exp_VassBenchmark_Ar_10xParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_10xParticles/')
#
NGP_RF_Benchmark_case1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1/')
# NGP_RF_Benchmark_case1_lowerRes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1_lowerRes/')
# NGP_RF_Benchmark_case1_10xParticles = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1_10xParticles/')
# NGP_VassBenchmark_Ar = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_VassBenchmark_Ar/')
# NGP_VassBenchmark_Ar_10xParticles = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_VassBenchmark_Ar_10xParticles/')
# NGP_RF_Benchmark_case1_restart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1_restart/')
# NGP_RFBenchmark_case1_retry = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_retry/')
# NGP_RFBenchmark_case1_Picard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_Picard/')
# NGP_RFBenchmark_case1_cluster2 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_cluster2/')
#
#
# Exp_RFBenchmark_EreminCase_resolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_EreminCase_resolved/')
# Exp_RFBenchmark_EreminCase_quadTime_p75Cells_eighthNumPart = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_EreminCase_quadTime_p75Cells_eighthNumPart/')
#
# NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p0005edge_sixteenthNumPart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p0005edge_sixteenthNumPart/')
# NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart/')
# NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard/')
# NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard/')
# NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_Picard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_Picard/')
# NGP_RFBenchmark_EreminCase_5xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_5xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard/')
# CIC_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard/')
#
# Exp_RFBenchmark_highDensity_resolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_resolved/')
# Exp_RFBenchmark_highDensity_MoreResolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_MoreResolved/')
# Exp_RFBenchmark_highDensity_EvenMoreResolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_EvenMoreResolved/')
# Exp_RFBenchmark_highDensity_SuperResolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_SuperResolved/')
# Exp_RFBenchmark_highDensity_p75Cells_doubleTime_halfParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_p75Cells_doubleTime_halfParticles/')
# Exp_RFBenchmark_highDensity_p75Cells_doubleTime_quarterParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_p75Cells_doubleTime_quarterParticles/')
#
# NGP_RFBenchmark_highDensity_sixteenthPart_129nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_129nodes/')
# NGP_RFBenchmark_highDensity_sixteenthPart_257nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_257nodes/')
# NGP_RFBenchmark_highDensity_sixteenthPart_257nodes_normTime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_257nodes_normTime/')
# NGP_RFBenchmark_highDensity_sixteenthPart_513nodes_normTime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_513nodes_normTime/')
# NGP_RFBenchmark_highDensity_200PPC_1024Cells = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_200PPC_1024Cells/')
#
# NGP_RFBenchmark_highDensity_fourthPart_1001nodes_10xtime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_fourthPart_1001nodes_10xtime/')
# NGP_RFBenchmark_highDensity_halfPart_1001nodes_10xtime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_halfPart_1001nodes_10xtime/')
# NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_halfPart_513nodes_10xtime/')
# NGP_RFBenchmark_highDensity_fourthPart_257nodes_10xtime_halfEvenHalfSin = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_fourthPart_257nodes_10xtime_halfEvenHalfSin/')
# NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime_restart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime_restart/')
# NGP_RFBenchmark_highDensity_halfPart_257nodes_10xtime_halfEvenHalfSin = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_257nodes_10xtime_halfEvenHalfSin/')
# NGP_RFBenchmark_highDensity_halfPart_2001_10xtime_even = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_2001_10xtime_even/')
#
CIC_RF_Benchmark_case1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1/')
# CIC_RF_Benchmark_case1_10xParticles = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1_10xParticles/')
# CIC_RF_Benchmark_case1_nonSmooth = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1_nonSmooth/')
# CIC_RF_Benchmark_case1_picard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1_picard/')
# CIC_RF_Benchmark_case1_picardHalf = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1_picardHalf/')
#
# CIC_RFBenchmark_highDensity_halfPart_1000nodes_10xtime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_halfPart_1000nodes_10xtime/')
# CIC_RFBenchmark_highDensity_200PPC_1024Cells = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_200PPC_1024Cells/')

# plt.figure()
# plt.errorbar(case1[:,0], case1[:,4], case1[:,6], linestyle = '-', marker = '.', elinewidth=2, label = 'Benchmark')
# plotAveDensity(Exp_RFBenchmark_case1, 'He+', 'Exp.')
# plotAveDensity(NGP_RF_Benchmark_case1, 'He+', 'NGP')
# plotAveDensity(CIC_RF_Benchmark_case1, 'He+', 'CIC')
# plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
# plt.xlabel('Distance (m)')
# plt.ylabel(r'He$^+$ Density (m$^{-3}$)')
# plt.legend(loc = 'best')
# plt.savefig('TurnerBenchmark/Case1_Exp_Imp.pdf')
# plt.close()
#
# plt.figure()
# plotAveDensity(Exp_RFBenchmark_case1, 'He+', 'Exp. Normal')
# plotAveDensity(Exp_RFBenchmark_case1_20xParticles, 'He+', 'Exp. 20x Particles')
# plotAveDensity(NGP_RF_Benchmark_case1_lowerRes, 'He+', r'NGP, $\varepsilon_a = 10^{-2}$')
# plotAveDensity(Exp_RFBenchmark_case1_halfTime_doubleCells_doubleParticles, 'He+', 'Exp. Double Res.')
# plotAveDensity(NGP_RF_Benchmark_case1_10xParticles, 'He+', 'NGP 10x Particles')
# plotAveDensity(CIC_RF_Benchmark_case1_10xParticles, 'He+', 'CIC 10x Particles')
# plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
# plt.xlabel('Distance (m)')
# plt.ylabel(r'He$^+$ Density (m$^{-3}$)')
# plt.legend(loc = 'best')
# plt.savefig('TurnerBenchmark/Case1_IncreasedRes.pdf')
# plt.close()
#
# plt.figure()
# plotAveDensity(Exp_VassBenchmark_Ar, 'e', 'Exp. Normal')
# plotAveDensity(Exp_VassBenchmark_Ar_10xParticles, 'e', 'Exp. 10x Particles')
# plotAveDensity(Exp_VassBenchmark_Ar_1024Cells_800PPC, 'e', 'Exp. 1024 Cells, 16x Particles')
# plotAveDensity(NGP_VassBenchmark_Ar, 'e', r'NGP Normal')
# plotAveDensity(NGP_VassBenchmark_Ar_10xParticles, 'e', r'NGP 10x Particles')
# plt.xlim(Exp_VassBenchmark_Ar.grid[0], Exp_VassBenchmark_Ar.grid[-1])
# plt.xlabel('Distance (m)')
# plt.ylabel(r'$n_e$ (m$^{-3}$)')
# plt.legend(loc = 'best')
# plt.savefig('TurnerBenchmark/Argon_particleComp.pdf')
# plt.close()
#
# plt.figure()
# plotAveDensity(Exp_VassBenchmark_Ar_noRam, 'e', 'Normal')
# plotAveDensity(Exp_VassBenchmark_Ar_noRam_10xParticles, 'e', '10x Particles')
# plt.xlim(Exp_VassBenchmark_Ar_noRam.grid[0], Exp_VassBenchmark_Ar_noRam.grid[-1])
# plt.xlabel('Distance (m)')
# plt.ylabel(r'$n_e$ (m$^{-3}$)')
# plt.legend(loc = 'best')
# plt.savefig('TurnerBenchmark/Argon_noRam_particleComp.pdf')
# plt.close()
#
# plt.figure()
# plotAveDensity(Exp_RFBenchmark_case2, 'He+', 'Benchmark')
# plotAveDensity(Exp_RFBenchmark_case2_10xParticles, 'He+', '10x Particles')
# plt.xlim(Exp_RFBenchmark_case2.grid[0], Exp_RFBenchmark_case2.grid[-1])
# plt.xlabel('Distance (m)')
# plt.ylabel(r'He$^+$ Density (m$^{-3}$)')
# plt.legend(loc = 'best')
# plt.savefig('TurnerBenchmark/Case2_particleComp.pdf')
# plt.close()
#
# plt.figure()
# plotAveDensity(Exp_RFBenchmark_case2_40amu, 'e', '64 PPC')
# plotAveDensity(Exp_RFBenchmark_case2_40amu_10xParticles, 'e', '640 PPC')
# plt.xlim(Exp_RFBenchmark_case2_40amu.grid[0], Exp_RFBenchmark_case2_40amu.grid[-1])
# plt.xlabel('Distance (m)')
# plt.ylabel(r'$n_e$ (m$^{-3}$)')
# plt.legend(loc = 'best')
# plt.savefig('TurnerBenchmark/Case2_40amu_particleComp.pdf')
# plt.close()
