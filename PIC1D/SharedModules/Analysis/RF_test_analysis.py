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

elasticAr = np.loadtxt('Y:/ImplicitPic1D/ImplicitPic1D/OpenMP/CollisionData/ElasticAr.txt')
elasticAr_mod = np.loadtxt('Y:/ImplicitPic1D/ImplicitPic1D/OpenMP/CollisionData/newElasticAr.txt')

#case 1
Exp_RFBenchmark_case1 = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1/')
# Exp_RFBenchmark_case1_0p05delT_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_0p05delT_10xParticles/')
# Exp_RFBenchmark_case1_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_10xParticles/')
Exp_RFBenchmark_case1_20xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_20xParticles/')
Exp_RFBenchmark_case1_halfTime_doubleCells_doubleParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_halfTime_doubleCells_doubleParticles/')
#Exp_RFBenchmark_case1_halfCells_doubleTime_halfParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_halfCells_doubleTime_halfParticles/')
# Exp_RFBenchmark_case1_34Cells_512PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_34Cells_512PPC/')
# Exp_RFBenchmark_case1_34Cells_1024PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_34Cells_1024PPC/')
# Exp_RFBenchmark_case1_34Cells_0p5delT_1024PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_34Cells_0p5delT_1024PPC/')
# Exp_RFBenchmark_case1_34Cells_0p5delT_512PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_34Cells_0p5delT_512PPC/')
#
Exp_RFBenchmark_case2 = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2/')
Exp_RFBenchmark_case2_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_10xParticles/')
Exp_RFBenchmark_case2_halfTime = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_halfTime/')
Exp_RFBenchmark_case2_quarterTime = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_quarterTime/')
# Exp_RFBenchmark_case2_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_10xParticles/')
#
# Exp_RFBenchmark_case2_40amu = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_40amu/')
# Exp_RFBenchmark_case2_40amu_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case2_40amu_10xParticles/')
#
Exp_VassBenchmark_Ar = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar/')
Exp_VassBenchmark_Ar_changeCross = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_changeCross/')
Exp_VassBenchmark_Ar_changeCross_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_changeCross_10xParticles/')
# Exp_VassBenchmark_Ar_noRam_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_noRam_10xParticles/')
# Exp_VassBenchmark_Ar_512Cells_200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_512Cells_200PPC/')
# Exp_VassBenchmark_Ar_1024Cells_200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_1024Cells_200PPC/')
# Exp_VassBenchmark_Ar_1024Cells_800PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_1024Cells_800PPC/')
Exp_VassBenchmark_Ar_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_10xParticles/')
Exp_VassBenchmark_Ar_lowerDensity = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_lowerDensity/')
Exp_VassBenchmark_Ar_lowerDensity_10xParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_VassBenchmark_Ar_lowerDensity_10xParticles/')
#
NGP_RF_Benchmark_case1 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1/')
NGP_RF_Benchmark_case2 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case2/')
NGP_RF_Benchmark_case2_noSmooth = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case2_noSmooth/')
NGP_RF_Benchmark_case2_noSmooth_changeMover = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case2_noSmooth_changeMover/')
NGP_RF_Benchmark_case2_halfTime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case2_halfTime/')
NGP_RF_Benchmark_case2_10xParticles = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case2_10xParticles/')
# NGP_RF_Benchmark_case1_lowerRes = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1_lowerRes/')
NGP_RF_Benchmark_case1_10xParticles = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1_10xParticles/')
# NGP_VassBenchmark_Ar = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_VassBenchmark_Ar/')
# NGP_VassBenchmark_Ar_10xParticles = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_VassBenchmark_Ar_10xParticles/')
# NGP_RF_Benchmark_case1_restart = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1_restart/')
# NGP_RFBenchmark_case1_retry = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_retry/')
# NGP_RFBenchmark_case1_Picard = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_Picard/')
# NGP_RFBenchmark_case1_cluster2 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_cluster2/')
#
#
# Exp_RFBenchmark_EreminCase_resolved = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_EreminCase_resolved/')
# Exp_RFBenchmark_EreminCase_quadTime_p75Cells_eighthNumPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_EreminCase_quadTime_p75Cells_eighthNumPart/')
#
# NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p0005edge_sixteenthNumPart = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p0005edge_sixteenthNumPart/')
# NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart/')
# NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard/')
# NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard/')
# NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_Picard = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_Picard/')
# NGP_RFBenchmark_EreminCase_5xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_5xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard/')
# CIC_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard/')
#
Exp_RFBenchmark_highDensity_lowRes = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_lowRes/')
Exp_RFBenchmark_highDensity_resolved = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_resolved/')
Exp_RFBenchmark_highDensity_MoreResolved = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_MoreResolved/')
Exp_RFBenchmark_highDensity_EvenMoreResolved = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_EvenMoreResolved/')
Exp_RFBenchmark_highDensity_EvenMoreResolved_2 = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_EvenMoreResolved_2/')
Exp_RFBenchmark_highDensity_SuperResolved = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_SuperResolved/')
Exp_RFBenchmark_highDensity_EvenMoreResolved_x3 = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_EvenMoreResolved_x3/')
Exp_RFBenchmark_highDensity_EvenMoreResolved_x6 = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_EvenMoreResolved_x6/')
Exp_RFBenchmark_highDensity_p75Cells_doubleTime_halfParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_p75Cells_doubleTime_halfParticles/')
Exp_RFBenchmark_highDensity_p75Cells_doubleTime_quarterParticles = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_p75Cells_doubleTime_quarterParticles/')

#
# NGP_RFBenchmark_highDensity_sixteenthPart_129nodes = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_129nodes/')
# NGP_RFBenchmark_highDensity_sixteenthPart_257nodes = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_257nodes/')
# NGP_RFBenchmark_highDensity_sixteenthPart_257nodes_normTime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_257nodes_normTime/')
# NGP_RFBenchmark_highDensity_sixteenthPart_513nodes_normTime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_513nodes_normTime/')
# NGP_RFBenchmark_highDensity_200PPC_1024Cells = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_200PPC_1024Cells/')
#

NGP_RFBenchmark_highDensity_halfPart_1001nodes = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_1001nodes/')

NGP_RFBenchmark_highDensity_eighthPart_1001nodes_epsNeg8 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_eighthPart_1001nodes_epsNeg8/')
NGP_RFBenchmark_highDensity_fourthPart_1001nodes_epsNeg8 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_fourthPart_1001nodes_epsNeg8/')
NGP_RFBenchmark_highDensity_fourthPart_1001nodes_epsNeg8_redo = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_fourthPart_1001nodes_epsNeg8_redo/')
NGP_RFBenchmark_highDensity_fourthPart_1001nodes_10xtime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_fourthPart_1001nodes_10xtime/')
NGP_RFBenchmark_highDensity_halfPart_1001nodes_10xtime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_halfPart_1001nodes_10xtime/')
NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_halfPart_513nodes_10xtime/')
NGP_RFBenchmark_highDensity_fourthPart_257nodes_10xtime_halfEvenHalfSin = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_fourthPart_257nodes_10xtime_halfEvenHalfSin/')
# NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime_restart = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime_restart/')
NGP_RFBenchmark_highDensity_halfPart_257nodes_10xtime_halfEvenHalfSin = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_257nodes_10xtime_halfEvenHalfSin/')
# NGP_RFBenchmark_highDensity_halfPart_2001_10xtime_even = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_2001_10xtime_even/')
NGP_RF_Benchmark_case4_201nodes_halfSinHalfEven = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case4_201nodes_halfSinHalfEven/')
NGP_RF_Benchmark_case4_201nodes_halfSinHalfEven_NoSmooth = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case4_201nodes_halfSinHalfEven_NoSmooth/')
NGP_RFBenchmark_highDensity_halfPart_257nodes_halfEvenHalfSin_NoSmooth = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_257nodes_halfEvenHalfSin_NoSmooth/')

CIC_RF_Benchmark_case1 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1/')
CIC_RF_Benchmark_case2 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case2/')
CIC_RF_Benchmark_case2_epsNeg8 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case2_epsNeg8/')
CIC_RF_Benchmark_case2_epsNeg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case2_epsNeg3/')
CIC_RF_Benchmark_case2_10xParticles = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case2_10xParticles/')
CIC_RF_Benchmark_case2_epsNeg2 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case2_epsNeg2/')
CIC_RF_Benchmark_case1_10xParticles = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1_10xParticles/')
# CIC_RF_Benchmark_case1_nonSmooth = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1_nonSmooth/')
# CIC_RF_Benchmark_case1_picard = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1_picard/')
# CIC_RF_Benchmark_case1_picardHalf = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1_picardHalf/')


CIC_RFBenchmark_highDensity_200nodes_320PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_200nodes_320PPC/')
CIC_RFBenchmark_highDensity_200nodes_96PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_200nodes_96PPC/')
CIC_RFBenchmark_highDensity_halfPart_1000nodes_epsNeg8 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_halfPart_1000nodes_epsNeg8/')
CIC_RFBenchmark_highDensity_fourthPart_1000nodes = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_fourthPart_1000nodes/')
CIC_RFBenchmark_highDensity_halfPart_1000nodes = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_halfPart_1000nodes/')
CIC_RFBenchmark_highDensity_threeQuartPart_1000nodes = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_threeQuartPart_1001nodes/')
CIC_RFBenchmark_highDensity_halfPart_256nodes_halfEvenHalfSin = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_halfPart_256nodes_halfEvenHalfSin/')
CIC_RFBenchmark_highDensity_halfPart_256nodes_halfEvenHalfSin_NoSmooth = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_halfPart_256nodes_halfEvenHalfSin_NoSmooth/')


CIC_RFBenchmark_highDensity_quarterPart_1000nodes_epsNeg5 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_quarterPart_1000nodes_epsNeg5/')
CIC_RFBenchmark_highDensity_eighthPart_1000nodes_epsNeg5 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_eighthPart_1000nodes_epsNeg5/')
CIC_RFBenchmark_highDensity_eighthPart_1000nodes_epsNeg8 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_eighthPart_1000nodes_epsNeg8/')
CIC_RFBenchmark_highDensity_halfPart_1000nodes_10xtime = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_halfPart_1000nodes_10xtime/')
# CIC_RFBenchmark_highDensity_200PPC_1024Cells = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_highDensity_200PPC_1024Cells/')

plt.figure()
plt.errorbar(case1[:,0], case1[:,4], case1[:,6], linestyle = '-', marker = '.', elinewidth=2, label = 'Benchmark')
plotAveDensity(Exp_RFBenchmark_case1, 'He+', 'Exp.')
plotAveDensity(NGP_RF_Benchmark_case1, 'He+', 'INGP')
plotAveDensity(CIC_RF_Benchmark_case1, 'He+', 'ICIC')
plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('TurnerBenchmark/Case1_Exp_Imp.pdf')
plt.close()

plt.figure()
plotAveVDF(Exp_RFBenchmark_case1, 'e', label = 'Exp.')
plotAveVDF(NGP_RF_Benchmark_case1, 'e', label = 'INGP')
plotAveVDF(CIC_RF_Benchmark_case1, 'e', label = 'ICIC')
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('TurnerBenchmark/Case1_Exp_Imp_VDF.pdf')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_case1, 'He+', 'Exp. Normal')
plotAveDensity(Exp_RFBenchmark_case1_20xParticles, 'He+', 'Exp. 20x Particles')
# plotAveDensity(NGP_RF_Benchmark_case1_lowerRes, 'He+', r'NGP, $\varepsilon_a = 10^{-2}$')
plotAveDensity(Exp_RFBenchmark_case1_halfTime_doubleCells_doubleParticles, 'He+', 'Exp. Double Res.')
plotAveDensity(NGP_RF_Benchmark_case1_10xParticles, 'He+', 'INGP 10x Particles')
plotAveDensity(CIC_RF_Benchmark_case1_10xParticles, 'He+', 'ICIC 10x Particles')
plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1]*1e2)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('TurnerBenchmark/Case1_IncreasedRes.pdf')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_case2, 'He+', 'Exp.')
plotAveDensity(Exp_RFBenchmark_case2_10xParticles, 'He+', 'Exp. 10x Particles')
plotAveDensity(NGP_RF_Benchmark_case2, 'He+', 'INGP')
plotAveDensity(NGP_RF_Benchmark_case2_10xParticles, 'He+', 'INGP 10x Particles')
plotAveDensity(CIC_RF_Benchmark_case2, 'He+', 'ICIC')
plotAveDensity(CIC_RF_Benchmark_case2_10xParticles, 'He+', 'ICIC 10x Particles')
plt.errorbar(case2[0:257:8,0], case2[0:257:8,4], case2[0:257:8,6], linestyle = '-', marker = '.', elinewidth=2, label = 'Benchmark')
plt.xlim(Exp_RFBenchmark_case1.grid[-1] * 0.25, Exp_RFBenchmark_case1.grid[-1] * 0.75)
plt.ylim(5e14, 9e14)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('TurnerBenchmark/Case2_Exp_Imp.pdf')
plt.close()

plt.figure()
plotAveDensity(Exp_VassBenchmark_Ar, 'e', 'Original')
plotAveDensity(Exp_VassBenchmark_Ar_10xParticles, 'e', 'Original 10x Particles')
plotAveDensity(Exp_VassBenchmark_Ar_changeCross, 'e', 'Altered')
plotAveDensity(Exp_VassBenchmark_Ar_changeCross_10xParticles, 'e', 'Altered 10x Particles')
plt.xlim(Exp_VassBenchmark_Ar.grid[0], Exp_VassBenchmark_Ar.grid[-1])
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'$n_e$ (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('TurnerBenchmark/Argon_particleComp.pdf')
plt.close()

plt.figure()
plt.plot(elasticAr[:,0], elasticAr[:, 1], '--o', markersize = 2, label = 'Original')
plt.plot(elasticAr_mod[:,0], elasticAr_mod[:,1], '--o', markersize = 2, label = 'Altered')
plt.xscale('symlog')
plt.yscale('log')
plt.xlabel('Center of Mass Energy (eV)', fontsize = 14)
plt.ylabel(r'Argon Elastic Cross Section (m$^2$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('TurnerBenchmark/Argon_elasticCross.pdf')
plt.close()

X_axis_label = []
length = len(Exp_VassBenchmark_Ar.binaryColl['e->Ar'])
freq = np.zeros((2, length))
for i in range(length):
    X_axis_label.append(Exp_VassBenchmark_Ar.binaryColl['e->Ar'][i]['type'])
    freq[0, i] = Exp_VassBenchmark_Ar.binaryColl['e->Ar'][i]['aveDiag']['ratio']/Exp_VassBenchmark_Ar.delT
    freq[1, i] = Exp_VassBenchmark_Ar_10xParticles.binaryColl['e->Ar'][i]['aveDiag']['ratio'] / Exp_VassBenchmark_Ar.delT
X_axis = np.arange(len(Exp_VassBenchmark_Ar.binaryColl['e->Ar']))
plt.figure()
plt.bar(X_axis, freq[1, :]/freq[0,:], 0.4, label = 'Original')
#plt.bar(X_axis + 0.2, freq[1,:], 0.4, label = 'Altered')
plt.xticks(X_axis, X_axis_label)
plt.ylabel(r'Ratio of Collision Frequencies')
plt.savefig('TurnerBenchmark/Argon_ratioPPCFreq.pdf')
plt.close()

plt.figure()
plotAveDensity(Exp_RFBenchmark_highDensity_lowRes, 'He+', 'Exp. #1')
plotAveDensity(Exp_RFBenchmark_highDensity_MoreResolved, 'He+', 'Exp. #2')
plotAveDensity(Exp_RFBenchmark_highDensity_EvenMoreResolved_2, 'He+', 'Exp. #3')
plotAveDensity(Exp_RFBenchmark_highDensity_SuperResolved, 'He+', 'Exp. #4')
plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
plt.ylim(0, None)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('HighDensityRF/Exp_diffCases.pdf')
plt.close()

# plt.figure()
# plotAveDensity(Exp_RFBenchmark_highDensity_SuperResolved, 'He+', 'Exp. #4')
# plotAveDensity(CIC_RFBenchmark_highDensity_fourthPart_1000nodes, 'He+', 'ICIC. #1')
# plotAveDensity(CIC_RFBenchmark_highDensity_halfPart_1000nodes, 'He+', 'ICIC. #2')
# plotAveDensity(CIC_RFBenchmark_highDensity_threeQuartPart_1000nodes, 'He+', 'ICIC. #3')
# plotAveDensity(CIC_RFBenchmark_highDensity_halfPart_256nodes_halfEvenHalfSin, 'He+', 'ICIC. #4')
# plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
# plt.ylim(0, None)
# plt.xlabel('Distance (m)', fontsize = 14)
# plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.legend(loc = 'lower center', fontsize = 12)
# plt.tight_layout()
# plt.savefig('HighDensityRF/CIC_diffCases.pdf')
# plt.close()

# plt.figure()
# plotAveDensity(Exp_RFBenchmark_highDensity_SuperResolved, 'He+', 'Exp. #4')
# plotAveDensity(Exp_RFBenchmark_highDensity_lowRes, 'He+', 'Exp. #1')
# plotAveDensity(CIC_RFBenchmark_highDensity_200nodes_96PPC, 'He+', 'ICIC. #5')
# plotAveDensity(CIC_RFBenchmark_highDensity_200nodes_320PPC, 'He+', 'ICIC. #6')
# plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
# plt.ylim(0, None)
# plt.xlabel('Distance (m)', fontsize = 14)
# plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.legend(loc = 'lower center', fontsize = 12)
# plt.tight_layout()
# plt.savefig('HighDensityRF/lowRes_diffCases.pdf')
# plt.close()


plt.figure()
plotAvePhi(Exp_RFBenchmark_highDensity_EvenMoreResolved_2, label = 'Exp. #3')
plotAvePhi(CIC_RFBenchmark_highDensity_threeQuartPart_1000nodes, label = 'ICIC #3')
plotAvePhi(CIC_RFBenchmark_highDensity_halfPart_256nodes_halfEvenHalfSin, label = 'ICIC #4')
plotAvePhi(CIC_RFBenchmark_highDensity_200nodes_320PPC, label = 'ICIC #6')
plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
plt.ylim(700, 860)
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower center', fontsize = 12)
plt.tight_layout()
plt.savefig('HighDensityRF/Phi_comp.pdf')
plt.close()

plt.figure()
plotAveEPF(Exp_RFBenchmark_highDensity_EvenMoreResolved_2, 'e', label = 'Exp. #3')
plotAveEPF(CIC_RFBenchmark_highDensity_threeQuartPart_1000nodes, 'e', label = 'ICIC #3')
plotAveEPF(CIC_RFBenchmark_highDensity_halfPart_256nodes_halfEvenHalfSin, 'e', label = 'ICIC #4')
plotAveEPF(CIC_RFBenchmark_highDensity_200nodes_320PPC, 'e', label = 'ICIC #6')
plt.ylim(0,None)
plt.xlim(0.1,150)
plt.xlabel('Energy (eV)', fontsize = 14)
plt.ylabel(r'EPDF (eV$^{-3/2}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('HighDensityRF/EPF_comp.pdf')
plt.close()

dataList = [CIC_RFBenchmark_highDensity_fourthPart_1000nodes, \
            CIC_RFBenchmark_highDensity_halfPart_1000nodes, \
            CIC_RFBenchmark_highDensity_threeQuartPart_1000nodes, \
            CIC_RFBenchmark_highDensity_halfPart_256nodes_halfEvenHalfSin]
labelList = ['ICIC #1', 'ICIC #2', 'ICIC #3', 'ICIC #4']

plt.figure()
plotAveDensityVsRef_percent(Exp_RFBenchmark_highDensity_SuperResolved, dataList, 'He+', labelList = labelList, marker = 'o', linestyle = '--')
plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'$\Delta n_{He+} / n_{He+}$', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('HighDensityRF/CIC_diffCases.pdf')
plt.close()


dataList = [Exp_RFBenchmark_highDensity_lowRes, \
            CIC_RFBenchmark_highDensity_200nodes_96PPC, \
            CIC_RFBenchmark_highDensity_200nodes_320PPC]
labelList = ['Exp #1', 'ICIC #5', 'ICIC #6']

plt.figure()
plotAveDensityVsRef_percent(Exp_RFBenchmark_highDensity_SuperResolved, dataList, 'He+', labelList = labelList, marker = 'o', linestyle = '--')
plt.xlim(Exp_RFBenchmark_case1.grid[0], Exp_RFBenchmark_case1.grid[-1])
plt.xlabel('Distance (m)', fontsize = 14)
plt.ylabel(r'$\Delta n_{He+} / n_{He+}$', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best', fontsize = 12)
plt.tight_layout()
plt.savefig('HighDensityRF/lowRes_diffCases.pdf')
plt.close()
