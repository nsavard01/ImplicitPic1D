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

#--------------------- Explicit ---------------------------------------------
Exp_test_RF = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_test_RF/')
Exp_test_typical = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_test_typical/')
Exp_RF_Helium_Animation = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RF_Helium_Animation/')

#case 1
Exp_RFBenchmark_case1 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1/')
Exp_RFBenchmark_case1_animation = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_case1_animation/')

NGP_RFBenchmark_case1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1/')
NGP_RF_Benchmark_case1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1/')
NGP_RF_Benchmark_case1_restart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RF_Benchmark_case1_restart/')
NGP_RFBenchmark_case1_retry = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_retry/')
NGP_RFBenchmark_case1_Picard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_Picard/')
NGP_RFBenchmark_case1_cluster2 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_case1_cluster2/')

CIC_RFBenchmark_case1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_case1/')

Exp_RFBenchmark_EreminCase_resolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_EreminCase_resolved/')
Exp_RFBenchmark_EreminCase_quadTime_p75Cells_eighthNumPart = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_EreminCase_quadTime_p75Cells_eighthNumPart/')

NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p0005edge_sixteenthNumPart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p0005edge_sixteenthNumPart/')
NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart/')
NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard/')
NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard/')
NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_Picard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_Picard/')
NGP_RFBenchmark_EreminCase_5xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_EreminCase_5xTime_65NodesSinusoid0p00025edge_sixteenthNumPart_nonPicard/')
CIC_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RFBenchmark_EreminCase_15xTime_65NodesEven_sixteenthNumPart_nonPicard/')

Exp_RFBenchmark_highDensity_resolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_resolved/')
Exp_RFBenchmark_highDensity_MoreResolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_MoreResolved/')
Exp_RFBenchmark_highDensity_EvenMoreResolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_EvenMoreResolved/')
Exp_RFBenchmark_highDensity_SuperResolved = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_SuperResolved/')
Exp_RFBenchmark_highDensity_p75Cells_doubleTime_halfParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_p75Cells_doubleTime_halfParticles/')
Exp_RFBenchmark_highDensity_p75Cells_doubleTime_quarterParticles = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_highDensity_p75Cells_doubleTime_quarterParticles/')


NGP_RFBenchmark_highDensity_sixteenthPart_129nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_129nodes/')
NGP_RFBenchmark_highDensity_sixteenthPart_257nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_257nodes/')
NGP_RFBenchmark_highDensity_sixteenthPart_257nodes_normTime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_257nodes_normTime/')
NGP_RFBenchmark_highDensity_sixteenthPart_513nodes_normTime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_sixteenthPart_513nodes_normTime/')
NGP_RFBenchmark_highDensity_fourthPart_1001nodes_10xtime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_fourthPart_1001nodes_10xtime/')
NGP_RFBenchmark_highDensity_halfPart_1001nodes_10xtime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_halfPart_1001nodes_10xtime/')
NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDenity_halfPart_513nodes_10xtime/')
NGP_RFBenchmark_highDensity_fourthPart_257nodes_10xtime_halfEvenHalfSin = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_fourthPart_257nodes_10xtime_halfEvenHalfSin/')
NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime_restart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_RFBenchmark_highDensity_halfPart_513nodes_10xtime_restart/')


CIC_RF_Benchmark_case1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_RF_Benchmark_case1/')
#vac convergenceDataExp = [Exp_128PPC_0p5Deb_0p2delT_2p5e16, NGP_2048PPC_64nodes_2p0delT_2p5e16]
# convergenceLabelExp = [r'Explicit', 'Implicit']

# compareModelToDatas(convergenceDataExp, convergenceLabelExp, modelOther, 'ICIS2023/converge')

# convergenceDataChangePPC = [Exp_32PPC_0p5Deb_0p2delT_2p5e16, Exp_128PPC_0p5Deb_0p2delT_2p5e16, NGP_128PPC_64nodes_2p0delT_2p5e16, NGP_2048PPC_64nodes_2p0delT_2p5e16]
# convergenceLabelChangePPC = [r'Exp., 32 PPC', 'Exp., 128 PPC', 'Imp., 128 PPC', 'Imp., 2048 PPC']

# compareModelToDatasRes(convergenceDataChangePPC, convergenceLabelChangePPC, modelOther)


# bitchPlease = [Exp_148PPD_0p5Deb_0p2delT_2p5e16, Exp_128PPC_0p5Deb_0p2delT_2p5e16_eps40]
# labelPlease = [r'Exp., 5 PPC, $\epsilon_r$ = $\epsilon_0$', r'Exp. 128 PPC, $\epsilon_r$ = 40$\epsilon_0$']
# compareModelToDatasRes(bitchPlease, labelPlease, modelOther)

# comparePPCData = [Exp_64PPC_0p5Deb_0p2delT_2p5e16, NGP_2048PPC_64nodes_2p0delT_2p5e16, NGP_512PPC_64nodes_2p0delT_2p5e16, NGP_64PPC_64nodes_2p0delT_2p5e16]
# comparePPCLabel = ['Exp., 64 PPC', 'Imp., 2048 PPC', 'Imp., 512 PPC', 'Imp., 64 PPC']

# compareRefToDatasRes(comparePPCData, comparePPCLabel, Exp_128PPC_0p15Deb_0p05delT_2p5e16)

# comparePPDData = [Exp_64PPC_0p5Deb_0p2delT_2p5e16, NGP_1902PPD_64nodes_2p0delT_2p5e16, NGP_951PPD_64nodes_2p0delT_2p5e16, 
#                   NGP_475PPD_64nodes_2p0delT_2p5e16]
# comparePPDLabel = ['Exp., 60864 NSP', 'Imp., 60864 NSP', 'Imp. 30432 NSP', 'Imp. 15200 NSP']
# compareRefToDatasAbsRes(comparePPDData, comparePPDLabel, Exp_128PPC_0p15Deb_0p05delT_2p5e16)
# compareRefToDatasRes(comparePPDData, comparePPDLabel, Exp_128PPC_0p15Deb_0p05delT_2p5e16, 'ICIS2023/refResNSP')

# comparePPDData = [Exp_148PPD_0p5Deb_0p2delT_2p5e16, NGP_100PPD_32nodes_2p0delT_2p5e16, NGP_75PPD_32nodes_2p0delT_2p5e16,
#                   NGP_50PPD_32nodes_2p0delT_2p5e16, CIC_50PPD_32nodes_2p0delT_2p5e16]
# comparePPDLabel = ['Exp., 4736 NSP', 'Imp., 3200 NSP', 'Imp., 2400 NSP', 'Imp., 1600 NSP', 'CIC, 1600 NSP']

# compareRefToDatasAbsRes(comparePPDData, comparePPDLabel, Exp_128PPC_0p15Deb_0p05delT_2p5e16)

# compareRefToDatasRes(comparePPDData, comparePPDLabel, Exp_128PPC_0p15Deb_0p05delT_2p5e16, 'ICIS2023/lowRes2p5e16')

# compareTimeData = [NGP_475PPD_64nodes_2p0delT_2p5e16, NGP_475PPD_64nodes_1p0delT_2p5e16, NGP_475PPD_64nodes_0p5delT_2p5e16]
# compareTimeLabel = [r'$\Delta t$ = $\frac{2}{\omega_p}$', r'$\Delta t$ = $\frac{1}{\omega_p}$', r'$\Delta t$ = $\frac{0.5}{\omega_p}$']

# compareRefToDatasRes(compareTimeData, compareTimeLabel, Exp_128PPC_0p15Deb_0p05delT_2p5e16)

# compareGwenaelData = [Exp_64PPC_0p5Deb_0p2delT_2p5e16, NGP_64PPC_64nodes_2p0delT_2p5e16, 
#                       NGP_951PPD_64nodes_2p0delT_2p5e16, NGP_475PPD_64nodes_2p0delT_2p5e16]
# compareGwenaelLabel = ['Exp., 60864 PPD',  'Imp., 64 PPC', 
#                       'Imp., 30432 PPD', 'Imp., 15200 PPD']
# compareRefToDatasRes(compareGwenaelData, compareGwenaelLabel, Exp_128PPC_0p15Deb_0p05delT_2p5e16)
# compareModelToDatas(compareGwenaelData, compareGwenaelLabel, model)