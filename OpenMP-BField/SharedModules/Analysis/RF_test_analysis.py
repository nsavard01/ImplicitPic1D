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
Exp_RFBenchmark_test = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test/')
Exp_RFBenchmark_test_PnullExp = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_PnullExp/')
Exp_RFBenchmark_test_PnullExp_noRepeats = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_PnullExp_noRepeats/')
Exp_RFBenchmark_test_PnullExp_noRepeats_changeSeed = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_PnullExp_noRepeats_changeSeed/')
Exp_RFBenchmark_test_PnullNoExp_noRepeats = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_PnullNoExp_noRepeats/')
Exp_RFBenchmark_test_PnullExp_noRepeats_randNew = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_PnullExp_noRepeats_randNew/')
Exp_RFBenchmark_test_PnullExp_noRepeats_PCG = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_PnullExp_noRepeats_PCG/')
Exp_RFBenchmark_test_totalColl = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_totalColl/')
Exp_RFBenchmark_test_noRepeats = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_noRepeats/')
Exp_RFBenchmark_test_serial = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_serial/')
Exp_RFBenchmark_test_serial_higheraverage = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_serial_higheraverage/')
Exp_RFBenchmark_test_serial_doublenumax = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_serial_doublenumax/')
Exp_RFBenchmark_test_PnullChange_doublenumax = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_PnullChange/')
Exp_RFBenchmark_test_PnullChange = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_PnullChange_Regnumax/')

# case 2
Exp_RFBenchmark_test2 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test2/')
Exp_RFBenchmark_test2_PnullExp_noRepeats = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test2_PnullExp_noRepeats/')
Exp_RFBenchmark_test2_PnullExp_noRepeats_PCG = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test2_PnullExp_noRepeats_PCG/')
Exp_RFBenchmark_test2_totalColl = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test2_totalColl/')
Exp_RFBenchmark_test2_MT = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test2_MT/')
Exp_RFBenchmark_test3 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test3/')
Exp_RFBenchmark_test_MT = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test_MT/')

Exp_RFBenchmark_test3_PnullExp_noRepeats = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test3_PnullExp_noRepeats/')
Exp_RFBenchmark_test3_PnullExp_noRepeats_randNew = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test3_PnullExp_noRepeats_randNew/')

Exp_RFBenchmark_test4_PnullExp_noRepeats_randNew = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test4_PnullExp_noRepeats_randNew/')

NGP_test_RF = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_test_RF/')
NGP_test_typical = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_test_typical/')

CIC_test_RF = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/CIC_test_RF/')
CIC_test_typical = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/CIC_test_typical/')

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