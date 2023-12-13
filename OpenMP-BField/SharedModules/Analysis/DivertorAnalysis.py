# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateDivertorExample import *
print('Loaded modules')
    
#%%

print('Loading data')

test_noField_uniformDomain = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_noField_uniformDomain/')
testExp_0deg = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test/')
testExp_0deg_Mach2 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_0deg_Mach2/')
testExp_45deg = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_45deg/')
testExp_45deg_eps100 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_45deg_eps100/')
testExp_75deg = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_75deg/')
testExp_extend_0p1delT_75deg = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_75deg_extend_0p1delT/')
testExp_extend_0p2delT_75deg = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_75deg_extend_0p2delT/')
testExp_extend_0p04delT_75deg = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_75deg_extend_0p04delT/')
test_45deg_fourTimesIonSpeed = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_45deg_fourTimesIonSpeed/')
testExp_norm_0p05T = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/test_0p05delT/')
test_45deg_uniformDomain = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/test_45deg_uniformDomain/')

testNGP = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP/')
testNGP_noField_uniformDomain_32nodes_2p0T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_noField_uniformDomain_32nodes_2p0T/')
testNGP_45deg_uniformDomain_32nodes_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_45deg_uniformDomain_32nodes_0p2T/')
testNGP_45deg_uniformDomain_400nodes_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_45deg_uniformDomain_400nodes_0p2T/')
testNGP_45deg_uniformDomain = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_45deg_uniformDomain/')
test_0deg_picard = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_0deg_picard/')
test_0deg = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/test_0deg_other/')
test_0deg_vacuumStart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/test_0deg_try/')
test_0deg_Mach2 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_0deg_Mach2/')
test_0deg_Mach2_sin = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_0deg_Mach2_sin/')
test_0deg_vacuumStart_noExactBoundary = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/test_0deg_vacuumStart_notExactBoundary/')
test_45deg_vacuumStart = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/test_45deg_vacuumStart/')
test_45deg_vacuumStart_noExactBoundary = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/test_45deg_vacuumStart_notExactBoundary/')
vacuumStart_45deg_PI_32Thread = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_45deg_PI_32thread/')
vacuumStart_45deg_AA_32Thread = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_45deg_AA_32thread/')
test_45deg_vacuumStart_noExactBoundary_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_45deg_0p2T/')
vacuumStart_45deg_noExactBoundary_400nodes_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_45deg_notExactBoundary_400nodes_0p2T/')
vacuumStart_75deg_noExactBoundary_400nodes_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_75deg_400nodes_0p2T/')
vacuumStart_75deg_noExactBoundary_433nodes_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_75deg_433nodes_0p2T/')


#CIC
testCIC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testCIC/')
testCIC_45deg_uniformDomain = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testCIC_45deg_uniformDomain/')
vacuumStart_45deg_CIC_2p0T_31nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_45deg_CIC_2p0T_31nodes/')
testCIC_45deg_CIC_0p2T_31nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testCIC_45deg_0p2delT/')
vacuumStart_45deg_CIC_0p2T_400nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/vacuumStart_45deg_CIC_0p2T_400nodes/')
testCIC_45deg_fourTimesIonSpeed = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testCIC_45deg_fourTimesIonSpeed/')
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