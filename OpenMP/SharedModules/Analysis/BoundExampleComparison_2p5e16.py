# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

    
    
#%%


model = np.load('BoundModels/modelNonTruncSparse_2p5e+16nave_5Te_1Ti.npy')
modelOther = np.load('BoundModels/modelNonTrunc_2p5e+16nave_5Te_1Ti.npy')

NGP_BoundModel_test = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_BoundModel_test/')

# Exp_1000PPD_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1000PPD_0p5Deb_0p2delT_2p5e16/')
# Exp_500PPD_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_500PPD_0p5Deb_0p2delT/')
# Exp_300PPD_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_300PPD_0p5Deb_0p2delT/')
# Exp_200PPD_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_200PPD_0p5Deb_0p2delT/')
# Exp_148PPD_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_148PPD_0p5Deb_0p2delT/')
# Exp_128PPC_0p5Deb_0p2delT_2p5e16_eps40 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_600PPD_0p5Deb_0p2delT_eps40/')
# Exp_100PPD_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_100PPD_0p5Deb_0p2delT/')
# Exp_128PPC_0p75Deb_0p2delT_2p5e16_CG = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p75Deb_0p2delT_CG/')
# Exp_128PPC_1Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_1Deb_0p2delT/')
# Exp_128PPC_2Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_2Deb_0p2delT/')
# Exp_64PPC_2Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_64PPC_2Deb_0p2delT/')
# Exp_128PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p75Deb_0p2delT/')
# Exp_64PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_64PPC_0p75Deb_0p2delT/')
# Exp_32PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_32PPC_0p75Deb_0p2delT/')
# Exp_128PPC_0p75Deb_0p2delT_2p5e16_fullAve = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p75Deb_0p2delT_fullAve/')
# Exp_256PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_256PPC_0p75Deb_0p2delT/')
# Exp_512PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_512PPC_0p75Deb_0p2delT/')
# Exp_128PPC_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p5Deb_0p2delT/')
# Exp_128PPC_0p25Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p25Deb_0p2delT/')
# Exp_128PPC_0p1Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p1Deb_0p2delT/')
# Exp_256PPC_0p1Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_256PPC_0p1Deb_0p2delT/')
# Exp_128PPC_2Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_2p0Deb_0p2delT/')
# Exp_128PPC_0p5Deb_0p1delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p5Deb_0p1delT/')
# Exp_128PPC_0p25Deb_0p05delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p25Deb_0p05delT/')
# Exp_128PPC_0p15Deb_0p05delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p15Deb_0p05delT/')
# Exp_64PPC_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_64PPC_0p5Deb_0p2delT/')
# Exp_32PPC_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_32PPC_0p5Deb_0p2delT/')

#NGP 2p5e16
# testAgain_boundExample = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/testAgain_boundExample/')
# NGP_1000PPD_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1000PPD_2p0delT_32nodes/')
# NGP_500PPD_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_500PPD_2p0delT_32nodes/')
# NGP_100PPD_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_100PPD_2p0delT_32nodes/')
# NGP_75PPD_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_75PPD_2p0delT_32nodes/')
# NGP_50PPD_32nodes_2p0delT_2p5e16_1e3eps = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_50PPD_2p0delT_32nodes_1e3eps/')
# NGP_50PPD_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_50PPD_2p0delT_32nodes/')
# NGP_1902PPD_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1902PPD_2p0delT_64nodes/')
# NGP_1902PPD_64nodes_2p0delT_2p5e16_saveField = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1902PPD_64nodes_2p0delT_2p5e16_saveEField/')
# NGP_1902PPD_64nodes_2p0delT_2p5e16_even = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1902PPD_2p0delT_64nodes_even/')
# NGP_951PPD_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_951PPD_2p0delT_64nodes/')
# NGP_475PPD_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_475PPD_2p0delT_64nodes/')
# NGP_475PPD_64nodes_1p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_475PPD_1p0delT_64nodes/')
# NGP_475PPD_64nodes_0p5delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_475PPD_0p5delT_64nodes/')
# NGP_128PPC_100nodes_2p0delT_2p5e16_even = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_100nodes_evenGrid/')
# NGP_128PPC_200nodes_2p0delT_2p5e16_even = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_200nodes_evenGrid/')
# NGP_128PPC_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_64nodes/')
# NGP_256PPC_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_256PPC_2p0delT_64nodes/')
# NGP_512PPC_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_512PPC_2p0delT_64nodes/')
# NGP_1024PPC_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_64nodes/')
# NGP_2048PPC_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_2048PPC_2p0delT_64nodes/')
# NGP_1024PPC_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_32nodes/')
# NGP_2048PPC_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_2048PPC_2p0delT_32nodes/')
# NGP_1024PPC_32nodes_3p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_3p0delT_32nodes/')
# NGP_1024PPC_32nodes_4p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_4p0delT_32nodes/')
# NGP_512PPC_32nodes_4p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_512PPC_4p0delT_32nodes/')
# NGP_128PPC_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_32nodes/')
# NGP_128PPC_32nodes_4p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_4p0delT_32nodes/')
# NGP_2048PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_2048PPC_2p0delT_16nodes/')
# NGP_4096PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_4096PPC_2p0delT_16nodes/')
# NGP_1024PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_16nodes/')
# NGP_1024PPC_16nodes_2p0delT_2p5e16_1e5eps = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_16nodes_1e-5eps/')
# NGP_1024PPC_16nodes_2p0delT_2p5e16_1e4eps = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_16nodes_1e-4eps/')
# NGP_1024PPC_16nodes_2p0delT_2p5e16_1e3eps = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_16nodes_1e-3eps/')
# NGP_512PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_512PPC_2p0delT_16nodes/')
# NGP_256PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_256PPC_2p0delT_16nodes/')
# NGP_128PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_16nodes/')
# NGP_256PPC_11nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_256PPC_2p0delT_11nodes/')
# NGP_64PPC_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_64PPC_2p0delT_32nodes/')
# NGP_64PPC_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_64PPC_2p0delT_64nodes/')
# NGP_1024PPC_32nodes_2p0delT_2p5e16_CG = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_32nodes_CG/')
# NGP_1902PPD_32nodes_2p0delT_2p5e16_noFieldTest = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_noField_boundExample_32nodes_2p0T/')
# NGP_1902PPD_64nodes_2p0delT_2p5e16_noFieldTest = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_noField_boundExample_64nodes_2p0T/')
# testNGP_noField_boundExample_64nodes_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_noField_boundExample_64nodes_0p2T/')
# testNGP_noField_boundExample_32nodes_2p0T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_noField_boundExample_32nodes_2p0T/')
# testNGP_noField_boundExample_32nodes_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_noField_boundExample_32nodes_0p2T/')
# testNGP_noFieldAlg_boundExample_32nodes_2p0T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testNGP_noFieldAlg_boundExample_32nodes_2p0T/')

#CIC 2p5e16
# CIC_475PPD_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_2p5e16_475PPD_2p0delT_64nodes/')
# CIC_50PPD_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_2p5e16_50PPD_2p0delT_32nodes/')
# testCIC_noField_boundExample_64nodes_2p0T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testCIC_noField_boundExample_64nodes_2p0T/')
# testCIC_noField_boundExample_64nodes_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/testCIC_noField_boundExample_64nodes_0p2T/')
# convergenceDataExp = [Exp_128PPC_0p5Deb_0p2delT_2p5e16, NGP_2048PPC_64nodes_2p0delT_2p5e16]
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