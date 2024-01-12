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
model = np.load('BoundModels/modelNonTruncSparse_2p5e+16nave_5Te_1Ti.npy')

#--------------------- Explicit ---------------------------------------------
Exp_noField = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField/')
Exp_45deg_0p05T = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField_45deg_0p05T/')
Exp_noField_keepBoris = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField_keepBoris/')
Exp_45deg_0p05T_3000PPD = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_45deg_0p05T_3000PPD/')
Exp_45deg_0p05T_3000PPD_invertV = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_45deg_0p05T_3000PPD_invertV/')
Exp_45deg_0p05T_reverseV = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_45deg_0p05T_reverseV/')

Exp_0p05T_fullDomain = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_0p05T_fullDomain/')
Exp_noField_fullDomain = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField_fullDomain/')
Exp_noField_halfDomain = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField_halfDomain/')
Exp_noField_fullDomain_injFlux = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField_fullDomain_injFlux/')
Exp_noField_fullDomain_heatingMax1 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField_fullDomain_heatingMax1/')
Exp_noField_fullDomain_heatingMax100 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField_fullDomain_heatingMax100/')
Exp_noField_fullDomain_heatingMax100_0p1delT = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_noField_fullDomain_heatingMax100_0p1delT/')

#---------------------- NGP -------------------------------------------------
NGP_noField_halfDomain_retest = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_halfDomain_retest/')
NGP_noField = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField/')
NGP_noFieldVStop = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noFieldVStop/')
NGP_45deg_0p05T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T/')
NGP_45deg_0p05T_reverseV_32nodes_0p2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_reverseV_32nodes_0p2delT/')
NGP_45deg_0p05T_reverseV_32nodes_2p0delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_reverseV_32nodes_2p0delT/')
NGP_45deg_0p05T_reverseV_32nodes_2p0delT_3000PPD = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_reverseV_32nodes_2p0delT_3000PPD/')
NGP_45deg_0p05T_reverseV_32nodes_2p0delT_step0p1w_c = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_reverseV_32nodes_2p0delT_step0p1w_c/')
NGP_45deg_0p05T_reverseV_32nodes_2p0delT_0p1w_cSubSteps = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_reverseV_32nodes_2p0delT_0p1w_cSubSteps/')
NGP_45deg_0p05T_0p2T = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_0p2T/')
NGP_45deg_0p05T_2p0T_minDelTau = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_2p0T_minDelTau/')
NGP_45deg_0p05T_2p0delT_secantMethod = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_2p0delT_secantMethod/')
NGP_45deg_0p05T_2p0delT_secantMethod2 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_2p0delT_secantMethod2/')
NGP_45deg_0p05T_2p0delT_secantMethod3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_2p0delT_secantMethod3/')

NGP_noField_fullDomain_heatingMax100_0p5delT_400nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_heatingMax100_0p5delT_400nodes/')
NGP_noField_fullDomain_200nodes_heatingMax100_2p0delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_200nodes_heatingMax100_2p0delT/')
NGP_noField_fullDomain_heatingMax100_2p0delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_heatingMax100_2p0delT/')
NGP_noField_fullDomain_heatingMax10_0p2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_heatingMax10_0p2delT/')
NGP_noField_fullDomain_heatingMax10 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_heatingMax10/')
NGP_noField_fullDomain_injFlux_2p0delT_VStop = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_injFlux_2p0delT_VStop/')
NGP_noField_fullDomain_injFlux_2p0delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_injFlux_2p0delT/')
NGP_noField_fullDomain_injFlux_0p2delT = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_injFlux_0p2delT/')
NGP_noField_fullDomain_injFlux = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_fullDomain_injFlux/')
NGP_noField_halfDomain = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_noField_halfDomain/')
NGP_45deg_noField_2p0delT_fullDomain = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_noField_2p0delT_fullDomain/')
NGP_45deg_0p05T_2p0delT_secantMethod_fullDomain = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData-BField/NGP_45deg_0p05T_2p0delT_secantMethod_fullDomain/')
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