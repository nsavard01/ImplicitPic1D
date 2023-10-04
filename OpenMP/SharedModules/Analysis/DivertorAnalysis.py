# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

    
    
#%%


test = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/test/')


#CIC 2p5e1

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