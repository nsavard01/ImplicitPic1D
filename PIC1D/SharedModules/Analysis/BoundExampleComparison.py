# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

    
    
#%%

model = np.load('BoundModels/modelNonTrunc_1p0e+18nave_5Te_1Ti.npy')

Exp_128PPC_0p75Deb_0p2delT_1e18_CG = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_128PPC_0p75Deb_0p2delT_1e18_CG/')
Exp_128PPC_0p75Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_128PPC_0p75Deb_0p2delT/')
Exp_128PPC_0p5Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_128PPC_0p5Deb_0p2delT/')
Exp_128PPC_0p1Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_128PPC_0p1Deb_0p2delT/')
Exp_16PPC_0p1Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_16PPC_0p1Deb_0p2delT/')
NGP_256PPC_120nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_256PPC_2p0delT_120nodes/')
NGP_1024PPC_64nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_1024PPC_2p0delT_64nodes/')
NGP_2048PPC_64nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_2048PPC_2p0delT_64nodes/')
NGP_1024PPC_16nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_1024PPC_2p0delT_16nodes/')
NGP_1024PPC_32nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_1024PPC_2p0delT_32nodes/')
NGP_1000PPD_32nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_1000PPD_2p0delT_32nodes/')
NGP_500PPD_32nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_500PPD_2p0delT_32nodes/')
NGP_250PPD_32nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_250PPD_2p0delT_32nodes/')
NGP_250PPD_32nodes_2p0delT_1e18_eps1e4 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_250PPD_2p0delT_32nodes_eps1e-4/')
NGP_125PPD_32nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_125PPD_2p0delT_32nodes/')
NGP_2048PPC_32nodes_2p0delT_1e18_CG = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_2048PPC_2p0delT_32nodes_CG/')

CIC_250PPD_32nodes_2p0delT_1e18_eps1e4 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_1e18_250PPD_2p0delT_32nodes_eps1e-4/')
CIC_250PPD_32nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/CIC_1e18_250PPD_2p0delT_32nodes/')

Exp_1000PPD_0p5Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_1000PPD_0p5Deb_0p2delT/')
Exp_1500PPD_0p5Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_1500PPD_0p5Deb_0p2delT/')
Exp_2000PPD_0p5Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_2000PPD_0p5Deb_0p2delT/')

# compareData = [Exp_2000PPD_0p5Deb_0p2delT_1e18, NGP_500PPD_32nodes_2p0delT_1e18, 
#                 NGP_250PPD_32nodes_2p0delT_1e18, NGP_250PPD_32nodes_2p0delT_1e18_eps1e4]
# compareLabel = ['Exp., 64000 NSP' , 'Imp., 16000 NSP', 'Imp., 8000 NSP', 'Imp., 8000 NSP, \n low res.']

# compareRefToDatasAbsRes(compareData, compareLabel, Exp_128PPC_0p5Deb_0p2delT_1e18)
#compareRefToDatasRes(compareData, compareLabel, Exp_128PPC_0p5Deb_0p2delT_1e18, 'ICIS2023/lowRes1e18')

compareDataCIC = [NGP_500PPD_32nodes_2p0delT_1e18, 
                NGP_250PPD_32nodes_2p0delT_1e18, CIC_250PPD_32nodes_2p0delT_1e18]
compareLabelCIC = ['NGP, 16000 NSP', 'NGP, 8000 NSP', 'CIC, 8000 NSP']

#compareRefToDatasAbsRes(compareData, compareLabel, Exp_128PPC_0p5Deb_0p2delT_1e18)
compareRefToDatasRes(compareDataCIC, compareLabelCIC, Exp_128PPC_0p5Deb_0p2delT_1e18, 'ICIS2023/lowResCIC1e18')

# compareDataGood = [NGP_2048PPC_32nodes_2p0delT_1e18_CG, Exp_128PPC_0p75Deb_0p2delT_1e18]
# compareLabelGood = ['Imp., 63488 PPD', 'Exp., 513408 PPD']
# compareRefToDatasRes(compareDataGood, compareLabelGood, Exp_128PPC_0p5Deb_0p2delT_1e18)


