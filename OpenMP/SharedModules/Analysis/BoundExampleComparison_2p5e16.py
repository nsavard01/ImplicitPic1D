# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

    
    
#%%


model = np.load('BoundModels/model_2p5e+16nave_5Te_1Ti.npy')

Exp_128PPC_0p75Deb_0p2delT_2p5e16_CG = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p75Deb_0p2delT_CG/')
Exp_128PPC_1Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_1Deb_0p2delT/')
Exp_128PPC_2Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_2Deb_0p2delT/')
Exp_64PPC_2Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_64PPC_2Deb_0p2delT/')
Exp_128PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p75Deb_0p2delT/')
Exp_64PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_64PPC_0p75Deb_0p2delT/')
Exp_32PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_32PPC_0p75Deb_0p2delT/')
Exp_128PPC_0p75Deb_0p2delT_2p5e16_fullAve = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p75Deb_0p2delT_fullAve/')
Exp_256PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_256PPC_0p75Deb_0p2delT/')
Exp_512PPC_0p75Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_512PPC_0p75Deb_0p2delT/')
Exp_128PPC_0p5Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p5Deb_0p2delT/')
Exp_128PPC_0p25Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p25Deb_0p2delT/')
Exp_128PPC_0p1Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p1Deb_0p2delT/')
Exp_256PPC_0p1Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_256PPC_0p1Deb_0p2delT/')
Exp_128PPC_2Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_2p0Deb_0p2delT/')
Exp_128PPC_0p5Deb_0p1delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p5Deb_0p1delT/')
Exp_128PPC_0p25Deb_0p05delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p25Deb_0p05delT/')
Exp_128PPC_0p15Deb_0p05delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_0p15Deb_0p05delT/')

#NGP 2p5e16
testNGP = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/test/')
NGP_128PPC_100nodes_2p0delT_2p5e16_even = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_100nodes_evenGrid/')
NGP_128PPC_200nodes_2p0delT_2p5e16_even = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_200nodes_evenGrid/')
NGP_128PPC_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_64nodes/')
NGP_256PPC_64nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_256PPC_2p0delT_64nodes/')
NGP_1024PPC_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_32nodes/')
NGP_2048PPC_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_2048PPC_2p0delT_32nodes/')
NGP_1024PPC_32nodes_3p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_3p0delT_32nodes/')
NGP_1024PPC_32nodes_4p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_4p0delT_32nodes/')
NGP_512PPC_32nodes_4p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_512PPC_4p0delT_32nodes/')
NGP_128PPC_32nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_32nodes/')
NGP_128PPC_32nodes_4p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_4p0delT_32nodes/')
NGP_2048PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_2048PPC_2p0delT_16nodes/')
NGP_1024PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_16nodes/')
NGP_1024PPC_16nodes_2p0delT_2p5e16_1e5eps = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_16nodes_1e-5eps/')
NGP_1024PPC_16nodes_2p0delT_2p5e16_1e4eps = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_16nodes_1e-4eps/')
NGP_1024PPC_16nodes_2p0delT_2p5e16_1e3eps = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_16nodes_1e-3eps/')
NGP_512PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_512PPC_2p0delT_16nodes/')
NGP_256PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_256PPC_2p0delT_16nodes/')
NGP_128PPC_16nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_128PPC_2p0delT_16nodes/')
NGP_256PPC_11nodes_2p0delT_2p5e16 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_256PPC_2p0delT_11nodes/')
NGP_1024PPC_32nodes_2p0delT_2p5e16_CG = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_2p5e16_1024PPC_2p0delT_32nodes_CG/')




