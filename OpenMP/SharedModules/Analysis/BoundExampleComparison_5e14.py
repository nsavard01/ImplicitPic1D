# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

    
    
#%%

#5e14
model5e14 = np.load('BoundModels/model_5e+14nave_5Te_1Ti.npy')

test = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/test/')
Exp_128PPC_0p75Deb_0p2delT_5e14 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_5e14_128PPC_0p75Deb_0p2delT/')
Exp_1024PPC_0p75Deb_0p2delT_5e14 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_5e14_1024PPC_0p75Deb_0p2delT/')
Exp_1024PPC_0p5Deb_0p2delT_5e14 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_5e14_1024PPC_0p5Deb_0p2delT/')
Exp_1024PPC_0p1Deb_0p025delT_5e14 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_5e14_1024PPC_0p1Deb_0p025delT/')
test = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/test/')

NGP_1024PPC_32nodes_2p0delT_5e14 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_5e14_1024PPC_2p0delT_32nodes/')
NGP_1024PPC_16nodes_2p0delT_5e14 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_5e14_1024PPC_2p0delT_16nodes/')
NGP_2048PPC_16nodes_2p0delT_5e14 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_5e14_2048PPC_2p0delT_16nodes/')
NGP_1024PPC_32nodes_0p2delT_5e14 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_5e14_1024PPC_0p2delT_32nodes/')
NGP_128PPC_64nodes_2p0delT_5e14 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_5e14_128PPC_2p0delT_64nodes/')
NGP_128PPC_128nodes_2p0delT_5e14 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_5e14_128PPC_2p0delT_128nodes/')
testNGP = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/test/')

