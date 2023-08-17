# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

    
    
#%%

#5e14
model = getBoundPlasmaSolutions(0.05, 500, 1e18, 5, 1, m_p)

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
NGP_2048PPC_32nodes_2p0delT_1e18_CG = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_2048PPC_2p0delT_32nodes_CG/')



