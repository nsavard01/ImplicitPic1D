# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

def compareModelToData(data):
    for name in data.particles.keys():
        if name != 'e':
            ion = name
            break
    deb = debye_length(data.T_e, data.n_ave)
    M = data.particles[ion]['mass']
    model = getBoundPlasmaSolutions(data.grid[-1] - data.grid[0], 300, data.n_ave, data.T_e, data.T_i, M)
    plt.figure()
    plotAveDensity(data, name = 'e', label = 'PIC')
    plt.plot(model[0], model[2], label = 'Model')
    plt.legend(loc = 'best')
    
    plt.figure()
    plotAveDensity(data, name = ion, label = 'PIC')
    plt.plot(model[0], model[3], label = 'Model')
    plt.legend(loc = 'best')
    
    plt.figure()
    plotAvePhi(data)
    plt.plot(model[0], model[1], label = 'Model')
    plt.legend(loc = 'best')
    
#%%

#5e14
Exp_128PPC_0p75Deb_0p2delT_5e14 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_5e14_128PPC_0p75Deb_0p2delT/')
Exp_1024PPC_0p75Deb_0p2delT_5e14 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_5e14_1024PPC_0p75Deb_0p2delT/')
Exp_1024PPC_0p5Deb_0p2delT_5e14 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_5e14_1024PPC_0p5Deb_0p2delT/')
Exp_1024PPC_0p1Deb_0p025delT_5e14 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_5e14_1024PPC_0p1Deb_0p025delT/')

NGP_1024PPC_32nodes_2p0delT_5e14 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_5e14_1024PPC_2p0delT_32nodes/')
NGP_1024PPC_32nodes_0p2delT_5e14 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_5e14_1024PPC_0p2delT_32nodes/')

Exp_128PPC_1Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_1Deb_0p2delT/')
Exp_128PPC_2Deb_0p2delT_2p5e16 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_2p5e16_128PPC_2Deb_0p2delT/')
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

Exp_128PPC_0p75Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_128PPC_0p75Deb_0p2delT/')
Exp_128PPC_0p5Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_128PPC_0p5Deb_0p2delT/')
Exp_128PPC_0p1Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_128PPC_0p1Deb_0p2delT/')
Exp_16PPC_0p1Deb_0p2delT_1e18 = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/Explicit_1e18_16PPC_0p1Deb_0p2delT/')
NGP_256PPC_120nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_256PPC_2p0delT_120nodes/')
NGP_1024PPC_64nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_1024PPC_2p0delT_64nodes/')
NGP_2048PPC_64nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_2048PPC_2p0delT_64nodes/')
NGP_1024PPC_16nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_1024PPC_2p0delT_16nodes/')
NGP_1024PPC_32nodes_2p0delT_1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitData/NGP_1e18_1024PPC_2p0delT_32nodes/')
test = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ExplicitData/test/')
expListPPC = [Exp_128PPC_0p75Deb_0p2delT_2p5e16, Exp_256PPC_0p75Deb_0p2delT_2p5e16,
                  Exp_512PPC_0p75Deb_0p2delT_2p5e16]
labelListPPC = [r'128 PPC', '256 PPC', '512 PPC']


dataList2p5e16 = [Exp_128PPC_0p75Deb_0p2delT_2p5e16, Exp_256PPC_0p75Deb_0p2delT_2p5e16,
                  Exp_512PPC_0p75Deb_0p2delT_2p5e16]
labelList2p5e16 = [r'128 PPC', '256 PPC', '512 PPC']

list1 = [Exp_128PPC_0p75Deb_0p2delT_2p5e16, Exp_512PPC_0p75Deb_0p2delT_2p5e16,
         Exp_128PPC_0p5Deb_0p2delT_2p5e16, Exp_128PPC_0p15Deb_0p05delT_2p5e16]
list2 = ['128 PPC, 0.2 $\Delta$ t, 0.75 $\lambda_{De}$', '512 PPC, 0.2 $\Delta$ t, 0.75 $\lambda_{De}$', 
         r'128 PPC, 0.2 $\Delta$ t, 0.5 $\lambda_{De}$', r'128 PPC, 0.05 $\Delta$ t, 0.15 $\lambda_{De}$']

NGPList = [Exp_128PPC_0p25Deb_0p05delT_2p5e16, NGP_128PPC_32nodes_2p0delT_2p5e16, 
           NGP_128PPC_32nodes_4p0delT_2p5e16, NGP_1024PPC_32nodes_2p0delT_2p5e16]
NGPlabel = ['Exp. Ref.', r'128 PPC, 2.0 $\Delta$t', r'128 PPC, 4.0 $\Delta$t'
           , r'1024 PPC, 2.0 $\Delta$t']

dataList1e18 = [Exp_128PPC_0p75Deb_0p2delT_1e18, Exp_128PPC_0p5Deb_0p2delT_1e18,
                  NGP_1024PPC_64nodes_2p0delT_1e18, NGP_2048PPC_64nodes_2p0delT_1e18]
labelList1e18 = [r'Exp., 0.2 $\Delta$t, 128PPC, 0.75 $\lambda_{De}$', r'Exp., 0.2 $\Delta$t, 128PPC, 0.5 $\lambda_{De}$',
                   r'NGP, 1024 PPC', r'NGP, 2048 PPC']

NGPListHighPPC = [NGP_1024PPC_32nodes_2p0delT_2p5e16, NGP_2048PPC_32nodes_2p0delT_2p5e16]
NGPListHighLabel = ['1024 PPC', '2048 PPC']

LowPPC = [NGP_128PPC_16nodes_2p0delT_2p5e16, NGP_256PPC_16nodes_2p0delT_2p5e16, NGP_512PPC_16nodes_2p0delT_2p5e16,
          Exp_64PPC_0p75Deb_0p2delT_2p5e16, Exp_32PPC_0p75Deb_0p2delT_2p5e16]
LowPPCLabel = ['NGP, 128 PPC', 'NGP, 256 PPC', 'NGP, 512 PPC', 'Exp., 64 PPC', 'Exp., 32 PPC']


