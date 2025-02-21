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


# Exp_ICP_test_lowDensity = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_lowDensity/')
# Exp_ICP_test_lowDensity_implicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_lowDensity_implicit/')
# Exp_ICP_test_lowDensity_fullyImplicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_lowDensity_fullyImplicit/')
#
# Exp_ICP_test_Meige = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_Meige/')
# Exp_ICP_test_Meige_quarterHeat = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_Meige_quarterHeat/')
#
# Exp_ICP_test_highDensity = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_implicit/')
# Exp_ICP_test_highDensity_regRF = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_regRF/')
# Exp_ICP_test_highDensity_regRF_6400delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_regRF_6400delT/')
# Exp_ICP_test_highDensity_regRF_6400delT_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_regRF_6400delT_doublePart/')
#
# Exp_ICP_test_highDensity_6400delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT/')
# Exp_ICP_test_highDensity_6400delT_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_doublePart/')
# Exp_ICP_test_highDensity_6400delT_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_quadPart/')
# Exp_ICP_test_highDensity_6400delT_x8Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_x8Part/')
# Exp_ICP_test_highDensity_6400delT_x16Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_x16Part/')
# Exp_ICP_test_highDensity_6400delT_x32Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_6400delT_x32Part/')

Exp_ICP_highDensity_highPressure_standRF_1300cells_4000delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_1300cells_4000delT_64PPC/')

Exp_ICP_highDensity_highPressure_standRF_6400delT_100PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_6400delT_100PPC/')
Exp_ICP_highDensity_highPressure_standRF_6400delT_200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_6400delT_200PPC/')
Exp_ICP_highDensity_highPressure_standRF_6400delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_6400delT_400PPC/')
Exp_ICP_highDensity_highPressure_standRF_6400delT_800PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_6400delT_800PPC/')
Exp_ICP_highDensity_highPressure_standRF_6400delT_1600PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_6400delT_1600PPC/')
Exp_ICP_highDensity_highPressure_standRF_6400delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_6400delT_3200PPC/')
Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC/')

NGP_ICP_highPressure_200Cells_uniform_2500delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_100PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_200PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_400PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_800PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_1600PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_3200PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_6400PPC/')

NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4/')

NGP_ICP_highPressure_128Cells_sinusoid_2500delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_100PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_200PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_800PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_12800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_12800PPC/')

NGP_ICP_highPressure_128Cells_sinusoid_2500delT_12800PPC_try2 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_12800PPC_try2/')

Exp_ICP_highDensity_highPressure_standRF_1500cells_5000delT_128PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_1500cells_5000delT_128PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC_try2 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC_try2/')

CIC_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4/')
# Exp_ICP_test_highDensity_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_doublePart/')
# Exp_ICP_test_highDensity_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_quadPart/')
#
# Exp_ICP_test_highDensity_lowPressure_6400delT_8mmHeating = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_8mmHeating/')
# Exp_ICP_test_highDensity_lowPressure_6400delT_8mmHeating_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_8mmHeating_doublePart/')
#
# Exp_ICP_test_highDensity_lowPressure_6400delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT/')
# Exp_ICP_test_highDensity_lowPressure_12800delT = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_12800delT/')
# Exp_ICP_test_highDensity_lowPressure_6400delT_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_doublePart/')
# Exp_ICP_test_highDensity_lowPressure_6400delT_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_quadPart/')
# Exp_ICP_test_highDensity_lowPressure_6400delT_x8Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_x8Part/')
# Exp_ICP_test_highDensity_lowPressure_6400delT_x16Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_6400delT_x16Part/')
#
#
# Exp_ICP_test_highDensity_lowPressure_x16Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_x16Part/')
# Exp_ICP_test_highDensity_lowPressure_x8Part = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_x8Part/')
# Exp_ICP_test_highDensity_lowPressure_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_quadPart/')
# Exp_ICP_test_highDensity_lowPressure_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure_doublePart/')
# Exp_ICP_test_highDensity_lowPressure = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_test_highDensity_lowPressure/')


plt.figure()
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_200PPC, 'He+', label = '200 PPC')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_400PPC, 'He+', label = '400 PPC')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_800PPC, 'He+', label = '800 PPC')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_1600PPC, 'He+', label = '1600 PPC')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_3200PPC, 'He+', label = '3200 PPC')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, 'He+', label = '3200 PPC, more res.')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/Exp_ICP_highPress.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_200PPC, label = '200 PPC')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_400PPC, label = '400 PPC')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_800PPC, label = '800 PPC')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_1600PPC, label = '1600 PPC')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_3200PPC, label = '3200 PPC')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, label = '3200 PPC, more res.')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/Exp_ICP_highPress_phi.png')

plt.figure()
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, 'He+', label = 'Exp High Res.')
plotAveDensity(NGP_ICP_highPressure_200Cells_uniform_2500delT_200PPC, 'He+', label = 'NGP 200 PPC')
plotAveDensity(NGP_ICP_highPressure_200Cells_uniform_2500delT_400PPC, 'He+', label = 'NGP 400 PPC')
plotAveDensity(NGP_ICP_highPressure_200Cells_uniform_2500delT_800PPC, 'He+', label = 'NGP 800 PPC')
plotAveDensity(NGP_ICP_highPressure_200Cells_uniform_2500delT_1600PPC, 'He+', label = 'NGP 1600 PPC')
plotAveDensity(NGP_ICP_highPressure_200Cells_uniform_2500delT_3200PPC, 'He+', label = 'NGP 3200 PPC')
plotAveDensity(NGP_ICP_highPressure_200Cells_uniform_2500delT_6400PPC, 'He+', label = 'NGP 6400 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_highPress_uniform.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, label = 'Exp High Res.')
plotAvePhi(NGP_ICP_highPressure_200Cells_uniform_2500delT_200PPC, label = 'NGP 200 PPC')
plotAvePhi(NGP_ICP_highPressure_200Cells_uniform_2500delT_400PPC, label = 'NGP 400 PPC')
plotAvePhi(NGP_ICP_highPressure_200Cells_uniform_2500delT_800PPC, label = 'NGP 800 PPC')
plotAvePhi(NGP_ICP_highPressure_200Cells_uniform_2500delT_1600PPC, label = 'NGP 1600 PPC')
plotAvePhi(NGP_ICP_highPressure_200Cells_uniform_2500delT_3200PPC, label = 'NGP 3200 PPC')
plotAvePhi(NGP_ICP_highPressure_200Cells_uniform_2500delT_6400PPC, label = 'NGP 6400 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_highPress_phi_uniform.png')


plt.figure()
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, label = 'Exp High Res.')
# plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_1300cells_4000delT_64PPC, label = 'Exp Low Res.')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC, label = 'NGP 400 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_800PPC, label = 'NGP 800 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC, label = 'NGP 1600 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC, label = 'NGP 3200 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC, label = 'NGP 6400 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_12800PPC, label = 'NGP 12800 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_highPress_phi.png')

plt.figure()
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, 'He+', label = 'Exp High Res.')
# plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_1300cells_4000delT_64PPC, 'He+', label = 'Exp Low Res.')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC, 'He+',label = 'NGP 400 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_800PPC, 'He+',label = 'NGP 800 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC, 'He+',label = 'NGP 1600 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC, 'He+',label = 'NGP 3200 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC, 'He+',label = 'NGP 6400 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_12800PPC, 'He+',label = 'NGP 12800 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_highPress.png')


# ----------------------------- low pressure ---------------------------

Exp_ICP_highDensity_lowPressure_standRF_7000delT_100PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_7000delT_100PPC/')
Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC/')
Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC/')
Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC/')
Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC/')
Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC/')
Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC/')

Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC/')
Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC/')
Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC_implicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC_implicit/')
# Exp_ICP_highDensity_lowPressure = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure/')
#
#
# NGP_ICP_test_lowDensity_fullyImplicit = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_test_lowDensity_fullyImplicit/')

NGP_ICP_lowPressure_2000Cells_sinusoid_420delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2000Cells_sinusoid_420delT_3200PPC/')
NGP_ICP_lowPressure_2000Cells_uniform_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2000Cells_uniform_640delT_3200PPC/')

NGP_ICP_lowPressure_128Cells_DenisMap_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_DenisMap_640delT_100PPC/')

NGP_ICP_lowPressure_200Cells_uniform_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_100PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_200PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_400PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_800PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_1600PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_3200PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_6400PPC/')

NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC/')
NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4/')

NGP_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC/')

NGP_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC_smoothing/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC_smoothing/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC_smoothing/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC_smoothing/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC_smoothing/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC_smoothing/')

CIC_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4/')
CIC_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC/')

CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC/')


plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC, 'He+', label = '200 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC, 'He+', label = '400 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC, 'He+', label = '800 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC, 'He+', label = '1600 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC, 'He+', label = '3200 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = '3200 PPC, more res.')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/Exp_ICP_lowPress.png')


plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC, label = '200 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC, label = '400 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC, label = '800 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC, label = '1600 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC, label = '3200 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = '3200 PPC, more res.')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/Exp_ICP_lowPress_phi.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = 'Exp High Res.')
# plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC, label = 'Exp Low Res.')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC, label = 'NGP 200 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, label = 'NGP 400 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC, label = 'NGP 800 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC, label = 'NGP 1600 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC, label = 'NGP 3200 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, label = 'NGP 6400 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_lowPress_phi.png')

plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = 'Exp High Res.')
# plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC, 'He+', label = 'Exp Low Res.')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC, 'He+',label = 'NGP 200 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, 'He+',label = 'NGP 400 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC, 'He+',label = 'NGP 800 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC, 'He+',label = 'NGP 1600 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC, 'He+',label = 'NGP 3200 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, 'He+',label = 'NGP 6400 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_lowPress.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = 'Exp High Res.')
plotAvePhi(NGP_ICP_lowPressure_200Cells_uniform_640delT_200PPC, label = 'NGP 200 PPC')
plotAvePhi(NGP_ICP_lowPressure_200Cells_uniform_640delT_400PPC, label = 'NGP 400 PPC')
plotAvePhi(NGP_ICP_lowPressure_200Cells_uniform_640delT_800PPC, label = 'NGP 800 PPC')
plotAvePhi(NGP_ICP_lowPressure_200Cells_uniform_640delT_1600PPC, label = 'NGP 1600 PPC')
plotAvePhi(NGP_ICP_lowPressure_200Cells_uniform_640delT_3200PPC, label = 'NGP 3200 PPC')
plotAvePhi(NGP_ICP_lowPressure_200Cells_uniform_640delT_6400PPC, label = 'NGP 6400 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_lowPress_phi_uniform.png')

plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = 'Exp High Res.')
plotAveDensity(NGP_ICP_lowPressure_200Cells_uniform_640delT_200PPC, 'He+',label = 'NGP 200 PPC')
plotAveDensity(NGP_ICP_lowPressure_200Cells_uniform_640delT_400PPC, 'He+',label = 'NGP 400 PPC')
plotAveDensity(NGP_ICP_lowPressure_200Cells_uniform_640delT_800PPC, 'He+',label = 'NGP 800 PPC')
plotAveDensity(NGP_ICP_lowPressure_200Cells_uniform_640delT_1600PPC, 'He+',label = 'NGP 1600 PPC')
plotAveDensity(NGP_ICP_lowPressure_200Cells_uniform_640delT_3200PPC, 'He+',label = 'NGP 3200 PPC')
plotAveDensity(NGP_ICP_lowPressure_200Cells_uniform_640delT_6400PPC, 'He+',label = 'NGP 6400 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_lowPress_density_uniform.png')


# -------------------- Trends between pressures  --------------------------------------------
plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = 'Exp High Res. low press')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, 'He+', label = 'Exp High Res. high press')
plotAveDensity(NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4, 'He+', label = 'NGP low press 64 cell 64 PPC')
plotAveDensity(NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4, 'He+', label = 'NGP high press 64 cell 64 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, 'He+', label = 'NGP low press 128 cell 6400 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/ICP_pressComp_density.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = 'Exp Res. low press')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, label = 'Exp Res. high press')
plotAvePhi(NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4, label = 'NGP low press 64 PPC')
plotAvePhi(NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4, label = 'NGP high press 64 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, label = 'NGP low press 6400 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/ICP_pressComp_phi.png')

data_points = np.array([200,400,800,1600,3200])

dataset_list_lowPress_Exp = [Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC, Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC,
                             Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC, Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC,
                             Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC]

dataset_list_lowPress_NGP_128cells = [NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC, NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC,
                             NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC, NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC,
                             NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC]

dataset_list_highPress_Exp = [Exp_ICP_highDensity_highPressure_standRF_6400delT_200PPC, Exp_ICP_highDensity_highPressure_standRF_6400delT_400PPC,
                             Exp_ICP_highDensity_highPressure_standRF_6400delT_800PPC, Exp_ICP_highDensity_highPressure_standRF_6400delT_1600PPC,
                             Exp_ICP_highDensity_highPressure_standRF_6400delT_3200PPC]


dataset_list_highPress_NGP_128cells = [NGP_ICP_highPressure_128Cells_sinusoid_2500delT_200PPC, NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC,
                             NGP_ICP_highPressure_128Cells_sinusoid_2500delT_800PPC, NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC,
                             NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC]


data_list_labels = ['Exp. low press', 'NGP low press 128 cells', 'Exp. high press', 'NGP high press 128 cells']


plt.figure()
plotWallCurrent(data_points, 'He+', dataset_list_lowPress_Exp, label = 'Exp. low press')
plotWallCurrent(data_points, 'He+', dataset_list_highPress_Exp, label = 'Exp. high press')
plotWallCurrent(data_points, 'He+', dataset_list_lowPress_NGP_128cells, label = 'NGP low press 128 cells')
plotWallCurrent(data_points, 'He+', dataset_list_highPress_NGP_128cells, label = 'NGP. high press 128 cells')
plotWallCurrent(np.array([64]), 'He+', [NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4], label = 'NGP. low press 64 cells')
plotWallCurrent(np.array([64]), 'He+', [NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4], label = 'NGP. high press 64 cells')
plt.legend(loc = 'best')
plt.ylabel(r'He$^+$ wall current ($\frac{A}{m^2}$)')
plt.xlabel(r'Particle-per-cell')
plt.savefig('RF_ICP/ICP_pressComp_current.png')

# ------------------------------------- 2/3 domain


Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_5000delT_128PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_5000delT_128PPC/')
Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_3200PPC/')
Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_1600PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_1600PPC/')
Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_800PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_800PPC/')
Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_400PPC/')

plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_400PPC, 'He+', label = '400 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_800PPC, 'He+', label = '800 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_1600PPC, 'He+', label = '1600 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_3200PPC, 'He+', label = '3200 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/Exp_ICP_lowPress_small.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_400PPC,  label = '400 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_800PPC,  label = '800 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_1600PPC, label = '1600 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_3200PPC,  label = '3200 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/Exp_ICP_lowPress_phi_small.png')

NGP_ICP_lowPressure_2thirdsDomain_128Cells_sinusoid_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_128Cells_sinusoid_640delT_100PPC/')
NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_100PPC_explicit = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_100PPC_explicit/')

NGP_ICP_lowPressure_2thirdsDomain_64Cells_sinusoid_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_64Cells_sinusoid_640delT_100PPC/')
NGP_ICP_lowPressure_2thirdsDomain_32Cells_sinusoid_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_32Cells_sinusoid_640delT_100PPC/')
NGP_ICP_lowPressure_2thirdsDomain_32Cells_sinusoid_640delT_50PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_32Cells_sinusoid_640delT_50PPC/')

NGP_ICP_lowPressure_2thirdsDomain_48Cells_sinusoid_640delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_48Cells_sinusoid_640delT_64PPC_epsneg4/')

NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_100PPC/')
NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_200PPC/')
NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_400PPC/')
NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_800PPC/')
NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_1600PPC/')
NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_3200PPC/')
NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_6400PPC/')
NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_12800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_12800PPC/')

NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_100PPC/')
NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_200PPC/')
NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_400PPC/')
NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_800PPC/')
NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_1600PPC/')
NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_3200PPC/')
NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_6400PPC/')
NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_12800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_12800PPC/')
NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_25600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2thirdsDomain_150Cells_uniform_640delT_25600PPC/')

CIC_ICP_lowPressure_2thirdsDomain_48Cells_sinusoid_640delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_2thirdsDomain_48Cells_sinusoid_640delT_64PPC_epsneg4/')







# -------------------- Trends between domain sizes  --------------------------------------------

plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = 'Exp Res. full length')
plotAveDensity(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_3200PPC, 'He+', label = 'Exp Res. 2/3 length')
plotAveDensity(NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4, 'He+', label = 'NGP full length 64 PPC')
plotAveDensity(NGP_ICP_lowPressure_2thirdsDomain_48Cells_sinusoid_640delT_64PPC_epsneg4, 'He+', label = 'NGP 2/3 length 64 PPC')
plotAveDensity(NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_6400PPC, 'He+', label = 'NGP 2/3 6400 PPC')
plt.legend(loc = 'lower center')
plt.xlim([0, Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC.grid[-1] * 1e2])
plt.savefig('RF_ICP/ICP_domainComp_density.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = 'Exp Res. full length')
plotAvePhi(Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_3200PPC, label = 'Exp Res. 3/4 length')
plotAvePhi(NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4,label = 'NGP full length 64 cells 64 PPC')
plotAvePhi(NGP_ICP_lowPressure_2thirdsDomain_48Cells_sinusoid_640delT_64PPC_epsneg4, label = 'NGP 3/4 length 64 cells 64 PPC')
plotAvePhi(NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_6400PPC, label = 'NGP 3/4 128 cells 6400 PPC')
plt.legend(loc = 'lower center')
plt.ylim([0,None])
plt.xlim([0, Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC.grid[-1] * 1e2])
plt.savefig('RF_ICP/ICP_domainComp_phi.png')


dataset_list_full_Exp = [Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC, Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC,
                             Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC, Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC,
                             Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC]

dataset_list_full_NGP_128cells = [NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC, NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC,
                             NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC, NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC,
                             NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC]

dataset_list_2Thirds_Exp = [Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_400PPC, Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_800PPC,
                             Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_1600PPC, Exp_ICP_highDensity_lowPressure_2thirdsDomain_standRF_1500cells_8000delT_3200PPC]


dataset_list_2Thirds_NGP_96cells = [NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_400PPC, NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_800PPC,
                             NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_1600PPC, NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_3200PPC,
                             NGP_ICP_lowPressure_2thirdsDomain_96Cells_sinusoid_640delT_6400PPC]


data_list_labels = ['Exp. low press', 'NGP low press 128 cells', 'Exp. high press', 'NGP high press 128 cells']


plt.figure()
plotWallCurrent(np.array([200,400,800,1600,3200]), 'He+', dataset_list_full_Exp, label = 'Exp. full length')
plotWallCurrent(np.array([400,800,1600,3200]), 'He+', dataset_list_2Thirds_Exp, label = 'Exp. 2/3 length')
plotWallCurrent(np.array([200,400,800,1600,3200]), 'He+', dataset_list_full_NGP_128cells, label = 'NGP full length 128 cells')
plotWallCurrent(np.array([400,800,1600,3200, 6400]), 'He+', dataset_list_2Thirds_NGP_96cells, label = 'NGP 2/3 length 96 cells')
plotWallCurrent(np.array([64]), 'He+', [NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4], label = 'NGP. full length 64 cells')
plotWallCurrent(np.array([64]), 'He+', [NGP_ICP_lowPressure_2thirdsDomain_48Cells_sinusoid_640delT_64PPC_epsneg4], label = 'NGP. 2/3 length 48 cells')
plt.legend(loc = 'best')
plt.ylabel(r'He$^+$ wall current ($\frac{A}{m^2}$)')
plt.xlabel(r'Particle-per-cell')
plt.savefig('RF_ICP/ICP_domainComp_current.png')



