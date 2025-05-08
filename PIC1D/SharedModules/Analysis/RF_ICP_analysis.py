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

Exp_ICP_highDensity_highPressure_standRF_640delT_400PPC_eps100 = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_highPressure_standRF_640delT_400PPC_eps100/')

ECExp_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC/')
ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_400PPC/')
ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_800PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_800PPC/')
ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_1600PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_1600PPC/')
ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_3200PPC/')
ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_6400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_6400PPC/')

ECExp_ICP_highPressure_64Cells_sinusoid_6400delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highPressure_64Cells_sinusoid_6400delT_64PPC/')

NGP_ICP_highPressure_200Cells_uniform_2500delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_100PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_200PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_400PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_800PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_1600PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_3200PPC/')
NGP_ICP_highPressure_200Cells_uniform_2500delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_200Cells_uniform_2500delT_6400PPC/')

NGP_ICP_highPressure_1000Cells_uniform_2500delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_1000Cells_uniform_2500delT_1600PPC/')
NGP_ICP_highPressure_1000Cells_uniform_2500delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_1000Cells_uniform_2500delT_3200PPC/')

NGP_ICP_highPressure_2000Cells_uniform_6400delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_2000Cells_uniform_6400delT_1600PPC/')
NGP_ICP_highPressure_2000Cells_uniform_6400delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_2000Cells_uniform_6400delT_3200PPC/')

NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4/')
NGP_ICP_highPressure_64Cells_sinusoid_2500delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_64Cells_sinusoid_2500delT_200PPC/')

NGP_ICP_highPressure_128Cells_sinusoid_2500delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_100PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_200PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_800PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_12800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_12800PPC/')

NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC_smoothing/')
NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC_smoothing/')

NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC_epsneg8 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC_epsneg8/')
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
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_200PPC, 'He+', label = '200 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_400PPC,'He+', label = '400 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_800PPC,'He+', label = '800 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_1600PPC,'He+', label = '1600 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_6400delT_3200PPC,'He+', label = '3200 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, 'He+', label = r'3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plt.xlim(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[0],Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/Exp_ICP_highPress_density.png')
plt.savefig('RF_ICP/Exp_ICP_highPress_density.pdf')
plt.close()

plt.figure()
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_200PPC, label = '200 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_400PPC, label = '400 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_800PPC, label = '800 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_1600PPC, label = '1600 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_6400delT_3200PPC, label = '3200 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, label = r'3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plt.xlim(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[0], Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/Exp_ICP_highPress_phi.png')
plt.savefig('RF_ICP/Exp_ICP_highPress_phi.pdf')
plt.close()


plt.figure()
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC, 'He+', label = 'INGP, 400 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC, 'He+', label = 'INGP, 1600 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC, 'He+', label = 'INGP, 6400 PPC')
plotAveDensity(ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_400PPC, 'He+', label = 'EC-PIC, 400 PPC')
plotAveDensity(ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_1600PPC, 'He+', label = 'EC-PIC, 1600 PPC')
plotAveDensity(ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_6400PPC, 'He+', label = 'EC-PIC, 6400 PPC')
plotAveDensity(NGP_ICP_highPressure_2000Cells_uniform_6400delT_3200PPC, 'He+', label = r'INGP, 3200 PPC, 2000 cells',  marker = '')
plt.xlim(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[0], Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_highPress_uneven_density.png')
plt.savefig('RF_ICP/ICP_highPress_uneven_density.pdf')
plt.close()

plt.figure()
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC, label = 'INGP, 400 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC, label = 'INGP, 1600 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC, label = 'INGP, 6400 PPC')
plotAvePhi(ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_400PPC,  label = 'EC-PIC, 400 PPC')
plotAvePhi(ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_1600PPC, label = 'EC-PIC, 1600 PPC')
plotAvePhi(ECExp_ICP_highPressure_128Cells_sinusoid_6400delT_6400PPC, label = 'EC-PIC, 6400 PPC')
plotAvePhi(NGP_ICP_highPressure_2000Cells_uniform_6400delT_3200PPC, label = r'INGP, 3200 PPC, 2000 cells', marker = '')
plt.xlim(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[0], Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_highPress_uneven_phi.png')
plt.savefig('RF_ICP/ICP_highPress_uneven_phi.pdf')
plt.close()

plt.figure()
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC, 'He+', label = '400 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC, 'He+', label = '3200 PPC')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC_smoothing, 'He+', label = '400 PPC, smoothing')
plotAveDensity(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC_smoothing, 'He+', label = '3200 PPC, smoothing')
plotAveDensity(NGP_ICP_highPressure_2000Cells_uniform_6400delT_3200PPC, 'He+', label = r'INGP, 3200 PPC, 2000 cells', marker = '')
plt.xlim(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[0], Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_highPress_uneven_smooth_density.png')
plt.savefig('RF_ICP/ICP_highPress_uneven_smooth_density.pdf')
plt.close()

plt.figure()
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC, label = '400 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC, label = '3200 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_400PPC_smoothing, label = '400 PPC, smoothing')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_3200PPC_smoothing, label = '3200 PPC, smoothing')
plotAvePhi(NGP_ICP_highPressure_2000Cells_uniform_6400delT_3200PPC, label = r'INGP, 3200 PPC, 2000 cells', marker = '')
plt.xlim(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[0], Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_highPress_uneven_smooth_phi.png')
plt.savefig('RF_ICP/ICP_highPress_uneven_smooth_phi.pdf')
plt.close()

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
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, 'He+', label = 'Exp High Res.')
plotAveDensity(NGP_ICP_highPressure_1000Cells_uniform_2500delT_1600PPC, 'He+', label = 'NGP 1600 PPC, 1000 Cells')
plotAveDensity(NGP_ICP_highPressure_1000Cells_uniform_2500delT_3200PPC, 'He+', label = 'NGP 3200 PPC, 1000 Cells')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_highPress_density_uniform_1000cells.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, label = 'Exp High Res.')
plotAvePhi(NGP_ICP_highPressure_1000Cells_uniform_2500delT_1600PPC, label = 'NGP 1600 PPC, 1000 Cells')
plotAvePhi(NGP_ICP_highPressure_1000Cells_uniform_2500delT_3200PPC, label = 'NGP 3200 PPC, 1000 Cells')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_highPress_phi_uniform_1000cells.png')

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

ECExp_ICP_lowPressure_2500Cells_uniform_7000delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_2500Cells_uniform_7000delT_3200PPC/')

Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC/')
Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC/')
Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC_subcycle = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC_subcycle/')

Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC_implicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_1500cells_5000delT_128PPC_implicit/')

Exp_ICP_highDensity_lowPressure_standRF_700delT_400PPC_eps100 = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_700delT_400PPC_eps100/')

ECExp_ICP_lowPressure_128Cells_sinusoid_10000delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_10000delT_400PPC/')
ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_400PPC/')
ECExp_ICP_lowPressure_128Cells_sinusoid_3500delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_3500delT_400PPC/')
ECExp_ICP_lowPressure_128Cells_sinusoid_1750delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_1750delT_400PPC/')
ECExp_ICP_lowPressure_128Cells_sinusoid_1750delT_400PPC_implicit = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_1750delT_400PPC_implicit/')

ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_400PPC_subcycle = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_400PPC_subcycle/')
ECExp_ICP_lowPressure_128Cells_sinusoidSmaller_7000delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoidSmaller_7000delT_400PPC/')
ECExp_ICP_lowPressure_128Cells_30CellEven_7000delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_30CellEven_7000delT_400PPC/')

ECExp_ICP_lowPressure_128Cells_20CellEven_7000delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_20CellEven_7000delT_400PPC/')

ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_800PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_800PPC/')
ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_1600PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_1600PPC/')
ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_3200PPC/')

ECExp_ICP_lowPressure_64Cells_sinusoid_7000delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_64Cells_sinusoid_7000delT_64PPC/')

NGP_ICP_lowPressure_2000Cells_sinusoid_420delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2000Cells_sinusoid_420delT_3200PPC/')
NGP_ICP_lowPressure_2000Cells_uniform_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2000Cells_uniform_640delT_3200PPC/')
NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC/')

NGP_ICP_lowPressure_1000Cells_uniform_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_1000Cells_uniform_640delT_1600PPC/')

NGP_ICP_lowPressure_128Cells_DenisMap_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_DenisMap_640delT_100PPC/')


NGP_ICP_lowPressure_128Cells_uniform_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_uniform_640delT_400PPC/')
NGP_ICP_lowPressure_128Cells_uniform_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_uniform_640delT_3200PPC/')
NGP_ICP_lowPressure_128Cells_uniform_640delT_3200PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_uniform_640delT_3200PPC_smoothing/')
NGP_ICP_lowPressure_128Cells_uniform_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_uniform_640delT_6400PPC/')

NGP_ICP_lowPressure_200Cells_uniform_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_100PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_200PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_400PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_800PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_1600PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_3200PPC/')
NGP_ICP_lowPressure_200Cells_uniform_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_200Cells_uniform_640delT_6400PPC/')

NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC/')
NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4/')
NGP_ICP_lowPressure_64Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_64Cells_sinusoid_640delT_200PPC/')

NGP_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC/')
NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC/')

NGP_ICP_lowPressure_128Cells_20CellEven_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_20CellEven_640delT_400PPC/')
NGP_ICP_lowPressure_128Cells_20CellEven_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_20CellEven_640delT_3200PPC/')
NGP_ICP_lowPressure_128Cells_20CellEven_640delT_3200PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_128Cells_20CellEven_640delT_3200PPC_smoothing/')

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


CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing_epsneg3 = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing_epsneg3/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC_smoothing/')
CIC_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC_smoothing = dataSet('Y:/ImplicitPic1D/ImplicitData/CIC_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC_smoothing/')


plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC, 'He+', label = '200 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC, 'He+', label = '400 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC, 'He+', label = '800 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC, 'He+', label = '1600 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC, 'He+', label = '3200 PPC', marker = '')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = r'3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/Exp_ICP_lowPress.pdf')
plt.savefig('RF_ICP/Exp_ICP_lowPress.png')
plt.close()

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC, label = '200 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC, label = '400 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC, label = '800 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC, label = '1600 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC, label = '3200 PPC', marker = '')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = r'3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/Exp_ICP_lowPress_phi.pdf')
plt.savefig('RF_ICP/Exp_ICP_lowPress_phi.png')
plt.close()

plt.figure()
plotAveEPF(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC, 'e', label = '200 PPC', marker = '')
plotAveEPF(Exp_ICP_highDensity_lowPressure_standRF_7000delT_400PPC, 'e', label = '400 PPC', marker = '')
plotAveEPF(Exp_ICP_highDensity_lowPressure_standRF_7000delT_800PPC, 'e', label = '800 PPC', marker = '')
plotAveEPF(Exp_ICP_highDensity_lowPressure_standRF_7000delT_1600PPC, 'e', label = '1600 PPC', marker = '')
plotAveEPF(Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC, 'e', label = '3200 PPC', marker = '')
plotAveEPF(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'e', label = r'3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plt.xlim(0, 40)
plt.ylim(2e-4, 5e-2)
plt.xlabel('Electron Energy (eV)', fontsize = 14)
plt.ylabel(r'Electron Energy Prob. Dist. (eV$^{-3/2}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower left', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/Exp_ICP_lowPress_EPF.pdf')
plt.savefig('RF_ICP/Exp_ICP_lowPress_EPF.png')
plt.close()


plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC, 'He+', label = r'MC-PIC, 2000 cells', marker = '')
plotAveDensity(ECExp_ICP_lowPressure_2500Cells_uniform_7000delT_3200PPC, 'He+', label = r'EC-PIC, 2500 cells', marker = '')
plotAveDensity(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, 'He+', label = r'INGP, 2000 cells', marker = '')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_lowPress_comp_density.pdf')
plt.savefig('RF_ICP/ICP_lowPress_comp_density.png')
plt.close()

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC, label = r'MC-PIC, 2000 cells', marker = '')
plotAvePhi(ECExp_ICP_lowPressure_2500Cells_uniform_7000delT_3200PPC, label = r'EC-PIC, 2500 cells', marker = '')
plotAvePhi(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, label = r'INGP, 2000 cells', marker = '')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_lowPress_comp_phi.pdf')
plt.savefig('RF_ICP/ICP_lowPress_comp_phi.png')
plt.close()

plt.figure()
plotAveEPF(Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC, 'e', label = r'MC-PIC, 2000 cells', marker = '')
plotAveEPF(ECExp_ICP_lowPressure_2500Cells_uniform_7000delT_3200PPC, 'e', label = r'EC-PIC, 2500 cells', marker = '')
plotAveEPF(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, 'e', label = r'INGP, 2000 cells', marker = '')
plt.xlim(0, 40)
plt.ylim(2e-4, 5e-2)
plt.xlabel('Electron Energy (eV)', fontsize = 14)
plt.ylabel(r'Electron Energy Prob. Dist. (eV$^{-3/2}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower left', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_lowPress_comp_EPF.pdf')
plt.savefig('RF_ICP/ICP_lowPress_comp_EPF.png')
plt.close()

plt.figure()
xi = np.arange(0, NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC.Nx-1) + 1
plt.plot(xi, np.diff(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC.grid) * 1e2, 'o', markersize = 4)
plt.xlim(xi[0], xi[-1])
plt.xlabel(r'Cell Number', fontsize = 14)
plt.ylabel(r'Cell Size (cm)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('RF_ICP/lowPress_maps.pdf')
plt.savefig('RF_ICP/lowPress_maps.png')
plt.close()

plt.figure()
plotAveDensity(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, 'He+', label = r'INGP 3200 PPC, 2000 cells, $\Delta t = \frac{\tau_{RF}}{7000}$', marker = '')
# plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC, 'He+', label = 'Exp Low Res.')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, 'He+',label = 'INGP 400 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC, 'He+',label = 'INGP 800 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC, 'He+',label = 'INGP 1600 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC, 'He+',label = 'INGP 3200 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, 'He+',label = 'INGP 6400 PPC')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/NGP_ICP_lowPress_noneven_density.pdf')
plt.savefig('RF_ICP/NGP_ICP_lowPress_noneven_density.png')
plt.close()

plt.figure()
plotAvePhi(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, label = r'INGP 3200 PPC, 2000 cells, $\Delta t = \frac{\tau_{RF}}{7000}$', marker = '')
# plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC, 'He+', label = 'Exp Low Res.')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC,label = 'INGP 400 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_800PPC,label = 'INGP 800 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC, label = 'INGP 1600 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC,label = 'INGP 3200 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, label = 'INGP 6400 PPC')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/NGP_ICP_lowPress_noneven_phi.png')
plt.savefig('RF_ICP/NGP_ICP_lowPress_noneven_phi.pdf')
plt.close()

plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = r'MC-PIC 3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC, 'He+', label = r'MC-PIC 64 PPC, 1100 cells, $\Delta t = \frac{\tau_{RF}}{4000}$', marker = '')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, 'He+',label = 'INGP 400 PPC')
plotAveDensity(ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_400PPC, 'He+',label = 'EC-PIC 400 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC, 'He+',label = 'INGP 3200 PPC')
plotAveDensity(ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_3200PPC, 'He+',label = 'EC-PIC 3200 PPC')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/NGP_ECExp_lowPress_noneven_density.pdf')
plt.savefig('RF_ICP/NGP_ECExp_lowPress_noneven_density.png')
plt.close()

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = r'MC-PIC 3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_1100cells_4000delT_64PPC, label = r'MC-PIC 64 PPC, 1100 cells, $\Delta t = \frac{\tau_{RF}}{4000}$', marker = '')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, label = 'INGP 400 PPC')
plotAvePhi(ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_400PPC, label = 'EC-PIC 400 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC, label = 'INGP 3200 PPC')
plotAvePhi(ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_3200PPC, label = 'EC-PIC 3200 PPC')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/NGP_ECExp_lowPress_noneven_phi.pdf')
plt.savefig('RF_ICP/NGP_ECExp_lowPress_noneven_phi.png')
plt.close()


plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = r'MC-PIC 3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plotAveDensity(ECExp_ICP_lowPressure_128Cells_sinusoid_10000delT_400PPC, 'He+',label = r'EC-PIC, $\frac{\tau_{RF}}{\Delta t} = 10000$')
plotAveDensity(ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_400PPC, 'He+',label = r'EC-PIC, $\frac{\tau_{RF}}{\Delta t} = 7000$')
plotAveDensity(ECExp_ICP_lowPressure_128Cells_sinusoid_3500delT_400PPC, 'He+',label = r'EC-PIC, $\frac{\tau_{RF}}{\Delta t} = 3500$')
plotAveDensity(ECExp_ICP_lowPressure_128Cells_sinusoid_1750delT_400PPC, 'He+',label = r'EC-PIC, $\frac{\tau_{RF}}{\Delta t} = 1750$')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ECExp_timeChange_density.pdf')
plt.savefig('RF_ICP/ECExp_timeChange_density.png')
plt.close()

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = r'MC-PIC 3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plotAvePhi(ECExp_ICP_lowPressure_128Cells_sinusoid_10000delT_400PPC, label = r'EC-PIC, $\frac{\tau_{RF}}{\Delta t} = 10000$')
plotAvePhi(ECExp_ICP_lowPressure_128Cells_sinusoid_7000delT_400PPC, label = r'EC-PIC, $\frac{\tau_{RF}}{\Delta t} = 7000$')
plotAvePhi(ECExp_ICP_lowPressure_128Cells_sinusoid_3500delT_400PPC, label = r'EC-PIC, $\frac{\tau_{RF}}{\Delta t} = 3500$')
plotAvePhi(ECExp_ICP_lowPressure_128Cells_sinusoid_1750delT_400PPC, label = r'EC-PIC, $\frac{\tau_{RF}}{\Delta t} = 1750$')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ECExp_timeChange_phi.pdf')
plt.savefig('RF_ICP/ECExp_timeChange_phi.png')
plt.close()


plt.figure()
plotAveDensity(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, 'He+', label =  r'2000 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, 'He+',label = 'INGP 400 PPC, Map 1')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_400PPC, 'He+',label = 'INGP 400 PPC, Map 2')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_uniform_640delT_400PPC, 'He+',label = 'INGP 400 PPC, uniform')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC, 'He+',label = '128 cells, Map 1')
plotAveDensity(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_3200PPC, 'He+',label = '128 cells, Map 2')
plotAveDensity(NGP_ICP_lowPressure_128Cells_uniform_640delT_3200PPC, 'He+',label = '128 cells, uniform')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/NGP_changeSheath_lowPress_noneven_density.pdf')
plt.savefig('RF_ICP/NGP_changeSheath_lowPress_noneven_density.png')
plt.close()

plt.figure()
plotAvePhi(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, label =  r'2000 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, 'He+',label = 'INGP 400 PPC, Map 1')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_400PPC, 'He+',label = 'INGP 400 PPC, Map 2')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_uniform_640delT_400PPC, 'He+',label = 'INGP 400 PPC, uniform')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC, label = '128 cells, Map 1')
plotAvePhi(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_3200PPC, label = '128 cells, Map 2')
plotAvePhi(NGP_ICP_lowPressure_128Cells_uniform_640delT_3200PPC, label = '128 cells, uniform')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/NGP_changeSheath_lowPress_noneven_phi.pdf')
plt.savefig('RF_ICP/NGP_changeSheath_lowPress_noneven_phi.png')
plt.close()

plt.figure()
xi = np.arange(0, NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC.Nx-1) + 1
plt.plot(xi, np.diff(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC.grid) * 1e2, 'o', markersize = 4,  label = 'Map 1')
plt.plot(xi, np.diff(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_400PPC.grid) * 1e2, 'o', markersize = 4,label = 'Map 2')
plt.plot(xi, np.diff(NGP_ICP_lowPressure_128Cells_uniform_640delT_400PPC.grid) * 1e2, 'o', markersize = 4,label = 'Uniform')
plt.xlim(xi[0], xi[-1])
plt.xlabel(r'Cell Number', fontsize = 14)
plt.ylabel(r'Cell Size (cm)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('RF_ICP/lowPress_maps_changeSheath.pdf')
plt.savefig('RF_ICP/lowPress_maps_changeSheath.png')
plt.close()

plt.figure()
plotAveDensity(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, 'He+', label =  r'2000 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, 'He+',label = 'INGP 400 PPC, Map 1')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_400PPC, 'He+',label = 'INGP 400 PPC, Map 2')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_uniform_640delT_400PPC, 'He+',label = 'INGP 400 PPC, uniform')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC_smoothing, 'He+',label = '128 cells, Map 1')
plotAveDensity(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_3200PPC_smoothing, 'He+',label = '128 cells, Map 2')
plotAveDensity(NGP_ICP_lowPressure_128Cells_uniform_640delT_3200PPC_smoothing, 'He+',label = '128 cells, uniform')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/NGP_changeSheath_smooth_lowPress_noneven_density.pdf')
plt.savefig('RF_ICP/NGP_changeSheath_smooth_lowPress_noneven_density.png')
plt.close()

plt.figure()
plotAvePhi(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, label =  r'2000 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC, 'He+',label = 'INGP 400 PPC, Map 1')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_400PPC, 'He+',label = 'INGP 400 PPC, Map 2')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_uniform_640delT_400PPC, 'He+',label = 'INGP 400 PPC, uniform')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_3200PPC_smoothing, label = '128 cells, Map 1')
plotAvePhi(NGP_ICP_lowPressure_128Cells_20CellEven_640delT_3200PPC_smoothing, label = '128 cells, Map 2')
plotAvePhi(NGP_ICP_lowPressure_128Cells_uniform_640delT_3200PPC_smoothing, label = '128 cells, uniform')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/NGP_changeSheath_smooth_lowPress_noneven_phi.pdf')
plt.savefig('RF_ICP/NGP_changeSheath_smooth_lowPress_noneven_phi.png')
plt.close()


plt.figure()
plotAveDensity(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, 'He+', label =  r'INGP, 2000 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
# plotAveDensity(CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing_epsneg3, 'He+',label = 'ICIC, 100 PPC epsneg3')
plotAveDensity(CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing, 'He+',label = 'ICIC, 100 PPC')
plotAveDensity(CIC_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC_smoothing, 'He+',label = 'ICIC, 200 PPC')
plotAveDensity(CIC_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC_smoothing, 'He+',label = 'ICIC, 400 PPC')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICIC_smooth_lowPress_noneven_density.pdf')
plt.savefig('RF_ICP/ICIC_smooth_lowPress_noneven_density.png')
plt.close()

plt.figure()
plotAvePhi(NGP_ICP_lowPressure_2000Cells_uniform_7000delT_3200PPC, label =  r'INGP, 2000 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
# plotAvePhi(CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing_epsneg3, label = 'ICIC, 100 PPC, epsneg3')
plotAvePhi(CIC_ICP_lowPressure_128Cells_sinusoid_640delT_100PPC_smoothing, label = 'ICIC, 100 PPC')
plotAvePhi(CIC_ICP_lowPressure_128Cells_sinusoid_640delT_200PPC_smoothing, label = 'ICIC, 200 PPC')
plotAvePhi(CIC_ICP_lowPressure_128Cells_sinusoid_640delT_400PPC_smoothing, label = 'ICIC, 400 PPC')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'Electric Potential (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICIC_smooth_lowPress_noneven_phi.pdf')
plt.savefig('RF_ICP/ICIC_smooth_lowPress_noneven_phi.png')
plt.close()

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

plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = 'Exp High Res.')
plotAveDensity(NGP_ICP_lowPressure_1000Cells_uniform_640delT_1600PPC, 'He+',label = 'NGP 1600 PPC, 1000 cells')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/NGP_ICP_lowPress_density_uniform1000cells.png')


# -------------------- Trends between pressures  --------------------------------------------
plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = r'MC-PIC, low press., 3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plotAveDensity(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, 'He+', label = r'MC-PIC, high press., 3200 PPC, 2500 cells, $\Delta t = \frac{\tau_{RF}}{10^4}$', marker = '')
plotAveDensity(NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4, 'He+', label = 'INGP, high press, 64 cell, 64 PPC')
plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, 'He+', label = 'INGP, low press, 128 cells, 6400 PPC')
plt.xlim(Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[0], Exp_ICP_highDensity_lowPressure_standRF_7000delT_200PPC.grid[-1]*1e2/2)
plt.ylim(0, None)
plt.xlabel('Distance (cm)', fontsize = 14)
plt.ylabel(r'He$^+$ Density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'lower right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_pressComp_density.png')
plt.savefig('RF_ICP/ICP_pressComp_density.pdf')
plt.close()

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = 'Exp Res. low press')
plotAvePhi(Exp_ICP_highDensity_highPressure_standRF_2500cells_10000delT_3200PPC, label = 'Exp Res. high press')
plotAvePhi(NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4, label = 'NGP low press 64 PPC')
plotAvePhi(NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4, label = 'NGP high press 64 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, label = 'NGP low press 6400 PPC')
plotAvePhi(NGP_ICP_highPressure_128Cells_sinusoid_2500delT_6400PPC, label = 'NGP high press 6400 PPC')
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


#-------------------------------- low pressure cut J in 3/4 -------------------------

Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_100PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_100PPC/')
Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_400PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_400PPC/')
Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_800PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_800PPC/')
Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_1600PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_1600PPC/')
Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_3200PPC/')

plt.figure()
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_400PPC, 'He+', label = '400 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_800PPC, 'He+', label = '800 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_1600PPC, 'He+', label = '1600 PPC')
plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_3200PPC, 'He+', label = '3200 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/Exp_ICP_lowPress_cutJ.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_400PPC,  label = '400 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_800PPC,  label = '800 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_1600PPC, label = '1600 PPC')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_3200PPC,  label = '3200 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/Exp_ICP_lowPress_phi_cutJ.png')

NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_100PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_100PPC/')
NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_200PPC/')
NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_400PPC/')
NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_800PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_800PPC/')
NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_1600PPC/')
NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_3200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_3200PPC/')
NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_6400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_6400PPC/')

NGP_ICP_lowPressure_cutJ_64Cells_sinusoid_640delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_64Cells_sinusoid_640delT_64PPC_epsneg4/')
NGP_ICP_lowPressure_cutJ_64Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_cutJ_64Cells_sinusoid_640delT_200PPC/')

ECExp_ICP_lowPressure_cutJ_64Cells_sinusoid_6400delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_cutJ_64Cells_sinusoid_6400delT_64PPC/')

# plt.figure()
# plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, 'He+', label = 'Exp Res. J_0')
# plotAveDensity(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_3200PPC, 'He+', label = 'Exp Res. 3/4 J_0')
# plotAveDensity(NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4, 'He+', label = 'NGP J_0 64 PPC')
# plotAveDensity(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, 'He+', label = 'NGP J_0 6400 PPC')
# plotAveDensity(NGP_ICP_lowPressure_cutJ_64Cells_sinusoid_640delT_64PPC_epsneg4, 'He+', label = 'NGP 3/4 J_0 64 PPC')
# plotAveDensity(NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_6400PPC, 'He+', label = 'NGP 3/4 J_0 6400 PPC')
# plt.legend(loc = 'lower center')
# plt.savefig('RF_ICP/ICP_domainComp_density.png')

plt.figure()
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_10000delT_2500cells_3200PPC, label = 'Exp Res. J_0')
plotAvePhi(Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_3200PPC, label = 'Exp Res. 3/4 J_0')
plotAvePhi(NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4, label = 'NGP J_0 64 cells 64 PPC')
plotAvePhi(NGP_ICP_lowPressure_128Cells_sinusoid_640delT_6400PPC, label = 'NGP J_0 128 cells 6400 PPC')
plotAvePhi(NGP_ICP_lowPressure_cutJ_64Cells_sinusoid_640delT_64PPC_epsneg4, label = 'NGP 3/4 J_0 64 cells 64 PPC')
plt.legend(loc = 'lower center')
plt.savefig('RF_ICP/ICP_JComp_phi.png')


# ------------------ low pressure cut J 1/2 ----------------

Exp_ICP_highDensity_lowPressure_standRF_halfJ_5000delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure_standRF_halfJ_5000delT_3200PPC/')
NGP_ICP_lowPressure_halfJ_64Cells_sinusoid_640delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_halfJ_64Cells_sinusoid_640delT_64PPC_epsneg4/')
NGP_ICP_lowPressure_halfJ_64Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_halfJ_64Cells_sinusoid_640delT_200PPC/')
NGP_ICP_lowPressure_halfJ_128Cells_sinusoid_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_halfJ_128Cells_sinusoid_640delT_400PPC/')
NGP_ICP_lowPressure_halfJ_128Cells_sinusoid_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_lowPressure_halfJ_128Cells_sinusoid_640delT_1600PPC/')

ECExp_ICP_lowPressure_halfJ_64Cells_sinusoid_6400delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_lowPressure_halfJ_64Cells_sinusoid_6400delT_64PPC/')

# -------------------- med pressure ----------------------------

Exp_ICP_highDensity_MedPressure_standRF_6400delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_MedPressure_standRF_6400delT_3200PPC/')
NGP_ICP_highDensity_MedPressure_64Cells_sinusoid_1200delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highDensity_MedPressure_64Cells_sinusoid_1200delT_64PPC_epsneg4/')
NGP_ICP_highDensity_MedPressure_64Cells_sinusoid_1200delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highDensity_MedPressure_64Cells_sinusoid_1200delT_200PPC/')
NGP_ICP_highDensity_MedPressure_128Cells_sinusoid_1200delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highDensity_MedPressure_128Cells_sinusoid_1200delT_400PPC/')
NGP_ICP_highDensity_MedPressure_128Cells_sinusoid_1200delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highDensity_MedPressure_128Cells_sinusoid_1200delT_1600PPC/')

ECExp_ICP_highDensity_MedPressure_64Cells_sinusoid_6400delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highDensity_MedPressure_64Cells_sinusoid_6400delT_64PPC/')

# ------------------- low pressure ------------

Exp_ICP_highDensity_lowPressure8_standRF_6400delT_3200PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_ICP_highDensity_lowPressure8_standRF_6400delT_3200PPC/')
NGP_ICP_highDensity_lowPressure8_64Cells_sinusoid_640delT_64PPC_epsneg4 = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highDensity_lowPressure8_64Cells_sinusoid_640delT_64PPC_epsneg4/')
NGP_ICP_highDensity_lowPressure8_64Cells_sinusoid_640delT_200PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highDensity_lowPressure8_64Cells_sinusoid_640delT_200PPC/')
NGP_ICP_highDensity_lowPressure8_128Cells_sinusoid_640delT_400PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highDensity_lowPressure8_128Cells_sinusoid_640delT_400PPC/')
NGP_ICP_highDensity_lowPressure8_128Cells_sinusoid_640delT_1600PPC = dataSet('Y:/ImplicitPic1D/ImplicitData/NGP_ICP_highDensity_lowPressure8_128Cells_sinusoid_640delT_1600PPC/')

ECExp_ICP_highDensity_lowPressure8_64Cells_sinusoid_6400delT_64PPC = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/ECExp_ICP_highDensity_lowPressure8_64Cells_sinusoid_6400delT_64PPC/')

list_exp = [Exp_ICP_highDensity_lowPressure_standRF_halfJ_5000delT_3200PPC,
        Exp_ICP_highDensity_lowPressure_standRF_cutJ_7000delT_3200PPC,
        Exp_ICP_highDensity_lowPressure_standRF_7000delT_3200PPC]

list_ECexp = [ECExp_ICP_lowPressure_halfJ_64Cells_sinusoid_6400delT_64PPC, ECExp_ICP_lowPressure_cutJ_64Cells_sinusoid_6400delT_64PPC,
                ECExp_ICP_lowPressure_64Cells_sinusoid_7000delT_64PPC]
list_NGP_low = [NGP_ICP_lowPressure_halfJ_64Cells_sinusoid_640delT_64PPC_epsneg4, NGP_ICP_lowPressure_cutJ_64Cells_sinusoid_640delT_64PPC_epsneg4,
                NGP_ICP_lowPressure_64Cells_sinusoid_640delT_64PPC_epsneg4]
list_NGP_med = [NGP_ICP_lowPressure_halfJ_64Cells_sinusoid_640delT_200PPC, NGP_ICP_lowPressure_cutJ_64Cells_sinusoid_640delT_200PPC,
                NGP_ICP_lowPressure_64Cells_sinusoid_640delT_200PPC]
list_NGP_high = [NGP_ICP_lowPressure_halfJ_128Cells_sinusoid_640delT_1600PPC, NGP_ICP_lowPressure_cutJ_128Cells_sinusoid_640delT_1600PPC,
                 NGP_ICP_lowPressure_128Cells_sinusoid_640delT_1600PPC]
J_list = [4000, 6000, 8000]

plt.figure()
aveQuantity_vs_parameter('density', [list_exp, list_ECexp, list_NGP_low, list_NGP_high], J_list, ['MC-PIC resolved', r'EC-PIC, 64 Cells, 64 PPC',r'INGP 64 Cells, 64 PPC','INGP 128 Cells, 1600 PPC'], 'He+')
plt.xlabel(r'Current Density (A/m$^2$)', fontsize = 14)
plt.ylabel(r'Average He$^+$ density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'upper left', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_Jcomp_runs_density.png')
plt.savefig('RF_ICP/ICP_Jcomp_runs_density.pdf')
plt.close()


plt.figure()
aveQuantity_vs_parameter('phi', [list_exp, list_ECexp, list_NGP_low, list_NGP_high], J_list, ['MC-PIC resolved', r'EC-PIC, 64 Cells, 64 PPC',r'INGP 64 Cells, 64 PPC','INGP 128 Cells, 1600 PPC'], 'He+')
plt.xlabel(r'Current Density (A/m$^2$)', fontsize = 14)
plt.ylabel(r'Average $\phi$ (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'center right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_Jcomp_runs_phi.png')
plt.savefig('RF_ICP/ICP_Jcomp_runs_phi.pdf')
plt.close()

plt.figure()
aveQuantity_vs_parameter('current', [list_exp, list_ECexp, list_NGP_low, list_NGP_high], J_list, ['MC-PIC resolved', r'EC-PIC, 64 Cells, 64 PPC',r'INGP 64 Cells, 64 PPC','INGP 128 Cells, 1600 PPC'], 'He+')
plt.xlabel(r'Current Density (A/m$^2$)', fontsize = 14)
plt.ylabel(r'Average He$^+$ wall current (A/m$^2$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'upper left', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_Jcomp_runs_current.png')
plt.savefig('RF_ICP/ICP_Jcomp_runs_current.pdf')
plt.close()

plt.figure()
aveQuantity_vs_parameter('temperature', [list_exp, list_ECexp, list_NGP_low, list_NGP_high], J_list, ['MC-PIC resolved', r'EC-PIC, 64 Cells, 64 PPC',r'INGP 64 Cells, 64 PPC','INGP 128 Cells, 1600 PPC'], 'e')
plt.xlabel(r'Current Density (A/m$^2$)', fontsize = 14)
plt.ylabel(r'Electron Temperature (eV)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'center', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_Jcomp_runs_Te.png')
plt.savefig('RF_ICP/ICP_Jcomp_runs_Te.pdf')
plt.close()

list_exp = [Exp_ICP_highDensity_lowPressure8_standRF_6400delT_3200PPC,
        Exp_ICP_highDensity_MedPressure_standRF_6400delT_3200PPC,
        Exp_ICP_highDensity_highPressure_standRF_6400delT_3200PPC]
list_ECexp = [ECExp_ICP_highDensity_lowPressure8_64Cells_sinusoid_6400delT_64PPC, ECExp_ICP_highDensity_MedPressure_64Cells_sinusoid_6400delT_64PPC,
                ECExp_ICP_highPressure_64Cells_sinusoid_6400delT_64PPC]
list_NGP_low = [NGP_ICP_highDensity_lowPressure8_64Cells_sinusoid_640delT_64PPC_epsneg4, NGP_ICP_highDensity_MedPressure_64Cells_sinusoid_1200delT_64PPC_epsneg4,
                NGP_ICP_highPressure_64Cells_sinusoid_2500delT_64PPC_epsneg4]
# list_NGP_med = [NGP_ICP_highDensity_lowPressure8_64Cells_sinusoid_640delT_200PPC, NGP_ICP_highDensity_MedPressure_64Cells_sinusoid_1200delT_200PPC,
#                 NGP_ICP_highPressure_64Cells_sinusoid_2500delT_200PPC]
list_NGP_high = [NGP_ICP_highDensity_lowPressure8_128Cells_sinusoid_640delT_1600PPC, NGP_ICP_highDensity_MedPressure_128Cells_sinusoid_1200delT_1600PPC,
                 NGP_ICP_highPressure_128Cells_sinusoid_2500delT_1600PPC]
P_list = [8, 16, 32.1]

plt.figure()
aveQuantity_vs_parameter('density', [list_exp, list_ECexp, list_NGP_low, list_NGP_high], P_list, ['MC-PIC resolved', r'EC-PIC, 64 Cells, 64 PPC',r'INGP 64 Cells, 64 PPC','INGP 128 Cells, 1600 PPC'], 'He+')
plt.xlabel(r'Helium Density ($10^{20}$ m$^{-3}$)', fontsize = 14)
plt.ylabel(r'Average He$^+$ density (m$^{-3}$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'upper left', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_press_runs_density.png')
plt.savefig('RF_ICP/ICP_press_runs_density.pdf')
plt.close()


plt.figure()
aveQuantity_vs_parameter('phi', [list_exp, list_ECexp, list_NGP_low, list_NGP_high], P_list, ['MC-PIC resolved', r'EC-PIC, 64 Cells, 64 PPC',r'INGP 64 Cells, 64 PPC','INGP 128 Cells, 1600 PPC'], 'He+')
plt.xlabel(r'Helium Density ($10^{20}$ m$^{-3}$)', fontsize = 14)
plt.ylabel(r'Average $\phi$ (V)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'upper right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_press_runs_phi.png')
plt.savefig('RF_ICP/ICP_press_runs_phi.pdf')
plt.close()

plt.figure()
aveQuantity_vs_parameter('current', [list_exp, list_ECexp, list_NGP_low, list_NGP_high], P_list, ['MC-PIC resolved', r'EC-PIC, 64 Cells, 64 PPC',r'INGP 64 Cells, 64 PPC','INGP 128 Cells, 1600 PPC'], 'He+')
plt.xlabel(r'Helium Density ($10^{20}$ m$^{-3}$)', fontsize = 14)
plt.ylabel(r'Average He$^+$ wall current (A/m$^2$)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'upper right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_press_runs_current.png')
plt.savefig('RF_ICP/ICP_press_runs_current.pdf')
plt.close()

plt.figure()
aveQuantity_vs_parameter('temperature', [list_exp, list_ECexp, list_NGP_low, list_NGP_high], P_list, ['MC-PIC resolved', r'EC-PIC, 64 Cells, 64 PPC',r'INGP 64 Cells, 64 PPC','INGP 128 Cells, 1600 PPC'], 'e')
plt.xlabel(r'Helium Density ($10^{20}$ m$^{-3}$)', fontsize = 14)
plt.ylabel(r'Electron Temperature (eV)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc = 'upper right', fontsize = 12)
plt.tight_layout()
plt.savefig('RF_ICP/ICP_press_runs_Te.png')
plt.savefig('RF_ICP/ICP_press_runs_Te.pdf')
plt.close()