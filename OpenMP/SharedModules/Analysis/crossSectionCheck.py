# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from import_libraries_constants import *

path = 'Y:/scratch/nsavard/ImplicitPic1D/ExplicitData-BField/Exp_RFBenchmark_test/CrossSections/'
E_array = np.fromfile(path + 'IncidentPart_e_energy.dat')
sizeEarray = len(E_array)
sigmaArray = np.fromfile(path + 'IncidentPart_e_sigma.dat')
amountXC = len(sigmaArray)//sizeEarray
sigmaArray = sigmaArray.reshape(amountXC, sizeEarray)

pathOG = 'Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/OpenMP-BField/CollisionData/turner_benchmark_he_electron_table.dat'

elasticData = np.loadtxt(pathOG, skiprows = 11, max_rows = 182-12+1)

#elastic
for i in range(sizeEarray):
    testPoint = np.interp(E_array[i], elasticData[:,0], elasticData[:,1])
    otherPoint = sigmaArray[0,i]
    if (abs((testPoint - otherPoint)/testPoint) > 1e-8):
        print('Issue with testpoint!')

#ionization
ionizData = np.loadtxt(pathOG, skiprows = 627, max_rows = 828-627)
for i in range(sizeEarray):
    testPoint = np.interp(E_array[i], ionizData[:,0], ionizData[:,1])
    otherPoint = sigmaArray[1,i]
    if (abs((testPoint - otherPoint)/testPoint) > 1e-8):
        print('Issue with testpoint ionization!')

#excit 1
excitSingletData = np.loadtxt(pathOG, skiprows = 413, max_rows = 614-413)
for i in range(sizeEarray):
    testPoint = np.interp(E_array[i], excitSingletData[:,0], excitSingletData[:,1])
    otherPoint = sigmaArray[2,i]
    if (abs((testPoint - otherPoint)/testPoint) > 1e-8):
        print('Issue with testpoint excit singlet!')

# excit 2
excitTripletData = np.loadtxt(pathOG, skiprows=197, max_rows=398-197)
for i in range(sizeEarray):
    testPoint = np.interp(E_array[i], excitTripletData[:, 0], excitTripletData[:, 1])
    otherPoint = sigmaArray[3, i]
    if (abs((testPoint - otherPoint) / testPoint) > 1e-8):
        print('Issue with testpoint excit singlet!')

E_array_ion = np.fromfile(path + 'IncidentPart_He+_energy.dat')
sizeEarray_ion = len(E_array_ion)
sigmaArray_ion = np.fromfile(path + 'IncidentPart_He+_sigma.dat')
amountXC_ion = len(sigmaArray_ion)//sizeEarray_ion
sigmaArray_ion = sigmaArray_ion.reshape(amountXC_ion, sizeEarray_ion)

cexgData = np.loadtxt(pathOG, skiprows = 845, max_rows = 946-845)
for i in range(sizeEarray_ion):
    testPoint = np.interp(E_array_ion[i], cexgData[:, 0], cexgData[:, 1])
    otherPoint = sigmaArray_ion[0, i] * 1e20
    if (abs((testPoint - otherPoint) / testPoint) > 1e-8):
        print('Issue with testpoint ion Charge Exchange!')

ionElastData = np.loadtxt(pathOG, skiprows=960, max_rows=1061-960)
for i in range(sizeEarray_ion):
    testPoint = np.interp(E_array_ion[i], ionElastData[:, 0], ionElastData[:, 1])
    otherPoint = sigmaArray_ion[1, i] * 1e20
    if (abs((testPoint - otherPoint) / testPoint) > 1e-8):
        print('Issue with testpoint ion Charge Exchange!')

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