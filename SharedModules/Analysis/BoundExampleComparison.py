# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

test = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/ExplicitPIC/test/')
data = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_100PPC_2p0delT_2p5e16nave_16nodes/')
data1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_200PPC_2p0delT_2p5e16nave_16nodes/')
data31Nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_100PPC_2p0delT_2p5e16nave_31nodes/')
data2 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_300PPC_2p0delT_2p5e16nave_16nodes/')
data3 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_100PPC_0p2delT_2p5e16nave_16nodes/')
data4 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_500PPC_2p0delT_2p5e16nave_16nodes/')
data1000PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_1000PPC_2p0delT_2p5e16nave_16nodes/')
data2000PPC = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_2000PPC_2p0delT_2p5e16nave_16nodes/')
data1000PPC_32Nodes = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_1000PPC_2p0delT_2p5e16nave_32nodes/')
data1e18 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_100PPC_2p0delT_1e18nave_16nodes/')
dataExplicit100PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/ExplicitPIC/Explicit_100PPC_0p2delT_2p5e16nave/')
dataExplicit300PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/ExplicitPIC/Explicit_300PPC_0p2delT_2p5e16nave/')
dataExplicit500PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/ExplicitPIC/Explicit_500PPC_0p2delT_2p5e16nave/')
dataExplicit1000PPC = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/ExplicitPIC/Explicit_1000PPC_0p2delT_2p5e16nave/')
dataList = [dataExplicit300PPC, dataExplicit500PPC, dataExplicit1000PPC, data4, data1000PPC, data2000PPC, data1000PPC_32Nodes]
labelList = ['300PPC-Explicit', '500PPC-Explicit', '1000PPC-Explicit', '500PPC-NGP', '1000PPC-NGP', '2000PPC-NGP', '1000PPC-NGP-32nodes']
def compareModelToData(data):
    for name in data.particles.keys():
        if name != 'e':
            ion = name
            break
    M = data.particles[ion]['mass']
    model = getBoundPlasmaSolutions(data.grid[-1] - data.grid[0], 100, data.n_ave, data.T_e, data.T_i, M)
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
    
    
def compareModeltoDatas(dataList, labelList):
    if (len(dataList) != len(labelList)):
        raise Warning("List of data and labels not the same size!")
        
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    for name in dataList[0].particles.keys():
        if name != 'e':
            ion = name
            break
    M = dataList[0].particles[ion]['mass']
    T_e = dataList[0].T_e
    T_i = dataList[0].T_i
    n_ave = dataList[0].n_ave
    model = getBoundPlasmaSolutions(dataList[0].grid[-1] - dataList[0].grid[0], 100, n_ave, T_e, T_i, M)
    for i in range(len(dataList)):
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        ax1.plot(data.grid, n_e, linestyle = '-', marker = '.',label = labelList[i])
        
        n_i = data.getAveDensity(ion)
        ax2.plot(data.grid, n_i, linestyle = '-', marker = '.',label = labelList[i])
        
        phi = data.getAvePhi()
        ax3.plot(data.grid, phi, linestyle = '-', marker = '.',label = labelList[i])
     
    
    ax1.plot(model[0], model[2], label = 'Model')
    ax2.plot(model[0], model[3], label = 'Model')
    ax3.plot(model[0], model[1], label = 'Model')
    ax1.legend(loc = 'best')
    ax2.legend(loc = 'best')
    ax3.legend(loc = 'best')
    ax1.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel(r'$n_e$ (m$^{-3}$)')
    
    ax2.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel(r'$n_+$ (m$^{-3}$)')
    
    ax3.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax3.set_xlabel('Distance (m)')
    ax3.set_ylabel('Voltage (V)')
        