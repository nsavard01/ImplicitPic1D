# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

data = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_100PPC_2p0delT_2p5e16nave_16nodes/')
data1 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_200PPC_2p0delT_2p5e16nave_16nodes/')
data2 = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/NGP_300PPC_2p0delT_2p5e16nave_16nodes/')
dataExplicit = dataSetExplicit('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/ExplicitPIC/Explicit_100PPC_0p2delT_2p5e16nave/')
dataList = [data, data1, data2, dataExplicit]
labelList = ['100PPC', '200PPC', '300PPC', '100PPC - Exp.']
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
    for i in range(len(dataList)):
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        ax1.plot(data.grid, n_e, linestyle = '-', marker = '.',label = labelList[i])
        
        n_i = data.getAveDensity(ion)
        ax2.plot(data.grid, n_i, linestyle = '-', marker = '.',label = labelList[i])
        
        phi = data.getAvePhi()
        ax3.plot(data.grid, phi, linestyle = '-', marker = '.',label = labelList[i])
     
    model = getBoundPlasmaSolutions(dataList[0].grid[-1] - dataList[0].grid[0], 100, n_ave, T_e, T_i, M)
    ax1.plot(model[0], model[2], label = 'Model')
    ax2.plot(model[0], model[3], label = 'Model')
    ax3.plot(model[0], model[1], label = 'Model')
    ax1.legend(loc = 'best')
    ax2.legend(loc = 'best')
    ax3.legend(loc = 'best')
    ax1.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax2.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax3.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
        