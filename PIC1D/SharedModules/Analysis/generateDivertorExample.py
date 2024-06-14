# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:07:07 2023

@author: Nicolas
"""

from import_libraries_constants import *


def compareRefToDatasRes(dataList, labelList, ref, saveFile):
    if (len(dataList) != len(labelList)):
        raise Warning("List of data and labels not the same size!")
        
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    for name in dataList[0].particles.keys():
        if name != 'e':
            ion = name
            break
    M = dataList[0].particles[ion]['mass']
    T_e = dataList[0].T_e
    T_i = dataList[0].T_i
    n_ave = dataList[0].n_ave
    for i in range(len(dataList)):
        print('Statistics for label', labelList[i])
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        n_e_ref = ref.getAveDensity('e')
        n_e_ref = np.interp(data.grid, ref.grid, n_e_ref)
        ax1.plot(data.grid, 100 * (n_e-n_e_ref)/n_e_ref, linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        print(r'$n_e$ is ', 100 * np.sum((n_e[1:-1]-n_e_ref[1:-1])/n_e_ref[1:-1])/(dataList[i].Nx - 2))
        
        n_i = data.getAveDensity(ion)
        n_i_ref = ref.getAveDensity('H+')
        n_i_ref = np.interp(data.grid, ref.grid, n_i_ref)
        ax2.plot(data.grid, 100 * (n_i-n_i_ref)/n_i_ref, linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        print(r'$n_i$ is ', 100 * np.sum((n_i[1:-1]-n_i_ref[1:-1])/n_i_ref[1:-1])/(dataList[i].Nx - 2))
        
        phi = data.getAvePhi()
        phi_ref = ref.getAvePhi()
        phi_ref = np.interp(data.grid, ref.grid, phi_ref)
        ax3.plot(data.grid, 100 * (phi-phi_ref)/(phi_ref + 1e-15), linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        print(r'$\phi$ is ', 100 * np.sum((phi[1:-1]-phi_ref[1:-1])/phi_ref[1:-1])/(dataList[i].Nx - 2))
        
        ax4.bar(labelList[i], data.totPotTime)
        ax4.text(i, data.totPotTime, data.totPotTime, ha = 'center')

    ax1.legend(loc = 'upper right', fontsize = 12)
    ax2.legend(loc = 'lower right', fontsize = 12)
    ax3.legend(loc = 'lower right', fontsize = 12)
    ax1.set_xlim(ref.grid[0], ref.grid[-1])
    ax1.set_xlabel('Distance (m)', fontsize = 14)
    ax1.set_ylabel(r'100 $\times$ $(n_e - n_{e,ref})/n_{e,ref}$', fontsize = 14)
    ax1.tick_params(labelsize = 14)
    
    ax2.set_xlim(ref.grid[0], ref.grid[-1])
    ax2.set_xlabel('Distance (m)', fontsize = 14)
    ax2.set_ylabel(r'100 $\times$ $(n_i - n_{i,ref})/n_{i,ref}$', fontsize = 14)
    ax2.tick_params(labelsize = 14)
    
    ax3.set_xlim(ref.grid[0], ref.grid[-1])
    ax3.set_xlabel('Distance (m)', fontsize = 14)
    ax3.set_ylabel(r'100 $\times$ $(\phi - \phi_{ref})/\phi_{ref}$', fontsize = 14)
    ax3.tick_params(labelsize = 14)
    
    ax4.set_ylabel(r'Total Wall Time for Solver (s)', fontsize = 14)
    
    fig1.savefig(saveFile + '_eDensity.pdf', bbox_inches = 'tight')
    fig2.savefig(saveFile + '_ionDensity.pdf', bbox_inches = 'tight')
    fig3.savefig(saveFile + '_voltage.pdf', bbox_inches = 'tight')
    fig4.savefig(saveFile + '_time.pdf', bbox_inches = 'tight')
    
def compareRefToDatasAbsRes(dataList, labelList, ref):
    if (len(dataList) != len(labelList)):
        raise Warning("List of data and labels not the same size!")
        
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    for name in dataList[0].particles.keys():
        if name != 'e':
            ion = name
            break
    M = dataList[0].particles[ion]['mass']
    T_e = dataList[0].T_e
    T_i = dataList[0].T_i
    n_ave = dataList[0].n_ave
    for i in range(len(dataList)):
        print('Statistics for label', labelList[i])
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        n_e_ref = ref.getAveDensity('e')
        n_e_ref = np.interp(data.grid, ref.grid, n_e_ref)
        ax1.plot(data.grid, 100 * abs((n_e-n_e_ref)/n_e_ref), linestyle = '--', marker = '.',label = labelList[i])
        print(r'$n_e$ is ', 100 * np.sum((n_e[1:-1]-n_e_ref[1:-1])/n_e_ref[1:-1])/(dataList[i].Nx - 2))
        
        n_i = data.getAveDensity(ion)
        n_i_ref = ref.getAveDensity('H+')
        n_i_ref = np.interp(data.grid, ref.grid, n_i_ref)
        ax2.plot(data.grid, 100 * abs((n_i-n_i_ref)/n_i_ref), linestyle = '--', marker = '.',label = labelList[i])
        print(r'$n_i$ is ', 100 * np.sum((n_i[1:-1]-n_i_ref[1:-1])/n_i_ref[1:-1])/(dataList[i].Nx - 2))
        
        phi = data.getAvePhi()
        phi_ref = ref.getAvePhi()
        phi_ref = np.interp(data.grid, ref.grid, phi_ref)
        ax3.plot(data.grid, 100 * abs((phi-phi_ref)/(phi_ref + 1e-15)), linestyle = '--', marker = '.',label = labelList[i])
        print(r'$\phi$ is ', 100 * np.sum((phi[1:-1]-phi_ref[1:-1])/phi_ref[1:-1])/(dataList[i].Nx - 2))
        
        ax4.bar(labelList[i], data.totPotTime)
        ax4.text(i, data.totPotTime, data.totPotTime, ha = 'center')

    ax1.legend(loc = 'best')
    ax2.legend(loc = 'best')
    ax3.legend(loc = 'best')
    ax1.set_xlim(ref.grid[0], ref.grid[-1])
    ax1.set_xlabel('Distance (m)', fontsize = 12)
    ax1.set_ylabel(r'100 $\times$ $|(n_e - n_{e,ref})/n_{e,ref}|$', fontsize = 12)
    
    ax2.set_xlim(ref.grid[0], ref.grid[-1])
    ax2.set_xlabel('Distance (m)', fontsize = 12)
    ax2.set_ylabel(r'100 $\times$ $|(n_i - n_{i,ref})/n_{i,ref}|$', fontsize = 12)
    
    ax3.set_xlim(ref.grid[0], ref.grid[-1])
    ax3.set_xlabel('Distance (m)', fontsize = 12)
    ax3.set_ylabel(r'100 $\times$ $|(\phi - \phi_{ref})/\phi_{ref}|$', fontsize = 12)
    
    ax4.set_ylabel(r'Total Wall Time for Solver (s)', fontsize = 12)
    
def compareModelToDatas(dataList, labelList, model, saveFile):
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
    factor = int(math.log10(n_ave))
    for i in range(len(dataList)):
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        ax1.plot(data.grid, n_e/(10**factor), linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        
        n_i = data.getAveDensity(ion)
        ax2.plot(data.grid, n_i/(10**factor), linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        
        phi = data.getAvePhi()
        ax3.plot(data.grid, phi, linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
     
    
    ax1.plot(model[0, :], model[2, :]/(10**factor), linewidth = 2, label = 'Model')
    ax2.plot(model[0, :], model[3, :]/(10**factor), linewidth = 2, label = 'Model')
    ax3.plot(model[0, :], model[1, :], linewidth = 2, label = 'Model')
    ax1.legend(loc = 'lower right', fontsize = 12)
    ax2.legend(loc = 'lower right', fontsize = 12)
    ax3.legend(loc = 'lower right', fontsize = 12)
    ax1.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax1.set_xlabel('Distance (m)', fontsize = 14)
    ax1.set_ylabel(r'$n_e$ (10$^{' + str(factor) + r'}$ m$^{-3}$)', fontsize = 14)
    ax1.tick_params(labelsize = 14)
    
    ax2.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax2.set_xlabel('Distance (m)', fontsize = 14)
    ax2.set_ylabel(r'$n_i$ (10$^{' + str(factor) + r'}$ m$^{-3}$)', fontsize = 14)
    ax2.tick_params(labelsize = 14)
    
    ax3.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax3.set_xlabel('Distance (m)', fontsize = 14)
    ax3.set_ylabel('Voltage (V)', fontsize = 14)
    ax3.tick_params(labelsize = 14)
    
    fig1.savefig(saveFile + '_eDensity.pdf', bbox_inches = 'tight')
    fig2.savefig(saveFile + '_ionDensity.pdf', bbox_inches = 'tight')
    fig3.savefig(saveFile + '_voltage.pdf', bbox_inches = 'tight')
    
def comparePhiToLaplace(data, angle):

    
    fileName = 'Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/OpenMP-BField/SharedModules/Analysis/DivertorComparison/profiles_LePIC25D3V/profiles_LePIC25D3V_alpha_' + str(angle) + '_t_4000000.txt'
    extrData = np.loadtxt(fileName, skiprows = 2)
    
    plt.figure()
    phi = data.getAvePhi()
    plt.plot(data.grid * 1e3, phi, 'o', label = 'Nicolas')
    plt.plot(extrData[:,0], extrData[:, 1], label = 'Laplace')
    plt.ylabel(r'Voltage (V)')
    plt.legend(loc = 'best')
        
def compareModelToDatasRes(dataList, labelList, model):
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
        modelN_e = np.interp(data.grid, model[0,:], model[2,:])
        ax1.plot(data.grid, 100 * (n_e-modelN_e)/modelN_e, linestyle = '--', marker = '.',label = labelList[i])
        
        n_i = data.getAveDensity(ion)
        modelN_i = np.interp(data.grid, model[0,:], model[3,:])
        ax2.plot(data.grid, 100 * (n_i-modelN_i)/modelN_i, linestyle = '--', marker = '.',label = labelList[i])
        
        phi = data.getAvePhi()
        modelPhi = np.interp(data.grid, model[0,:], model[1,:])
        ax3.plot(data.grid, 100 * (phi-modelPhi)/modelPhi, linestyle = '--', marker = '.',label = labelList[i])
     
    
    ax1.legend(loc = 'best')
    ax2.legend(loc = 'best')
    ax3.legend(loc = 'best')
    ax1.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax1.set_xlabel('Distance (m)', fontsize = 12)
    ax1.set_ylabel(r'100 $\times$ $(n_e - n_{e,ref})/n_{e,ref}$', fontsize = 12)
    
    ax2.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax2.set_xlabel('Distance (m)', fontsize = 12)
    ax2.set_ylabel(r'100 $\times$ $(n_i - n_{i,ref})/n_{i,ref}$', fontsize = 12)
    
    ax3.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax3.set_xlabel('Distance (m)', fontsize = 12)
    ax3.set_ylabel(r'100 $\times$ $(\phi - \phi_{ref})/\phi_{ref}$', fontsize = 12)   