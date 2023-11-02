# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:39:21 2023

@author: Nicolas
"""
from import_libraries_constants import *
from dataSet import *
from dataSetExplicit import *


# ---------------------------- Averages -------------------------------------
def plotAveDensity(dataSet, name = "", label = ""):
    if name == "":
        colors = ['b', 'r', 'g', 'k', 'c', 'm', 'yAve']
        for i,name in enumerate(dataSet.particles.keys()):
            n = dataSet.getAveDensity(name)
            plt.plot(dataSet.grid, n,  linestyle = '-', marker = 'o', color = colors[i], label = r'$n_{' + name +  '}$')
        plt.xlabel('Distance (m)')
        plt.ylabel('Particle Density (1/m^3)')
        plt.xlim([0, dataSet.grid[-1]])
        plt.legend(loc = 'best')
    else:
        if name not in dataSet.particles.keys():
            raise Warning("For average density, particle", name, "does not exist in the dataSet!")
        else:
            n = dataSet.getAveDensity(name)
            plt.plot(dataSet.grid, n,  linestyle = '-', marker = 'o', label = label)
            plt.xlabel('Distance (m)')
            plt.ylabel(name + ' Density (1/m^3)')
            plt.xlim([0, dataSet.grid[-1]])
 
def plotAvePhi(dataSet, label = ''):
    phi = dataSet.getAvePhi()
    plt.plot(dataSet.grid, phi, 'o-', label = label)
    plt.xlabel('Distance (m)')
    plt.ylabel('Potential (V)')
    plt.xlim([0, dataSet.grid[-1]])
    plt.show()
    
def maxwellEDVF(x, T):
    return np.sqrt(m_e/2/np.pi / e/ T) * np.exp(- m_e * x**2 / 2 / e/ T)  
    
def plotAveEVDF(dataSet):
    Vbins, VHist = dataSet.getAveEVDF()
    dv = Vbins[1] - Vbins[0]
    Norm = np.sum(VHist * dv)
    VHist = VHist/Norm
    popt, pcov = opt.curve_fit(maxwellEDVF, Vbins, VHist, p0 = [dataSet.T_e])
    plt.plot(Vbins, VHist, 'o-', label = 'Global EDVF')
    plt.plot(Vbins, maxwellEDVF(Vbins, popt[0]), label = 'Best Fit T_e = ' + '{:2.2f}'.format(popt[0]))
    plt.legend(loc = 'best')
    

    
#-------------------------- Time Dependent ---------------------
    
def update_plot_Phi(i, dataSet, ax):
    ax.clear()
    phi = dataSet.getPhi(i)
    ax.plot(dataSet.grid, phi, 'o-')
        
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Potential (V)')
    ax.set_xlim([0, dataSet.grid[-1]])
    #ax.set_ylim([-20, 40])

    
def phiAnimation(dataSet, boolMakeAnimation = False, savePath = "Figures/BoundPlasmaPhi.gif"):
    if boolMakeAnimation:
        numframes = dataSet.numDiag
        fig, ax = plt.subplots()
        phi = dataSet.getPhi(0)
        ax.plot(dataSet.grid, phi, 'o-')
            
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Potential (V)')
        ax.set_xlim([0, dataSet.grid[-1]])
        #ax.set_ylim([-20, 40])
        ani = animation.FuncAnimation(fig, update_plot_Phi, frames=range(numframes), interval = 100,fargs=(dataSet, ax))
        
        ani.save(savePath)
        plt.show()
        
        
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(dataSet.numDiag):
        
            plt.cla()
            
            phi = dataSet.getPhi(y)
            plt.plot(dataSet.grid, phi, 'o-')
            plt.xlabel('Distance (m)')
            plt.ylabel('Potential (V)')
            plt.xlim([0, dataSet.grid[-1]])
            plt.pause(0.05)  

def update_plot_Density(i, dataSet, ax, nameList):
    ax.clear()
    for name in nameList:
        n = dataSet.getDensity(name, i)
        ax.plot(dataSet.grid, n, 'o-', label = name)
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Particle Density (1/m^3)')
    ax.set_xlim([0, dataSet.grid[-1]])
    plt.legend(loc = 'lower center')            
            
def densityAnimation(dataSet, nameList,boolMakeAnimation = False, savePath = "Figures/BoundPlasmaDensity.gif"):
    if boolMakeAnimation:
        numframes = dataSet.numDiag
        fig, ax = plt.subplots()
        for name in nameList:
            n = dataSet.getDensity(name, 0)
            ax.plot(dataSet.grid, n, 'o-', label = name)
            
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel(r'Density (m$^{-3}$)')
        ax.set_xlim([0, dataSet.grid[-1]])
        #ax.set_ylim([-20, 40])
        ani = animation.FuncAnimation(fig, update_plot_Density, frames=range(numframes), interval = 100,fargs=(dataSet, ax, nameList))
        
        ani.save(savePath)
        plt.show()
        
        
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(dataSet.numDiag):
        
            plt.cla()
            
            for name in nameList:
                n = dataSet.getDensity(name, y)
                plt.plot(dataSet.grid, n, 'o-', label = name)
            plt.xlabel('Distance (m)')
            plt.ylabel(r'Density (m$^{-3}$)')
            plt.xlim([0, dataSet.grid[-1]])
            plt.pause(0.05)  
            
            
def update_plot_PhaseSpace(i, scat, dataSet, name):
    phaseSpace = dataSet.getPhaseSpace(name, i)
    scat.set_offsets(phaseSpace[:,0:2])
    return scat,

def PhaseSpaceAnimation(dataSet, name, boolMakeAnimation = False, savePath = 'PhaseSpaceAnimation.gif'): 
    if name not in dataSet.particles.keys():
        raise Warning("Particle", name, "does not exist in the dataSet!")
    if boolMakeAnimation:
        numframes = dataSet.numDiag
        phaseSpace = dataSet.getPhaseSpace(name, 0)
        fig = plt.figure()
        scat = plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
        plt.xlabel('Distance (m)')
        plt.ylabel('Speed (m/s)')
        vmax = abs(np.max(phaseSpace[:,1]))*2.0
        plt.axis([dataSet.grid[0], dataSet.grid[-1], -vmax, vmax])
        ani = animation.FuncAnimation(fig, update_plot_PhaseSpace, frames=range(numframes), interval = 100,
                                      fargs=(scat, dataSet, name))
        ani.save(savePath)
        plt.show()
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(dataSet.numDiag):
            plt.cla()
            phaseSpace = dataSet.getPhaseSpace(name, y)
            plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
            plt.xlabel('Distance (m)')
            plt.ylabel('Particle velocity (m/s)')
            plt.xlim([0, dataSet.grid[-1]])
            plt.pause(0.1)
            
    
 
#------------ Bounded Model comparison ---------------------------            