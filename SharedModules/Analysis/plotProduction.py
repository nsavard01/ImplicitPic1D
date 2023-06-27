# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:39:21 2023

@author: Nicolas
"""
from import_libraries_constants import *
from dataSet import *

test = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/testDiag/')

# ---------------------------- Averages -------------------------------------
def plotAveDensity(dataSet, name = "", label = ""):
    if name == "":
        colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']
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
            plt.plot(dataSet.grid, n,  linestyle = '-', marker = 'o', label = "")
            plt.xlabel('Distance (m)')
            plt.ylabel('Particle Density (1/m^3)')
            plt.xlim([0, dataSet.grid[-1]])
 
def plotAvePhi(dataSet):
    phi = dataSet.getAvePhi()
    plt.plot(dataSet.grid, phi, 'o-')
    plt.xlabel('Distance (m)')
    plt.ylabel('Potential (V)')
    plt.xlim([0, dataSet.grid[-1]]) 
    
def plotAveEVDF(dataSet):
    Vbins, VHist = dataSet.getAveEVDF()
    dv = Vbins[1] - Vbins[0]
    Norm = np.sum(VHist * dv)
    plt.plot(Vbins, VHist/Norm, 'o-')
    
def plotMaxwellV(v_x, T, m):
    f= np.sqrt(m/2/np.pi / e/ T) * np.exp(- m * v_x**2 / 2 / e/ T)
    plt.plot(v_x, f)
    
#-------------------------- Time Dependent ---------------------
    
def update_plot_Phi(i, dataSet, ax):
    ax.clear()
    phi = dataSet.getPhi(i)
    ax.plot(dataSet.grid, phi, 'o-')
        
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Potential (V)')
    ax.set_xlim([0, dataSet.grid[-1]])
    ax.set_ylim([-20, 40])

    
def phiAnimation(dataSet, boolMakeAnimation = False, savePath = "BoundPlasmaPhi.gif"):
    if boolMakeAnimation:
        numframes = dataSet.numDiag
        fig, ax = plt.subplots()
        phi = dataSet.getPhi(0)
        ax.plot(dataSet.grid, phi, 'o-')
            
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Potential (V)')
        ax.set_xlim([0, dataSet.grid[-1]])
        ax.set_ylim([-20, 40])
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
            
            
def update_plot_PhaseSpace(i, scat, dataSet, name):
    phaseSpace = dataSet.getPhaseSpace(name, i)
    scat.set_offsets(phaseSpace[:,0:2])
    return scat,

def PhaseSpaceAnimation(dataSet, name, boolMakeAnimation = False, savePath = 'PhaseSpaceAnimation.gif'): 
    if boolMakeAnimation:
        numframes = numDiagnosticTimes+1
        filename = '../Data/PhaseSpace/phaseSpace_'+ ParticleName + '_0.dat'
        phaseSpace = extractPhaseSpace(filename, grid)
        fig = plt.figure()
        scat = plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
        plt.xlabel('Distance (m)')
        plt.ylabel('Speed (m/s)')
        plt.axis([0, grid[-1], -6500, 6500])
        ani = animation.FuncAnimation(fig, update_plot_PhaseSpace, frames=range(numframes), interval = 100,
                                      fargs=(scat, ParticleName))
        ani.save('PostProcessing/twoStreamAnimation.gif')
        plt.show()
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(numDiagnosticTimes+1):
            plt.cla()
            filename = '../Data/PhaseSpace/phaseSpace_'+ ParticleName + '_' + str(y) +'.dat'
            phaseSpace = extractPhaseSpace(filename, grid)
            plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
            plt.xlabel('Distance (m)')
            plt.ylabel('Particle velocity (m/s)')
            plt.xlim([0, grid[-1]])
            plt.pause(0.1)