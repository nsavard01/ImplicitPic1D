# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 15:07:18 2022

@author: Nicolas
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
import scipy.optimize as opt
from scipy.stats import cosine
import scipy.sparse as sp
import time
import pandas as pd
import glob, os
import matplotlib.animation as animation
import itertools
import shutil
import sys
import BoundAnalytical




eps_0 = scipy.constants.epsilon_0
c = scipy.constants.c
m_e = scipy.constants.m_e
m_p = scipy.constants.m_p
mu_0 = scipy.constants.mu_0
k_boltz = scipy.constants.k
e = scipy.constants.e


dataFolder = ''
InitialConditions = None
GlobalDiagnostics = None
ParticleProperties = None
GlobalDiagnosticsAveraged = None
boolAverageFile = None
grid = None
delX = None
numDiagnosticTimes = None


def extractPhaseSpace(filename):
    phaseSpace = np.fromfile(filename, dtype = 'float', offset = 4)
    phaseSpace = phaseSpace.reshape((int(phaseSpace.size/4), 4))
    phaseSpace[:,0] = (phaseSpace[:,0]-1) * delX
    return phaseSpace



#------------------- Animation -------------------

def update_plot_Density(i, ax, ParticleProperties, colors):
    ax.clear()
    for j,name in enumerate(ParticleProperties['name']):
        n = np.fromfile(dataFolder + 'Density/density_' + name + '_' + str(i) + '.dat', dtype = 'float', offset = 4)
        ax.plot(grid, n, 'o-', c=colors[j], label = r'$n_{' + name[1:-1] +  '}$')
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Particle Density (1/m^3)')
    ax.set_xlim([0, grid[-1]])
    ax.set_ylim([0, 9e14])
    plt.legend(loc = 'lower center')
    
def densityAnimation(ParticleProperties, boolMakeAnimation):
    colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']
    if boolMakeAnimation:
        numframes = numDiagnosticTimes + 1
        fig, ax = plt.subplots()
        for j,name in enumerate(ParticleProperties['name']):
            n = np.fromfile(dataFolder + 'Density/density_' + name + '_0.dat', dtype = 'float', offset = 4)
            ax.plot(grid, n, 'o-', c = colors[j], label = r'$n_{' + name[1:-1] +  '}$')
            
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Particle Density (1/m^3)')
        ax.set_xlim([0, grid[-1]])
        ax.set_ylim([0, 9e14])
        plt.legend(loc = 'lower center')
        ani = animation.FuncAnimation(fig, update_plot_Density, frames=range(numframes), interval = 100,
                                      fargs=(ax, ParticleProperties, colors))
        
        ani.save('PostProcessing/BoundPlasmaDensity.gif')
        plt.show()
        
        
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(numDiagnosticTimes+1):
        
            plt.cla()
            
            for i,name in enumerate(ParticleProperties['name']):
                n = np.fromfile(dataFolder + 'Density/density_' + name + '_' + str(y) + '.dat', dtype = 'float', offset = 4)
                plt.plot(grid, n, colors[i], label = r'$n_{' + name[1:-1] +  '}$')
            plt.xlabel('Distance (m)')
            plt.ylabel('Particle Density (1/m^3)')
            plt.xlim([0, grid[-1]])
            plt.legend(loc = 'best')
            plt.pause(0.05)
            
    
def plotAverageDensity(ParticleProperties):
    colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']
    for i,name in enumerate(ParticleProperties['name']):
        n = np.fromfile(dataFolder + 'Density/density_' + name + '_Average.dat', dtype = 'float', offset = 4)
        plt.plot(grid, n,  linestyle = '-', marker = 'o', color = colors[i], label = r'$n_{' + name[1:-1] +  '}$')
    plt.xlabel('Distance (m)')
    plt.ylabel('Particle Density (1/m^3)')
    plt.xlim([0, grid[-1]])
    plt.legend(loc = 'best')
    plt.pause(0.05)
    
def plotAverageDensityParticle(ParticleProperties, particleName, label):
    n = np.fromfile(dataFolder + 'Density/density_' + particleName + '_Average.dat', dtype = 'float', offset = 4)
    plt.plot(grid, n,  linestyle = '-', marker = 'o', label = label)
    plt.xlabel('Distance (m)')
    plt.ylabel('Particle Density ' + particleName[1:-1] + ' (1/m^3)')
    plt.xlim([0, grid[-1]])
    plt.legend(loc = 'lower center')
    plt.pause(0.05)
 
def update_plot_PhaseSpace(i, scat, ParticleName):
    filename = dataFolder + 'PhaseSpace/phaseSpace_'+ ParticleName + '_' + str(i) +'.dat'
    phaseSpace = extractPhaseSpace(filename, grid)
    scat.set_offsets(phaseSpace[:,0:2])
    return scat,

def PhaseSpaceAnimation(ParticleName, boolMakeAnimation): 
    if boolMakeAnimation:
        numframes = numDiagnosticTimes+1
        filename = dataFolder + 'PhaseSpace/phaseSpace_'+ ParticleName + '_0.dat'
        phaseSpace = extractPhaseSpace(filename)
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
            filename = dataFolder + 'PhaseSpace/phaseSpace_'+ ParticleName + '_' + str(y) +'.dat'
            phaseSpace = extractPhaseSpace(filename)
            plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
            plt.xlabel('Distance (m)')
            plt.ylabel('Particle velocity (m/s)')
            plt.xlim([0, grid[-1]])
            plt.pause(0.1)
        
def TwoStreamEnergyIncrease():
    
    
    PE = np.zeros(numDiagnosticTimes+1)
    t = np.zeros(numDiagnosticTimes+1)
    t[1::] = GlobalDiagnostics['time(s)']
    for y in range(numDiagnosticTimes+1):
        phi = np.fromfile(dataFolder + 'Phi/phi_' + str(y) + '.dat', dtype = 'float', offset = 4)
        PE[y] = 0.0 #0.5 * eps_0 * np.sum((np.diff(phi)**2)) / delX
        for i in range(InitialConditions['n_x'].values[0]-1):
            PE[y] = PE[y] + ((phi[i-1] + phi[i + 1])/2/delX)**2
        PE[y] = 0.5 * eps_0 * PE[y]
        
        
    plasmaPeriod = InitialConditions['del_t'].values/InitialConditions['FracFreq'].values 
    PE_growth = PE[0] * np.exp( 0.5 * (1/plasmaPeriod) * t)

    plt.figure()
    plt.plot(t/plasmaPeriod, PE, label = 'Simulation')
    plt.plot(t/plasmaPeriod, PE_growth, label = 'Analytical')
    plt.plot()
    plt.yscale('log')
    plt.xlabel(r'Time (normalized to $\frac{1}{\omega_p}$)')
    plt.ylabel('Potential Energy (J)')
    plt.title(r'Two Stream Instability with $\Delta t$ = ' + "{:.{}f}".format(InitialConditions['FracFreq'].values[0] , 1) + r"$(\frac{1}{\omega_p})$")
    plt.legend(loc = 'best')
    plt.savefig('PostProcessing/twoStreamGrowth.png')
        
    
def update_plot_Phi(i, ax, upperLimit):
    ax.clear()
    phi = np.fromfile(dataFolder + 'Phi/phi_'+str(i) + '.dat', dtype = 'float', offset = 4)
    ax.plot(grid, phi, 'o-')
        
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Potential (V)')
    ax.set_xlim([0, grid[-1]])
    ax.set_ylim([0, upperLimit])
    
def phiAnimation(boolMakeAnimation, upperLimit):
    if boolMakeAnimation:
        numframes = numDiagnosticTimes + 1
        fig, ax = plt.subplots()
        phi = np.fromfile(dataFolder + 'Phi/phi_0.dat', dtype = 'float', offset = 4)
        ax.plot(grid, phi, 'o-')
            
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Potential (V)')
        ax.set_xlim([0, grid[-1]])
        ax.set_ylim([0, upperLimit])
        ani = animation.FuncAnimation(fig, update_plot_Phi, frames=range(numframes), interval = 100,fargs=(ax,upperLimit))
        
        ani.save('PostProcessing/BoundPlasmaPhi.gif')
        plt.show()
        
        
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(numDiagnosticTimes+1):
        
            plt.cla()
            
            phi = np.fromfile(dataFolder + 'Phi/phi_' + str(y) + '.dat', dtype = 'float', offset = 4)
            plt.plot(grid, phi, 'o-')
            plt.xlabel('Distance (m)')
            plt.ylabel('Potential (V)')
            plt.xlim([0, grid[-1]])
            plt.pause(0.05)        
        
def temperatureAnimation(boolMakeAnimation):
    if boolMakeAnimation:
        numframes = numDiagnosticTimes + 1
        fig, ax = plt.subplots()
        phaseSpace = extractPhaseSpace(dataFolder + 'PhaseSpace/phaseSpace_[e]_' + str(i) +'.dat')
        KE = np.sum(phaseSpace[:, 1::]**2, axis = 1) * 0.5 * m_e / e
        E_binned = np.histogram(phaseSpace[:,0], bins = grid, weights = KE)[0]
        num_binned = np.clip(np.histogram(phaseSpace[:,0], bins = grid)[0], a_min = 1, a_max = None)
        ax.plot(halfGrid, E_binned*2.0/num_binned/3.0, 'o-')
            
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel(r'$T_e$ (eV)')
        ax.set_xlim([0, grid[-1]])
        ax.set_ylim([0, 7])
        ani = animation.FuncAnimation(fig, update_plot_Phi, frames=range(numframes), interval = 100,fargs=(ax,))
        
        ani.save('PostProcessing/BoundPlasmaTemp.gif')
        plt.show()
        
        
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(numDiagnosticTimes+1):
        
            plt.cla()
            
            phaseSpace = extractPhaseSpace(dataFolder + 'PhaseSpace/phaseSpace_[e]_' + str(y) +'.dat')
            KE = np.sum(phaseSpace[:, 1::]**2, axis = 1) * 0.5 * m_e / e
            E_binned = np.histogram(phaseSpace[:,0], bins = grid, weights = KE)[0]
            num_binned = np.clip(np.histogram(phaseSpace[:,0], bins = grid)[0], a_min = 1, a_max = None)
            plt.plot(grid[0:-1], E_binned*2.0/num_binned/3.0, 'o-')
            plt.xlabel('Distance (m)')
            plt.ylabel(r'$T_e$ (eV)')
            plt.xlim([0, grid[-1]])
            plt.pause(0.05)        

    
def plotAveragePhi():
    phi = np.fromfile(dataFolder + 'Phi/phi_Average.dat', dtype = 'float', offset = 4)
    plt.plot(grid, phi, 'o-')
    plt.xlabel('Distance (m)')
    plt.ylabel('Potential (V)')
    plt.xlim([0, grid[-1]])
    plt.pause(0.05)        

def getData(filename):
    global dataFolder
    global InitialConditions
    global GlobalDiagnostics
    global ParticleProperties
    global GlobalDiagnosticsAveraged
    global boolAverageFile
    global grid
    global delX
    global numDiagnosticTimes
    
    
    dataFolder = '../' + filename + '/'
    diagList = ['time(s)', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)']
    diagAverageList = ['steps', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)']
    InitialConditionsList = ['n_x', 't_final', 'del_t', 'FracFreq', 'delX', 'Power', 'heatSteps', 'nu_h']
    # list names in diagnostic file for each time steps
    InitialConditions =  pd.read_csv(dataFolder + 'InitialConditions.dat', skiprows = 1, delim_whitespace=True, names = InitialConditionsList)
    GlobalDiagnostics = pd.read_csv(dataFolder + 'GlobalDiagnosticData.dat', skiprows = 1, delim_whitespace=True, names = diagList)
    
    ParticleProperties = pd.read_csv(dataFolder + 'ParticleProperties.dat', skiprows = 1, names = ['name', 'mass', 'q', 'w_p'], delim_whitespace = True)
    boolAverageFile = os.path.isfile(dataFolder + 'GlobalDiagnosticDataAveraged.dat')
    GlobalDiagnosticsAveraged = pd.DataFrame()
    if (boolAverageFile):
        GlobalDiagnosticsAveraged = pd.read_csv(dataFolder + 'GlobalDiagnosticDataAveraged.dat', skiprows = 1, delim_whitespace=True, names = diagAverageList)
    
    grid = np.fromfile(dataFolder + 'domainGrid.dat', dtype = 'float', offset = 4)
    delX = InitialConditions['delX'].values[0]
    numDiagnosticTimes = GlobalDiagnostics.shape[0]

def saveData(saveFile):
    shutil.copytree('../Data', '../' + saveFile)
    
# ---------------------- Final Plots --------------------------

def finalElectronEEDFVsMaxwellian():
    phaseSpace = extractPhaseSpace(dataFolder + 'PhaseSpace/phaseSpace_[e]_' + str(numDiagnosticTimes) +'.dat')
    KE = np.sum(phaseSpace[:, 1::]**2, axis = 1) * 0.5 * m_e / e
    Ehist = np.histogram(KE, bins = 100, density = True)
    T_elec = np.mean(KE)* (2/3)
    Ebins = (Ehist[1][0:-1] + Ehist[1][1::])/2
    plt.figure()
    plt.title('Mean electron temperature' + '{:10.4f}'.format(T_elec) +  ' eV')
    plt.plot(Ebins, Ehist[0]/np.sqrt(Ebins), label = 'Binned Electrons')
    plt.plot(Ebins, 2/np.sqrt(np.pi) * (1/T_elec)**(3/2) * np.exp(-Ebins/T_elec), label = r'Maxwellian PDF')
    plt.yscale('log')
    plt.ylabel(r'PDF ($\frac{1}{eV^{3/2}}$)')
    plt.xlabel('Electron Energy (eV)')
    plt.legend(loc = 'best')


getData('Data')
if(boolAverageFile):
    print('Total power is:', GlobalDiagnosticsAveraged['Ploss(W/m^2)'] +GlobalDiagnosticsAveraged['P_wall(W/m^2)'] )









