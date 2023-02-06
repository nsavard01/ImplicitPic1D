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


eps_0 = scipy.constants.epsilon_0
c = scipy.constants.c
m_e = scipy.constants.m_e
m_p = scipy.constants.m_p
mu_0 = scipy.constants.mu_0
k_boltz = scipy.constants.k
e = scipy.constants.e

def extractPhaseSpace(filename, grid):
    phaseSpace = np.fromfile(filename, dtype = 'float', offset = 4)
    phaseSpace = phaseSpace.reshape((int(phaseSpace.size/4), 4))
    d = phaseSpace[:,0] - phaseSpace[:,0].astype(int)
    phaseSpace[:,0] = grid[phaseSpace[:,0].astype(int)-1] + d * (grid[phaseSpace[:,0].astype(int)] - grid[phaseSpace[:,0].astype(int)-1])
    return phaseSpace

# name for each initial condition simulation
n_x = int(np.loadtxt('Data/InitialConditions.dat', skiprows = 1)[0])
t_final = np.loadtxt('Data/InitialConditions.dat', skiprows = 1)[1]
# list names in diagnostic file for each time steps
diagList = ['time', 'Ploss', 'I_wall', 'P_wall']

GlobalDiagnostics = np.loadtxt('Data/GlobalDiagnosticData.dat', skiprows = 1)
numDiagnosticTimes = GlobalDiagnostics.shape[0]


grid = np.fromfile('Data/domainGrid.dat', dtype = 'float', offset = 4)
halfGrid = (grid[0:-1] + grid[1::])/2
ne = np.zeros(n_x)
ni = np.zeros(n_x)

test_adaptive = np.fromfile('test_adaptive.dat', dtype = 'float', offset = 4)
test_normal = np.fromfile('test_normal.dat', dtype = 'float', offset = 4)

#------------------- Animation -------------------

plt.figure(figsize = (5,4), dpi = 80)
for y in range(numDiagnosticTimes+1):

    plt.cla()
    ne = np.fromfile('Data/Density/density_e_' + str(y) + '.dat', dtype = 'float', offset = 4)
    ni = np.fromfile('Data/Density/density_H+_' + str(y) + '.dat', dtype = 'float', offset = 4)
    plt.plot(grid, ne, 'b', label = r'$n_e$')
    plt.plot(grid, ni, 'r', label = r'$n_i$')
    plt.xlabel('Distance (m)')
    plt.ylabel('Particle Density (1/m^3)')
    plt.xlim([0, grid[-1]])
    plt.legend(loc = 'best')
    plt.pause(0.05)
 
    
# plt.figure(figsize = (5,4), dpi = 80)
# for y in range(numDiagnosticTimes+1):
#     plt.cla()
#     filename = 'Data/PhaseSpace/phaseSpace_e_' + str(y) +'.dat'
#     phaseSpace = extractPhaseSpace(filename, grid)
#     plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
#     plt.xlabel('Distance (m)')
#     plt.ylabel('Particle velocity (m/s)')
#     plt.xlim([0, grid[-1]])
#     plt.pause(0.1)
    
plt.figure(figsize = (5,4), dpi = 80)
for y in range(numDiagnosticTimes+1):
    
    plt.cla()
    phi = np.fromfile('Data/Phi/phi_' + str(y) + '.dat', dtype = 'float', offset = 4)
    plt.plot(grid, phi, 'o-')
    plt.xlabel('Distance (m)')
    plt.ylabel('Potential (V)')
    plt.xlim([0, grid[-1]])
    plt.pause(0.1)
    
# plt.figure(figsize = (5,4), dpi = 80)
# for y in range(numDiagnosticTimes+1):
    
#     plt.cla()
#     T_e = np.fromfile('Data/ElectronTemperature/eTemp_' + str(y) + '.dat', dtype = 'float', offset = 4)
#     plt.plot(halfGrid, T_e,'o-')
#     plt.xlabel('Distance (m)')
#     plt.ylabel('Electron Temperature (eV)')
#     plt.xlim([0, grid[-1]])
#     plt.pause(0.05)


    
    
# ---------------------- Final Plots --------------------------
phaseSpace = extractPhaseSpace('Data/PhaseSpace/phaseSpace_e_' + str(numDiagnosticTimes) +'.dat', grid)
KE = np.sum(phaseSpace[:, 1::]**2, axis = 1) * 0.5 * m_e / e
Ehist = np.histogram(KE, bins = 100, density = True)
T_elec = np.mean(KE)* (2/3)
Ebins = (Ehist[1][0:-1] + Ehist[1][1::])/2
plt.figure()
plt.title('Mean electron temperature' + '{:10.4f}'.format(T_elec) +  ' eV')
plt.plot(Ebins, Ehist[0]/np.sqrt(Ebins), label = 'Binned Electrons')
plt.plot(Ebins, 2/np.sqrt(np.pi) * (1/T_elec)**(3/2) * np.exp(-Ebins/T_elec), label = r'Maxwellian PDF')
plt.yscale('log')
plt.legend(loc = 'best')







