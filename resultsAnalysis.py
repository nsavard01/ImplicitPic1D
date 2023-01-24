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

# name for each initial condition simulation
n_x = int(np.loadtxt('Data/InitialConditions.dat', skiprows = 1)[0])
t_final = np.loadtxt('Data/InitialConditions.dat', skiprows = 1)[1]
# list names in diagnostic file for each time steps
diagList = ['time', 'Eloss', 'I_wall', 'P_wall']

ScalarDiagnostics = np.loadtxt('Data/ScalarDiagnosticData.dat', skiprows = 1)
numDiagnosticTimes = ScalarDiagnostics.shape[0]

grid = np.fromfile('Data/domainGrid.dat', dtype = 'float', offset = 4)
ne = np.zeros(n_x)
ni = np.zeros(n_x)


#------------------- Animation -------------------

plt.figure(figsize = (5,4), dpi = 80)
for y in range(numDiagnosticTimes+1):
    print(y)
    plt.cla()
    ne = np.fromfile('/mnt/SuperComputer/home/nsavard/ImplicitPic1D/ImplicitPic1D/Data/density_e_' + str(y) + '.dat', dtype = 'float', offset = 4)
    ni = np.fromfile('/mnt/SuperComputer/home/nsavard/ImplicitPic1D/ImplicitPic1D/Data/density_H+_' + str(y) + '.dat', dtype = 'float', offset = 4)
    plt.plot(grid, ne, 'b', label = r'$n_e$')
    plt.plot(grid, ni, 'r', label = r'$n_i$')
    plt.xlabel('Distance (m)')
    plt.ylabel('Particle Density (1/m^3)')
    plt.xlim([0, grid[-1]])
    plt.legend(loc = 'best')
    plt.pause(0.1)






