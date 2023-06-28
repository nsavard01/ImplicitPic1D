# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *

data = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/test_1e18_1000PartPerCell/')

def compareModelToData(data):
    for name in data.particles.keys():
        if name != 'e':
            ion = name
            break
    M = data.particles[ion]['mass']
    model = getBoundPlasmaSolutions(data.grid[-1] - data.grid[0], 100, data.n_ave, data.T_e, data.T_i, M)
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