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

eps_0 = scipy.constants.epsilon_0
c = scipy.constants.c
m_e = scipy.constants.m_e
m_p = scipy.constants.m_p
mu_0 = scipy.constants.mu_0
k_boltz = scipy.constants.k
e = scipy.constants.e

# l_p = np.fromfile('record_particlePosition.dat', dtype = float)
# v_p = np.fromfile('record_particleVelocity.dat', dtype = float)
# v_p = v_p.reshape((3, int(v_p.size/3)))
rho = np.fromfile('record_Rho.dat', dtype = float)
grid = np.fromfile('record_Grid.dat', dtype = float)

rho_e = np.trapz(rho, grid)/grid[-1]
phi_anal = (rho_e/2/eps_0) * grid * (grid[-1] - grid)
plt.plot(grid, phi_anal)

