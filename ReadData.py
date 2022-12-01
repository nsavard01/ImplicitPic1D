#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:42:42 2022

@author: nicolas
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

def maxBoltzDist(v, T, m): # T in eV
    return (m/2/np.pi/e/T)**(1/2) * np.exp(-m * v**2 /2 / e / T)

def maxBoltzDistE(E, T):
    return 2 * np.sqrt(E/np.pi) * (1/T)**(1.5) * np.exp(-E/T)

position = np.fromfile('src/record_particlePosition.dat', dtype=float)
velocity = np.fromfile('src/record_particleVelocity.dat', dtype=float)
velocity = velocity.reshape((3, int(velocity.size/3)))

plt.figure()
plt.hist(velocity[2, :], bins = 100, density = True)

T_e = 5
v_test = np.linspace(-np.sqrt(T_e * e/m_e) * 5, np.sqrt(T_e * e/m_e) * 5, 10000)
plt.plot(v_test, maxBoltzDist(v_test, T_e, m_e))

plt.figure()
E = np.sum(velocity**2, axis = 0) * m_e * 0.5 / e
E_test = np.linspace(0, T_e * 3 *25/2, 10000)
plt.hist(E, bins = 100, density = True)
plt.plot(E_test, maxBoltzDistE(E_test, T_e))
