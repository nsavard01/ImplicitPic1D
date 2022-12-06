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

l_p = np.fromfile('record_particlePosition.dat', dtype = float)
v_p = np.fromfile('record_particleVelocity.dat', dtype = float)
v_p = v_p.reshape((3, int(v_p.size/3)))
rho = np.fromfile('record_Rho.dat', dtype = float)
grid = np.fromfile('record_Grid.dat', dtype = float)
