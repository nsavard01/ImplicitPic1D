import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import scipy.constants as constants
import scipy.optimize as opt
import pandas as pd
import scipy
import glob, os
import matplotlib.animation as animation
import math

phi = np.fromfile('phi_0.dat', dtype = 'float', offset = 4)
grid = np.fromfile('domainGrid.dat', dtype='float', offset = 4)
boundary_conditions = np.fromfile('domainBoundaryConditions.dat', dtype='int', offset = 4)
density = np.fromfile('density_He+_0.dat', dtype = 'float', offset = 4)
temp = np.fromfile('Temp_e_0.dat', dtype = 'float', offset = 4)