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

phi = np.fromfile('test.dat', dtype = 'float', offset = 4)