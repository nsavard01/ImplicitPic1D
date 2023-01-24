# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 16:50:43 2023

@author: Nicolas
"""
import numpy as np

data = np.fromfile('dataNormal.dat', dtype='float', offset = 4)

dataObject = np.fromfile('dataObject.dat', dtype='float', offset=4)