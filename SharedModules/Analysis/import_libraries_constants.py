# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:18:53 2023

@author: Nicolas
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constants
import pandas as pd
import glob, os
import matplotlib.animation as animation

eps_0 = constants.epsilon_0
c = constants.c
m_e = constants.m_e
m_p = constants.m_p
mu_0 = constants.mu_0
k_boltz = constants.k
e = constants.e