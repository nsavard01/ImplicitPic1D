# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:20:31 2023

@author: Nicolas
"""

from import_libraries_constants import *

import tkinter
from tkinter import filedialog



class dataSet:
    
    def __init__(self, filename = ""):
        if filename == "":
            tkinter.Tk().withdraw()
            path = filedialog.askdirectory()
            self.path = path + '/'
        else:
            self.path = filename
        
        initialCond = np.loadtxt(self.path + 'InitialConditions.dat', skiprows = 1)
        if int(initialCond[0])== 0:
            self.scheme = 'NGP'
        elif int(initialCond[0])== 1:
            self.scheme = 'CIC'
        else:
            self.scheme = 'newCIC'
        self.Nx = int(initialCond[1])
        self.T_e = initialCond[2]
        self.T_i = initialCond[3]
        self.n_ave = initialCond[4]
        self.simTimeTotal = initialCond[5]
        self.delT = initialCond[6]
        self.fracTime = initialCond[7]
        self.numDiag = int(initialCond[11])
        ParticleProperties = pd.read_csv(self.path + 'ParticleProperties.dat', skiprows = 1, names = ['name', 'mass', 'q', 'w_p'], delim_whitespace = True)
        self.particles = {}
        for i in range(len(ParticleProperties)):
            self.particles[ParticleProperties.iloc[i]['name']] = {}
            self.particles[ParticleProperties.iloc[i]['name']]['mass'] = ParticleProperties.iloc[i]['mass']
            self.particles[ParticleProperties.iloc[i]['name']]['q'] = ParticleProperties.iloc[i]['q']
            self.particles[ParticleProperties.iloc[i]['name']]['w_p'] = ParticleProperties.iloc[i]['w_p']
            
        self.grid = np.fromfile(self.path + 'domainGrid.dat', dtype = 'float', offset = 4)
        self.dx_dl = np.fromfile(self.path + 'domainDxDl.dat', dtype = 'float', offset = 4)
        endDiag = np.loadtxt(self.path + 'SimulationFinalData.dat', skiprows=1)
        self.totPotTime = endDiag[0]
        self.totCollTime = endDiag[1]
        self.totTimeSteps = int(endDiag[2])
        self.totSplitSteps = int(endDiag[3])
        
    def getAveDiag(self):
        boolAverageFile = os.path.isfile(self.path + 'GlobalDiagnosticDataAveraged.dat')
        if (boolAverageFile):
            diagAverageList = ['steps', 'time(s)', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)']
            GlobalDiagnosticsAveraged = pd.read_csv(self.path + 'GlobalDiagnosticDataAveraged.dat', skiprows = 1, delim_whitespace=True, names = diagAverageList)
        else:
            GlobalDiagnosticsAveraged = 0
            print("No averaging done for this simulation!")
        return GlobalDiagnosticsAveraged

    
    
        
            