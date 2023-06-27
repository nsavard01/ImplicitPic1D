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
        self.numDiag = int(initialCond[11]) + 1
        ParticleProperties = pd.read_csv(self.path + 'ParticleProperties.dat', skiprows = 1, names = ['name', 'mass', 'q', 'w_p'], delim_whitespace = True)
        self.particles = {}
        partDiag = ['time', 'leftCurrLoss', 'rightCurrLoss', 'leftPowerLoss', 'rightPowerLoss']
        for i in range(len(ParticleProperties)):
            name = ParticleProperties.iloc[i]['name']
            self.particles[name] = {}
            self.particles[name]['mass'] = ParticleProperties.iloc[i]['mass']
            self.particles[name]['q'] = ParticleProperties.iloc[i]['q']
            self.particles[name]['w_p'] = ParticleProperties.iloc[i]['w_p']
            self.particles[name]['diag'] = pd.read_csv(self.path + 'ParticleDiagnostic_' + name + '.dat', skiprows = 1, delim_whitespace=True, names = partDiag)
            
        self.grid = np.fromfile(self.path + 'domainGrid.dat', dtype = 'float', offset = 4)
        self.dx_dl = np.fromfile(self.path + 'domainDxDl.dat', dtype = 'float', offset = 4)
        endDiag = np.loadtxt(self.path + 'SimulationFinalData.dat', skiprows=1)
        self.totTime = endDiag[0]
        self.totPotTime = endDiag[1]
        self.totCollTime = endDiag[2]
        self.totTimeSteps = int(endDiag[3])
        self.totSplitSteps = int(endDiag[4])
        diagList = ['time(s)', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)', 'TotalEnergy(J/m^2)', 'chargeError', 'energyError', 'numPicardIter']
        self.globDiag = pd.read_csv(self.path + 'GlobalDiagnosticData.dat', skiprows = 1, delim_whitespace=True, names = diagList)
        self.boolAverageFile = os.path.isfile(self.path + 'GlobalDiagnosticDataAveraged.dat')
        if (self.boolAverageFile):
            diagAverageList = ['steps', 'time(s)', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)']
            self.aveGlobDiag = pd.read_csv(self.path + 'GlobalDiagnosticDataAveraged.dat', skiprows = 1, delim_whitespace=True, names = diagAverageList)
        else:
            self.aveGlobDiag = None
            print("No averaging done for this simulation!")
        solver = np.loadtxt(self.path + 'SolverState.dat', skiprows = 1)
        if (int(solver[0]) == 0):
            self.solverType = 'AAc'
        else:
            self.solverType = 'JFNK'
        self.eps_r = solver[1]
        self.m_And = int(solver[2])
        self.beta_k = solver[3]
        self.maxIter = int(solver[4])
            
            
    def getPhaseSpace(self, name, i):
        if (name not in self.particles,keys()):
            raise Warning('No such particle', name, 'in simulation!')
        if (i <= self.numDiag - 1):
            phaseSpace = np.fromfile(self.path + 'PhaseSpace/phaseSpace_'+ name + '_' + str(i) + '.dat', dtype = 'float', offset = 4)
            phaseSpace = phaseSpace.reshape((int(phaseSpace.size/4), 4))
            d = phaseSpace[:,0] - phaseSpace[:,0].astype(int)
            phaseSpace[:,0] = self.grid[phaseSpace[:,0].astype(int)-1] + d * (self.grid[phaseSpace[:,0].astype(int)] - self.grid[phaseSpace[:,0].astype(int)-1])
            return phaseSpace
        else:
            raise Warning("No such i diagnostic!")
    
    def getPhi(self, i):
        if (i <= self.numDiag - 1):
            phi = np.fromfile(self.path + 'Phi/phi_' + str(i) + '.dat', dtype = 'float', offset = 4)
            return phi
        else:
            raise Warning("No such i diagnostic!")
            
    def getDensity(self, name, i):
        if (name not in self.particles.keys()):
            raise Warning('No such particle', name, 'in simulation!')
        if (i <= self.numDiag - 1):
            n = np.fromfile(self.path + 'Density/density_' + name + '_' + str(i) + '.dat', dtype = 'float', offset = 4)
            return n
        else:
            raise Warning("No such i diagnostic!")
            
    def getETemp(self, i):
        if (i <= self.numDiag - 1):
            temp = np.fromfile(self.path + 'ElectronTemperature/eTemp_' + str(i) + '.dat', dtype = 'float', offset = 4)
            return temp
        else:
            raise Warning("No such i diagnostic!")
            
    def getAvePhi(self):
        if (self.boolAverageFile):
            phi = np.fromfile(self.path + 'Phi/phi_Average.dat', dtype = 'float', offset = 4)
        else:
            raise Warning('No averaging done!')
        return phi
            
    def getAveDensity(self, name):
        if (name not in self.particles.keys()):
            raise Warning('No such particle', name, 'in simulation!')
        if (self.boolAverageFile):
            n = np.fromfile(self.path + 'Density/density_' + name + '_Average.dat', dtype = 'float', offset = 4)
        else:
            raise Warning('No averaging done!')
        return n
    
    def getAveEVDF(self):
        if (self.boolAverageFile):
            VTot = np.fromfile(self.path + 'ElectronTemperature/EVDF_average.dat', dtype = 'float', offset = 4)
            VHist = VTot[0:-1]
            VMax = VTot[-1]
            size = VHist.size/2
            Vedge = np.arange(-size,size+1) * VMax/size
            Vbins = (Vedge[0:-1] + Vedge[1:])/2
        else:
            raise Warning('No averaging done!')
        return Vbins, VHist
            
        

    
    
        
            