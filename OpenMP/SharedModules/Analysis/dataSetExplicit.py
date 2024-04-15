from import_libraries_constants import *

import tkinter
from tkinter import filedialog



class dataSetExplicit:
    
    def __init__(self, filename = ""):
        if filename == "":
            tkinter.Tk().withdraw()
            path = filedialog.askdirectory()
            self.path = path + '/'
        else:
            self.path = filename
        
        initialCond = np.loadtxt(self.path + 'InitialConditions.dat', skiprows = 1)
        self.scheme = 'Explicit'
        self.Nx = int(initialCond[0])
        self.simTimeTotal = initialCond[1]
        self.delT = initialCond[2]
        self.fracTime = initialCond[3]
        self.delX = initialCond[4]
        self.n_ave = initialCond[5]
        self.T_e = initialCond[6]
        self.T_i = initialCond[7]
        self.numDiag = int(initialCond[8]) + 1
        self.numChargedParticles = int(initialCond[9])
        self.numThreads = int(initialCond[10])
        self.RF_rad_frequency = initialCond[11]
        self.RF_half_amplitude = initialCond[12]
        ParticleProperties = pd.read_csv(self.path + 'ParticleProperties.dat', skiprows = 1, names = ['name', 'mass', 'q', 'w_p', 'maxIdx'], delim_whitespace = True)
        self.particles = {}
        partDiag = ['time', 'leftCurrLoss', 'rightCurrLoss', 'leftPowerLoss', 'rightPowerLoss', 'N_p', 'Temp']
        for i in range(len(ParticleProperties)):
            name = ParticleProperties.iloc[i]['name']
            self.particles[name] = {}
            self.particles[name]['mass'] = ParticleProperties.iloc[i]['mass']
            self.particles[name]['q'] = ParticleProperties.iloc[i]['q']
            self.particles[name]['w_p'] = ParticleProperties.iloc[i]['w_p']
            self.particles[name]['maxIdx'] = ParticleProperties.iloc[i]['maxIdx']
            self.particles[name]['diag'] = pd.read_csv(self.path + 'ParticleDiagnostic_' + name + '.dat', skiprows = 1, delim_whitespace=True, names = partDiag)

        boundConditions = np.fromfile(self.path + 'domainBoundaryConditions.dat', dtype = np.int32, offset = 4)[0:-1]
        if (boundConditions[0] == 1):
            self.leftBoundary = 'Dirichlet'
        elif (boundConditions[0] == 2):
            self.leftBoundary = 'Neumann'
        elif (boundConditions[0] == 3):
            self.leftBoundary = 'Periodic'
        elif (boundConditions[0] == 4):
            self.leftBoundary = 'RF-Dirichlet'
        else:
            raise Warning("Left boundary not defined!")

        if (boundConditions[-1] == 1):
            self.rightBoundary = 'Dirichlet'
        elif (boundConditions[-1] == 2):
            self.rightBoundary = 'Neumann'
        elif (boundConditions[-1] == 3):
            self.rightBoundary = 'Periodic'
        elif (boundConditions[-1] == 4):
            self.rightBoundary = 'RF-Dirichlet'
        else:
            raise Warning("Right boundary not defined!")
        self.grid = np.fromfile(self.path + 'domainGrid.dat', offset = 4)
        self.x_min = self.grid[0]
        self.x_max = self.grid[-1]
        endDiag = np.loadtxt(self.path + 'SimulationFinalData.dat', skiprows=1)
        self.totTime = endDiag[0]
        self.totPotTime = endDiag[1]
        self.totCollTime = endDiag[2]
        self.totTimeSteps = int(endDiag[3])
        diagList = ['time(s)', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)']
        self.globDiag = pd.read_csv(self.path + 'GlobalDiagnosticData.dat', skiprows = 1, delim_whitespace=True, names = diagList)
        self.boolAverageFile = os.path.isfile(self.path + 'GlobalDiagnosticDataAveraged.dat')
        if (self.boolAverageFile):
            diagAverageList = ['steps', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)']
            self.aveGlobDiag = pd.read_csv(self.path + 'GlobalDiagnosticDataAveraged.dat', skiprows = 1, delim_whitespace=True, names = diagAverageList)
        else:
            self.aveGlobDiag = None
            print("No averaging done for this simulation!")
            
            
    def getPhaseSpace(self, name):
        if (name not in self.particles.keys()):
            raise Warning('No such particle', name, 'in simulation!')

        phaseSpace = np.fromfile(self.path + 'PhaseSpace/phaseSpace_'+ name + '.dat', dtype = 'float', offset = 4)
        phaseSpace = phaseSpace.reshape((int(phaseSpace.size/4), 4))
        d = phaseSpace[:,0] - phaseSpace[:,0].astype(int)
        phaseSpace[:,0] = self.grid[phaseSpace[:,0].astype(int)-1] + d * (self.grid[phaseSpace[:,0].astype(int)] - self.grid[phaseSpace[:,0].astype(int)-1])
        return phaseSpace

    
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
            
    def getTemp(self, name, i):
        if (name not in self.particles.keys()):
            raise Warning('No such particle', name, 'in simulation!')
        if (i <= self.numDiag - 1):
            temp = np.fromfile(self.path + 'Temperature/Temp_' + name + '_' + str(i) + '.dat', dtype = 'float', offset = 4)
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
    
    def getBoundaryConditions(self):
        cond = np.fromfile(self.path + 'domainBoundaryConditions.dat', dtype = np.int32, offset = 4)[0:-1]
        return cond