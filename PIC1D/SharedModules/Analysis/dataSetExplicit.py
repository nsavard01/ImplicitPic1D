from import_libraries_constants import *

import tkinter
from tkinter import filedialog
import re


class dataSetExplicit:
    
    def __init__(self, filename = ""):
        if filename == "":
            tkinter.Tk().withdraw()
            path = filedialog.askdirectory()
            self.path = path + '/'
        else:
            self.path = filename
        dateTimeBool = os.path.isfile(self.path + 'DateTime.dat')
        if (dateTimeBool):
            dateTime = np.loadtxt(self.path + 'DateTime.dat', skiprows=1)
            self.dateTime = str(dateTime[0])[0:4] + '-' + str(dateTime[0])[4:6] + '-' + str(dateTime[0])[6:8] + ' ' + \
                        str(int(str(dateTime[1])[0:2]) - 7) + ':' + str(dateTime[1])[2:4] + ':' + str(dateTime[1])[4:6]
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
        ParticleProperties = pd.read_csv(self.path + 'ParticleProperties.dat', skiprows = 1, names = ['name', 'mass', 'q', 'w_p', 'maxIdx'], sep='\s+')
        self.particles = {}
        partDiag = ['time', 'leftCurrLoss', 'rightCurrLoss', 'leftPowerLoss', 'rightPowerLoss', 'N_p', 'Temp']
        for i in range(len(ParticleProperties)):
            name = ParticleProperties.iloc[i]['name']
            self.particles[name] = {}
            self.particles[name]['mass'] = ParticleProperties.iloc[i]['mass']
            self.particles[name]['q'] = ParticleProperties.iloc[i]['q']
            self.particles[name]['w_p'] = ParticleProperties.iloc[i]['w_p']
            self.particles[name]['maxIdx'] = ParticleProperties.iloc[i]['maxIdx']
            self.particles[name]['diag'] = pd.read_csv(self.path + 'ParticleDiagnostic_' + name + '.dat', skiprows = 1, sep='\s+', names = partDiag)
            if (os.path.isfile(self.path + 'ParticleAveDiagnostic_' + name + '.dat')):
                avePartDiag = ['leftCurrLoss', 'rightCurrLoss', 'leftPowerLoss', 'rightPowerLoss']
                self.particles[name]['aveDiag'] = pd.read_csv(self.path + 'ParticleAveDiagnostic_' + name + '.dat', skiprows = 1, sep='\s+', names = avePartDiag)

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
        if (dateTimeBool):
            self.totMoveTime = endDiag[2]
            self.totCollTime = endDiag[3]
            self.totTimeSteps = int(endDiag[4])
        else:
            self.totCollTime = endDiag[2]
            self.totTimeSteps = int(endDiag[3])
        diagList = ['time(s)', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)', 'TotalMomentum(kg/m/s)', 'TotalEnergy(J/m^2)']
        self.globDiag = pd.read_csv(self.path + 'GlobalDiagnosticData.dat', skiprows = 1, sep='\s+', names = diagList)
        timeDataBool = os.path.isfile(self.path + 'SimulationTimeData.dat')
        if (timeDataBool):
            diagList = ['Time', 'PotTime', 'MoverTime', 'CollTime', 'Step']
            self.timeDiag = pd.read_csv(self.path + 'SimulationTimeData.dat', skiprows = 1, sep='\s+', names = diagList)
            if (self.globDiag.shape[0] > self.numDiag):
                self.numDiag = self.globDiag.shape[0]
        self.boolAverageFile = os.path.isfile(self.path + 'GlobalDiagnosticDataAveraged.dat')
        if (self.boolAverageFile):
            diagAverageList = ['steps', 'Ploss(W/m^2)', 'I_wall(A/m^2)', 'P_wall(W/m^2)']
            self.aveGlobDiag = pd.read_csv(self.path + 'GlobalDiagnosticDataAveraged.dat', skiprows = 1, sep='\s+', names = diagAverageList)
        else:
            self.aveGlobDiag = None
            print("No averaging done for this simulation!")
        if (os.path.isdir(self.path + 'BinaryCollisions')):
            self.binaryColl = {}
            pattern = '(.*)_on_(.*)'
            for filename in os.listdir(self.path + 'BinaryCollisions'):
                test = re.search(pattern, filename)
                primary = test.group(1)
                target = test.group(2)
                react = primary + '->' + target
                prop = np.loadtxt(self.path + 'BinaryCollisions/' + filename + '/CollisionProperties.dat', skiprows = 1)
                self.binaryColl[react] = list(range(prop.shape[0]))
                collDiag = ['ratio', 'aveEnergyLoss', 'aveIncEnergy', 'P_loss(W/m^2)', 'aveFreq(Hz/m^2)']
                for i in range(prop.shape[0]):
                    self.binaryColl[react][i] = {}
                    type = int(prop[i, 1])
                    if (type == 1):
                        self.binaryColl[react][i]['type'] = 'Elastic'
                    elif (type == 2):
                        self.binaryColl[react][i]['type'] = 'Ionization'
                    elif (type == 3):
                        self.binaryColl[react][i]['type'] = 'Excitation'
                    elif (type == 4):
                        self.binaryColl[react][i]['type'] = 'ChargeExchange'
                    self.binaryColl[react][i]['E_thres'] = prop[i,2]
                    self.binaryColl[react][i]['maxSigma'] = prop[i, 3]
                    self.binaryColl[react][i]['EatMaxSigma'] = prop[i, 4]
                    self.binaryColl[react][i]['diag'] = pd.read_csv(self.path + 'BinaryCollisions/' + filename + '/CollisionDiag_' + str(i+1) + '.dat', skiprows = 1, sep='\s+', names = collDiag)
                if (self.boolAverageFile):
                    prop = np.loadtxt(self.path + 'BinaryCollisions/' + filename + '/AveCollisionDiag.dat', skiprows = 1)
                    for i,coll in enumerate(self.binaryColl[react]):
                        coll['aveDiag'] = {}
                        coll['aveDiag']['ratio'] = prop[i,1]
                        coll['aveDiag']['aveEnergyLoss'] = prop[i,2]
                        coll['aveDiag']['aveIncEnergy'] = prop[i, 3]
                        coll['aveDiag']['P_loss(W/m^2)'] = prop[i, 4]
                        coll['aveDiag']['aveFreq(Hz/m^2)'] = prop[i, 5]







    def getPhaseSpace(self, name):
        if (name not in self.particles.keys()):
            raise Warning('No such particle', name, 'in simulation!')
        phaseSpace = []
        for i in range(self.numThreads):
            temp = np.fromfile(self.path + 'PhaseSpace/phaseSpace_' + name + '_thread' + str(i+1) + '.dat', dtype='float',
                               offset=0)
            temp = temp.reshape((int(temp.size / 4), 4))
            d = temp[:, 0] - temp[:, 0].astype(int)
            temp[:, 0] = self.grid[temp[:, 0].astype(int) - 1] + d * (
                        self.grid[temp[:, 0].astype(int)] - self.grid[temp[:, 0].astype(int) - 1])
            phaseSpace.append(temp)

        phaseSpace = np.concatenate(phaseSpace)
        return phaseSpace

    
    def getPhi(self, i):
        if (i <= self.numDiag - 1):
            phi = np.fromfile(self.path + 'Phi/phi_' + str(i) + '.dat', dtype = 'float', offset = 4)
            return phi
        else:
            raise Warning("No such i diagnostic!")

    def getEField(self, i):
        if (i <= self.numDiag - 1):
            phi = np.fromfile(self.path + 'Phi/phi_' + str(i) + '.dat', dtype='float', offset=4)
            EField = np.zeros(self.Nx)
            EField[1:-1] = (phi[0:-2] - phi[2::])/(self.grid[2::] - self.grid[0:-2])
            EField[0] = (phi[-1] - phi[1])/(2 * self.grid[1])
            EField[-1] = EField[0]
            return (self.grid, EField)
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

    def getAveTemp(self, name):
        if (self.boolAverageFile):
            EHist, Ebin = self.getAveEDF(name)
            Norm = np.sum(EHist*Ebin)/ np.sum(EHist)
            T = Norm * 2 / 3
        else:
            raise Warning('No averaging done!')
        return T
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

    def getAveVDF(self, name):
        if (self.boolAverageFile):
            VTot = np.fromfile(self.path + 'Temperature/Temp_' + str(name) + '_average.dat', dtype='float', offset=4)
            VHist = VTot[0:-1]
            VMax = VTot[-1]

            Vedge = np.linspace(-VMax, VMax, VHist.size)
        else:
            raise Warning('No averaging done!')
        return VHist, Vedge

    def getAveEDF(self, name):
        if (self.boolAverageFile):
            ETot = np.fromfile(self.path + 'Temperature/TempEnergy_' + str(name) + '_average.dat', dtype='float', offset=4)
            EHist = ETot[0:int(ETot.size/2)]
            Ebin = ETot[int(ETot.size/2)::]
        else:
            raise Warning('No averaging done!')
        return EHist, Ebin
    
    def getBoundaryConditions(self):
        cond = np.fromfile(self.path + 'domainBoundaryConditions.dat', dtype = np.int32, offset = 4)[0:-1]
        return cond