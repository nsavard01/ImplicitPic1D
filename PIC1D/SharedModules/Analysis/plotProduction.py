# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:39:21 2023

@author: Nicolas
"""
from import_libraries_constants import *
from dataSet import *
from dataSetExplicit import *


# ---------------------------- Averages -------------------------------------
def plotAveDensity(dataSet, name = "", label = "", marker = 'o', linestyle = '--'):
    if name == "":
        colors = ['b', 'r', 'g', 'k', 'c', 'm', 'yAve']
        for i,name in enumerate(dataSet.particles.keys()):
            n = dataSet.getAveDensity(name)
            plt.plot(dataSet.grid, n,  linestyle = linestyle, marker = marker, markersize = 2,color = colors[i], label = r'$n_{' + name +  '}$')
        plt.xlabel('Distance (m)')
        plt.ylabel('Particle Density (1/m^3)')
        plt.xlim([dataSet.x_min, dataSet.x_max])
        plt.legend(loc = 'best')
    else:
        if name not in dataSet.particles.keys():
            raise Warning("For average density, particle", name, "does not exist in the dataSet!")
        else:
            n = dataSet.getAveDensity(name)
            plt.plot(dataSet.grid, n,  linestyle = linestyle, marker = marker, markersize = 2, label = label)
            plt.xlabel('Distance (m)')
            plt.ylabel(name + ' Density (1/m^3)')
            plt.xlim([dataSet.x_min, dataSet.x_max])
 
def plotAvePhi(dataSet, label = '', marker = 'o', linestyle = '--'):

    phi = dataSet.getAvePhi()
    plt.plot(dataSet.grid, phi, marker = marker, linestyle = linestyle, markersize = 2, label = label)
    plt.xlabel('Distance (m)')
    plt.ylabel('Potential (V)')
    plt.xlim([dataSet.x_min, dataSet.x_max])
    plt.show()
    
def maxwellEDVF(x, T):
    return np.sqrt(m_e/2/np.pi / e/ T) * np.exp(- m_e * x**2 / 2 / e/ T)  
    
def plotAveVDF(dataSet, name, CurveFit = False, marker = '', linestyle = '--', label = ''):
    VHist, Vedge = dataSet.getAveVDF(name)
    dv = Vedge[1] - Vedge[0]
    Norm = np.sum(VHist[1:-1] * dv) + 0.5 * (VHist[0] + VHist[-1]) * dv
    VHist = VHist/Norm
    plt.plot(Vedge, VHist, marker = marker, linestyle = linestyle, linewidth = 2, label=label)
    plt.ylabel('VDF (s/m)')
    plt.xlabel('Velocity (m/s)')
    if CurveFit:
        popt, pcov = opt.curve_fit(maxwellEDVF, Vedge, VHist, p0 = [dataSet.T_e])
        plt.plot(Vedge, maxwellEDVF(Vedge, popt[0]), label = label + 'Fit T_e = ' + '{:2.2f}'.format(popt[0]))
        plt.legend(loc = 'best')

def plotAveEDF(dataSet, name, CurveFit = False, label = ''):
    EHist, Ebin = dataSet.getAveEDF(name)
    EHist = EHist/Ebin
    Norm = np.trapz(EHist, Ebin)
    EHist = EHist/Norm

    if CurveFit:
        f_log = np.log(EHist/np.sqrt(Ebin))
        fit = np.polyfit(Ebin, f_log, 1)
        plt.plot(Ebin, EHist/np.sqrt(Ebin), '--', label=label)
        dist = np.poly1d(fit)
        plt.plot(Ebin, np.exp(dist(Ebin)), label = label + 'Fit T_e = ' + '{:2.2f}'.format(-1/fit[0]))
        plt.ylabel(r'EPDF eV$^{-3/2}$')
        plt.yscale('log')
        plt.xlabel('Particle Energy (eV)')
        plt.legend(loc = 'best')
    else:
        plt.plot(Ebin, EHist, 'o-', label=label)
        plt.ylabel('EDF (1/eV)')
        plt.xlabel('Particle Energy (eV)')


def plotAveEPF(dataSet, name, lower_limit = 0.2, label = '', marker = 'o', linestyle = '--'):
    EHist, Ebin = dataSet.getAveEDF(name)
    EHist = EHist/Ebin
    Norm = np.trapz(EHist, Ebin)
    EHist = EHist/Norm
    lower_indx = np.where(Ebin < lower_limit)[0][-1]
    plt.plot(Ebin[lower_indx::], EHist[lower_indx::]/np.sqrt(Ebin[lower_indx::]), linestyle = linestyle, marker = marker, markersize = 2, label=label)
    plt.ylabel(r'EPDF eV$^{-3/2}$')
    plt.yscale('log')
    plt.xlabel('Particle Energy (eV)')

def getAveDensityFiles(dataName, name = 'e', diagNumber = 0):
    # Given name of file which will have _1, _2, etc appended
    if ('Explicit' in dataName):
        data = dataSetExplicit(dataName + '/')
    else:
        data = dataSet(dataName + '/')
    density = data.getDensity(name, diagNumber)
    fileList = glob.glob(dataName + '_[0-9]')
    number = len(fileList)
    for filename in fileList:
        if ('Explicit' in dataName):
            data = dataSetExplicit(filename + '/')
        else:
            data = dataSet(filename + '/')
        density = density + data.getDensity(name, diagNumber)
    print('Finished averaging for', dataName)
    print('Contained', number+1, 'files')
    return density/(number+1)

def plotAveDensityVsRef(refData, dataList, name, labelList = [''], marker = 'o', linestyle = '--'):
    #Plotting average density vs a reference density
    refDensity = refData.getAveDensity(name)
    refGrid = refData.grid
    for indx,data in enumerate(dataList):
        density = data.getAveDensity(name)
        interpDensity = np.interp(data.grid, refGrid, refDensity)
        relDensity = (density - interpDensity)
        if (labelList[0] == ''):
            labelStr = ''
        else:
            labelStr = labelList[indx]
        plt.plot(data.grid, relDensity, linestyle=linestyle, marker=marker, markersize=2,
             label= labelStr)

def plotAveDensityVsRef_percent(refData, dataList, name, labelList = [''], marker = 'o', linestyle = '--'):
    #Plotting average density vs a reference density
    refDensity = refData.getAveDensity(name)
    refGrid = refData.grid
    for indx,data in enumerate(dataList):
        density = data.getAveDensity(name)
        interpDensity = np.interp(data.grid, refGrid, refDensity)
        relDensity = (density - interpDensity)/interpDensity
        if (labelList[0] == ''):
            labelStr = ''
        else:
            labelStr = labelList[indx]
        plt.plot(data.grid, relDensity, linestyle=linestyle, marker=marker, markersize=2,
             label= labelStr)
    
#-------------------------- Time Dependent ---------------------
    
def update_plot_Phi(i, dataSet, ax, ylim):

    ax.clear()
    phi = dataSet.getPhi(i)
    ax.plot(dataSet.grid, phi, 'o-')
        
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Potential (V)')
    ax.set_xlim([dataSet.x_min, dataSet.x_max])
    if (ylim[0] != 0 and ylim[1] != 0):
        ax.set_ylim(ylim)

    
def phiAnimation(dataSet, boolMakeAnimation = False, savePath = "Figures/BoundPlasmaPhi.gif", pauseTime = 0.05, ylim = [0,0]):

    if boolMakeAnimation:
        numframes = dataSet.numDiag
        fig, ax = plt.subplots()
        phi = dataSet.getPhi(0)
        ax.plot(dataSet.grid, phi, 'o-')
            
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Potential (V)')
        ax.set_xlim([dataSet.x_min, dataSet.x_max])
        if (ylim[0] != 0 and ylim[1] != 0):
            ax.set_ylim(ylim)
        ani = animation.FuncAnimation(fig, update_plot_Phi, frames=range(numframes), interval = 100,fargs=(dataSet, ax, ylim))
        
        ani.save(savePath)
        plt.show()
        
        
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(dataSet.numDiag):
        
            plt.cla()
            
            phi = dataSet.getPhi(y)
            plt.plot(dataSet.grid, phi, 'o-')
            plt.xlabel('Distance (m)')
            plt.ylabel('Potential (V)')
            plt.xlim([dataSet.x_min, dataSet.x_max])
            if (ylim[0] != 0 and ylim[1] != 0):
                plt.ylim(ylim)
            plt.pause(pauseTime)

def update_plot_Density(i, dataSet, ax, nameList):

    ax.clear()
    for name in nameList:
        n = dataSet.getDensity(name, i)
        ax.plot(dataSet.grid, n, 'o-', label = name)
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Particle Density (1/m^3)')
    ax.set_xlim([dataSet.x_min, dataSet.x_max])
    plt.legend(loc = 'lower center')            
            
def densityAnimation(dataSet, nameList,boolMakeAnimation = False, savePath = "Figures/BoundPlasmaDensity.gif", pauseTime = 0.05):

    if boolMakeAnimation:
        numframes = dataSet.numDiag
        fig, ax = plt.subplots()
        for name in nameList:
            n = dataSet.getDensity(name, 0)
            ax.plot(dataSet.grid, n, 'o-', label = name)
            
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel(r'Density (m$^{-3}$)')
        ax.set_xlim([dataSet.x_min, dataSet.x_max])
        #ax.set_ylim([-20, 40])
        ani = animation.FuncAnimation(fig, update_plot_Density, frames=range(numframes), interval = 100,fargs=(dataSet, ax, nameList))
        
        ani.save(savePath)
        plt.show()
        
        
    else:
        plt.figure(figsize = (5,4), dpi = 80)
        for y in range(dataSet.numDiag):
        
            plt.cla()
            
            for name in nameList:
                n = dataSet.getDensity(name, y)
                plt.plot(dataSet.grid, n, 'o-', label = name)
            plt.xlabel('Distance (m)')
            plt.ylabel(r'Density (m$^{-3}$)')
            plt.xlim([dataSet.x_min, dataSet.x_max])
            plt.legend(loc='lower center')
            plt.pause(pauseTime)
            
def plotPhaseSpace(dataSet, name):
    phaseSpace = dataSet.getPhaseSpace(name)
    plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
# def update_plot_PhaseSpace(i, scat, dataSet, name):
#     phaseSpace = dataSet.getPhaseSpace(name, i)
#     scat.set_offsets(phaseSpace[:,0:2])
#     return scat,
#
# def PhaseSpaceAnimation(dataSet, name, boolMakeAnimation = False, savePath = 'PhaseSpaceAnimation.gif'):
#     if name not in dataSet.particles.keys():
#         raise Warning("Particle", name, "does not exist in the dataSet!")
#     if boolMakeAnimation:
#         numframes = dataSet.numDiag
#         phaseSpace = dataSet.getPhaseSpace(name, 0)
#         fig = plt.figure()
#         scat = plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
#         plt.xlabel('Distance (m)')
#         plt.ylabel('Speed (m/s)')
#         vmax = abs(np.max(phaseSpace[:,1]))*2.0
#         plt.axis([dataSet.grid[0], dataSet.grid[-1], -vmax, vmax])
#         ani = animation.FuncAnimation(fig, update_plot_PhaseSpace, frames=range(numframes), interval = 100,
#                                       fargs=(scat, dataSet, name))
#         ani.save(savePath)
#         plt.show()
#     else:
#         plt.figure(figsize = (5,4), dpi = 80)
#         for y in range(dataSet.numDiag):
#             plt.cla()
#             phaseSpace = dataSet.getPhaseSpace(name, y)
#             plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
#             plt.xlabel('Distance (m)')
#             plt.ylabel('Particle velocity (m/s)')
#             plt.xlim([0, dataSet.grid[-1]])
#             plt.pause(0.1)
            
    
 
#------------ Calculations ---------------------------

def get_thermalization_parameters(dataSet, reaction):
    nu = 0
    del_t = dataSet.delT
    for obj in dataSet.binaryColl[reaction]:
        nu += obj['aveDiag']['ratio']
        if (obj['type'] == 'Ionization'):
            nu_ion = obj['aveDiag']['ratio']/del_t
    nu = nu / del_t
    n_e = dataSet.getAveDensity('e').mean()

    EHist, Ebin = dataSet.getAveEDF('e')
    EHist = EHist/Ebin
    Norm = np.trapz(EHist, Ebin)
    EHist = EHist/Norm
    T_e = np.trapz(EHist*Ebin, Ebin) * 2/3
    deb = debye_length(T_e, n_e)
    L = dataSet.grid[-1] - dataSet.grid[0]
    N_p = dataSet.particles['e']['diag']['N_p'].values[-1]
    N_D = N_p/(L/deb)
    print('N_D is', N_D)
    plas_freq = plasmaFreq(n_e)
    denom = (N_D**(-2)) + 28 * (1/N_D) * (nu/plas_freq)
    tau_R = 34.4 / denom / plas_freq

    nu_coulomb = crossSectionCoulomb(T_e, n_e)
    nu_coulomb = nu_coulomb * n_e * np.sqrt(2 * e * T_e / np.pi / m_e)
    return tau_R, 1/nu, 1/nu_coulomb, 1/nu_ion


