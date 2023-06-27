# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:39:21 2023

@author: Nicolas
"""
from import_libraries_constants import *
from dataSet import *

test = dataSet('Y:/scratch/nsavard/ImplicitPic1D/ImplicitPic1D/NGP/Data/')


# def PhaseSpaceAnimation(dataSet, boolAnimation = False): 
#     if boolMakeAnimation:
#         numframes = numDiagnosticTimes+1
#         filename = '../Data/PhaseSpace/phaseSpace_'+ ParticleName + '_0.dat'
#         phaseSpace = extractPhaseSpace(filename, grid)
#         fig = plt.figure()
#         scat = plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
#         plt.xlabel('Distance (m)')
#         plt.ylabel('Speed (m/s)')
#         plt.axis([0, grid[-1], -6500, 6500])
#         ani = animation.FuncAnimation(fig, update_plot_PhaseSpace, frames=range(numframes), interval = 100,
#                                       fargs=(scat, ParticleName))
#         ani.save('PostProcessing/twoStreamAnimation.gif')
#         plt.show()
#     else:
#         plt.figure(figsize = (5,4), dpi = 80)
#         for y in range(numDiagnosticTimes+1):
#             plt.cla()
#             filename = '../Data/PhaseSpace/phaseSpace_'+ ParticleName + '_' + str(y) +'.dat'
#             phaseSpace = extractPhaseSpace(filename, grid)
#             plt.scatter(phaseSpace[:,0], phaseSpace[:,1])
#             plt.xlabel('Distance (m)')
#             plt.ylabel('Particle velocity (m/s)')
#             plt.xlim([0, grid[-1]])
#             plt.pause(0.1)