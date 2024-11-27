# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:42 2023

@author: Nicolas
"""

from plotProduction import *
from generateBoundExample import *
print('Loaded modules')
    
#%%

print('Loading data')

def readBenchmarkResults():
    firstIndx = 0
    data = np.loadtxt('TurnerBenchmark/turner_benchmark_results.dat')
    case1 = data[0:129,:]
    firstIndx = firstIndx + 129
    case2 = data[firstIndx:firstIndx+257,:]
    firstIndx = firstIndx + 257
    case3 = data[firstIndx:firstIndx+513, :]
    firstIndx = firstIndx+513
    case4 = data[firstIndx:firstIndx+513]
    return case1,case2,case3,case4

def readBenchmarkResultsRefined():
    data = np.loadtxt('TurnerBenchmark/turner_benchmark_refined_results.dat')
    indx = np.where(data[:,0] == 0)[0]
    case1 = data[indx[0]:indx[1], :]
    case2 = data[indx[1]:indx[2], :]
    case3 = data[indx[2]:indx[3], :]
    case4 = data[indx[3]:,:]
    return case1, case2, case3, case4

case1, case2, case3, case4 = readBenchmarkResults()
case1_refined, case2_refined, case3_refined, case4_refined = readBenchmarkResultsRefined()

Exp_RFBenchmark_Denis = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis/')
Exp_RFBenchmark_Denis_moreCells = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells/')
Exp_RFBenchmark_Denis_moreCells_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_doublePart/')
Exp_RFBenchmark_Denis_moreCells_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_quadPart/')
Exp_RFBenchmark_Denis_moreCells_40xPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_40xPart/')

Exp_RFBenchmark_Denis_moreCells_lowerDen = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_doublePart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_doublePart/')
Exp_RFBenchmark_Denis_moreCells_lowerDen_quadPart = dataSetExplicit('Y:/ImplicitPic1D/ExplicitData/Exp_RFBenchmark_Denis_moreCells_lowerDen_quadPart/')

