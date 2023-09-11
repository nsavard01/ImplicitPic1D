# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:07:07 2023

@author: Nicolas
"""

from import_libraries_constants import *

def otherFunc(phi, tau):
    coeff_1 = 2.0 * np.exp(-tau *  phi**2)/np.sqrt(np.pi * tau)
    return coeff_1 * scipy.special.dawsn(phi) + scipy.special.erf(np.sqrt(tau) * phi) - 1

def funcPhi_s(phi, B, s):
    return (2/np.pi/B) * scipy.special.dawsn(phi) - s

def funcNonLinear(eta, A, B, tau, I_G, t, dah):
    res = np.zeros(eta.size)
    dt = t[1] - t[0]
    laplace = eta[0:-2] - 2*eta[1:-1] + eta[2:]
    firstDer = eta[2:] - eta[0:-2]
    bracket = laplace * (1 - t[1:-1])**(2-2*dah) / dah**2 / dt**2 + firstDer*(1-t[1:-1])**(1-2*dah) * (dah-1) / dah**2 / dt
    res[1:-1] = A * bracket * np.exp(-eta[0]) + np.exp(-eta[1:-1]) - B * np.exp(eta[1:-1] * tau) * I_G[1:-1]
    res[0] = np.exp(-eta[0]) - B * np.exp(eta[0] * tau) * I_G[0] #A * (-2 * eta[0] + 2 * eta[1]) * np.exp(-eta[0]) +
    res[-1] = eta[-1] - 0
    return res

def funcNonLinearTrunc(eta, A, B, tau, I_G, truncFact):
    res = np.zeros(eta.size)
    laplace = eta[0:-2] - 2*eta[1:-1] + eta[2:]
    firstDer = eta[2:] - eta[0:-2]
    bracket = laplace * (1 - t[1:-1])**(2-2*dah) / dah**2 / dt**2 + firstDer*(1-t[1:-1])**(1-2*dah) * (dah-1) / dah**2 / dt
    res[1:-1] = A * bracket * np.exp(-eta[0]) + np.exp(-eta[1:-1]) * truncFact[1:-1] - B * np.exp(eta[1:-1] * tau) * I_G[1:-1]
    res[0] = np.exp(-eta[0]) - B * np.exp(eta[0] * tau) * I_G[0] #A * (-2 * eta[0] + 2 * eta[1]) * np.exp(-eta[0]) +
    res[-1] = eta[-1] - 0
    return res

def getBoundPlasmaSolutions(L, numNodes, n_ave, T_e, T_i, M):
    fracDebye = 0.1
    lambDebye = debye_length(T_e, n_ave)
    t = np.linspace(0, 1, numNodes)
    dah = np.log(fracDebye * lambDebye/L)/np.log(1/(numNodes-1))
    s = 1 - (1-t)**dah
    dt = t[1] - t[0]
    phi = np.zeros(numNodes)
    tau = T_e/T_i
    M_m = M/m_e
    if (tau < 36):
        phi_1 = opt.fsolve(otherFunc, 0.5, args = (tau,))**2
    else:
        phi_1 = 0.854

    coeff_1 = np.pi * np.sqrt(M_m / 4 / np.pi) / (1 + T_i/T_e) / 2 / scipy.special.dawsn(np.sqrt(phi_1))
    phi_w = -np.log(coeff_1)
    B = 0.5 * np.sqrt(M_m/np.pi) * (tau / (1 + tau)) * np.exp(phi_w)

    phi[0] = - phi_w * T_e
    phi[-1] = 0.0
    for i in range(1, numNodes-1):
        phi[i] = phi[0] - (opt.fsolve(funcPhi_s, 0.2, args = (B, s[i]))**2) * T_e
        

    eta = (- phi)/T_e
    eta[-1] = 0.0
    n_0 = n_ave / np.exp(eta[0])/ np.trapz(np.exp(-eta) * dah * (1-t)**(dah-1), dx = dt)
    B = 0.5 * np.sqrt(M_m * tau)
    A = eps_0 * T_e / n_0 / e / L**2 

    I_G = np.zeros(numNodes)

    for i in range(numNodes):
        I_int = np.where(eta[i] < eta, (1 - t)**(dah-1) * np.exp(-eta * tau), (1 - t)**(dah-1) * np.exp(-eta*tau) * scipy.special.erfc(np.sqrt((eta[i] - eta)*tau)))
        I_G[i] = np.trapz(I_int, dx = dt)
    I_G = I_G * dah

    n_e = n_0 * np.exp(eta[0] - eta)
    n_i = n_0 * np.exp(eta*tau + eta[0]) * B * I_G

    initialR = np.sum(funcNonLinear(eta, A, B, tau, I_G, t, dah)**2)/numNodes
    eta = opt.fsolve(funcNonLinear, eta, args = (A, B, tau, I_G, t, dah))
    for i in range(1000):
        n_0 = n_ave / np.exp(eta[0])/ np.trapz(np.exp(-eta) * dah * (1-t)**(dah-1), dx = dt)
        A = eps_0 * T_e / n_0 / e / L**2
        for k in range(numNodes):
            I_int = np.where(eta[k] < eta, (1 - t)**(dah-1) * np.exp(-eta * tau), (1 - t)**(dah-1) * np.exp(-eta*tau) * scipy.special.erfc(np.sqrt((eta[k] - eta)*tau)))
            I_G[k] = np.trapz(I_int, dx = dt)
        I_G = I_G * dah
        pastEta = eta
        #currRes = np.sum(funcNonLinear(eta, A, B, tau, I_G, t, dah)**2)/numNodes
        eta = opt.fsolve(funcNonLinear, pastEta, args = (A, B, tau, I_G, t, dah))
        currRes = np.sqrt(np.sum((eta - pastEta)**2) / numNodes)
        if (currRes < 1e-8):
            print('Took', i, 'iterations for convergence!')
            break
        if (i == 999):
            raise Warning('Iteration has not converged for bounded plasma model')
            
    
    n_e = n_0 * np.exp(eta[0] - eta)
    n_i = n_0 * np.exp(eta*tau + eta[0]) * B * I_G
    phi = -eta * T_e
    sol = np.zeros((4, numNodes))
    sol[0,:] = L-s*L
    sol[1,:] = phi
    sol[2, :] = n_e
    sol[3, :] = n_i
    return np.flip(sol, axis = 1)

def compareRefToDatasRes(dataList, labelList, ref, saveFile):
    if (len(dataList) != len(labelList)):
        raise Warning("List of data and labels not the same size!")
        
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    for name in dataList[0].particles.keys():
        if name != 'e':
            ion = name
            break
    M = dataList[0].particles[ion]['mass']
    T_e = dataList[0].T_e
    T_i = dataList[0].T_i
    n_ave = dataList[0].n_ave
    for i in range(len(dataList)):
        print('Statistics for label', labelList[i])
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        n_e_ref = ref.getAveDensity('e')
        n_e_ref = np.interp(data.grid, ref.grid, n_e_ref)
        ax1.plot(data.grid, 100 * (n_e-n_e_ref)/n_e_ref, linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        print(r'$n_e$ is ', 100 * np.sum((n_e[1:-1]-n_e_ref[1:-1])/n_e_ref[1:-1])/(dataList[i].Nx - 2))
        
        n_i = data.getAveDensity(ion)
        n_i_ref = ref.getAveDensity('H+')
        n_i_ref = np.interp(data.grid, ref.grid, n_i_ref)
        ax2.plot(data.grid, 100 * (n_i-n_i_ref)/n_i_ref, linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        print(r'$n_i$ is ', 100 * np.sum((n_i[1:-1]-n_i_ref[1:-1])/n_i_ref[1:-1])/(dataList[i].Nx - 2))
        
        phi = data.getAvePhi()
        phi_ref = ref.getAvePhi()
        phi_ref = np.interp(data.grid, ref.grid, phi_ref)
        ax3.plot(data.grid, 100 * (phi-phi_ref)/(phi_ref + 1e-15), linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        print(r'$\phi$ is ', 100 * np.sum((phi[1:-1]-phi_ref[1:-1])/phi_ref[1:-1])/(dataList[i].Nx - 2))
        
        ax4.bar(labelList[i], data.totPotTime)
        ax4.text(i, data.totPotTime, data.totPotTime, ha = 'center')

    ax1.legend(loc = 'upper right', fontsize = 12)
    ax2.legend(loc = 'lower right', fontsize = 12)
    ax3.legend(loc = 'lower right', fontsize = 12)
    ax1.set_xlim(ref.grid[0], ref.grid[-1])
    ax1.set_xlabel('Distance (m)', fontsize = 14)
    ax1.set_ylabel(r'100 $\times$ $(n_e - n_{e,ref})/n_{e,ref}$', fontsize = 14)
    ax1.tick_params(labelsize = 14)
    
    ax2.set_xlim(ref.grid[0], ref.grid[-1])
    ax2.set_xlabel('Distance (m)', fontsize = 14)
    ax2.set_ylabel(r'100 $\times$ $(n_i - n_{i,ref})/n_{i,ref}$', fontsize = 14)
    ax2.tick_params(labelsize = 14)
    
    ax3.set_xlim(ref.grid[0], ref.grid[-1])
    ax3.set_xlabel('Distance (m)', fontsize = 14)
    ax3.set_ylabel(r'100 $\times$ $(\phi - \phi_{ref})/\phi_{ref}$', fontsize = 14)
    ax3.tick_params(labelsize = 14)
    
    ax4.set_ylabel(r'Total Wall Time for Solver (s)', fontsize = 14)
    
    fig1.savefig(saveFile + '_eDensity.pdf', bbox_inches = 'tight')
    fig2.savefig(saveFile + '_ionDensity.pdf', bbox_inches = 'tight')
    fig3.savefig(saveFile + '_voltage.pdf', bbox_inches = 'tight')
    fig4.savefig(saveFile + '_time.pdf', bbox_inches = 'tight')
    
def compareRefToDatasAbsRes(dataList, labelList, ref):
    if (len(dataList) != len(labelList)):
        raise Warning("List of data and labels not the same size!")
        
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    for name in dataList[0].particles.keys():
        if name != 'e':
            ion = name
            break
    M = dataList[0].particles[ion]['mass']
    T_e = dataList[0].T_e
    T_i = dataList[0].T_i
    n_ave = dataList[0].n_ave
    for i in range(len(dataList)):
        print('Statistics for label', labelList[i])
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        n_e_ref = ref.getAveDensity('e')
        n_e_ref = np.interp(data.grid, ref.grid, n_e_ref)
        ax1.plot(data.grid, 100 * abs((n_e-n_e_ref)/n_e_ref), linestyle = '--', marker = '.',label = labelList[i])
        print(r'$n_e$ is ', 100 * np.sum((n_e[1:-1]-n_e_ref[1:-1])/n_e_ref[1:-1])/(dataList[i].Nx - 2))
        
        n_i = data.getAveDensity(ion)
        n_i_ref = ref.getAveDensity('H+')
        n_i_ref = np.interp(data.grid, ref.grid, n_i_ref)
        ax2.plot(data.grid, 100 * abs((n_i-n_i_ref)/n_i_ref), linestyle = '--', marker = '.',label = labelList[i])
        print(r'$n_i$ is ', 100 * np.sum((n_i[1:-1]-n_i_ref[1:-1])/n_i_ref[1:-1])/(dataList[i].Nx - 2))
        
        phi = data.getAvePhi()
        phi_ref = ref.getAvePhi()
        phi_ref = np.interp(data.grid, ref.grid, phi_ref)
        ax3.plot(data.grid, 100 * abs((phi-phi_ref)/(phi_ref + 1e-15)), linestyle = '--', marker = '.',label = labelList[i])
        print(r'$\phi$ is ', 100 * np.sum((phi[1:-1]-phi_ref[1:-1])/phi_ref[1:-1])/(dataList[i].Nx - 2))
        
        ax4.bar(labelList[i], data.totPotTime)
        ax4.text(i, data.totPotTime, data.totPotTime, ha = 'center')

    ax1.legend(loc = 'best')
    ax2.legend(loc = 'best')
    ax3.legend(loc = 'best')
    ax1.set_xlim(ref.grid[0], ref.grid[-1])
    ax1.set_xlabel('Distance (m)', fontsize = 12)
    ax1.set_ylabel(r'100 $\times$ $|(n_e - n_{e,ref})/n_{e,ref}|$', fontsize = 12)
    
    ax2.set_xlim(ref.grid[0], ref.grid[-1])
    ax2.set_xlabel('Distance (m)', fontsize = 12)
    ax2.set_ylabel(r'100 $\times$ $|(n_i - n_{i,ref})/n_{i,ref}|$', fontsize = 12)
    
    ax3.set_xlim(ref.grid[0], ref.grid[-1])
    ax3.set_xlabel('Distance (m)', fontsize = 12)
    ax3.set_ylabel(r'100 $\times$ $|(\phi - \phi_{ref})/\phi_{ref}|$', fontsize = 12)
    
    ax4.set_ylabel(r'Total Wall Time for Solver (s)', fontsize = 12)
    
def compareModelToDatas(dataList, labelList, model, saveFile):
    if (len(dataList) != len(labelList)):
        raise Warning("List of data and labels not the same size!")
        
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    for name in dataList[0].particles.keys():
        if name != 'e':
            ion = name
            break
    M = dataList[0].particles[ion]['mass']
    T_e = dataList[0].T_e
    T_i = dataList[0].T_i
    n_ave = dataList[0].n_ave
    factor = int(math.log10(n_ave))
    for i in range(len(dataList)):
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        ax1.plot(data.grid, n_e/(10**factor), linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        
        n_i = data.getAveDensity(ion)
        ax2.plot(data.grid, n_i/(10**factor), linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
        
        phi = data.getAvePhi()
        ax3.plot(data.grid, phi, linestyle = '--', linewidth = 2, marker = '.',label = labelList[i])
     
    
    ax1.plot(model[0, :], model[2, :]/(10**factor), linewidth = 2, label = 'Model')
    ax2.plot(model[0, :], model[3, :]/(10**factor), linewidth = 2, label = 'Model')
    ax3.plot(model[0, :], model[1, :], linewidth = 2, label = 'Model')
    ax1.legend(loc = 'lower right', fontsize = 12)
    ax2.legend(loc = 'lower right', fontsize = 12)
    ax3.legend(loc = 'lower right', fontsize = 12)
    ax1.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax1.set_xlabel('Distance (m)', fontsize = 14)
    ax1.set_ylabel(r'$n_e$ (10$^{' + str(factor) + r'}$ m$^{-3}$)', fontsize = 14)
    ax1.tick_params(labelsize = 14)
    
    ax2.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax2.set_xlabel('Distance (m)', fontsize = 14)
    ax2.set_ylabel(r'$n_i$ (10$^{' + str(factor) + r'}$ m$^{-3}$)', fontsize = 14)
    ax2.tick_params(labelsize = 14)
    
    ax3.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax3.set_xlabel('Distance (m)', fontsize = 14)
    ax3.set_ylabel('Voltage (V)', fontsize = 14)
    ax3.tick_params(labelsize = 14)
    
    fig1.savefig(saveFile + '_eDensity.pdf', bbox_inches = 'tight')
    fig2.savefig(saveFile + '_ionDensity.pdf', bbox_inches = 'tight')
    fig3.savefig(saveFile + '_voltage.pdf', bbox_inches = 'tight')
    
def compareModelToData(data, model):
    for name in data.particles.keys():
        if name != 'e':
            ion = name
            break
    deb = debye_length(data.T_e, data.n_ave)
    M = data.particles[ion]['mass']
    plt.figure()
    n_e = data.getAveDensity('e')
    plt.plot(data.grid, n_e, 'o', label = 'PIC')
    plt.plot(model[0, :], model[2, :], label = 'Model')
    plt.ylabel(r'$n_e$ (m$^-3$)')
    plt.legend(loc = 'best')
    modelN_e = np.interp(data.grid[1:-1], model[0,:], model[2,:])
    res = np.sum(abs((modelN_e - n_e[1:-1])/modelN_e)/(data.Nx - 2))
    print('Percent diff. root mean square in n_e is:', 100*res)
    
    plt.figure()
    n_i = data.getAveDensity(ion)
    plt.plot(data.grid, n_i, 'o',  label = 'PIC')
    plt.plot(model[0, :], model[3, :], label = 'Model')
    plt.ylabel(r'$n_i$ (m$^-3$)')
    plt.legend(loc = 'best')
    modelN_i = np.interp(data.grid[1:-1], model[0,:], model[3,:])
    res = np.sum(abs((modelN_i - n_i[1:-1])/modelN_i))/(data.Nx - 2)
    print('Percent diff. root mean square in n_i is:', 100*res)
    
    plt.figure()
    phi = data.getAvePhi()
    plt.plot(data.grid, phi, 'o', label = 'PIC')
    plt.plot(model[0, :], model[1, :], label = 'Model')
    plt.ylabel(r'Voltage (V)')
    plt.legend(loc = 'best')
    modelPhi = np.interp(data.grid[1:-1], model[0,:], model[1,:])
    res = np.sum(abs((modelPhi - phi[1:-1])/modelPhi)/(data.Nx - 2))
    print('Percent diff. root mean square in phi is:', 100 * res)
        
def compareModelToDatasRes(dataList, labelList, model):
    if (len(dataList) != len(labelList)):
        raise Warning("List of data and labels not the same size!")
        
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    for name in dataList[0].particles.keys():
        if name != 'e':
            ion = name
            break
    M = dataList[0].particles[ion]['mass']
    T_e = dataList[0].T_e
    T_i = dataList[0].T_i
    n_ave = dataList[0].n_ave
    for i in range(len(dataList)):
        data = dataList[i]
        
        n_e = data.getAveDensity('e')
        modelN_e = np.interp(data.grid, model[0,:], model[2,:])
        ax1.plot(data.grid, 100 * (n_e-modelN_e)/modelN_e, linestyle = '--', marker = '.',label = labelList[i])
        
        n_i = data.getAveDensity(ion)
        modelN_i = np.interp(data.grid, model[0,:], model[3,:])
        ax2.plot(data.grid, 100 * (n_i-modelN_i)/modelN_i, linestyle = '--', marker = '.',label = labelList[i])
        
        phi = data.getAvePhi()
        modelPhi = np.interp(data.grid, model[0,:], model[1,:])
        ax3.plot(data.grid, 100 * (phi-modelPhi)/modelPhi, linestyle = '--', marker = '.',label = labelList[i])
     
    
    ax1.legend(loc = 'best')
    ax2.legend(loc = 'best')
    ax3.legend(loc = 'best')
    ax1.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax1.set_xlabel('Distance (m)', fontsize = 12)
    ax1.set_ylabel(r'100 $\times$ $(n_e - n_{e,ref})/n_{e,ref}$', fontsize = 12)
    
    ax2.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax2.set_xlabel('Distance (m)', fontsize = 12)
    ax2.set_ylabel(r'100 $\times$ $(n_i - n_{i,ref})/n_{i,ref}$', fontsize = 12)
    
    ax3.set_xlim(dataList[0].grid[0], dataList[0].grid[-1])
    ax3.set_xlabel('Distance (m)', fontsize = 12)
    ax3.set_ylabel(r'100 $\times$ $(\phi - \phi_{ref})/\phi_{ref}$', fontsize = 12)   