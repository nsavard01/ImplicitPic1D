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

def getBoundPlasmaSolutions(L, numNodes, n_ave, T_e, T_i, M):
    fracDebye = 0.2
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
    for i in range(500):
        n_0 = n_ave / np.exp(eta[0])/ np.trapz(np.exp(-eta) * dah * (1-t)**(dah-1), dx = dt)
        A = eps_0 * T_e / n_0 / e / L**2
        for i in range(numNodes):
            I_int = np.where(eta[i] < eta, (1 - t)**(dah-1) * np.exp(-eta * tau), (1 - t)**(dah-1) * np.exp(-eta*tau) * scipy.special.erfc(np.sqrt((eta[i] - eta)*tau)))
            I_G[i] = np.trapz(I_int, dx = dt)
        I_G = I_G * dah
        currRes = np.sum(funcNonLinear(eta, A, B, tau, I_G, t, dah)**2)/numNodes
        eta = opt.fsolve(funcNonLinear, eta, args = (A, B, tau, I_G, t, dah))
        if (currRes < initialR*1e-8):
            break
        if (i == 499):
            raise Warning('Iteration has not converged for bounded plasma model')
    n_e = n_0 * np.exp(eta[0] - eta)
    n_i = n_0 * np.exp(eta*tau + eta[0]) * B * I_G
    phi = -eta * T_e
    return L-s*L, phi, n_e, n_i
    