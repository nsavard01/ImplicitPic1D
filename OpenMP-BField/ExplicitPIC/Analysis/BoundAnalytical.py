# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:50:57 2023

@author: Nicolas
"""
import numpy as np
import scipy
import scipy.optimize as opt







def otherFunc(phi, tau):
    coeff_1 = 2.0 * np.exp(-tau *  phi**2)/np.sqrt(np.pi * tau)
    return coeff_1 * scipy.special.dawsn(phi) + scipy.special.erf(np.sqrt(tau) * phi) - 1

def funcPhi_s(phi, B, s):
    return (2/np.pi/B) * scipy.special.dawsn(phi) - s
  

# phi_1 = getPhi1(tau)
# phi_w = getPhiWall(phi_1)
# B = 0.5 * np.sqrt(m_p/m_e/np.pi) * (tau / (1 + tau)) * np.exp(phi_w)

# phi[0] = - phi_w * T_e
# phi[-1] = 0.0
# for i in range(1, numNodes-1):
#     phi[i] = phi[0] - (opt.fsolve(funcPhi_s, 0.2, args = (B, s[i]))**2) * T_e
    
# n_e_0 = n_ave * L / np.trapz(np.exp((phi - phi[0])/T_e), dx = dz)
# n_e = n_e_0 * np.exp((phi - phi[0])/T_e)
# S_0 = n_e_0 * np.sqrt(e*T_e/2/np.pi/m_e) * np.exp(-phi[0]/T_e)/L
# n_i = np.zeros(numNodes)
# n_i_coeff = S_0 * np.sqrt(2 * np.pi * m_p * e * T_i) / 2 / e / T_i

# n_i[0] = np.trapz(np.exp((phi-phi[0])/T_i), dx = dz)
# for i in range(1, numNodes):
#     eq_less = np.exp((phi[:(i+1)]-phi[i])/T_i) * scipy.special.erfc(np.sqrt((phi[:(i+1)]-phi[i])/T_i))
#     eq_more = np.exp((phi[i:]-phi[i])/T_i)
#     n_i[i] = np.trapz(eq_less, dx = dz) + np.trapz(eq_more, dx=dz)
# n_i = n_i * n_i_coeff

# res = np.zeros(numNodes-2)

# laplace = (phi[0:-2] - 2 * phi[1:-1] + phi[2::])/(dz**2)
# res = (e / eps_0) * (n_e[1:-1] - n_i[1:-1]) - laplace

# #unknowns
# X_0 = phi[:-1]

def getPhiNeAnalBound(L, n_ave, T_e, T_i, M_m):
    tau = T_e/T_i
    phi_1 = 0
    if (tau < 36):
        phi_1 = opt.fsolve(otherFunc, 0.5, args = (tau,))**2
    else:
        phi_1 = 0.854
    
    coeff_1 = np.pi * np.sqrt(M_m / 4 / np.pi) / (1 + T_i/T_e) / 2 / scipy.special.dawsn(np.sqrt(phi_1))
    phi_w = -np.log(coeff_1)
    B = 0.5 * np.sqrt(M_m/np.pi) * (tau / (1 + tau)) * np.exp(phi_w)
    
    s = np.linspace(0, 1, 1000)
    dz = (s[1] - s[0])*L
    phi = np.zeros(s.size)
    phi[0] = - phi_w * T_e
    phi[-1] = 0.0
    for i in range(1, 999):
        phi[i] = phi[0] - (opt.fsolve(funcPhi_s, 0.2, args = (B, s[i]))**2) * T_e
       
    
    n_e_0 = n_ave * L / np.trapz(np.exp((phi - phi[0])/T_e), dx = dz)
    n_e = n_e_0 * np.exp((phi - phi[0])/T_e)
    return s*L, phi, n_e
    


# trapCoeff = np.trapz(np.exp((X_0-X_0[0])/T_e), dx = dz) + dz * (np.exp((X_0[-1] - X_0[0])/T_e) + np.exp(-X_0[0]/T_e)) * 0.5
# n_e_0 = n_ave * L / trapCoeff
# S_0 = n_e_0 * np.sqrt(e*T_e/2/np.pi/m_e) * np.exp(-X_0[0]/T_e)/L
# n_i_coeff = S_0 * np.sqrt(2 * np.pi * m_p * e * T_i) / 2 / e / T_i
# n_i_test = np.zeros(numNodes-1)
# lastStep = dz * (np.exp((X_0[-1] - X_0[0])/T_i) + np.exp(-X_0[0]/T_i)) * 0.5
# n_i_test[0] = np.trapz(np.exp((X_0-X_0[0])/T_i), dx = dz) + lastStep
# for i in range(1, numNodes-1):
#     eq_less = np.exp((X_0[:(i+1)]-X_0[i])/T_i) * scipy.special.erfc(np.sqrt((X_0[:(i+1)]-X_0[i])/T_i))
#     eq_more = np.exp((X_0[i:]-X_0[i])/T_i)
#     lastStep = dz * (np.exp((X_0[-1] - X_0[i])/T_i) + np.exp(-X_0[i]/T_i)) * 0.5
#     n_i_test[i] = np.trapz(eq_less, dx = dz) + np.trapz(eq_more, dx=dz) + lastStep
# n_i_test = n_i_test * n_i_coeff
# n_e_test = np.exp((X_0 - X_0[0])/T_e) * n_e_0


# def func(X, n_e, n_i):
#     sol = np.zeros(X.size)
    

#     laplace_test = (X[0:-2] - 2 * X[1:-1] + X[2::])/(dz**2)
#     sol =(e/eps_0) * ( n_e - n_i)
#     sol[1:-1] = sol[1:-1] - laplace_test
#     sol[0] = sol[0] - (-2 * X[0] + 2 * X[1]) / (dz**2)
#     sol[-1] = sol[-1] - (X[-2] - 2 *X[-1]) / (dz**2)
#     return sol
    
# for i in range(2):
#     trapCoeff = np.trapz(np.exp((X_0-X_0[0])/T_e), dx = dz) + dz * (np.exp((X_0[-1] - X_0[0])/T_e) + np.exp(-X_0[0]/T_e)) * 0.5
#     n_e_0 = n_ave * L / trapCoeff
#     S_0 = n_e_0 * np.sqrt(e*T_e/2/np.pi/m_e) * np.exp(-X_0[0]/T_e)/L
#     n_i_coeff = S_0 * np.sqrt(2 * np.pi * m_p * e * T_i) / 2 / e / T_i
#     n_i_test = np.zeros(numNodes-1)
#     lastStep = dz * (np.exp((X_0[-1] - X_0[0])/T_i) + np.exp(-X_0[0]/T_i)) * 0.5
#     n_i_test[0] = np.trapz(np.exp((X_0-X_0[0])/T_i), dx = dz) + lastStep
#     for i in range(1, numNodes-1):
#         eq_less = np.exp((X_0[:(i+1)]-X_0[i])/T_i) * scipy.special.erfc(np.sqrt((X_0[:(i+1)]-X_0[i])/T_i))
#         eq_more = np.exp((X_0[i:]-X_0[i])/T_i)
#         lastStep = dz * (np.exp((X_0[-1] - X_0[i])/T_i) + np.exp(-X_0[i]/T_i)) * 0.5
#         n_i_test[i] = np.trapz(eq_less, dx = dz) + np.trapz(eq_more, dx=dz) + lastStep
#     n_i_test = n_i_test * n_i_coeff
#     n_e_test = np.exp((X_0 - X_0[0])/T_e) * n_e_0
    
#     final = opt.fsolve(func, X_0, args = (n_e_test, n_i_test))
#     X_0 = final
    


