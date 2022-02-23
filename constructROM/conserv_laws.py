# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 15:03:53 2022

@author: mlols
"""

import sympy as sym
import numpy as np

def conserv_laws(p_modes, t_modes, rhs_psi, rhs_theta):
    """
    Compute selected quantities that should be conserved in the 
    dissipationless limit nu, kappa -> 0.

    Parameters
    ----------
    p_modes : list
        List of psi modes, represented as tuples.
        Each tuple contains the horizontal and vertical wavenumbers.
    t_modes : list
        List of theta modes, represented as tuples.
    rhs_psi : list
        Expressions on the right-hand side of each psi ODE.
    rhs_theta : list
        Expressions on the right-hand side of each theta ODE.

    Returns
    -------
    None.

    """
    k = sym.Symbol('k')
    psi = np.array([sym.symbols('psi_%i_%i' % mode) for mode in p_modes])
    theta = np.array([sym.symbols('theta_%i_%i' % mode) for mode in t_modes])
    rho_p = [(mode[0]*k)**2 + mode[1]**2 for mode in p_modes]

    K = 0
    for i in range(len(p_modes)):
        mode = p_modes[i]
        aux = 1/4*rho_p[i]*psi[i]*rhs_psi[i]
        if mode[0] == 0:
            aux *= 2
        K += aux
            
    U = 0
    for i in range(len(t_modes)):
        mode = t_modes[i]
        if mode[0] == 0:
            n = mode[1]
            U += (-1)**(n+1)/n*rhs_theta[i]
            
    T = 0
    for i in range(len(t_modes)):
        mode = t_modes[i]
        if mode[0] == 0 and (mode[1]+1) % 2 == 0:
            n = mode[1]
            T += rhs_theta[i]/n
    
    T2 = 0 
    for i in range(len(t_modes)):
        mode = t_modes[i]
        if mode[0] == 0:
            T2 += 1/2*theta[i]*rhs_theta[i]
        else:
            T2 += 1/4*theta[i]*rhs_theta[i]
    
    W = 0
    for idx,mode in enumerate(p_modes):
        if (mode[0] == 0) and (mode[1] % 2 == 1):
            W += -2*rho_p[idx]/(mode[1])*rhs_psi[idx]
            
    print('dE/dt = ' + str(sym.simplify(K+U)))
    print('dT/dt = ' + str(sym.simplify(T)))
    print('dT^2/dt = ' + str(sym.simplify(T2 + U)))
    print('dW/dt = ' + str(sym.simplify(W)))
        
    
    return 