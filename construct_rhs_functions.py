import sympy as sym
import numpy as np


def construct_psi(p_modes, t_modes, a=None, dissip_limit=False):
    """
    Construct right-hand side of ODEs for psi modes

    Parameters
    ----------
    p_modes : list
        List of psi modes, represented as tuples.
        Each tuple contains the horizontal and vertical wavenumbers.
    t_modes : list
        List of theta modes, represented as tuples.
    a : array, optional
        Scale factors for model rescaling. Can be numbers or symbols.
        To produce unscaled model, pass in integer array of ones.
    dissip_limit : bool, optional
        If True, computes model in the dissipationless limit nu, kappa -> 0.
        Used for checking conservation properties.

    Returns
    -------
    rhs_psi : list
        Expressions on the right-hand side of each psi ODE.
    """    
    psi = np.array([sym.symbols('psi_%i_%i' % mode) for mode in p_modes])
    theta = np.array([sym.symbols('theta_%i_%i' % mode) for mode in t_modes])
    k = sym.Symbol('k')
    if dissip_limit:
        sigma = 1
    else:
        sigma = sym.Symbol('sigma')
    R = sym.Symbol('r')
    
    if a is None:
        a = np.ones(len(p_modes) + len(t_modes), dtype=int)
    
    rho_p = [(mode[0]*k)**2 + mode[1]**2 for mode in p_modes]
    rhs_psi = []
    
    #Loop over all psi modes
    for idx,mode in enumerate(p_modes):
        
        #Compute linear terms
        rhs_term1 = -sigma*rho_p[idx]*psi[idx]/R
        sgn = sum(mode)
        if mode in t_modes:
            t_index = t_modes.index(mode)
            coeff1 = (-1)**(sgn)*sigma*mode[0]*k
            coeff2 = a[t_index+len(p_modes)]/a[idx]
            rhs_term2 = coeff1*coeff2*theta[t_index]/rho_p[idx]
        else:
            rhs_term2 = 0
            
        
        rho_p = [(mode[0]*k)**2 + mode[1]**2 for mode in p_modes]
        
        #Compute nonlinear terms from poisson bracket
        m,n = mode
        nonlin = 0
        idx = p_modes.index(mode)
        for idx1,t1 in enumerate(p_modes): #Loop over psi modes
            r,s = t1
            #Find compatible triples, avoiding duplicates
            compat = [(p,q) for (p,q) in p_modes[idx1:] if 
                      (m == (p+r) or m == p-r) and 
                      (n == (q+s) or n == abs(q-s))]
            for (p,q) in compat:
                idx2 = p_modes.index((p,q))
                mu_1 = -1
                if (r+s) % 2 == 1: #Some cases for sign flips
                    mu_1 *= B(p,m,r)
                    if (m+n) % 2 == 0:
                        mu_1 *= -1
                if r == 0: #Compute divisor
                    d = 2
                else:
                    d = 4
                coeff1 = (B(p,m,r)*B(s,n,q)*p*s - B(q,n,s)*q*r)
                coeff2 = a[idx1]*a[idx2]/a[idx]
                coeff3 = k/d*mu_1*(rho_p[idx2] - rho_p[idx1])
                nonlin += coeff1*coeff2*coeff3*psi[idx1]*psi[idx2]
        
        rhs_term3 = nonlin/rho_p[idx]
        if dissip_limit:
             rhs_psi.append(rhs_term2+rhs_term3)
        else:
             rhs_psi.append(rhs_term1+rhs_term2+rhs_term3)
        
    return rhs_psi

def construct_theta(p_modes, t_modes, a=None, dissip_limit=False):
    """
    Construct right-hand side of ODEs for theta modes

    Parameters
    ----------
    p_modes : list
        List of psi modes, represented as tuples.
        Each tuple contains the horizontal and vertical wavenumbers.
    t_modes : list
        List of theta modes, represented as tuples.
    a : array, optional
        Scale factors for model rescaling. Can be numbers or symbols.
        To produce unscaled model, pass in integer array of ones.
    dissip_limit : bool, optional
        If True, computes model in the dissipationless limit nu, kappa -> 0.
        Used for checking conservation properties.

    Returns
    -------
    rhs_theta : list
        Expressions on the right-hand side of each theta ODE.
    """    
    psi = np.array([sym.symbols('psi_%i_%i' % mode) for mode in p_modes])
    theta = np.array([sym.symbols('theta_%i_%i' % mode) for mode in t_modes])
    k = sym.Symbol('k')
    R = sym.Symbol('r')
    if a is None:
        a = np.ones(len(p_modes) + len(t_modes), dtype=int)
    
    rho_t = [(mode[0]*k)**2 + mode[1]**2 for mode in t_modes]
    
    rhs_theta = []
    #Loop over all theta modes
    for idx,mode in enumerate(t_modes):
        
        #Compute linear terms
        rhs_term1 = -rho_t[idx]*theta[idx]/R
        sgn = sum(mode)
        if mode in p_modes:
            p_index = p_modes.index(mode)
            coeff1 = (-1)**(sgn)*mode[0]*k
            coeff2 = a[p_index]/a[idx+len(p_modes)]
            rhs_term2 = coeff1*coeff2*psi[p_index]
        else:
            rhs_term2 = 0
            
        m, n = mode 
        
        nonlin = 0
        for idx1,t1 in enumerate(t_modes):
            r,s = t1
            #Find compatible triples, avoiding duplicates
            compat = [(p,q) for (p,q) in p_modes if 
                  (m == (p+r) or m == abs(p-r)) and 
                  (n == (q+s) or n == abs(q-s))]
            for (p,q) in compat:
                idx2 = p_modes.index((p,q))
                mu_2 = 1  
                if (p+q) % 2 == 1: #Cases for sign flips
                    mu_2 *= B(r,p,m)
                    if (m+n)% 2 == 1:
                        mu_2 *= -1
                if (r+s) % 2 == 0:
                    mu_2 *= B(p,m,r)
                if m == 0 and (p+q)%2 == 1:
                    mu_2 *= -1
                if r == 0 or p == 0: #Compute Divisor
                    d = 2
                else:
                    d = 4
                if m == 0:
                    mu_3 = -1
                else:
                    mu_3 = 1
                coeff1 = (B(p,m,r)*B(s,n,q)*p*s - mu_3*B(q,n,s)*B(r,p,m)*q*r)
                coeff2 = a[idx1+len(p_modes)]*a[idx2]/a[idx+len(p_modes)]
                coeff3= mu_2*k/d
                nonlin += coeff1*coeff2*coeff3*psi[idx2]*theta[idx1]   
            
        rhs_term3 = nonlin
        
        if dissip_limit:
            rhs_theta.append(rhs_term2+rhs_term3)
        else:
            rhs_theta.append(rhs_term1+rhs_term2+rhs_term3)
    return rhs_theta


def B(i,j,k):
    """
    Tensor B used in constructing ROMs. 

    Parameters
    ----------
    i : int
    j : int
    k : int
        Indices in the tensor.

    Returns
    -------
    int
        Tensor output.
    """
    if i == j + k:
        return -1
    elif j == i + k or k == i + j:
        return 1
    else:
        msg = "Possible Error: Indices ({},{},{})".format(i,j,k)
        print(msg)
        return 0


    


