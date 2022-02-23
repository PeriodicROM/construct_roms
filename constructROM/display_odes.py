def latex_print(rhs):
    """
    Wrapper to print ROM to console.

    Parameters
    ----------
    func : symbolic
        Expression to be printed (rhs_psi or rhs_theta).

    Returns
    -------
    None.

    """
    from IPython.display import display
    import sympy as sym
    sym.init_printing(forecolor='White')
    [display(eqn) for eqn in rhs]
    sym.init_printing(pretty_print=False)
    
    return 

def matlab_print(rhs, psi, theta, file):
    """
    Wrapper to print ROM to a matlab ODE file.

    Parameters
    ----------
    func : symbolic
        Expression to be printed (rhs_psi or rhs_theta).
    psi : symbolic array
        Array of symbols for psi modes
    theta : symbolic array
        Array of symbols for theta modes
    file : text file
        File to print on.

    Returns
    -------
    None
    """
    old_vars = [str(p) for p in psi] + [str(t) for t in theta]
    
    new_vars = ['x({})'.format(i) for i in range(1,len(old_vars)+1) ]
    
    new_string = str(rhs)[1:-1]
    N = len(old_vars)
    
    for i in range(len(old_vars)):
        new_string = new_string.replace(old_vars[N-i-1],new_vars[N-i-1])
    
    new_string = new_string.replace('sigma','s')
    new_string = new_string.replace('**','^')
    new_string = new_string.replace(',',',...\n')

    file.write(new_string)
    
    return 

def python_print(rhs, psi, theta, file):
    """
    Wrapper to print ROM to a python ODE file.

    Parameters
    ----------
    func : symbolic
        Expression to be printed (rhs_psi or rhs_theta).
    psi : symbolic array
        Array of symbols for psi modes
    theta : symbolic array
        Array of symbols for theta modes
    file : text file
        File to print on.

    Returns
    -------
    None
    """
    new_string = str(rhs)[1:-1]
    
    old_vars_x = [str(p) for p in psi] + [str(t) for t in theta]
    N = len(old_vars_x)
    new_vars_x = ['x[{}]'.format(i) for i in range(N)]
    
    old_vars_a = ['a({})'.format(i+1) for i in range(N)]
    new_vars_a = ['a[{}]'.format(i) for i in range(N)]    
    
    for i in range(N):
        new_string = new_string.replace(old_vars_x[N-i-1],new_vars_x[N-i-1])
        new_string = new_string.replace(old_vars_a[N-i-1],new_vars_a[N-i-1])
        
    new_string = new_string.replace('sigma','s')
    new_string = new_string.replace(',',',\n\t\t')
    
    file.write(new_string)
    
    return 

