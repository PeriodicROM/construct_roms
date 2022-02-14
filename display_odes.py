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

def matlab_print(func, psi, theta, file):
    """
    Wrapper to print ROM to a matlab file.

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
    rhs = func
    new_string = convert_rhs(rhs, psi, theta)
    file.write(new_string)
    return 

        
# def convert_rhs(rhs,psi,theta):
#     old_vars = [str(p) for p in psi] + [str(t) for t in theta]
    
#     new_vars = ['x({})'.format(i) for i in range(1,len(old_vars)+1) ]
    
#     new_string = str(rhs)[1:-1]
    
#     for i in range(len(old_vars)):
#         new_string = new_string.replace(old_vars[i],new_vars[i])
      
#     new_string = new_string.replace('sig','s')
#     new_string = new_string.replace('**','^')
#     new_string = new_string.replace(',',',...\n')
    
    
#     return new_string

def convert_rhs(rhs, psi, theta):
    """
    Convert equations to matlab format and clean up.

    Parameters
    ----------
    rhs : string
        Original equations (python sym form)
    psi : symbolic array
        Array of psi modes.
    theta : symbolic array
        Array of theta modes.

    Returns
    -------
    new_string : string
        Equations in matlab format.

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

    
    return new_string


        

