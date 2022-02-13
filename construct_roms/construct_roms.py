import sympy as sym
import numpy as np
from timeit import default_timer as timer
from os import path, mkdir, chdir
import sys
from display_odes import latex_print, matlab_print
import construct_rhs_functions as RHS
from hk_modes import hk_modes


def construct_roms(mode_sel = 'hk', p_modes = [(1,1)],
                   t_modes = [(1,1), (0,2)], 
                   hier_num=1, 
                   print_matlab=True, print_console=False, 
                   scale_factors=True, check_conserv=False, 
                   out_dir='Matlabfiles'):
    """
    Construct ROM for the Rayleigh--Benard system. Each
    model satisfies stress-free, isothermal boundaries in a rectangular
    domain.

    Parameters
    ----------
    mode_sel : string, optional
        Method of mode selection. Options:
            'hk' : select model from hk hierarchy (model number = hier_num)
            'input' : input mode list manually (p_modes, t_modes)
    p_modes : list, optional
        List of psi modes, represented as tuples.
        Each tuple contains the horizontal and vertical wavenumbers.
        Only necessary if mode_type = 'input'. 
        Default (Lorenz): [(1,1)].
    t_modes : list, optional 
        List of theta modes, represented as tuples.
        Only necessary if mode_type = 'input'.
        Default (Lorenz): [(0,2), (1,1)]
    hier_num : int, optional
        Number in the HK hierarchy (hier_num = n means the nth model).
    print_matlab : bool, optional
        If True, print ODE equations to a matlab file. 
        Modes will be formatted as x(i) in order [p_modes, t_modes].
        The default is True.
    print_console : bool, optional
        If True, print equations to the console. 
        The default is False.
    scale_factors : bool, optional
        If True, includes scale factors for each variable, 
        i.e. x(i) -> a(i)*x(i). Useful for numerical scaling for SDP.
        The default is True.
    check_conserv : bool, optional
        If True, verifies that certain conservation laws are satisfied.
        The default is False.
    out_dir : string, optional
        Name of directory to output matlab files. Must be of form dir1/dir2.
        If print_matlab is False, this argument does nothing.    
        The default is 'Matlabfiles'

    Returns
    -------
    None.
    
    Examples
    --------
    construct_roms()
    """   
    start = timer()
    
    if type(hier_num) is not int:
        print('Error: hier_num must be an integer. Input the model number' + 
              'in the HK hierarchy.')
        return
    if type(print_matlab) is not bool:
        print('Error: print_matlab must be True or False.')
    if type(print_console) is not bool:
        print('Error: print_console must be True or False.')
    if type(check_conserv) is not bool:
        print('Error: check_conserv must be True or False.')
    if type(scale_factors) is not bool:
        print('Error: scale_factors must be True or False.')
        
    #Make necessary folders
    if print_matlab and not path.isdir(out_dir):
        dir_names = out_dir.split(sep='/')
        for name in dir_names:
            mkdir(name)
            chdir(name)
        chdir('../'*len(dir_names)) #Back to parent directory
            
    if mode_sel == 'hk':
        p_modes, t_modes = hk_modes(hier_num)
    elif mode_sel == 'input':
        #Make sure modes are in correct order
        p_modes.sort()
        t_modes.sort()
    
    num_modes = len(p_modes) + len(t_modes)
    
    if scale_factors:
        a = np.array([sym.symbols('a(%i)' % (n)) for n in range(1,num_modes+1)])
    else:
        a = np.ones(num_modes, dtype=int)
        
    rhs_psi = RHS.construct_psi(p_modes, t_modes, a=a)
    rhs_theta = RHS.construct_theta(p_modes, t_modes, a=a)
    
    psi = np.array([sym.symbols('psi_%i_%i' % mode) for mode in p_modes])
    theta = np.array([sym.symbols('theta_%i_%i' % mode) for mode in t_modes])
        
    if print_console:
        latex_print(rhs_psi)
        latex_print(rhs_theta)
        
    #Construct eqns in dissipationless limit and check quantities are conserved
    #Currently checking energy, total temperature, and vorticity
    if check_conserv:
        from conserv_laws import conserv_laws
        a = np.ones(num_modes, dtype=int)
        rhs_psi = RHS.construct_psi(p_modes, t_modes, a=a, dissip_limit=True)
        rhs_theta = RHS.construct_theta(p_modes, t_modes, a=a, 
                                        dissip_limit=True)
        conserv_laws(p_modes, t_modes, rhs_psi, rhs_theta)

        
    if print_matlab:
        fname = 'Matlabfiles/hk' + str(num_modes) + '_rhs.m'
        if path.exists(fname):
            print('System already constructed')
            overwrite = input('Overwrite file (y/n)? ')
            if overwrite != 'y':
                end = timer()
                print('Time = ' + str(end-start))
                sys.exit()
        file = open(fname,'w')
        file.write('function f = hk' + str(num_modes) + '_rhs(x,a,s,k,R) \n')
        pstr = ' '.join(['p%i,%i' % mode for mode in p_modes]) 
        tstr = ' '.join(['t%i,%i' % mode for mode in t_modes])
        file.write('% Vars: x = [' + pstr + ' ' + tstr + '] \n \n')
        file.write('r = R^(1/2);\n\n')
        file.write('f = [')
        matlab_print(rhs_psi, psi, theta, file)
        file.write(',...\n')
        matlab_print(rhs_theta, psi, theta, file)
        file.write('];')
        file.close()
    
    return
