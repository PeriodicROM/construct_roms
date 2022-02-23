import sympy as sym
import numpy as np
from timeit import default_timer as timer
from os.path import isdir, exists
from .construct_rhs_functions import construct_psi, construct_theta
from .hk_modes import hk_modes
from .display_odes import latex_print, matlab_print, python_print
from .conserv_laws import conserv_laws

def construct_roms(mode_sel = 'hk', p_modes = [(1,1)],
                   t_modes = [(1,1), (0,2)], hier_num=1, 
                   print_to='console', out_dir='ODE_files', 
                   scaling='unit cube',
                   scale_factors=True, check_conserv=False):
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
    print_to : string, optional
        Specify where to print output data. 
            'console' : prints ODEs to console in LaTeX format.
            'matlab' : creates Matlab ODE file (.m)
            'python' : creates Python ODE file (.py)
            'none' : don't print ROMs as output
        Default: 'console'
    out_dir : string, optional
        Name of directory to output matlab files. Must be of form dir1/dir2.
        If print_matlab is False, this argument does nothing.    
        The default is 'ODE_files'
    scaling : string, optional
        Variable scaling used in model
        "unit cube" : psi -> R^(1/2)*psi, t -> R^(-1/2)*t, useful for
            SDP computations (Default)
        "standard" : standard scaling of Rayleigh-Benard
        "normalize theta" : theta -> theta/R, used in some prior ROM studies
    scale_factors : bool, optional
        If True, includes scale factors for each variable, 
        i.e. x(i) -> a(i)*x(i). Useful for numerical scaling for SDP.
        The default is True.
    check_conserv : bool, optional
        If True, verifies that certain conservation laws are satisfied.
        The default is False.

    Returns
    -------
    None.
    
    Examples
    --------
    construct_roms(mode_sel='hk', hier_num=2, out_dir='HK8_model')
        Generates Matlab ODE file for the HK8 model with scale factors in the
        directory 'HK8_model'.
    construct_roms(mode_sel='hk', hier_num=1, print_to='console', 
                   scale_factors=False)
        Prints ODEs for the HK4 model to the console without scale factors.
    construct_roms(model_sel='input', p_modes=[(1,1)], t_modes=[(0,2), (1,1)],
                   print_to='none', check_conserv=True)
        Checks energy and vorticity balance laws for the Lorenz equations.
    
    """   
    start = timer()
    
    #%% Check for errors, initialize
    if type(hier_num) is not int:
        print('Error: hier_num must be an integer. Input the model number' + 
              'in the HK hierarchy.')
        return
    
    if type(print_to) is not str:
        print('Error: print_to must be a string')
        print('Choose from: "console", "matlab", "python", "none"')
        return
    
    print_to = print_to.lower()
    
    if print_to not in ['console', 'matlab', 'python', 'none']:
        print('Error: Invalid option for print_to.')
        print('Choose from: "console", "matlab", "python", "none"')
        return 
    
    if type(out_dir) is not str:
        print('Error: out_dir must be a string')
        return
    
    if print_to in ['matlab','python'] and not isdir(out_dir):
        print('Error: Directory ' + out_dir + ' does not exist.')
        return
    
    if type(mode_sel) is not str:
        print('Error: mode_sel must be a string')
        print('Choose from: "hk", "input"')
        return
    
    mode_sel = mode_sel.lower()
    
    if mode_sel not in ['hk','input']:
        print('Error: Invalid option for mode_sel')
        print('Choose from: "hk", "input"')
        return
    
    if type(scaling) is not str:
        print('Error: scaling must be a string')
        print('Choose from: "unit cube", "standard", "normalize theta"')
        return
    
    scaling = scaling.lower()
    
    if scaling not in ['unit cube','standard', 'normalize theta']:
        print('Error: Invalid option for mode_sel')
        print('Choose from: "unit cube", "standard", "normalize theta"')
        return
    
    if type(scale_factors) is not bool:
        print('Error: scale_factors must be bool')
        return
    
    if type(check_conserv) is not bool:
        print('Error: check_conserv must be bool')
        return
            
    if mode_sel == 'hk':
        p_modes, t_modes = hk_modes(hier_num)
    elif mode_sel == 'input':
        #Make sure modes are in correct order
        p_modes.sort()
        t_modes.sort()
        
    logic_1 = type(p_modes) is list 
    logic_2 = all(isinstance(n, tuple) and len(n) == 2 for n in p_modes)
    if not (logic_1 and logic_2):
        print('Incorrect format for p_modes. See documentation')
        return
    logic_1 = type(t_modes) is list 
    logic_2 = all(isinstance(n, tuple) and len(n) == 2 for n in t_modes)
    if not (logic_1 and logic_2):
        print('Incorrect format for t_modes. See documentation')
        return
    
    #%% Construt ROMs
    num_modes = len(p_modes) + len(t_modes)
    
    if scale_factors:
        a = np.array([sym.symbols('a(%i)' % n) for n in range(1,num_modes+1)])
    else:
        a = np.ones(num_modes, dtype=int)
        
    rhs_psi = construct_psi(p_modes, t_modes, a=a, scaling=scaling)
    rhs_theta = construct_theta(p_modes, t_modes, a=a, scaling=scaling)
    
    psi = np.array([sym.symbols('psi_%i_%i' % mode) for mode in p_modes])
    theta = np.array([sym.symbols('theta_%i_%i' % mode) for mode in t_modes])
        
    #%% Print output
    
    if print_to == 'console':
        if num_modes > 20:
            print('Console output will be very large.')
            opt = input('Are you sure you sure you want to continue (y/n)? ')
            if opt != 'y':
                return
        latex_print(rhs_psi)
        latex_print(rhs_theta)
        
    #Construct eqns in dissipationless limit and check quantities are conserved
    #Currently checking energy, total temperature, and vorticity
    if check_conserv:
        a = np.ones(num_modes, dtype=int)
        rhs_psi = construct_psi(p_modes, t_modes, a=a, dissip_limit=True)
        rhs_theta = construct_theta(p_modes, t_modes, a=a, 
                                        dissip_limit=True)
        conserv_laws(p_modes, t_modes, rhs_psi, rhs_theta)

        
    if print_to in ['matlab', 'python']:
        if out_dir[-1] != '/':
            out_dir += '/'
        if print_to == 'matlab':
            ext = '.m'
        else:
            ext = '.py'
            
        fname = out_dir + 'hk' + str(num_modes) + ext
        if exists(fname):
            print('System already constructed')
            overwrite = input('Overwrite file (y/n)? ')
            if overwrite != 'y':
                end = timer()
                print('Time = ' + str(end-start))
                return
        file = open(fname,'w')
        
        pstr = ' '.join(['p%i,%i' % mode for mode in p_modes]) 
        tstr = ' '.join(['t%i,%i' % mode for mode in t_modes])
        var_names = 'Vars: x = [' + pstr + ' ' + tstr + ']'
        
        if print_to == 'matlab':
            
            file.write('function f = hk' + str(num_modes) + 
                       '(x, a, s, k, R) \n')
 
            file.write('% ' + var_names + '\n \n')
            if scaling == 'unit cube':
                file.write('r = R^(1/2);\n\n')
            file.write('f = [')
            matlab_print(rhs_psi, psi, theta, file)
            file.write(',...\n')
            matlab_print(rhs_theta, psi, theta, file)
            file.write('];')
        else:
            file.write('def hk' + str(num_modes) + 
                       '(x, a, s, k, R): \n')
            file.write('# ' + var_names + '\n \n')
            if scaling == 'unit cube':
                file.write('\tr = R**(1/2) \n\n')
            file.write('\tf = [')
            python_print(rhs_psi, psi, theta, file)
            file.write(',\n\t\t')
            python_print(rhs_theta, psi, theta, file)
            file.write(']\n\n\t')
            file.write('return f')
            
        file.close()
    
    return
