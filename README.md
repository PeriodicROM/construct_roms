# Construct ROMs
Construct ROMs is a Python package that builds distinguished reduced-order models (ROMs) that approximate the dynamics of 2D Rayleigh&ndash;Bénard convection. The model domain is a horizontally periodic domain with stress-free, isothermal walls at top and bottom. ROMs are built using Galerkin expansion with Fourier basis functions.
## Model equations:
![Boussinesq equations](https://latex.codecogs.com/svg.image?\begin{align*}&space;\partial_t&space;\nabla^2&space;\psi&space;-&space;\\{\psi,&space;\nabla^2&space;\psi\\}&space;&=&space;\sigma&space;\nabla^4&space;\psi&space;&plus;&space;\sigma&space;\mathcal{R}&space;\partial_x&space;\theta,&space;\\\\\partial_t&space;\theta&space;-&space;\\{\psi,&space;\theta\\}&space;&=&space;\nabla^2&space;\theta&space;&plus;&space;\partial_x&space;\psi,\end{align*})

## Features:
- ROMs can be printed to the consoles in LaTeX format or printed to a MATLAB or Python ODE file for further analysis
- Models can be selected from the HK hierarchy or the user can input them manually
- Include scale factors for easier SDP scaling
- Can check certain conservation properties (energy and vorticity balance)

# Package Requirements
- numpy (Version 1.5 or later)
- sympy (Version 1.6 or later)

# Installation
To install the package, either:
- download or clone this repository and use `from constructROM import construct_roms` from the directory containing the package, OR
- install package directly using `pip install construct-roms` then `from constructROM import contruct_roms`

# Instructions
To construct a system of ROMs, use the command `construct_roms(*args)`. Options can be passed as function arguments as detailed below.

## Options:
- `mode_sel` : Method of mode selection. Current options:
  - "hk" : select model from hk hierarchy
  - "input" : input list of modes manually (as p_modes, t_modes arguments)
- `p_modes` : List of psi modes, represented as a list of tuples. Only necessary if `mode_type="input"` 
- `t_modes` : List of theta modes, formatted as a list of tuples. Same format as `p_modes`
- `hier_num` : Number in the HK hierarchy (e.g. `hier_num = 1` denotes the first model)
- `print_to` : Specify where to print output data. Options:
  - "console" : prints ODEs to console in LaTeX format
  - "matlab" : creates Matlab ODE file (.m)
  - "python" : creates Python ODE file (.py)
  - "none:" : don't print ROMs as output
- `out_dir` : Name of directory to print output
- `scaling` :   Variable scaling used in model. Options:
  - "unit cube" : psi -> R^(1/2)*psi, t -> R^(-1/2)*t, useful for SDP computations
  - "standard" : standard scaling of Rayleigh-Benard
  - "normalize theta" : theta -> theta/R, used in some prior ROM studies
- `scale_factors` : If True, includes individual scale factors for each variable (True/False)
- `check_conserv` : If True, verifies that certain conservation laws are satisfied

## Examples:
`construct_roms(mode_sel="hk", hier_num=1, print_to="console" scale factors=False)`  
  Generates the first model in the HK hierarchy (Lorenz model plus one shear mode) without scale factors and prints output to the console
  
`construct_roms(mode_sel="hk", hier_num=2, print_to="matlab", out_dir="HK8_model")`  
  Generates Matlab ODE file for the HK8 model with scale factors to folder "HK8_model"
    
`construct_roms(model_sel="input", p_modes=[(1,1)], t_modes=[(0,2), (1,1)], print_to="none", check_conserv=True)`  
  Checks energy and vorticity balance laws for the Lorenz equations
