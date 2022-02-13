# Construct ROMs
Construct ROMs is a Python package that builds distinguished reduced-order models that approximate the dynamics of 2D Rayleigh-Benard convection. Features:
- ROMs can be printed in LaTeX format or printed to a MATLAB ODE file for further analysis
- ROMs built from Fourier-truncated Galerkin expansions applied to Rayleigh-Benard PDE
- Boundary conditions: horizontally periodic rectangular domain with free-slip, isothermal boundaries above and below
- Models can be selected from an existing hierarchy (HK) or the user can input them manually
- Can check certain conservation properties (energy and vorticity balance)
