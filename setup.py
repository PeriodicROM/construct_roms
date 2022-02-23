from setuptools import setup, find_packages 

with open("README.md", "r") as readme_file:
    readme = readme_file.read()
    
requirements = ["ipython>=6", "numpy", "sympy"]

description = ("A Python package to build reduced order models for " + 
              "Rayleigh-Benard convection")

setup(
      name="construct_roms",
      version="0.0.5",
      author="Matt Olson",
      author_email="mlolson@umich.edu",
      description=description,
      long_description=readme,
      long_description_content_type="text/markdown",
      url="https://github.com/PeriodicROM/construct_roms/",
      packages=["construct_roms", "construct_rhs_functions", "display_odes", "conserv_laws", "hk_modes"],
      install_requires=requirements,
      classifiers=[
          "Programming Language :: Python :: 3.7",
          "License :: OSI Approved :: MIT License",
      ],
) 