from setuptools import setup, find_packages 

with open("README.md", "r") as readme_file:
    readme = readme_file.read()
    
requirements = ["ipython>=6", "numpy>=1", "sympy>=1.5"]

description = ("A Python package to build reduced order models for " + 
              "Rayleigh-Benard convection")

setup(
      name="construct_roms",
      version="0.0.1",
      author="Matt Olson",
      author_email="mlolson@umich.edu",
      description=description,
      long_description=readme,
      long_description_content_type="text/markdown",
      url="https://github.com/PeriodicROM/construct_roms/",
      packages=find_packages(),
      install_requires=requirements,
      classifiers=[
          "Programming Language :: Python :: 3.7",
          "License :: OSI Approved :: MIT License",
      ],
  )