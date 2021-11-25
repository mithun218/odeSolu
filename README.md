--------------------------------------------
## odeSolu                                v1.0.0                                                                       #
### odeSolu is a Python package that determines power series solution of nonlinear ODE.  #
___________
-----------
 
## External packages: 
    odeSolu version 1.0.0 uses the following external packages:
     1. NumPy>=1.14.5
     2. SymPy>=1.0.0
     3. Cython 
     4. matplotlib (for plottings)
    
	
    
## Installation:

     **Currently odeSolu only supports Python >= 3.8 **

    To install odeSolu from a local copy of the project open a terminal
    in this directory and type: 

 	python setup.py build_ext install


## Examples and Testings:
    The "odeSolu-source/examples" folder contains some examples from which we can know the usage of 
    the odeSolu package to solve an ODE. This folder contains two examples, example1.py and example2.py. 
    One can test the package using these example files by running the commands:

	python example1.py
	python example2.py

    To run the example1.py we only need the package odeSolu. example2.py shows how to plot output results from odeSolu, 
    and so to run this Python file, we also require the plotting library matplotlib. The output results obtained from 
    these two example Python files are also provided in this folder.


## Brief Description of all files:

   1. In "odeSolu-source" folder

    README.txt          : This file.
    LICENSE.txt         : License file.
    setup.py            : odeSolu Python installation script.

   2. In "odeSolu-source/examples" folder

    example1.py         : Some examples to solve nonlinear ODE using odeSolu (without plottings).
    example2.py         : Some examples to solve nonlinear ODE using odeSolu (with plottings).

    Besides these two example files, there present some other 22 files in this folder. These files show the output results produced from
    the example files.
    
   3. In "odeSolu-source/odeSolu" folder

    __init__.py         : Initialize seriesSolu module. 
    seriesSolu.py       : Implementation of the module seriesSolu. It calculates the series solution using adomianMat module.
    adomianMat.pyx      : Implementation of the module adomianMat. This Cython source code finds Adomian polynomials.
