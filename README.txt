#######################################################################################
# odeSolu README file                                                                  #
# Version: 1.0.0                                                                       #
# odeSolu is a Python package that determines power series solution of nonlinear ODE.  #
#######################################################################################

 
(1) External packages: 
    odeSolu version 1.0.0 uses the following external packages:
     (a) NumPy>=1.14.5
     (b) SymPy>=1.0.0
     (c) Cython 
     (d) matplotlib (for plottings)
    
	
    
(2) Installation:
     =================================================
     * Currently odeSolu only supports Python >= 3.8 * 
     ================================================= 

    To install odeSolu from a local copy of the project open a terminal
    in this directory and type: 

 	python setup.py build_ext install


(6) Examples and Testings:
    The "odeSolu-source/examples" folder contains some examples from which we can know the usage of 
    the odeSolu package to solve an ODE. This folder contains two examples, example1.py and example2.py. 
    One can test the package using these example files by running the commands:

	python example1.py
	python example2.py

    To run the example1.py we only need the package odeSolu. example2.py shows how to plot output results from odeSolu, 
    and so to run this Python file, we also require the plotting library matplotlib. The output results obtained from 
    these two example Python files are also provided in this folder.


(9) Brief Description of all files:

   (a) In "odeSolu-source" folder

    README.txt          : This file.
    LICENSE.txt         : License file.
    setup.py            : odeSolu Python installation script.

   (b) In "odeSolu-source/examples" folder

    example1.py         : Some examples to solve nonlinear ODE using odeSolu (without plottings).
    example2.py         : Some examples to solve nonlinear ODE using odeSolu (with plottings).

    Besides these two example files, there present some other 22 files in this folder. These files show the output results produced from
    the example files.
    
   (h) In "odeSolu-source/odeSolu" folder

    __init__.py         : Initialize seriesSolu module. 
    seriesSolu.py       : Implementation of the module seriesSolu. It calculates the series solution using adomianMat module.
    adomianMat.pyx      : Implementation of the module adomianMat. This Cython source code finds Adomian polynomials.
