# This is a collections of examples to demonstrate how to use odeSolu

import odeSolu as os
from sympy import symbols,Function,diff
import numpy as np


######## Example 1 (Bernoulli equation) ##################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,1)+(0.001)*y(x)-(0.1)*y(x)**3
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),300)

file = open('Bernoulli.txt','w')
file.write(np.array_str(solu))
file.close()
##########################################################

######## Example 2 (Riccati equation) ####################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,1)-(2*x**3+3)-(5*x**2+2*x+1)*y(x)-(x+7)*y(x)**2
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),300)

file = open('Riccati.txt','w')
file.write(np.array_str(solu))
file.close()
##########################################################

######## Example 3 (Abel's differential equation of the first kind) ####################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,1)+4+(5)*y(x)+(0.1*x)*y(x)**2+(0.2*x**2)*y(x)**3
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),300)

file = open('Abel.txt','w')
file.write(np.array_str(solu))
file.close()
#######################################################################################

######## Example 4 (A second-order ODE with quartic nonlinearity) ####################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,2)+4+(10**(-1))*diff(y(x),x,1)+y(x)**4
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 1)),500)

file = open('Secondorder.txt','w')
file.write(np.array_str(solu))
file.close()
#######################################################################################

######## Example 5 (De Boer-Ludford equation) ####################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,2)+x**2*y(x)-y(x)**4
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0)),500)

file = open('DeBoer.txt','w')
file.write(np.array_str(solu))
file.close()
##################################################################

######## Example 6 (Van der Pol equation) ####################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,2)-0.05*(1-y(x)**2)*diff(y(x),x,1)+y(x)
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5)),500)

file = open('VanderPol.txt','w')
file.write(np.array_str(solu))
file.close()
##############################################################################

######## Example 7 (Painleve-Ince equation) ##################################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,2)+3*y(x)*diff(y(x),x,1)+y(x)**3
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5)),500)

file = open('Painleve.txt','w')
file.write(np.array_str(solu))
file.close()
##############################################################################

######## Example 8 (Blasius equation) ########################################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,3)+y(x)*diff(y(x),x,2)
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 1),(diff(y(x), x, 2), 0)),500)

file = open('Blasius.txt','w')
file.write(np.array_str(solu))
file.close()
##############################################################################

######## Example 9 (Falkner–Skan equation) ####################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,3)+y(x)*diff(y(x),x,2)+2*(1-diff(y(x),x,1)**2)
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0.5),(diff(y(x), x, 2), 1)),500)

file = open('FalknerSkan.txt','w')
file.write(np.array_str(solu))
file.close()
##############################################################################

######## Example 10 (fourth-order ODE) #####################################################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,4)-x**2*diff(y(x),x,3)+3*x*y(x)*diff(y(x),x,2)-6*(diff(y(x),x,1))+2*x**2+x
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5),(diff(y(x), x, 2), 1),(diff(y(x), x, 3), 1)),1000)

file = open('Fourthorder.txt','w')
file.write(np.array_str(solu))
file.close()
############################################################################################

######## Example 11 (5th-order) #####################################################
x = symbols('x')
y = Function('y')
ode = diff(y(x),x,5)-0.001*y(x)**2*diff(y(x),x,4)-2*x*y(x)*diff(y(x),x,3)**2+0.5*x*y(x)*diff(y(x),x,2)**4\
    -diff(y(x),x,1)+y(x)**3*x**2
solu=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0),
                                  (diff(y(x), x, 2), 1),(diff(y(x), x, 3), 1),(diff(y(x), x, 4), 0.5)),500)

file = open('Fifthorder.txt','w')
file.write(np.array_str(solu))
file.close()
############################################################################################