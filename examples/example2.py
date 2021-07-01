# This is a collections of examples to demonstrate how to use odeSolu and plot the series solution

import odeSolu as os
from sympy import symbols,Function,diff


###### For plottings import these libraries ##############
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as nppoly
###########################################################


######## Example 1 (Bernoulli equation) ##################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,1)+(10**(-2))*y(x)-(1/10)*y(x)**3
T2=50 
T3=300
U=5

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),T3)

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=300')

plt.legend(loc='upper left')

plt.savefig('exa1.png')
plt.clf()
##########################################################


######## Example 2 (Riccati equation) ####################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,1)-(2*x**3+3)-(5*x**2+2*x+1)*y(x)-(x+7)*y(x)**2
T2=50 
T3=300
U=0.116

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),T3)

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=300')

plt.legend(loc='upper left')

plt.savefig('exa2.png')
plt.clf()
##########################################################


######## Example 3 (Abel's differential equation of the first kind) ####################
x=symbols('x')
y=Function('y')
ode=diff(y(x),x,1)+4+(5)*y(x)+(0.1*x)*y(x)**2+(0.2*x**2)*y(x)**3
T2=50
T3=300
U=0.44

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 1),),T3)

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=300')

plt.legend(loc='upper right')

plt.savefig('exa3.png')
plt.clf()
#######################################################################################


######## Example 4 (A second-order ODE with quartic nonlinearity) ####################
x=symbols('x')
y=Function('y')
ode=diff(y(x),x,2)+4+(10**(-1))*diff(y(x),x,1)+y(x)**4
T2=50 
T3=500
U=1.064

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 1)),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 1)),T3)

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=500')

plt.legend(loc='upper right')

plt.savefig('exa4.png')
plt.clf()
#######################################################################################


######## Example 5 (De Boer-Ludford equation) ####################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,2)+x**2*y(x)-y(x)**4

T2=50 
T3=500
U=1.5

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0)),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0)),T3)

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=500')

plt.legend(loc='upper left')

plt.savefig('exa5.png')
plt.clf()
##################################################################


######## Example 6 (Van der Pol equation) ####################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,2)-0.05*(1-y(x)**2)*diff(y(x),x,1)+y(x)
T2=50
T3=500
U=3.6

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5)),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5)),T3)

x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=500')

plt.legend(loc=(0.25,.6))

plt.savefig('exa6.png')
plt.clf()
##############################################################################


######## Example 7 (Painleve-Ince equation) ##################################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,2)+3*y(x)*diff(y(x),x,1)+y(x)**3
T2=50
T3=500
U=1.99

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5)),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5)),T3)

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=500')

plt.legend(loc='upper left')

plt.savefig('exa7.png')
plt.clf()
##############################################################################


######## Example 8 (Blasius equation) ########################################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,3)+y(x)*diff(y(x),x,2)
T2=50
T3=500
U=50

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 1),(diff(y(x), x, 2), 0)),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 1),(diff(y(x), x, 2), 0)),T3)

x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=500')

plt.legend(loc='upper left')

plt.savefig('exa8_1.png')
plt.clf()
##############################################################################


######## Example 9 (Falkner–Skan equation) ####################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,3)+y(x)*diff(y(x),x,2)+2*(1-diff(y(x),x,1)**2)
T2=50 #50,100,500
T3=500
U=2.3

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0.5),(diff(y(x), x, 2), 1)),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0.5),(diff(y(x), x, 2), 1)),T3)

x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=500')

plt.legend(loc='upper right')

plt.savefig('exa8.png')
plt.clf()
##############################################################################


######## Example 10 (fourth-order ODE) #####################################################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,4)-x**2*diff(y(x),x,3)+3*x*y(x)*diff(y(x),x,2)-6*(diff(y(x),x,1))+2*x**2+x
T2=100
T3=1000
U=2.585

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5),(diff(y(x), x, 2), 1),(diff(y(x), x, 3), 1)),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 0),(diff(y(x), x, 1), 0.5),(diff(y(x), x, 2), 1),(diff(y(x), x, 3), 1)),T3)

x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=100')
plt.plot(x3,y3,'orange',linestyle='--', label='N=1000')

plt.legend(loc='center')

plt.savefig('exa9.png')
plt.clf()
############################################################################################


######## Example 11 (5th-order) #####################################################
x=symbols('x')
y=Function('y')

ode=diff(y(x),x,5)-0.001*y(x)**2*diff(y(x),x,4)-2*x*y(x)*diff(y(x),x,3)**2+0.5*x*y(x)*diff(y(x),x,2)**4\
    -diff(y(x),x,1)+y(x)**3*x**2
T2=50
T3=500
U=1.6

solu2=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0),
                                  (diff(y(x), x, 2), 1),(diff(y(x), x, 3), 1),(diff(y(x), x, 4), 0.5)),T2)
solu3=os.seriesSolu(ode,(x,y(x)),((y(x), 1),(diff(y(x), x, 1), 0),
                                  (diff(y(x), x, 2), 1),(diff(y(x), x, 3), 1),(diff(y(x), x, 4), 0.5)),T3)

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("y")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='N=50')
plt.plot(x3,y3,'orange',linestyle='--', label='N=500')

plt.legend(loc='upper left')

plt.savefig('exa10.png')
plt.clf()
############################################################################################