# This is a collection of examples to demonstrate how to use odeSolu to plot the series solution
# with convergence test.

import odeSolu as os
from sympy import symbols,Function,diff


###### For plottings import these libraries ##############
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as nppoly
###########################################################


######## Example 1 (Bernoulli equation) ##################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,1)+(10**(-2))*u(x)-(1/10)*u(x)**3
T2=50 
T3=300
U=5

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T3,(0,U))

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=300')

plt.legend(loc='upper left')

plt.savefig('Fig1.png')
plt.clf()
##########################################################


######## Example 2 (Riccati equation) ####################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,1)-(2*x**3+3)-(5*x**2+2*x+1)*u(x)-(x+7)*u(x)**2
T2=50 
T3=300
U=0.1

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T3,(0,U))

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=300')

plt.legend(loc='upper left')

plt.savefig('Fig2.png')
plt.clf()
##########################################################


######## Example 3 (Abel's differential equation of the first kind) ####################
x=symbols('x')
u=Function('u')
ode=diff(u(x),x,1)+4+(5)*u(x)+(0.1*x)*u(x)**2+(0.2*x**2)*u(x)**3
T2=50
T3=300
U=0.42

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T3,(0,U))

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=300')

plt.legend(loc='upper right')

plt.savefig('Fig3.png')
plt.clf()
#######################################################################################


######## Example 4 (A second-order ODE with quartic nonlinearity) ####################
x=symbols('x')
u=Function('u')
ode=diff(u(x),x,2)+4+(10**(-1))*diff(u(x),x,1)+u(x)**4
T2=50 
T3=500
U=1.0

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 1)),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 1)),T3,(0,U))

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.legend(loc='upper right')

plt.savefig('Fig4.png')
plt.clf()
#######################################################################################


######## Example 5 (De Boer-Ludford equation) ####################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,2)+x**2*u(x)-u(x)**4

T2=50 
T3=500
U=1.36

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0)),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0)),T3,(0,U))

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.legend(loc='upper left')

plt.savefig('Fig5.png')
plt.clf()
##################################################################


######## Example 6 (Van der Pol equation) ####################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,2)-0.05*(1-u(x)**2)*diff(u(x),x,1)+u(x)
T2=50
T3=500
U=3.55

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5)),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5)),T3,(0,U))

x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.legend(loc=(0.25,.6))

plt.savefig('Fig6.png')
plt.clf()
##############################################################################


######## Example 7 (Painleve-Ince equation) ##################################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,2)+3*u(x)*diff(u(x),x,1)+u(x)**3
T2=50
T3=500
U=1.92

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5)),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5)),T3,(0,U))

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.legend(loc='upper left')

plt.savefig('Fig7.png')
plt.clf()
##############################################################################


######## Example 8 (Blasius equation) ########################################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,3)+u(x)*diff(u(x),x,2)
T2=50
T3=500
U=50

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 1),(diff(u(x), x, 2), 0)),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 1),(diff(u(x), x, 2), 0)),T3,(0,U))

x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.legend(loc='upper left')

plt.savefig('Fig8.png')
plt.clf()
##############################################################################


######## Example 9 (Falkner–Skan equation) ####################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,3)+u(x)*diff(u(x),x,2)+2*(1-diff(u(x),x,1)**2)
T2=50 #50,100,500
T3=500
U=2.25

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0.5),(diff(u(x), x, 2), 1)),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0.5),(diff(u(x), x, 2), 1)),T3,(0,U))

x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.legend(loc='upper right')

plt.savefig('Fig9.png')
plt.clf()
##############################################################################


######## Example 10 (fourth-order ODE) #####################################################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,4)-x**2*diff(u(x),x,3)+3*x*u(x)*diff(u(x),x,2)-6*(diff(u(x),x,1))+2*x**2+x
T2=100
T3=1000
U=2.0

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5),(diff(u(x), x, 2), 1),(diff(u(x), x, 3), 1)),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5),(diff(u(x), x, 2), 1),(diff(u(x), x, 3), 1)),T3,(0,U))

x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=100')
plt.plot(x3,y3,'orange',linestyle='--', label='n=1000')

plt.legend(loc='center')

plt.savefig('Fig10.png')
plt.clf()
############################################################################################


######## Example 11 (5th-order) #####################################################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,5)-0.001*u(x)**2*diff(u(x),x,4)-2*x*u(x)*diff(u(x),x,3)**2+0.5*x*u(x)*diff(u(x),x,2)**4\
    -diff(u(x),x,1)+u(x)**3*x**2
T2=50
T3=500
U=1.4

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0),
                                  (diff(u(x), x, 2), 1),(diff(u(x), x, 3), 1),(diff(u(x), x, 4), 0.5)),T2,(0,U))
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0),
                                  (diff(u(x), x, 2), 1),(diff(u(x), x, 3), 1),(diff(u(x), x, 4), 0.5)),T3,(0,U))

x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2[0])
y3=nppoly.polyval(x3,solu3[0])

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.legend(loc='upper left')

plt.savefig('Fig11.png')
plt.clf()
############################################################################################