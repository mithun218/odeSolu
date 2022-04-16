# This is a collections of examples that have been used in the papers to compare the series solutions
# with numerical solutions.

import odeSolu as os
from sympy import symbols,Function,diff


###### For plottings import these libraries ##############
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as nppoly
###########################################################

### The numerical solutions of the ODEs obtained using this SciPy's module ########
from scipy.integrate import odeint
##############################################################################


######## Example 1 (Bernoulli equation) ##################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,1)+(10**(-2))*u(x)-(1/10)*u(x)**3
T3=300
U=5

solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T3, (0,5))

x3=np.linspace(0,U,T3)


plt.xlabel("x")
plt.ylabel("u")

y3=nppoly.polyval(x3,solu3[0])

plt.plot(x3,y3,'orange', label='n=300')


plt.legend(loc='upper left')

plt.savefig('exa1.pdf')
plt.clf()
##########################################################



######## Example 2 (Abel's differential equation of the first kind) ####################
x=symbols('x')
u=Function('u')
ode=diff(u(x),x,1)+4+(5)*u(x)+(0.1*x)*u(x)**2+(0.2*x**2)*u(x)**3
T2=50
T3=300
U=0.44

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T2)
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),),T3)

def odeintFunc(u,x):
    return -(4+(5)*u+(0.1*x)*u**2+(0.2*x**2)*u**3)

x1=np.linspace(0,U,500)
x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)


odeintSol=odeint(odeintFunc,[1],x1)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=300')

plt.plot(x1,odeintSol[:,0],'b', alpha=0.4, label='N.S.')

plt.legend(loc='upper right')

plt.savefig('exa2.pdf')
plt.clf()
#######################################################################################


######## Example 3 (A second-order ODE with quartic nonlinearity) ####################
x=symbols('x')
u=Function('u')
ode=diff(u(x),x,2)+4+(10**(-1))*diff(u(x),x,1)+u(x)**4
T2=50 
T3=500
U=1.064

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 1)),T2)
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 1)),T3)

def odeintFunc(u,x):   
    return u[1],-(4+(10**(-1))*u[1]+u[0]**4)

x1=np.linspace(0,U,500)
x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)


odeintSol=odeint(odeintFunc,[0,1],x1)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.plot(x1,odeintSol[:,0],'b', alpha=0.4, label='N.S.')

plt.legend(loc='upper right')

plt.savefig('exa3.pdf')
plt.clf()
#######################################################################################


######## Example 4 (De Boer-Ludford equation) ####################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,2)+x**2*u(x)-u(x)**4

T2=50 
T3=500
U=1.5

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0)),T2)
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0)),T3)

def odeintFunc(u,x):
    return [u[1],-((x)**2*u[0]-u[0]**4)]

x1=np.linspace(0,U,500)
x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)

odeintSol=odeint(odeintFunc,[1,0],x1)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.plot(x1,odeintSol[:,0],'b', alpha=0.4, label='N.S.')

plt.legend(loc='upper left')

plt.savefig('exa4.pdf')
plt.clf()
##################################################################


######## Example 5 (Van der Pol equation) ####################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,2)-0.05*(1-u(x)**2)*diff(u(x),x,1)+u(x)
T2=50
T3=500
U=3.6

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5)),T2)
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5)),T3)

def odeintFunc(u,x):
    return [u[1],-(-0.05*(1-u[0]**2)*u[1]+u[0])]

x11=np.linspace(0,-U,500)
x12=np.linspace(0,U,500)
x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)


odeintSol1=odeint(odeintFunc,[0,0.5],x11)
odeintSol2=odeint(odeintFunc,[0,0.5],x12)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.plot(x11,odeintSol1[:,0],'b', alpha=0.4, label='N.S.')
plt.plot(x12,odeintSol2[:,0],'b', alpha=0.4)

plt.legend(loc=(0.25,.6))

plt.savefig('exa5.pdf')
plt.clf()
##############################################################################


######## Example 6 (Painleve-Ince equation) ##################################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,2)+3*u(x)*diff(u(x),x,1)+u(x)**3
T2=50
T3=500
U=1.99

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5)),T2)
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5)),T3)

def odeintFunc(u,x):
    return [u[1],-(3*u[0]*u[1]+u[0]**3)]

x1=np.linspace(0,U,500)
x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)


odeintSol=odeint(odeintFunc,[0,0.5],x1)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.plot(x1,odeintSol[:,0],'b', alpha=0.4, label='N.S.')

plt.legend(loc='upper left')

plt.savefig('exa6.pdf')
plt.clf()
##############################################################################


######## Example 7 (Falkner–Skan equation) ####################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,3)+u(x)*diff(u(x),x,2)+2*(1-diff(u(x),x,1)**2)
T2=50 #50,100,500
T3=500
U=2.3

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0.5),(diff(u(x), x, 2), 1)),T2)
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0.5),(diff(u(x), x, 2), 1)),T3)

def odeintFunc(u,x):
    return [u[1],u[2],-(u[0]*u[2]+2*(1-u[1]**2))]

x11=np.linspace(0,-U,500)
x12=np.linspace(0,U,500)
x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)


odeintSol1=odeint(odeintFunc,[1,0.5,1],x11)
odeintSol2=odeint(odeintFunc,[1,0.5,1],x12)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.plot(x11,odeintSol1[:,0],'b', alpha=0.4, label='N.S.')
plt.plot(x12,odeintSol2[:,0],'b', alpha=0.4)

plt.legend(loc='upper right')

plt.savefig('exa7.pdf')
plt.clf()
##############################################################################


######## Example 8 (fourth-order ODE) #####################################################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,4)-x**2*diff(u(x),x,3)+3*x*u(x)*diff(u(x),x,2)-6*(diff(u(x),x,1))+2*x**2+x
T2=100
T3=1000
U=2.585

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5),(diff(u(x), x, 2), 1),(diff(u(x), x, 3), 1)),T2)
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 0),(diff(u(x), x, 1), 0.5),(diff(u(x), x, 2), 1),(diff(u(x), x, 3), 1)),T3)

def odeintFunc(u,x):
    return [u[1],u[2],u[3],-(-x**2*u[3]+3*x*u[0]*u[2]-6*(u[1])+2*x**2+x)]

x11=np.linspace(0,-U,500)
x12=np.linspace(0,U,500)
x2=np.linspace(-U,U,T2)
x3=np.linspace(-U,U,T3)


odeintSol1=odeint(odeintFunc,[0,0.5,1,1],x11)
odeintSol2=odeint(odeintFunc,[0,0.5,1,1],x12)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='n=100')
plt.plot(x3,y3,'orange',linestyle='--', label='n=1000')

plt.plot(x11,odeintSol1[:,0],'b', alpha=0.4, label='N.S.')
plt.plot(x12,odeintSol2[:,0],'b', alpha=0.4)

plt.legend(loc='center')

plt.savefig('exa8.pdf')
plt.clf()
############################################################################################


######## Example 9 (5th-order) #####################################################
x=symbols('x')
u=Function('u')

ode=diff(u(x),x,5)-0.001*u(x)**2*diff(u(x),x,4)-2*x*u(x)*diff(u(x),x,3)**2+0.5*x*u(x)*diff(u(x),x,2)**4\
    -diff(u(x),x,1)+u(x)**3*x**2
T2=50
T3=500
U=1.6

solu2=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0),
                                  (diff(u(x), x, 2), 1),(diff(u(x), x, 3), 1),(diff(u(x), x, 4), 0.5)),T2)
solu3=os.seriesSolu(ode,(x,u(x)),((u(x), 1),(diff(u(x), x, 1), 0),
                                  (diff(u(x), x, 2), 1),(diff(u(x), x, 3), 1),(diff(u(x), x, 4), 0.5)),T3)

def odeintFunc(u,x):
    return [u[1],u[2],u[3],u[4],
            -0.001*u[0]**2*u[4]-2*x*u[0]*u[3]**2+0.5*x*u[0]*u[2]**4-u[1]+u[0]**3*x**2]

x1=np.linspace(0,U,600)
x2=np.linspace(0,U,T2)
x3=np.linspace(0,U,T3)


odeintSol=odeint(odeintFunc,[1,0,1,1,0.5],x1)

plt.xlabel("x")
plt.ylabel("u")

y2=nppoly.polyval(x2,solu2)
y3=nppoly.polyval(x3,solu3)

plt.plot(x2,y2,'r',linestyle='--', label='n=50')
plt.plot(x3,y3,'orange',linestyle='--', label='n=500')

plt.plot(x1,odeintSol[:,0],'b', alpha=0.4, label='N.S.')

plt.legend(loc='upper left')

plt.savefig('exa9.pdf')
plt.clf()
############################################################################################