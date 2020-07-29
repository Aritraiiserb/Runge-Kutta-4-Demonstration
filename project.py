# This is the project that might change my life. Good luck!!
''' I am interested in solving a system of quantum of Harmonic Oscillator. I should be solving the Schrodinger
Equation. This is the following plan to solve the problem.
1. The concerned differential equation is the Time independent Schrodinger equation: H*psi = E*psi,
    where H = Cd2/dx2 where C is a constant derived from physical constants, E-Energy eigenstates and psi is the Wavefunction of the system.
The idea is to solve this ODE using Runge Kutta 4th order by parameterising dpsi/dx = y giving two coupled differential eqns
                dpsi/dx = y :  dy/dx = E*psi
2. After Solving it, comparing it with Analytical solution. If successfuly completed, extend it to three dimension.
3. After that generalise it by finding the time dependent solution by calculating the Cn's.
4. Animating the wavefunction psi. '''
import sys
#import vpython as vp
from math import *
sys.path.append('/home/chakri/python/PrOjEcT/Modules')
from RK4CD_module import * # Module for solving the coupled differential equation
psi0 = 0.0
phi0 = 1.0
def f2(x,phi,psi):
    return phi
def f1(x,phi,psi):
    global E
    return (x**2-E)*psi
def solver(e):
    global E
    E = e
    A = 0.0
    psi,phi,x = ODE_runge_4th_CD(f1,f2,phi0,psi0,-5,5)
    for i in range(len(psi)):
        A += psi[i]**2
    for i in range(len(psi)):
        psi[i] = psi[i]/sqrt(1.0*A)
    return psi[-1]
def solver_fullfunc(e):
    global E
    E = e
    A = 0.0
    psi,phi,x = ODE_runge_4th_CD(f1,f2,phi0,psi0,-5,5)
    for i in range(len(psi)):
        A += psi[i]**2
    for i in range(len(psi)):
        psi[i] = psi[i]/sqrt(1.0*A)
    return psi,x
E_min = 0.0
E_max = 15.0
dE = 0.1
Edash = linspace(E_min,E_max,(E_max-E_min)/dE+1)
boundary = []
for e in Edash:
    boundary.append(solver(e))
plot(Edash,boundary)
ylabel(r'$\psi$(a)')
xlabel('E\'')
show()
#----------------------------------------------
e1 = []
e2 = []
wf1 = []
wf2 = []
PSIS = []
EIGEN = []
PSI_Im = []
PSI_Re = []
#gr = vp.graph(xmin = -6.0, xmax = 6.0, ymin = -0.1, ymax = 0.1)
#psireal = vp.gcurve(color = vp.color.red)
#psiimag = vp.gcurve(color = vp.color.green)
for i in range(len(Edash)-1):
    if boundary[i]*boundary[i+1]<=0:
        e1.append(Edash[i]);e2.append(Edash[i+1])
        wf1.append(boundary[i]);wf2.append(boundary[i+1])
for i in range(len(e1)):
    x0 = e1[i]
    x1 = e2[i]
    dx = x1 - x0
    while(abs(dx)>0.00000001):
        d = solver(x1) - solver(x0)
        x2 = x1 - (solver(x1)*dx)/d
        x0 = x1
        x1 = x2
        dx = x1 - x0
        psi,x = solver_fullfunc(x2)
    PSIS.append(psi)
    plot(x,psi)
    EIGEN.append(x2)
    print("Eigenvalue no %lf:%lf"%(i+1,x2))
show()
t_max = 10.0
t = 0.0
dt = 0.01
while t <= t_max:
    #vp.rate(30)
    append1 = []
    append2 = []
    for i in range(len(PSIS[0])):
        append1.append(-1*PSIS[0][i]*sin(t*EIGEN[0]*0.5))
        append2.append(PSIS[0][i]*cos(t*EIGEN[0]*0.5))
    PSI_Im.append(append1)
    PSI_Re.append(append2)
    #psireal.plot(pos = (x,PSI_Re[j]))
    t += dt
t=0.0
#for j in range(len(PSI_Im)):
    #gr = vp.graph()
    #psiimag.delete()
    #psireal.delete()
for i in range(len(PSI_Im[0])):
    plot(x,PSI_Im[i])
    plot(x,PSI_Re[i])
    #vp.rate(1500)
show()
