import sys
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
    plot(x,psi)
    print("Eigenvalue no %lf:%lf"%(i+1,x2))
show()
