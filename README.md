# Runge-Kutta-4-Demonstration
 The repository demobstrates the application of the Runge Kutta 4 method to solve the Schrodinger equation for the Quantum Harmonic Oscillator potential.
 The following is the plan to solve the problem.
1. The concerned differential equation is the time independent Schrodinger equation: H*psi = E*psi,
    where H = Cd2/dx2 where C is a constant derived from physical constants, E-Energy eigenstates and psi is the Wavefunction of the system.
The idea is to solve this ODE using Runge Kutta 4th order by parameterising dpsi/dx = y giving two coupled differential eqns
                dpsi/dx = y :  dy/dx = E*psi
2. After Solving it, comparing it with Analytical solution. If successfuly completed, extend it to three dimension.
3. After that generalise it by finding the time dependent solution by calculating the Cn's.
4. Animating the wavefunction psi. 
