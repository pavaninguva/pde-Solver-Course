"""
This problem is based of exercise 3a. 

We seek to evaluate how the different convection schemes perform.

The following schemes are available within FiPY: 
1. CentralDifferenceConvectionTerm
2. ExponentialConvectionTerm
3. HybridConvectionTerm
4. PowerLawConvectionTerm
5. UpwindConvectionTerm
6. ExplicitUpwindConvectionTerm
7. VanLeerConvectionTerm
"""

from fipy import *
from fipy.tools import numerix

# Setting up a mesh of domain size 4L
nx = 400
dx = 0.01
L = 1.0
mesh = PeriodicGrid1D(nx=nx, dx=dx)

#define parameters
D = 1.0
U = 100.0
peclet = (U*L/D)
convCoeff= (1.0,)
sigma = 0.05*L

# We define four variables to test four convection schemes
c1 = CellVariable(mesh=mesh, name=r"$c1$")
c2 = CellVariable(mesh=mesh, name=r"$c2$")
c3 = CellVariable(mesh=mesh, name=r"$c3$")
c4 = CellVariable(mesh=mesh, name=r"$c4$")

# Setting initial conditions
x = mesh.cellCenters[0]
c1.value=0.0
c1.setValue(0.05*(1.0/(sigma*numerix.sqrt(2*numerix.pi)))*numerix.exp(-0.5*((x-1.0)/(sigma))**2.0))
c2.value=0.0
c2.setValue(0.05*(1.0/(sigma*numerix.sqrt(2*numerix.pi)))*numerix.exp(-0.5*((x-1.0)/(sigma))**2.0))
c3.value=0.0
c3.setValue(0.05*(1.0/(sigma*numerix.sqrt(2*numerix.pi)))*numerix.exp(-0.5*((x-1.0)/(sigma))**2.0))
c4.value=0.0
c4.setValue(0.05*(1.0/(sigma*numerix.sqrt(2*numerix.pi)))*numerix.exp(-0.5*((x-1.0)/(sigma))**2.0))

# Defining equations
eq1 = TransientTerm(var=c1) + VanLeerConvectionTerm(coeff=convCoeff, var=c1) == 0
eq2 = TransientTerm(var=c2) + HybridConvectionTerm(coeff=convCoeff, var=c2) == 0
eq3 = TransientTerm(var=c3) + PowerLawConvectionTerm(coeff=convCoeff, var=c3) == 0
eq4 = TransientTerm(var=c4) + ExplicitUpwindConvectionTerm(coeff=convCoeff, var=c4) == 0


# Time stepping
dt = 0.001
time_stride = 10
timestep = 0
run_time = 2.0
t = timestep * dt

if __name__ == "__main__":
    viewer = Viewer(vars=(c1, c2, c3, c4), datamin=0., datamax=1.)

while t < run_time:
    t += dt
    timestep += 1
    eq1.solve(var=c1, dt=dt)
    eq2.solve(var=c2, dt=dt)
    eq3.solve(var=c3, dt=dt)
    eq4.solve(var=c4, dt=dt)
    if (timestep % time_stride ==0):
        print ("Beep")
        if __name__ == '__main__':
            viewer.plot()
if __name__ == '__main__':
    input("Press <return> to proceed...")

