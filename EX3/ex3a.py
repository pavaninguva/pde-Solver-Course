from fipy import *
from fipy.tools import numerix

# Setting up a mesh of domain size 4L
nx = 400
dx = 0.01
L = 1.0

# The use of FiPy's PeriodicGrid1D / 2D gives us a mesh with the necessary BC implemented
mesh = PeriodicGrid1D(nx=nx, dx=dx)

#define parameters
D = 1.0
U = 100.0
peclet = (U*L/D)
# This is the syntax needed to specify U
convCoeff= (1.0,)

# intial condition width (standard deviation)
sigma = 0.05*L

#Defining the variable
c = CellVariable(mesh=mesh, name=r"$c$")

# Setting initial conditions
x = mesh.cellCenters[0]
c.value=0.0
c.setValue(0.05*(1.0/(sigma*numerix.sqrt(2*numerix.pi)))*numerix.exp(-0.5*((x-1.0)/(sigma))**2.0))

# defining the equation
# We write three equations: pure convection, pure diffusion and convection-diffusion
# eq = TransientTerm() + VanLeerConvectionTerm(coeff=convCoeff) == 0
# eq = TransientTerm() - DiffusionTerm(coeff=(1/peclet)) == 0
eq = TransientTerm() + VanLeerConvectionTerm(coeff=convCoeff) - DiffusionTerm(coeff=(1.0/peclet)) == 0


# The choice of dt and dx in a convection problem is a bit more complicated 
dt = 0.001

# We might not want to see the output from every single timestep, so we define a stride parameter
time_stride = 10
timestep = 0
run_time = 2.0
t = timestep * dt


if __name__ == "__main__":
    viewer = Viewer(vars=(c), datamin=0., datamax=1.)

while t < run_time:
    t += dt
    timestep += 1
    eq.solve(var=c, dt=dt)
    if (timestep % time_stride ==0):
        print ("Beep")
        if __name__ == '__main__':
            viewer.plot()
if __name__ == '__main__':
    input("Press <return> to proceed...")
