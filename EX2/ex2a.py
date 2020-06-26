# Necessary imports

from fipy import Variable, FaceVariable, CellVariable, Grid1D, ExplicitDiffusionTerm, TransientTerm, DiffusionTerm, Viewer
from fipy.tools import numerix


# Setting up a mesh
nx = 50
dx = 0.02
mesh = Grid1D(nx=nx, dx=dx)

# Defining the cell variable and initialising it
c = CellVariable(mesh= mesh, name=r"$c$", value = 0.0)

# Provide value for diffusion coefficient 
D = 1.0

# Setting boundary conditions
valueLeft = 1.0
valueRight = 0.0

c.constrain(valueLeft, mesh.facesLeft)
c.constrain(valueRight, mesh.facesRight)

# Specifying our alpha value
alpha = 0

# Defining the equation with alpha which gives us the flexibility to choose the time-stepping scheme
eq = (TransientTerm() == DiffusionTerm(coeff= alpha* D) + ExplicitDiffusionTerm(coeff = (1.0 - alpha)*D))


# dt = 0.00018 should be stable for forward Euler with the default simulation conditions
dt = 0.00018

# We might not want to see the output from every single timestep, so we define a stride parameter
time_stride = 5
timestep = 0
run_time = 1


# We need to implement the analytical solution. We can do two solutions: 1) Full solution, 2) S.S. solution
pi = numerix.pi
x = mesh.cellCenters[0]
t = Variable(0.0)

# This part of the code may be a little too much, dont worry about it, essentially we are recreating the analytical solution

# First we define a function that gives me an individual (nth) fourier mode: 
def one_fourier (n):
    mode = (-2.0/(n*pi))*numerix.sin(n*pi*x)*numerix.exp(-t*(pi**2)*(n**2))
    return mode

# Next we define a second function that gives me an truncated version of the analytical solution with k modes
def analytical_expression (k):
    Fourier = 0
    for mode in range(1, k+1):
        Fourier = Fourier + one_fourier(mode)
    Approx = 1.0 - x + Fourier
    return Approx

analytical_solution_transient = analytical_expression(10)
analytical_solution_transient.name = "Full Analytical Solution"

# analytical_solution_transient = 0.5 + x - 100*t

analytical_solution_steady = 1.0 - x

# analytical_transient = CellVariable(mesh=mesh, name="Full Analytical Solution", value = analytical_solution_transient, hasOld=1)


analytical_steady = CellVariable(mesh=mesh, name="Steady State Solution", value = analytical_solution_steady)


# Set up the time stepping

# First we got to initialise the viewer:
if __name__ == "__main__":
    viewer = Viewer(vars=(c, analytical_steady, analytical_solution_transient), datamin=0., datamax=1.)

while t < run_time:
    t.value = t.value + dt
    timestep += 1
    eq.solve(var=c, dt=dt)
    if (timestep % time_stride ==0):
        print ("Beep")
        if __name__ == '__main__':
            viewer.plot()
if __name__ == '__main__':
    input("Press <return> to proceed...")

