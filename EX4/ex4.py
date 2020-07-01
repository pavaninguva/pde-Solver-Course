
from fipy import CellVariable, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, VTKCellViewer, LinearLUSolver
from fipy.tools import numerix
import time
from fipy import PeriodicGrid2D, Viewer, UniformNoiseVariable

# Simulation parameters
# Mesh resolution and domain size
nx = ny = 50
dx = dy = 1.0
# Initial composition
a_0 = 0.77
# Noise magnitude
noise_mag = 0.03
# Flory-Huggins interaction parameter
chi_AB = 0.006
# Length of polymer chains
n_a = 1000
n_b = 1000

# Defining the mesh
mesh = PeriodicGrid2D(nx=nx, ny=ny, dx=dx, dy=dy)

# Defining the variables
a = CellVariable (name=r"$a$", mesh=mesh, hasOld=1)
mu_AB = CellVariable(name=r"$\mu_{AB}", mesh=mesh, hasOld=1)

# Setting the initial composition of the system with noise
# noise = GaussianNoiseVariable(mesh=mesh, mean=a_0, variance=noise_mag).value
noise = UniformNoiseVariable(mesh=mesh, minimum=(a_0-noise_mag), maximum=(a_0+noise_mag))
a[:] = noise

# differentiate g(a)
dgda = ((1.0/n_a) - (1.0/n_b)) + (1.0/n_a)*numerix.log(a) - (1.0/n_b)*numerix.log(1.0 - a) + chi_AB*(1.0 - 2*a)
d2gda2 = (1.0/(n_a*a)) + (1.0/(n_b*(1.0 - a))) - 2*chi_AB

# Evaluate kappa
kappa = (2.0/3.0)*chi_AB

# Defining the equations
eq1 = (TransientTerm(var=a)) == DiffusionTerm(coeff=a*(1.0-a), var=mu_AB)
eq2 = (ImplicitSourceTerm(coeff=1.0, var=mu_AB)) == dgda - DiffusionTerm(coeff=kappa, var=a)
# eq2 = (ImplicitSourceTerm(coeff=1.0, var=mu_AB)) == ImplicitSourceTerm(coeff=d2gda2, var=a) - d2gda2*a + dgda - DiffusionTerm(coeff=kappa, var=a)

# Coupling the equations
eq = eq1 & eq2

# Setting up the solver
solver = LinearLUSolver(tolerance=1e-9, iterations=50, precon="ilu")

# Set up time stepping
dt = 10.0
duration = 20000
time_stride = 100
timestep = 0
elapsed =0

# Initialising the viewer
if __name__ == "__main__":
    viewer = Viewer(vars=(a), datamin=0., datamax=1.)

# start the time
start = time.time()

# Time stepping
while elapsed < duration: 
    elapsed += dt
    timestep += 1
    a.updateOld()
    mu_AB.updateOld()
    res = 1e4
    while res > 1e-10:
        res = eq.sweep(dt=dt, solver=solver)
        print ("sweep")
        print (res)
    print (elapsed)
    end = time.time()
    print (end-start)
    if (timestep % time_stride ==0):
        print ("Beep")
        if __name__ == '__main__':
            viewer.plot()
if __name__ == '__main__':
    input("Press <return> to proceed...")

# Code for VTK ouput
# while elapsed < duration: 
#     if (timestep == 0):
#         vw = VTKCellViewer(vars=(a, mu_AB))
#         vw.plot(filename="0_output.vtk")
#     elapsed += dt
#     timestep += 1
#     a.updateOld()
#     mu_AB.updateOld()
#     res = 1e+10
#     while res > 1e-10:
#         res = eq.sweep(dt=dt, solver=solver)
#         print ("sweep!")
#     print (elapsed)
#     end = time.time()
#     print(end-start)
#     if (timestep % time_stride ==0):
#         vw = VTKCellViewer(vars=(a, mu_AB))
#         vw.plot(filename="%s_output.vtk" %(elapsed))


    