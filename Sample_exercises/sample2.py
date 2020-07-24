"""
This problem is based of exercise 2. 

We consider diffusion in a 1D slab  

The diffusion coefficient in the first half of the slab is double that in the second half
"""
from fipy import Variable, FaceVariable, CellVariable, Grid1D, TransientTerm, DiffusionTerm, Viewer
from fipy.tools import numerix

# Setting up a mesh
nx = 50
dx = 0.02
L = nx * dx
mesh = Grid1D(nx=nx, dx=dx)
c = CellVariable(mesh= mesh, name=r"$c$", value = 0.0)

# Setting up the spatially varying diffusion coefficient
# Due to the mechanics of the solver, the diffusion coefficient variable must be defined on the mesh face
D = FaceVariable(mesh=mesh, value=2.0)
x = mesh.faceCenters[0]
D.setValue(1.0, where=(x > L/2.0))

# Boundary conditions
valueLeft = 1.0
valueRight = 0.0
c.constrain(valueLeft, mesh.facesLeft)
c.constrain(valueRight, mesh.facesRight)

#Equation
eq = (TransientTerm() == DiffusionTerm(coeff= D))

dt = 0.001
t = Variable(0.0)
time_stride = 20
timestep = 0
run_time = 1

if __name__ == "__main__":
    viewer = Viewer(vars=(c), datamin=0., datamax=1.)

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




