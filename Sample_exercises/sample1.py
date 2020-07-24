"""
This sample exercise is based of exercise 1. 

We explore how different boundary conditions affect the solution namely: 
1. If more than one vertex has a Dirichlet BC of \phi = 1.0
2. If one or more vertex has a funky BC e.g. a sine function
3. Combination of Dirichlet and Neumann BC

A helpful way of thinking about the results for this simulation is to visualise the steady state solution for thermal diffusion.
"""

# Necessary imports
from fipy import *
from fipy.viewers.matplotlibViewer.matplotlibViewer import _ColorBar

# Setting up the mesh
nx = ny = 100
dx = dy = 0.01
mesh = Grid2D (dx=dx, dy=dy, nx=nx, ny=ny)

phi = CellVariable(mesh=mesh, name=r"$\phi$", value= 0.0)

#### BCs for the first case
# valueTop = valueLeft =  1.0
# valueBottom = valueRight = 0.0

# phi.constrain(valueTop, mesh.facesTop)
# phi.constrain(valueBottom, mesh.facesBottom)
# phi.constrain(valueLeft, mesh.facesLeft)
# phi.constrain(valueRight, mesh.facesRight)

#### BCs for the second case
# x = mesh.faceCenters[0]
# valueTop = numerix.sin(x*2.0*numerix.pi)
# valueBottom = -numerix.sin(x*2.0*numerix.pi)
# valueLeft = valueRight = 0.0

# phi.constrain(valueTop, mesh.facesTop)
# phi.constrain(valueBottom, mesh.facesBottom)
# phi.constrain(valueLeft, mesh.facesLeft)
# phi.constrain(valueRight, mesh.facesRight)

#### BCs for the third case. 
valueTop = 1.0
valueLeft = 0.0
valueRight = 0.0
valueBottom = 0.0

# phi.faceGrad.constrain(valueTop, where=mesh.facesTop)
phi.faceGrad.constrain(valueLeft, where=mesh.facesLeft)
phi.faceGrad.constrain(valueRight, where=mesh.facesRight)
# phi.faceGrad.constrain(valueBottom, where=mesh.facesBottom)

phi.constrain(valueTop, mesh.facesTop)
# phi.constrain(valueLeft, mesh.facesLeft)
# phi.constrain(valueRight, mesh.facesRight)
phi.constrain(valueBottom, mesh.facesBottom)

# Defining the equation
phi.equation = (DiffusionTerm(coeff = 1.0) == 0)

# Solving the equation
phi.equation.solve(var=phi)

from fipy import input
if __name__ == '__main__':
    vi = Viewer(phi, colorbar=None)
    # vi.colorbar = _ColorBar(viewer=vi, vmin=0.0, vmax=1.0)
    vi.plot()
    input("Press <return> to proceed")