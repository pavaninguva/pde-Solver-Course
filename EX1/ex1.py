# Necessary imports
from fipy import *
from fipy.viewers.matplotlibViewer.matplotlibViewer import _ColorBar
import sys


# Setting up the mesh

# We set up a unit square mesh of 100 by 100 cells

nx = ny = 100
dx = dy = 0.01
mesh = Grid2D (dx=dx, dy=dy, nx=nx, ny=ny)

phi = CellVariable(mesh=mesh, name=r"$\phi$", value= 3.2)

# Applying Dirichlet boundary conditions only
valueTop = 4.0
valueBottom = valueLeft = valueRight = 0.0

# Top BC
phi.constrain(valueTop, mesh.facesTop)

# Bottom BC
phi.constrain(valueBottom, mesh.facesBottom)

# Left BC
phi.constrain(valueLeft, mesh.facesLeft)

# Right BC
phi.constrain(valueRight, mesh.facesRight)

# Code to implement a Neumann BC
# phi.faceGrad.constrain(0.0, where=mesh.facesLeft)

# Defining the equation
phi.equation = (DiffusionTerm(coeff = 1.0) == 0)


# Setting up the analytical solution
pi = numerix.pi

x = mesh.cellCenters[0]
y = mesh.cellCenters[1]

analytical_solution = ( (4.0/(pi*numerix.sinh(pi)))*(numerix.sin(pi*x))*(numerix.sinh(pi*y)) + 
                    (4.0/(3.0*pi*numerix.sinh(3.0*pi)))*(numerix.sin(3.0*pi*x))*(numerix.sinh(3.0*pi*y)) +
                    (4.0/(5.0*pi*numerix.sinh(5.0*pi)))*(numerix.sin(5.0*pi*x))*(numerix.sinh(5.0*pi*y)) +
                    (4.0/(7.0*pi*numerix.sinh(7.0*pi)))*(numerix.sin(7.0*pi*x))*(numerix.sinh(7.0*pi*y)) +
                    (4.0/(9.0*pi*numerix.sinh(9.0*pi)))*(numerix.sin(9.0*pi*x))*(numerix.sinh(9.0*pi*y)) +
                    (4.0/(11.0*pi*numerix.sinh(11.0*pi)))*(numerix.sin(11.0*pi*x))*(numerix.sinh(11.0*pi*y)) +
                    (4.0/(13.0*pi*numerix.sinh(13.0*pi)))*(numerix.sin(13.0*pi*x))*(numerix.sinh(13.0*pi*y)) 
)

analytical = CellVariable(mesh=mesh, name=r"$\phi_{Analytical}$", value = analytical_solution )



phi.equation.solve(var=phi)

from fipy import input
if __name__ == '__main__':
    vi = Viewer(phi, colorbar=None)
    vi.colorbar = _ColorBar(viewer=vi, vmin=0.0, vmax=1.0)
    vi.plot()
    input("Press <return> to proceed")


# Setting up the second viewer for the analytical solution
if __name__ == '__main__':
    vi2 = Viewer(analytical, colorbar=None)
    vi2.colorbar = _ColorBar(viewer=vi, vmin=0.0, vmax=1.0)
    input("Press <return> to proceed")
