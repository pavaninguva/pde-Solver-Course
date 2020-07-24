
# Imports
from fipy import Variable, FaceVariable, CellVariable, Grid1D, ExplicitDiffusionTerm, TransientTerm, DiffusionTerm, Viewer, ImplicitSourceTerm, TSVViewer
from fipy.tools import numerix
import pandas as pd
import glob as glob
import re
import os

# Setting up a mesh
nx = 50
dx = 0.02
L = nx*dx
mesh = Grid1D(nx=nx, dx=dx)

# Defining the cell variable
c = CellVariable(mesh= mesh, name=r"$c$")

# Initialise c with a funky profile
x = mesh.cellCenters[0]
c.value=0.0
c.setValue((0.3*numerix.sin(2.0*x*numerix.pi)) + 0.5)

# Setting up the no-flux boundary conditions
c.faceGrad.constrain(0.0, where=mesh.facesLeft)
c.faceGrad.constrain(0.0, where=mesh.facesRight)


# Provide value for diffusion and reaction coefficient 
D = 1.0
k = -2.0

# Specifying our alpha value. We will only use backwards Euler going forward
alpha = 1

# Defining the equation with the reaction term. 
eq = TransientTerm() == DiffusionTerm(coeff=  D)  + ImplicitSourceTerm(coeff=k)


# We can use a much larger dt now since we are using an implicit method. 
dt = 0.001

# We might not want to see the output from every single timestep, so we define a stride parameter
time_stride = 100
timestep = 0
run_time = 0.5
t = timestep * dt

# We first initialise the default viewer to get things started:

# if __name__ == "__main__":
#     viewer = Viewer(vars=(c), datamin=0., datamax=1.)

# while t < run_time:
#     t += dt
#     timestep += 1
#     eq.solve(var=c, dt=dt)
#     if (timestep % time_stride ==0):
#         print ("Beep")
#         if __name__ == '__main__':
#             viewer.plot()
# if __name__ == '__main__':
#     input("Press <return> to proceed...")


if __name__ == "__main__":
    vw = TSVViewer(vars=(c))

while t < run_time:
    if (timestep == 0):
        print ("Beep Initial Profile Output")
        if __name__ == '__main__':
            vw.plot(filename="ex2c_%s.tsv"%(t))
    t += dt
    timestep += 1
    eq.solve(var=c, dt=dt)
    if (timestep % time_stride ==0):
        print ("Beep")
        if __name__ == '__main__':
            vw.plot(filename="ex2c_%s.tsv"%(t))
if __name__ == '__main__':
    input("Press <return> to proceed...")

# The next section of code is to convert the tsv files into csv which is more useable:

# We first generate a list of the tsv files in the current folder
tsv_list = glob.glob("./*.tsv")

# Next we use a for loop to go through the list
for tsv in tsv_list:
    csv_filename = tsv.replace("tsv", "csv")

    # This is a hack to remove the first line of the tsv file 
    with open(tsv, "r") as fin:
        data = fin.read().splitlines(True)
    with open(tsv, "w") as fout:
        fout.writelines(data[1:])

    # We read the tsv file into a dataframe
    csv_table = pd.read_table(tsv, sep="\t")
    csv_table.to_csv(csv_filename)

    # Remove the tsv files
    os.remove(tsv)
