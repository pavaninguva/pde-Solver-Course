from fipy import *
from fipy.tools import numerix
import pandas as pd
import glob as glob
import re
import os

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
# This is the syntax needed to specify u_x
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
eq = TransientTerm() + VanLeerConvectionTerm(coeff=convCoeff) == 0
# eq = TransientTerm() + ExponentialConvectionTerm(coeff=convCoeff) == 0
# eq = TransientTerm() - DiffusionTerm(coeff=(1/peclet)) == 0
# eq = TransientTerm() + VanLeerConvectionTerm(coeff=convCoeff) - DiffusionTerm(coeff=(1.0/peclet)) == 0


# The choice of dt and dx in a convection problem is a bit more complicated 
dt = 0.0001

# We might not want to see the output from every single timestep, so we define a stride parameter
time_stride = 200
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

# if __name__ == "__main__":
#     vw = TSVViewer(vars=(c))

# while t < run_time:
#     if (timestep == 0):
#         print ("Beep Initial Profile Output")
#         if __name__ == '__main__':
#             vw.plot(filename="ex3a_initial_%s.tsv"%(t))
#     t += dt
#     timestep += 1
#     eq.solve(var=c, dt=dt)
#     if (timestep % time_stride ==0):
#         print ("Beep")
#         if __name__ == '__main__':
#             vw.plot(filename="ex3a_advection.tsv"%(t))
# if __name__ == '__main__':
#     input("Press <return> to proceed...")

# # The next section of code is to convert the tsv files into csv which is more useable:

# # We first generate a list of the tsv files in the current folder
# tsv_list = glob.glob("./*.tsv")

# # Next we use a for loop to go through the list
# for tsv in tsv_list:
#     csv_filename = tsv.replace("tsv", "csv")

#     # This is a hack to remove the first line of the tsv file 
#     with open(tsv, "r") as fin:
#         data = fin.read().splitlines(True)
#     with open(tsv, "w") as fout:
#         fout.writelines(data[1:])

#     # We read the tsv file into a dataframe
#     csv_table = pd.read_table(tsv, sep="\t")
#     csv_table.to_csv(csv_filename)

#     # Remove the tsv files
#     os.remove(tsv)
