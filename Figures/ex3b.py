from fipy import *
from fipy.tools import numerix
from matplotlib import rc
import matplotlib.pyplot as plt

"""
Settings for matplotlib
"""
rc('font', family='serif')
rc('xtick', labelsize='xx-large')
rc('ytick', labelsize='xx-large')
rc("axes", titlesize="xx-large")
rc("axes", labelsize="xx-large")
rc("legend", fontsize="xx-large")
 
plt.rcParams["figure.figsize"] = (10,10)

# Setting up a mesh of domain size L = 4
nx = 400
dx = 0.01

mesh = Grid1D(nx=nx, dx=dx)

# Define the variables and boundary conditions: 
a = CellVariable(mesh=mesh, name=r"$a$", value = 0.0, hasOld=True)
b = CellVariable(mesh=mesh, name=r"$b$", value = 0.0, hasOld=True)
c = CellVariable(mesh=mesh, name=r"$c$", value = 0.0, hasOld=True)

# Inlet values for a, b, c
a_in = 1.0
b_in = 0.5
c_in = 0.0 

# Boundary conditions for inlet: 
a.constrain(a_in, mesh.facesLeft)
b.constrain(b_in, mesh.facesLeft)
c.constrain(c_in, mesh.facesLeft)

# Boundary conditions for the outlet
a.faceGrad.constrain([0], mesh.facesRight)
b.faceGrad.constrain([0], mesh.facesRight)
c.faceGrad.constrain([0], mesh.facesRight)

# Define the parameters: 
D = 1.0
k = 30.0
convCoeff = (3.0,)


# Defining the equations: 
eq_a = (TransientTerm(var=a)) + VanLeerConvectionTerm(coeff=convCoeff, var=a) - DiffusionTerm(coeff=D, var =a) + ImplicitSourceTerm(coeff = k*b, var=a) == 0
eq_b = (TransientTerm(var=b)) + VanLeerConvectionTerm(coeff=convCoeff, var=b) - DiffusionTerm(coeff=D, var =b) + ImplicitSourceTerm(coeff = k*b, var=a) == 0

eq_c = (TransientTerm(var=c)) + VanLeerConvectionTerm(coeff=convCoeff, var=c) - DiffusionTerm(coeff=D, var =c) + ImplicitSourceTerm(coeff = -k*b, var=a) == 0


# Setting up timestepping
dt = 0.01
time_stride = 10
timestep = 0
run_time = 5.0
t = timestep * dt

# coupling equations eq_a and eq_b 
eqn_main = eq_a & eq_b 

# Initialise the viewer

if __name__ == "__main__":
    viewer = Viewer((a, b, c),datamin=0., datamax=1.,  legend="upper right")
    viewer.axes.set_ylabel("Concentration")
    viewer.axes.set_xlabel("Distance")


while t < run_time:
    t += dt
    timestep += 1
    a.updateOld()
    b.updateOld()
    c.updateOld()
    res1 = 1e+10
    res2 = 1e+10
    while res1 > 1e-4 and res2 > 1e-4:    
        res1 = eqn_main.sweep(dt=dt)
        res2 =  eq_c.sweep(dt=dt)
        print("sweep")
    if (timestep % time_stride ==0):
        print ("Beep")
        if __name__ == '__main__':
            viewer.plot()
if __name__ == '__main__':
    input("Press <return> to proceed...")

