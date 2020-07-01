from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, Viewer
from fipy.tools import numerix


mesh = Grid2D(nx=30, ny=30, dx=0.25, dy=0.25)
phi = CellVariable(name=r"$\phi$", mesh=mesh)
psi = CellVariable(name=r"$\psi$", mesh=mesh)

noise = GaussianNoiseVariable(mesh=mesh,
                               mean=0.5,
                               variance=0.01).value

phi[:] = noise

if __name__ == "__main__":
    viewer = Viewer(vars=(phi,), datamin=0., datamax=1.)

D = a = epsilon = 1.

dfdphi = a**2 * phi * (1 - phi) * (1 - 2 * phi)
d2fdphi2 = a**2 * (1 - 6 * phi * (1 - phi))

eq1 = (TransientTerm(var=phi) == DiffusionTerm(coeff=D, var=psi))
eq2 = (ImplicitSourceTerm(coeff=1., var=psi) == ImplicitSourceTerm(coeff=d2fdphi2, var=phi) - d2fdphi2 * phi + dfdphi - DiffusionTerm(coeff=epsilon**2, var=phi))

eq = eq1 & eq2

dexp = -5
elapsed = 0.
if __name__ == "__main__":
    duration = 1000.
else:
    duration = .5e-1


while elapsed < duration:
    dt = min(100, numerix.exp(dexp))
    elapsed += dt
    dexp += 0.01
    eq.solve(dt=dt)
    if __name__ == "__main__":
        viewer.plot()

from fipy import input
if __name__ == '__main__':
    input("Coupled equations. Press <return> to proceed...")