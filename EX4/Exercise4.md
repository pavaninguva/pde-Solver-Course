# Exercise 4

We have reached a point where we have explored a variety of mathematical problems and features within `FiPy`. At this stage, it is possible to understand and implement a model that is still an area of active research from both an theoretical and application stand point. Whilst implementing and running the model is relatively straightforward, making progress on the problems which I will highlight is much more complicated. 

## Learning objectives

1. Understand how the Cahn-Hilliard equation is derived
2. Implement the model in `FiPy`
3. Understand the numerical issues associated with this equation system
4. (Not too important) Output files as `*.vtk` for visualization in Paraview

## Problem statement and derivation

The Cahn-Hilliard equation can model the demixing and precipitation of components from a mixture. If you recall Thermodynamics II, we can study the Gibbs energy of mixing $\Delta G _{mix}$ which gives us insight on whether the system will exhibit liquid-liquid equilibria (LLE).  We first start by considering the conventional diffusive flux expression that you should be familiar with for a species $A$: 
$$
\mathbf{j} = -D \nabla a
\\
\mathbf{j} = -D \frac{da}{dx} \ (\mathrm{in \ one \ dimension})
$$
where $\mathbf{j}$ is the flux, $D$ is the diffusion coefficient and $a$ is the mole fraction of species A. We can then apply the continuity / mass conversation equation: 
$$
\frac{\partial a}{\partial t} + \nabla \cdot\mathbf{j} = 0
$$
Which gives us the transient diffusion equation was have explored extensively before:
$$
\frac{\partial a}{\partial t} = D \nabla^{2}a
$$
When we model demixing, the mass transport is against the concentration gradient (we are moving in the direction of increasing concentration). This has been called "uphill" diffusion. We can see that the concentration gradient $\nabla a$ is no longer suitable as the driving force. We can actually write a more generalizable definition for diffusion which can model uphill diffusion. The driving force is gradients in the chemical potential $\mu$ rather than $a$. 

The flux expression for species $i$ , $\mathbf{j}_{i}$ is written as follows: 
$$
\mathbf{j}_{i} = -\sum_{j} L_{ij}\nabla \mu_{j}
$$
Where $L_{ij}$ is the mobility coefficient and $\mu_{i}$ is the chemical potential of species $i$. For species $1$ in a binary mixture: 
$$
\mathbf{j}_{1} = -(L_{11}\nabla \mu_{1} + L_{12}\nabla \mu_{2})
$$
The following constraints apply: 
$$
L_{ij} = L_{ji} \\
\sum_{i} L_{ij} = 0 \\
\sum_{i} \mathbf{j}_{i} = 0
$$
We can therefore write: 
$$
L_{11} + L_{12} = 0 
\\ L_{11} = - L_{12}
$$
When substituting this result back into $\mathbf{j}_{1}$: 
$$
\mathbf{j}_{1} = L_{12}(\nabla \mu_{1} - \nabla \mu_{2})
$$
Since the gradient operator is linear, we can further compactify the above expression: 
$$
\mathbf{j}_{i} = L_{12}\nabla\mu_{12}
$$
where $\mu_{12}  = \mu_{1} - \mu_{2}$. Working in differences for chemical potentials i.e. $\mu_{ij} = \mu_{i} - \mu_{j}$ is convenient for subsequent analysis. The mobility coefficient needs to be a function of composition, but the diffusion coefficient can be a constant: 
$$
L_{12} = D_{12} x_{1} (1-x_{1}) 
$$
where $D_{12}$ is the diffusion coefficient, $x_{1}$ is the mole fraction of species $i$.  

We can then apply the same continuity equation as we did for the vanilla diffusion equation which gives us the following equation: 
$$
\frac{\partial x_{1}}{\partial t} = \nabla \cdot (D_{12}x_{1}(1-x_{1}) \nabla \mu_{12})
$$
We now need to figure out an expression for the chemical potential To figure that out, we need to consider the total Gibbs energy of the system. We consider the Landau-Ginzburg free energy functional: 
$$
 G_{\text{system}}= \int_{V} g(x_{1},x_{2}...x_{N}) + \sum_{i}^{N-1}\frac{\kappa_{i}}{2}(\nabla x_{i})^2 +
	   \sum_{j>i}\sum_{i}^{N-1}\kappa_{ij}(\nabla x_{i})(\nabla x_{j}) \ dV
$$
For two components: 
$$
G_{\text{system}}= \int_{V} g(x_{1}) + \frac{\kappa}{2}(\nabla x_{1})^2 \ dV
$$
where $G_{system}$ is the total Gibbs energy of the system, $g(x_{1})$ is the homogenous free energy of mixing and $\kappa$ is the gradient energy parameter. When a mixture demixes, this is a spontaneous which means that the total Gibbs energy of the system decreases.  There are two forces that are considered: 

1. The system can reduces $G_{system}$ by firstly demixing which results in species that have unfavorable interactions concentrating in different phases 
2. However, demixing results in the formation of an interface betweeen the two phases. The interface itself has an energy associated with it and it is unfavourable. Hence the $\kappa$ term penalizes demixing. 

By definition, we know that the chemical potential $\mu_{i}$ for a conventional system can be written as follows: 
$$
\mu_{i} = \bigg( \frac{\partial G}{\partial x_{i}} \bigg)_{T, P}
$$
However, this definition is inadequate for inhomogeneous systems. We can apply the variational derivative to obtain a generalized expression for the chemical potential:
$$
\mu_{i} = \frac{\delta G_{system}}{\delta x_{i}} = \frac{\partial G}{\partial x_{i}} - \nabla \cdot \frac{\partial G}{\partial \nabla x_{i}}
$$
Therefore, $\mu_{i}$ can be written as follows:
$$
\mu_{1} = \frac{\partial g}{\partial x_{1}} - \nabla \cdot (\kappa \nabla x_{1})
$$
We assume that $\kappa$ is not dependent on composition, so: 
$$
\mu_{1} = \frac{\partial g}{\partial x_{1}} - \kappa \nabla^{2} x_{1}
$$
and: 
$$
\mu_{2} = 0
$$

$$
\therefore \mu_{12} = \frac{\partial g}{\partial x_{1}} - \kappa \nabla^{2}x_{1}
$$

We can observe that the system forms a 4th order PDE. However, solving a 4th order PDE imposes severe time-stepping requirements (referencing exercise 2). Hence, the Cahn-Hilliard equation is typically treated as a set of coupled 2nd order PDEs: 
$$
\frac{\partial x_{1}}{\partial t} = \nabla \cdot (D_{12}x_{1}(1-x_{1}) \nabla \mu_{12}) \\
\mu_{12} = \frac{\partial g}{\partial x_{1}} - \kappa \nabla^{2}x_{1}
$$
We now need an expression for $g(x_{1})$. For polymer blends and polymer-solvent systems, we can use the Flory-Huggins equation: 
$$
g(x_{1}) = \frac{x_{1}}{N_{1}}\ln{x_{1}} + \frac{(1-x_{1})}{N_{2}}\ln{(1-x_{1})} + \chi_{12}x_{1}(1-x_{1})
$$
where $N_{1}, N_{2}$ is the number of segments of the polymer chain of species $1$ and $2$ respectively and $\chi_{12}$ is the interaction parameter between species $1$ and $2$. From literature, the following expression can be used for $\kappa$: 
$$
\kappa = \frac{2}{3} R_{G}^{2} \chi_{12}
$$


Where $R_{G}$ is the radius of gyration of the polymer. To wrap things up, we scale our equations by introducing the following scaling relationships: 
$$
\tilde{t} = \frac{D_{12} t}{R_{G}^{2}} \\
\mathbf{x} = \tilde{\mathbf{x}} R_{G}
$$
We thus obtain the final equation system of interest for us: 
$$
\tilde{\mu}_{12} = \frac{\partial g}{\partial x_{1}} - \tilde{\kappa} \tilde{\nabla}^{2}x_{1} \\
\frac{\partial x_{1}}{\partial \tilde{t}} = \tilde{\nabla}\cdot (x_{1}(1-x_{1})\tilde{\nabla}\tilde{\mu}_{12})
$$

### Ideal limit

For sanity's sake, this part outlines how the modified Cahn-Hilliard equation reduces to the vanilla diffusion equation when the system is ideal. When the system is ideal, $\kappa =0$ since there are no enthalpic / residual interactions / contributions.  

For an ideal system, we know that: 
$$
g(x_{1}) = x_{1}\ln{x_{1}} + (1-x_{1}) \ln{(1-x_{1})} 
$$
Therefore: 
$$
\mu_{12} =  \ln{x_{1}} - \ln{(1-x_{1})} = \ln{\frac{x_{1}}{1-x_{1}}}
$$

$$
\nabla \mu_{12} = \frac{\partial \mu_{12}}{\partial x_{1}} \nabla x_{1} = \frac{1}{x_{1}(1-x_{1})} \nabla x_{1}
$$

Considering the first equation: 
$$
\frac{\partial x_{1}}{\partial t} = \nabla \cdot (D_{12}x_{1}(1-x_{1}) \nabla \mu_{12}) =  \nabla \cdot \bigg(D_{12}x_{1}(1-x_{1})  \frac{1}{x_{1}(1-x_{1})} \nabla x_{1} \bigg) = D_{12}\nabla^{2}x_{1}
$$
We recover our original equation for simple diffusion. 



## Numerical Implementation

The problem is defined on a 2D grid with periodic boundary conditions. For the initial conditions, we need to perturb the system with some noise (i.e. random spatial variation in the concentration field). We replace $1, 2$ with $_{a}, _{b}$ and $x_{1}$ with $a$ for clarity. 



The simulation parameters are defined as follows:

```python
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
```

Defining the mesh: 

```python
mesh = PeriodicGrid2D(nx=nx, ny=ny, dx=dx, dy=dy)
```

Defining the variables. We need to perform sweeps at each timestep due to the non-linear mobility coefficient: 

```python
# Defining the variables
a = CellVariable (name=r"$a$", mesh=mesh, hasOld=1)
mu_AB = CellVariable(name=r"$\mu_{AB}", mesh=mesh, hasOld=1)
```

Adding the initial noise:

```python
# Setting the initial composition of the system with noise
noise = UniformNoiseVariable(mesh=mesh, minimum=(a_0-noise_mag), maximum=(a_0+noise_mag))
a[:] = noise
```

Setting up the equations. We need to manually differentiate $g(x)$. The `FiPy` documentation introduces a trick which is also presented in the code, but is commented out currently, that you can explore. 

```python
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

```

We now introduce the idea of different solvers and solver settings. We can define this in `FiPy` quite simply. But do not worry about this. 

```python
# Setting up the solver
solver = LinearLUSolver(tolerance=1e-9, iterations=50, precon="ilu")
```

We now introduce the code needed to output `*.VTK` files that can be processed in software such as Paraview. That will be covered in a separate document. 

```python
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
```

## Results and area for study

