# Exercise 3

Thus far, we have solved equations that are comparatively simple. We now attempt to solve a more complicated problem of interest namely the Convection-Diffusion-Reaction (CDR) equation. This exercise is broken up into 2 parts. The first part demonstrates the difference between convective and diffusive transport. The second part outlines how reactions can be modelled. In particular, we explore coupling the transport equations for multiple species. 



## Part 1

### Problem statement

The model equation is written as follows: 
$$
\frac{\partial c}{\partial t} = D \frac{\partial ^{2}c}{\partial x^{2}} + u_{x}\frac{\partial c}{\partial x}
$$
where $D$ is the diffusion coefficient and $u_{x}$ is the velocity. $c$ is the mole fraction and is by definition, non-dimensional. Typically, the equations are scaled and solved in non-dimensional form which enables us to work in dimensionless groups instead. However, some solvers do not accommodate this  so readily e.g. `OpenFOAM`. 

We introduce the following scalings: 
$$
x = L \tilde{x}
$$

$$
u_{x} = U \tilde{u}_{x}
$$

$$
t = \frac{L^{2}}{D} \tilde{t}
$$

The model equation then becomes: 
$$
\frac{D}{L^{2}} \frac{\partial c}{\partial \tilde {t}} = \frac{D}{L^{2}} \frac{\partial ^{2}c}{\partial \tilde{x}^{2}} + \frac{U}{L} \tilde\frac{}{}
$$
