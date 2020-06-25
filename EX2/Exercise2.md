# Exercise 2 

In the first exercise, we solved a simple "steady-state" problem (i.e. there was no time derivative) in the problem. We will know look at one of the simplest transient problems that you should be familiar with from `Heat and Mass Transfer 1`. 

This exercise is broken into two parts as we want to explore different aspects of the same equation system. In the first part, we focus more on how transient problems are tackled and we consider the case where Dirichlet boundary conditions are applied. In the second part, we approach the problem more from a transport perspective and how a reaction in the system can be modelled. 

## Learning objectives

By the end of this exercise, we hope to achieve the following goals: 

1. Appreciate the difference between numerically solving transient and steady-state problems
2. Appreciate how different time-stepping methods can impact the solution
3. Visualize how diffusion takes place in different conditions
4. Perform basic post-processing to generate relevant plots 
5. Appreciate how reactions and similar processes are treating in the context of numerical solutions



## Part 1

The equation system of for the transient problem can be written as follows: 
$$
\frac{\partial c}{\partial t} = D \nabla ^{2}c,
$$
where $c$ is the concentration / mole fraction, $t$ is the time and $D$ is the diffusion coefficient. The treatment of the diffusion term is not too dissimilar from the previous exercise, but we need to think about how we handle the time derivative. 

### Time stepping

For this exercise, we will consider two simple ways of handling the time derivative:

1. Fully explicit forward Euler time-stepping 
2. Fully implicit backward Euler time-stepping

This part may get a bit technical, so I will attempt to flesh out every step and keep it as simple as possible. We are able to write an approximation for the full equation as follows for $x$ at a point $i$ as follows:
$$
\frac{c_{i}^{new} - c_{i}^{old}}{\Delta t} = \alpha D \frac{c^{new}_{i+1}- 2c^{new}_{i} + c^{new}_{i-1}}{\Delta x^{2}} + (1-\alpha) D \frac{c^{old}_{i+1}- 2c^{old}_{i} + c^{old}_{i-1}}{\Delta x^{2}}
$$
We introduce a new parameter $\alpha$ which readily enables us to pick a specific numerical scheme. The three values of $\alpha$ that are of particular interest are (We are only considering the first 2): 

1. $\alpha = 1$: Fully implicit backward Euler time-stepping. 
2. $\alpha=0$: Fully explicit forward Euler time-stepping
3. $\alpha = 0.5$: Crank-Nicolson time-stepping. 

When $\alpha = 1$: 
$$
\frac{c_{i}^{new} - c_{i}^{old}}{\Delta t} =  D \frac{c^{new}_{i+1}- 2c^{new}_{i} + c^{new}_{i-1}}{\Delta x^{2}}
$$
When $\alpha = 0$: 
$$
\frac{c_{i}^{new} - c_{i}^{old}}{\Delta t} =  D \frac{c^{old}_{i+1}- 2c^{old}_{i} + c^{old}_{i-1}}{\Delta x^{2}}
$$
When a fully explicit scheme such as the forward Euler is used, the new value (moving forward in time) can be evaluated explicitly from values from the old values from the previous timestep which makes it much more easy to implement than implicit methods. However, the trade off is that often, explicit time stepping requires excessively small timesteps for numerical stability while the implicit backwards Euler time-stepping is far more stable. 

#### Stability for forward Euler time stepping

For diffusion problems (We can treat heat and mass transfer as diffusion problems collectively), we can consider the following metric: 

We can write down the dimensionless mesh Fourier number $Fo_{Mesh}$: 
$$
Fo_{Mesh} = D \frac{\Delta t}{\Delta x ^{2}}
$$
For fullness, a not entirely rigorous derivation of the stability criteria for the explicit Forward Euler method is presented. Feel free to gloss over this section. 

So we start with the discretization scheme of the explicit forward Euler: 
$$
\frac{c_{i}^{new} - c_{i}^{old}}{\Delta t} =  D \frac{c^{old}_{i+1}- 2c^{old}_{i} + c^{old}_{i-1}}{\Delta x^{2}}
$$
Rejigging things: 
$$
c_{i}^{new} = c_{i}^{old} + D\frac{\Delta t}{\Delta x^{2}} \bigg( c^{old}_{i+1}- 2c^{old}_{i} + c^{old}_{i-1} \bigg)
$$
Introducing $Fo_{Mesh}$:
$$
c_{i}^{new} = c_{i}^{old} + Fo_{Mesh} \bigg( c^{old}_{i+1}- 2c^{old}_{i} + c^{old}_{i-1} \bigg)
$$
We can further rejig this to group all the $c_{i}^{old}$ terms together:
$$
c_{i}^{new} = (1- 2 Fo_{Mesh}) c_{i}^{old} + 2Fo_{Mesh} \bigg( \frac{c^{old}_{i+1} + c^{old}_{i-1} } {2} \bigg)
$$
The insight to be gained from rejigging the equation to the above form is that $x_{i}^{new}$ is the weighted average of the old value $x$ at that point $i$ and the average of the neighboring points' $i+1$  and $i-1$ old value of $c$. If $Fo_{Mesh} = 0$, then its easy to see that we are not advancing forward in time. When $Fo_{Mesh} = 0.5$: 
$$
c_{i}^{new} = \frac{1}{2} (c^{old}_{i+1} + c^{old}_{i-1})
$$
We see that the new value is the average of the old value of the two adjacent points. When we have $Fo_{Mesh} > 0.5$, we run into the problem where updated value unphysically overshoots the average of the neighboring values. So this gives us the following stability criteria which we can use: 
$$
\Delta t \leq \frac{\Delta x ^{2}}{2D}
$$
Another useful insight that can be gleaned from this relationship is that when we increase the mesh resolution (i.e. make $\Delta x$ smaller by a 2X), we have to decrease the timestep $\Delta t$ 4X

### Problem statement

We first intend to model a simple 1D slab of size $1$  where one end is held at a fixed concentration and the other end is held is zero concentration. There is no reaction and the concentration throughout the block is initially zero. The model equation of interest can be written as follows: 
$$
\frac{\partial c}{\partial t} = D \nabla^{2}c
$$
With the following boundary conditions: 
$$
c(x=0, t) = 1
\\
c(x=1, t) = 0
$$
And the initial conditions: 
$$
c(x, t=0) = 0
$$

### Analytical solution

Lets work through the problem. We quite quickly encounter the fact that due to the non-homogenous boundary conditions, a naïve application of separation of variables will not work. The boundary conditions require modification. 

Fortunately, the physical nature of the problems lends to the solution. By fixing the concentrations at the ends of the slab and not considering a flux condition, the system will eventually tend to a steady state solution. So we can employ this idea to perform a variable transformation that enables us to reformulate the boundary conditions as homogenous BCs. 

The steady state solution can be obtained as follows (We cancel the time derivative which simplifies the problem to a 1D ODE) : 
$$
\frac{d^{2}c}{dx^{2}} = 0
$$
The solution of this is given by: 
$$
c_{S.S.} (x) = c_{1}x + c_{2}
$$
Substituting the boundary conditions, we get: 
$$
c_{S.S.}(x) = 1 -x
$$
We now introduce the following variable transform: 
$$
\phi(x,t) = c(x,t) - c_{S.S.}(x)
$$
Evaluating the relevant derivatives: 
$$
\frac{\partial \phi}{\partial t} = \frac{\partial c}{\partial t}
$$

$$
\frac{\partial ^{2} \phi}{\partial x^{2}} = \frac{\partial^{2}c}{\partial x^{2}}
$$

This is a nice result! The equation does not change shape in anyway with the variable transform!! Lets consider how the BCs evolve: 
$$
\phi(0, t) = c(0,t) - c_{S.S.}(0) = 1 -1 = 0
$$

$$
\phi(1,t) = c(1,t) - c_{S.S.}(1) = 0 -0 = 0
$$

We get nice homogenous BCs which enable us to solve for $\phi(x,t)$ using separation of variables. We obtain the following result: 
$$
\phi (x,t) = \sum_{n=1}^{\infin} a_{n} \sin{(n \pi x)}e^{-Dn^{2}\pi^{2} t}
$$
To evaluate the value of the coefficients $a_{n}$: we need to consider the initial condition $\phi(x, 0)$:
$$
\phi(x,0) = c(x,0) - c_{S.S.}(x) = 0 - (1-x) = x -1
$$
$a_{n}$ is given as follows: 
$$
a_{n} = 2 \int_{0}^{1} (x-1) \sin{(n\pi x)} \ dx, \ \ \  n=1,2,3 \dots
$$
This is a chore to evaluate. We employ the power of `WolframAlpha`! 
$$
a_{n} = \frac{2}{\pi^{2}n^{2}}( \sin{(\pi n)} - \pi n) = \frac{-2}{n \pi}, \ \ \ n = 1,2,3 \dots
$$
Now that we have a nice clean expression for $a_{n}$, let us finally write down our full analytical solution to this problem. 
$$
c(x,t) = \phi(x,t) + c_{S.S.}(x) = 1 - x - \sum_{n=1}^{\infin}\frac{2}{n \pi} \sin{(n\pi x) e^{-Dn^{2}\pi^{2} t}}
$$

### Numerical solution

Setting up the problem for the most part is similar, the only difference being how the equation is defined and how time-stepping is performed. I will detail those two sections here for clarity. 

#### Defining the equation

We make use of the full time discretization scheme involving $\alpha$ as outlined above as this gives us the convenience of toggling between explicit and implicit time-stepping very easily

```python
eq = (TransientTerm() == DiffusionTerm(coeff= alpha* D) + ExplicitDiffusionTerm(coeff = (1.0 - alpha)*D))
```

We can now set $\alpha$ as a simulation parameter easily and not worry about rejigging the equation. 



#### Time-stepping

First we need to pick an $\alpha$ value. Suppose we pick $\alpha = 0$ for now, we can use our stability criteria as follows: 
$$
\Delta t = \frac{\Delta x ^{2}}{2D} \times (\% Margin) = \frac{4\times10^{-4}}{2.0} \times 0.9 = 0.00018
$$


Now that we have a time-step, the next question we need to figure out is how long to run the simulation for! Now that may seem like a straightforward question in that this should normally be pre-defined. But in this case, we dont immediately have a fixed value for the time taken for the system to reach steady state. So lets employ a bit of street-fighting for this:

- Our simulation domain is $1m$ long
- We set our diffusion coefficient to be $1 m^{2} / s$ 
- So an approximation for the time taken for a molecule to diffuse across the entire domain is $ t \sim \frac{D}{L^{2}}$ by considering simple dimensional analysis. 
- Therefore, the time would be approximately $1s$ 