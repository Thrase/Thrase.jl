# Numerical Methods
This page describes how we solve the benchmark problem 1 (BP1-QD) by using numerical methods.

## Computational Domain
The first consideration to make is that we must convert the semi-infinite domain problem from the original description into a finite problem that we can compute.

First, we assume a fixed y value, which is the plane in which displacement happens. 
This means that our fields depend only on the x and z values, creating a two dimensional problem.

We restrict the x-axis of the computational domain in the following ways:
* x<sub>1</sub> - x<sub>2</sub> becomes one side of the problem, where we assume mirrored behavior on the other half
* x<sub>1</sub> can thus always be 0
* x<sub>2</sub> can be any value, a good default is 100 (km)
* the fault is located along x<sub>1</sub>

We then restrict the z-axis of the computational domain as such:
* z<sub>1</sub>=0 can be assumed as the location where the ground meets the air 
* z<sub>2</sub> can be any depth value, a good default is 100 (km)
* The external forces along z<sub>1</sub> equal 0 because it is a "free surface" where external forces can be assumed to have no effect
* The external forces along z<sub>2</sub> also equal 0 because it is far enough down to also be considered a "free surface"

![](img/ThraseBP1.jpg)

## Governing Equations
We now discuss the governing equations for our problem setup. 

(1) We are solving a 2D poisson equation given by:
```math
0= \mu(\frac{\partial^2u}{\partial^2x}+\frac{\partial^2u}{\partial^2z})
```

(2) We consider the material displacement over time as u(x,z,t). The slow tectonic loading (input parameter plate rate = V<sub>p</sub>) at the far right boundary (x<sub>2</sub>) is given by:
```math
u(x=x_2, z, t) = \frac{V_pt}{2}
```

(3) The material displacement at the far left boundary (x<sub>1</sub>) is given by:
```math
u(x=x_1, z, t) = \delta(z,t)
```

(4) The "free surface" at z<sub>1</sub> is given by:
```math
\mu\frac{\partial u}{\partial z}(x, z = z_1, t) = 0
```

(5) Likewise, the "free surface" at z<sub>2</sub> is given by:
```math
\mu\frac{\partial u}{\partial z}(x, z = z_2, t) = 0
```

## Converting $\theta$ into $\psi$
In the benchmark description we denote the state variable as $\theta$, but for computation we prefer to use the equivalent (mathematically consistent) $\psi$ as the state variable. We do this because it ranges over a smaller order of magnitude and is thus quicker(?) to compute. We describe how we convert from $\theta$ into $\psi$ for computing this problem.

Since (TODO: where did we get this from?)
```math
\theta = \frac{D_c}{V_0}exp[\frac{\psi-f_0}{b}]
```
We can solve for $\psi$ and obtain: 
```math
\psi = f_0 + b * ln(\frac{V_0\theta}{D_c})
```

TODO: and then I get lost...and in the code we start with all zeros for $\psi\delta$ and I'm confused by that too.

## Frictional Fault Boundary Condition Details
Now, using $\psi$ instead of $\theta$ we can define the frictional strength at the fault (x = x<sub>1</sub>) as:
TODO (how did we establish this from the benchmark description? I can't understand the correlation)
```math
\tau - \eta V = F(V,\psi)
```
where $\tau$ is the fault shear stress...and $F(V,\psi)$ is the frictional strength, and also something abou the slip rate...:
```math
\tau = \tau(z,t)
```
```math
F(V,\psi) = \sigma_n * a * sinh^{-1}[\frac{V}{2V_0}e^{\psi/a}]
```

## Numerical Time-Stepping Method

We illustrate our timestepping method using Forward Euler (from t<sup>n</sup> -> t<sup>n+1</sup> in one step). However, please note that in the code we use the TSit5() function which actually utilizes a Runge-Kutta method. The Runge-Kutta method employs the same steps as Forward Euler but calculates intermediate values between t<sup>n</sup> and t<sup>n+1</sup> to arrive at a more accurate solution for t<sup>n+1</sup>, thus it is more computationally intense. 

Assuming we know everything at time t<sup>n</sup> we take the following steps to calculate values at t<sup>n+1</sup>:

(1) Integrate $\delta$ and $\psi$, 
```math
\delta^{n+1} = \delta^n + dt * V^n 
```
```math
\psi^{n+1} = \psi^n + dt * G(V^n, \psi^n)
```
(2) Solve the poisson equation using the Summation-By-Parts Simultaneous Approximation Term (SBP-SAT) finite difference method.
This means that we solve the linear system Au<sup>n+1</sup> = b<sup>n+1</sup> for u<sup>n+1</sup>. Where u<sup>n+1</sup> is equal to the displacement everywhere. In this step we obtain the answer for u(x,z,t<sup>n+1</sup>) by solving:

```math
0= \mu(\frac{\partial^2u}{\partial^2x}+\frac{\partial^2u}{\partial^2z})
```
with the following boundary conditions for t<sup>n+1</sup>:
```math
u(x=x_2, z, t^{n+1}) = \frac{V_pt^{n+1}}{2}
```
```math
u(x=x_1, z, t^{n+1}) = \frac{\delta^{n+1}}{2}
```
```math
\mu\frac{\partial u}{\partial z}(x, z = z_1, t^{n+1}) = 0
```
```math
\mu\frac{\partial u}{\partial z}(x, z = z_2, t^{n+1}) = 0
```

(3) Compute $\tau^{n+1}$, since:
```math
\tau^{n+1} = \mu\frac{\partial u^{n+1}}{\partial x}
```
We solve numerically at x=x1:
```math
\left.\mu\frac{\partial u^{n+1}}{\partial x}\right\vert_{x=x_1}
```
(4) Solve for the new slip rate V<sup>n+1</sup> by imposing friciton (F). Our nonlinear equation is given by the following where everything is known except for V<sup>n+1</sup>:
```math
\tau^{n+1} -\eta V^{n+1} = F(V^{n+1}, \psi^{n+1}) 
```
We use Newton's method to solve for V<sup>n+1</sup> with: (TODO: Do we need this equation?)
```math
\bar V = V(z,t)
```

(5) Return to step 1 for timestep t<sup>n+2</sup>

## Other Considerations in the Code

Here should we talk about $\delta\psi$ and the Index-1 DAE? Or work it in somewhere else?
