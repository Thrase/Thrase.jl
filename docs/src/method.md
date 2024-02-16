## Numerical Methods
This page describes how we solve the benchmark problem 1 by using numerical methods.

### Computational Domain
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


### Governing Equations
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
