# Thrase


[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://thrase.github.io/Thrase.jl/dev/)

<img src="kobul.jpg" alt="">

# About

Thrase is a GPU-enabled, high-order accurate SBP-SAT finite difference code on unstructured meshes for SEAS (Sequences of Earthquakes and Aseismic Slip) problems written entirely in Julia.  

# Features:
- high-order accurate finite difference spatial discretizations featuring provably stability
- direct and matrix-free iterative solvers for static problems via Julia
- non-stiff methods for fully-dynamic (wave propagation) problems involving rate-and-state friction
- high-order accurate, adaptive time-stepping via Julia
- unstructured hexahedral meshes
- non-planar boundaries and interfaces 

# Dependencies: 
- Thrase is written entirely in <a href="https://julialang.org">Julia</a>

# Installation: 

## As a julia registered package (recommended):
- Download the latest release of <a href="https://julialang.org">Julia</a>
- Coming soon!: add the Thrase package: from the julia repl:
```
] add Thrase
```
- Test the installation:
``` 
] test Thrase
```

## From source (in your terminal):
- Download the latest release of <a href="https://julialang.org">Julia</a>
- Clone the repository
```
git clone https://github.com/Thrase/Thrase.jl.git
```
- Activate the project and download and compile dependencies
```
cd Thrase.jl
julia --project=.
] activate .
```
- Install any dependencies (only necessary the first time)
```
] instantiate
```
- Test the installation 
```
include("test/runtests.jl");
```

# Getting Started:
To run the code for the <a href="https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf">first SEAS benchmark problem 

```
localARGS = ["examples/bp1-qd.dat"]
include("examples/stripped_qd_driver.jl");

```
To run the code for <a href="https://strike.scec.org/cvws/seas/download/SEAS_BP6.pdf">BP6

```
localARGS = ["examples/bp6.dat"]
include("examples/BP6_driver.jl");


```
Note that the parameters set in the input file "examples/bp1-qd.dat" are slightly modified from the benchmark description to allow for fast simulation on most personal computers. The driver file calls a stripped down version of the code (for training purposes).
 
# Documentation: 
<a href="https://thrase.github.io/Thrase.jl/dev/"> More to come soon...

# How to cite: 
- Erickson, B. A., Kozdon, J. E., and Harvey, T. (2022), <a href="https://link.springer.com/article/10.1007/s10915-022-01961-1">A non-stiff summation-by-parts finite difference method for the wave equation in second order form: Characteristic boundary conditions and nonlinear interfaces</a>, Journal of Scientific Computing, doi: 10.1007/s10915-022-01961-1. 
- Kozdon, J. E., Erickson, B. A., and Wilcox, L. C. (2020), <a href="https://link.springer.com/article/10.1007/s10915-021-01448-5">Hybridized summation-by-parts finite difference methods</a>, Journal of Scientific Computing, doi: 10.1007/s10915-021-01448-5.
- Erickson, B. A. and Dunham, E. M. (2014), <a href="https://ix.cs.uoregon.edu/~bae/resources/Erickson_Dunham_jgrb50593.pdf">
  An efficient numerical method for earthquake cycles in heterogeneous media: Alternating sub-basin and surface-rupturing events on faults crossing a sedimentary basin</a>, Journal of  Geophysical Research, doi:10.1002/2013JB010614.

# License: 
- Distributed under the MIT License. See LICENSE.txt for more information.

# Contributing: 
Thrase is an open-source project and welcomes:
- contributions via forks and pull requests
- questions, feature requests, or bug reports via issues

# Contact:
- Brittany A. Erickson (bae@uoregon.edu)


