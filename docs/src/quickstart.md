# Installation 

## As a julia registered package (recommended)
- Download the latest release of [Julia](https://julialang.org)
- Coming soon!: add the Thrase package: from the julia repl:
```
] add Thrase
```
- Test the installation:
``` 
] test Thrase
```

## From source (in your terminal)
- Download the latest release of [Julia](https://julialang.org)
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
