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

# Getting Started
To run the code for the first [SEAS benchmark problem](https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf)

```
localARGS = ["examples/bp1-qd.dat"]
include("examples/stripped_qd_driver.jl");
```
Note that the parameters set in the input file "examples/bp1-qd.dat" are slightly modified from the benchmark description to allow for fast simulation on most personal computers. The driver file calls a stripped down version of the code (for training purposes).
