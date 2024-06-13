
Once you have installed the Thrase project you should be able to run any of the below examples. If you want to change any input parameters you can do so by updating the .dat files for each benchmark problem.

# Benchmark Problem 1
To run the code for the first <a href="https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf">SEAS benchmark problem

```
localARGS = ["examples/bp1-qd.dat"]
include("examples/stripped_qd_driver.jl");
```

Note that the parameters set in the input file "examples/bp1-qd.dat" are slightly modified from the benchmark description to allow for fast simulation on most personal computers. The driver file calls a stripped down version of the code (for training purposes).


# Benchmark Problem 6
To run the example code for the <a href="https://strike.scec.org/cvws/seas/download/SEAS_BP6.pdf">sixth SEAS benchmark problem

```
localARGS = ["examples/BP6/bp6.dat"]
include("examples/BP6/BP6_driver.jl");
```

# Benchmark Problem 5 (GPU):
To run the example code for the <a href="https://strike.scec.org/cvws/seas/download/SEAS_BP5.pdf">fifth SEAS benchmark problem</a> you will need GPU resources to solve in 3 dimensions, without a GPU it will not be able to complete in a timely manner
```
localARGS = ["examples/BP5/bp5.dat"]
include("examples/BP5/BP5-QD.jl");
```
The setup of the domains are done in "examples/BP5/domain.jl" and "examples/BP5/domain_256.jl" for the different problem sizes respectively. In addition to updating the .dat file, you may need to change values in the domain file if needing to configure the domain size differently.