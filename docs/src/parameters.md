# Input Parameters

In the example input file for benchmark problem 1 (Thrase.jl/examples/bp1-qd.dat) there are multiple parameters that can be updated to change the problem.
They are desribed in the table below. 

Sometimes the [benchmark description paper](https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf) contains slightly different names than the input file, we indicate in the description in parantheses if this is the case.

| Input File Parameter | Units | Data Type | Description | 
| :-------- | :--------: | :--------: | :-------- | 
| pth | - | string | path to the directory that stores output data 
| generate_plots | - | boolean | option to generate example plots
| stride_space | - | Real <: Number | write-out every "stride_space" grid points along fault
| stride_time | - | Real <: Number | write-out every "stride_time" time steps
| Lx1 | km | Real <: Number | x-domain is (Lx1, Lx2)
| Lx2 | km | Real <: Number | x-domain is (Lx1, Lx2)
| Lz1 | km | Real <: Number | z-domain is (Lz1, Lz2)
| Lz2 | km | Real <: Number | z-domain is (Lz1, Lz2)
| Nx | - | Real <: Number | number of grid points in each dimension
| Nz | - | Real <: Number | number of grid points in each dimension
| sim_years | years | Real <: Number | how many years to simulate
| Vp | m/s | Real <: Number | plate rate
| ρ | kg/$m^3$ | Real <: Number | density
| cs | km/s | Real <: Number | shear wave speed
| σn | MPa | Real <: Number | effective normal stress on fault
| RSamin | - | Real <: Number | rate-and-state parameter ($a_0$)
| RSamax | - | Real <: Number | rate-and-state parameter ($a_max$)
| RSb | - | Real <: Number | rate-and-state parameter ($b_0$)
| RSDc | m | Real <: Number | critical slip distance ($D_c$)
| RSf0 | m/s | Real <: Number | reference friction 
| RSV0 | m/s | Real <: Number | reference slip rate initial slip rate ($V_0$)
| RSVinit | m/s | Real <: Number | initial slip rate ($V_{init}$)
| RSH1 | km | Real <: Number | depth extent of uniform VS region (H)
| RSH2 | km | Real <: Number | width of VS-VS transition zone (h)
| RSWf | km | Real <: Number | width of rate-and-state fault ($W_f$)
| SBPp | - | Real <: Number | SBP interior spatial order of accuracy
