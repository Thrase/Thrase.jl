# Input Parameters

In the example input file for benchmark problem 1 there are multiple parameters that can be updated to change the problem.
They are desribed in the table below. 

Additional information is provided link them if they have slightly different names in the paper: https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf

| Input File Parameter | Units | Data Type | Description | 
| :-------- | :--------: | :--------: | :-------- | 
| pth | - | string | path to the direcotory that stores output data 
| stride_space | - | positive integer | write-out every "stride_space" grid points along fault
| stride_time | - | positive integer | write-out every "stride_time" time steps
| x1 | km | positive integer | x-domain is (x1, x2)
| x2 | km | positive integer | x-domain is (x1, x2)
| z1 | km | positive integer | z-domain is (z1, z2)
| z2 | km | positive integer | z-domain is (z1, z2)
| Nx | - | positive integer | number of grid points in each dimension, make sure to resolve h* and l_b 
| Nz | - | positive integer | number of grid points in each dimension, make sure to resolve h* and l_b
| sim_years | years | positive integer | how many years to simulate
| Vp | m/s | floating point | plate rate
| ρ | kg/m<sup>3</sup> | floating point | density
| cs | km/s | floating point | shear wave speed
| σn | MPa | positive integer | effective normal stress on fault
| RSamin | - | floating point | rate-and-state parameter (a<sub>0</sub>)
| RSamax | - | floating point | rate-and-state parameter (a<sub>max</sub>)
| RSb | - | floating point | rate-and-state parameter (b<sub>0</sub>)
| RSDc | m | floating point | critical slip distance (D<sub>c</sub>)
| RSf0 | m/s | floating point | reference friction 
| RSV0 | m/s | floating point | reference slip rate initial slip rate (V<sub>0</sub>)
| RSVinit | m/s | floating point | initial slip rate (V<sub>init</sub>)
| RSH1 | km | positive integer | depth extent of uniform VS region (H)
| RSH2 | km | positive integer | width of VS-VS transition zone (h)
| RSWf | km | positive integer | width of rate-and-state fault (W<sub>f</sub>)
| SBPp | - | positive integer | SBP interior spatial order of accuracy
