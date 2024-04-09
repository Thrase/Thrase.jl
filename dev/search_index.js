var documenterSearchIndex = {"docs":
[{"location":"method/#Numerical-Methods","page":"Numerical Method","title":"Numerical Methods","text":"","category":"section"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"This page describes how we solve the benchmark problem 1 (BP1-QD).","category":"page"},{"location":"method/#Computational-Domain","page":"Numerical Method","title":"Computational Domain","text":"","category":"section"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"The first consideration to make is that we must convert the semi-infinite domain problem from the original description into a finite domain problem to which we can apply finite difference methods.","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"The anti-plane assumption means that fields depend only on the x and z values, creating a two dimensional problem.","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"We restrict our domain to (x z) in (0 L_x) x (0 L_z), with the fault at x = 0, and assume anti-symmetry for x in (-L_x 0). ","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"z = 0","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"refers to Earth's free surface, and at z = L_z we also assume a traction-free boundary. ","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(Image: BP schematic 1)","category":"page"},{"location":"method/#Governing-Equations","page":"Numerical Method","title":"Governing Equations","text":"","category":"section"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"As described in the benchmark description, the governing equations are","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(1) the 2D poisson equation given by:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"0= mu(fracpartial^2upartial x^2+fracpartial^2upartial z^2)","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(2) The material displacement is given by u(x,z,t). Slow tectonic loading (input parameter plate rate = V_p) is imposed at the far right boundary (x_2) is given by:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"u(x=x_2 z t) = fracV_pt2","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(3) The material displacement at the fault is given by:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"u(x=x_1 z t) = delta(zt)","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"where delta(z t) is the fault slip.","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(4) The free surface at z = 0 is given by:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"mufracpartial upartial z(x z = 0 t) = 0","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(5) Likewise, the free surface at z<sub>2</sub> is given by:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"mufracpartial upartial z(x z = L_z t) = 0","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(Image: BP schematic 2)","category":"page"},{"location":"method/#Converting-\\theta-into-\\psi","page":"Numerical Method","title":"Converting theta into psi","text":"","category":"section"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"In the benchmark description the state variable is given in terms of theta, but we prefer to use the equivalent (mathematically consistent) psi as the state variable, defined by","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"psi = f_0 + b ln(fracV_0thetaD_c)","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"This defines the following aging law for psi:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"fracdpsidt = G(V psi)","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"where","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"G(V psi) = fracbV_0D_c e^fracf_0-psib-fracVV_0","category":"page"},{"location":"method/#Frictional-Fault-Boundary-Condition-Details","page":"Numerical Method","title":"Frictional Fault Boundary Condition Details","text":"","category":"section"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"Now, using psi instead of theta we can define the frictional strength at the fault (x = 0) from equation 5 in the benchmark as:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"tau = F(Vpsi)","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"where tau is the fault shear stress...and F(Vpsi) is the frictional strength:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"F(Vpsi) = sigma_n a sinh^-1fracV2V_0e^psia","category":"page"},{"location":"method/#Numerical-Time-Stepping-Method","page":"Numerical Method","title":"Numerical Time-Stepping Method","text":"","category":"section"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"We formulate the governing equations as an Index-1 differential algebraic equation (DAE) where slip and state evolve in time and a nonlinear equation for slip rate must be solved at each timestep.","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"We illustrate our timestepping method using Forward Euler (from t^n -> t^n+1 in one step). However, please note that in the code we use Julia's TSit5() function (a 4/5-order adaptive Runge-Kutta method). ","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"Assuming we know all fields at time t^n we take the following steps to calculate values at t^n+1:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(1) Integrate delta and psi, ","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"delta^n+1 = delta^n + dt V^n ","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"psi^n+1 = psi^n + dt G(V^n psi^n)","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(2) Solve the poisson equation using the Summation-By-Parts Simultaneous Approximation Term (SBP-SAT) finite difference method. This amounts to solving the linear system Au^n+1 = b^n+1 for u^n+1. Where u^n+1 is the displacement within the entire domain at time t^n+1. In this step we assume the following boundary conditions at t^n+1:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"u(x=L_x z t^n+1) = fracV_pt^n+12","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"u(x=0 z t^n+1) = fracdelta^n+12","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"mufracpartial upartial z(x z = 0 t^n+1) = 0","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"mufracpartial upartial z(x z = L_z t^n+1) = 0","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(3) In the rate-and-state portion of the fault we next compute the change in shear stress due to quasi-static deformation Delta tau^n+1 via:","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"Deltatau^n+1 = leftmufracpartial u^n+1partial xrightvert_x=0","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(4) and then solve for the new slip rate V^n+1 by imposing friction. This yields a nonlinear equation (where everything is known except for V^n+1):","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"tau_0 + Deltatau^n+1 - eta V^n+1 = F(V^n+1 psi^n+1) ","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"which is solved usign a bracketed Newton method (i.e. safe-guarded with bisection). ","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"For fault depths below the rate-and-state region, V^n+1 is set to the plate rate V_p as specified in the benchmark. ","category":"page"},{"location":"method/","page":"Numerical Method","title":"Numerical Method","text":"(5) Return to step 1 for timestep t^n+1","category":"page"},{"location":"quickstart/#Installation:","page":"Install and Quick Start","title":"Installation:","text":"","category":"section"},{"location":"quickstart/#As-a-julia-registered-package-(recommended):","page":"Install and Quick Start","title":"As a julia registered package (recommended):","text":"","category":"section"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"Download the latest release of Julia\nComing soon!: add the Thrase package: from the julia repl:","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"] add Thrase","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"Test the installation:","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"] test Thrase","category":"page"},{"location":"quickstart/#From-source-(in-your-terminal):","page":"Install and Quick Start","title":"From source (in your terminal):","text":"","category":"section"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"Download the latest release of Julia\nClone the repository","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"git clone https://github.com/Thrase/Thrase.jl.git","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"Activate the project and download and compile dependencies","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"cd Thrase.jl\njulia --project=.\n] activate .","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"Install any dependencies (only necessary the first time)","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"] instantiate","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"Test the installation ","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"include(\"test/runtests.jl\");","category":"page"},{"location":"quickstart/#Getting-Started:","page":"Install and Quick Start","title":"Getting Started:","text":"","category":"section"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"To run the code for the first SEAS benchmark problem","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"localARGS = [\"examples/bp1-qd.dat\"]\ninclude(\"examples/stripped_qd_driver.jl\");","category":"page"},{"location":"quickstart/","page":"Install and Quick Start","title":"Install and Quick Start","text":"Note that the parameters set in the input file \"examples/bp1-qd.dat\" are slightly modified from the benchmark description to allow for fast simulation on most personal computers. The driver file calls a stripped down version of the code (for training purposes).","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Thrase","category":"page"},{"location":"#Thrase","page":"Home","title":"Thrase","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Thrase.","category":"page"},{"location":"#About","page":"Home","title":"About","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Thrase is a GPU-enabled, high-order accurate SBP-SAT finite difference code on unstructured meshes for SEAS (Sequences of Earthquakes and Aseismic Slip) problems written entirely in Julia.  ","category":"page"},{"location":"#Features:","page":"Home","title":"Features:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"high-order accurate finite difference spatial discretizations featuring provably stability\ndirect and matrix-free iterative solvers for static problems via Julia\nnon-stiff methods for fully-dynamic (wave propagation) problems involving rate-and-state friction\nhigh-order accurate, adaptive time-stepping via Julia\nunstructured hexahedral meshes\nnon-planar boundaries and interfaces ","category":"page"},{"location":"#Dependencies:","page":"Home","title":"Dependencies:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Thrase is written entirely in Julia","category":"page"},{"location":"#How-to-cite:","page":"Home","title":"How to cite:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Erickson, B. A., Kozdon, J. E., and Harvey, T. (2022), A non-stiff summation-by-parts finite difference method for the wave equation in second order form: Characteristic boundary conditions and nonlinear interfaces, Journal of Scientific Computing, doi: 10.1007/s10915-022-01961-1. \nKozdon, J. E., Erickson, B. A., and Wilcox, L. C. (2020), Hybridized summation-by-parts finite difference methods, Journal of Scientific Computing, doi: 10.1007/s10915-021-01448-5.\nErickson, B. A. and Dunham, E. M. (2014), An efficient numerical method for earthquake cycles in heterogeneous media: Alternating sub-basin and surface-rupturing events on faults crossing a sedimentary basin, Journal of  Geophysical Research, doi:10.1002/2013JB010614.","category":"page"},{"location":"#License:","page":"Home","title":"License:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Distributed under the MIT License. See LICENSE.txt for more information.","category":"page"},{"location":"#Contact:","page":"Home","title":"Contact:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Brittany A. Erickson (bae@uoregon.edu)","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Thrase]","category":"page"},{"location":"parameters/#Input-Parameters","page":"Input Parameters","title":"Input Parameters","text":"","category":"section"},{"location":"parameters/","page":"Input Parameters","title":"Input Parameters","text":"In the example input file for benchmark problem 1 (Thrase.jl/examples/bp1-qd.dat) there are multiple parameters that can be updated to change the problem. They are desribed in the table below. ","category":"page"},{"location":"parameters/","page":"Input Parameters","title":"Input Parameters","text":"Sometimes the benchmark description paper contains slightly different names than the input file, we indicate in the description in parantheses if this is the case.","category":"page"},{"location":"parameters/","page":"Input Parameters","title":"Input Parameters","text":"Input File Parameter Units Data Type Description\npth - string path to the directory that stores output data\ngenerate_plots - boolean option to generate example plots\nstride_space - Real <: Number write-out every \"stride_space\" grid points along fault\nstride_time - Real <: Number write-out every \"stride_time\" time steps\nx1 km Real <: Number x-domain is (x1, x2)\nx2 km Real <: Number x-domain is (x1, x2)\nz1 km Real <: Number z-domain is (z1, z2)\nz2 km Real <: Number z-domain is (z1, z2)\nNx - Real <: Number number of grid points in each dimension\nNz - Real <: Number number of grid points in each dimension\nsim_years years Real <: Number how many years to simulate\nVp m/s Real <: Number plate rate\nρ kg/m^3 Real <: Number density\ncs km/s Real <: Number shear wave speed\nσn MPa Real <: Number effective normal stress on fault\nRSamin - Real <: Number rate-and-state parameter (a<sub>0</sub>)\nRSamax - Real <: Number rate-and-state parameter (a<sub>max</sub>)\nRSb - Real <: Number rate-and-state parameter (b<sub>0</sub>)\nRSDc m Real <: Number critical slip distance (D<sub>c</sub>)\nRSf0 m/s Real <: Number reference friction\nRSV0 m/s Real <: Number reference slip rate initial slip rate (V<sub>0</sub>)\nRSVinit m/s Real <: Number initial slip rate (V<sub>init</sub>)\nRSH1 km Real <: Number depth extent of uniform VS region (H)\nRSH2 km Real <: Number width of VS-VS transition zone (h)\nRSWf km Real <: Number width of rate-and-state fault (W<sub>f</sub>)\nSBPp - Real <: Number SBP interior spatial order of accuracy","category":"page"}]
}
