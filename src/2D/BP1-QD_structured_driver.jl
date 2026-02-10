# Include a coordinate tranformation.
# On a structured mesh of quads:
# Solves the SEAS benchmark problem BP1-QD:
# https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf

using Thrase
using LinearAlgebra
using DifferentialEquations
using DiffEqCallbacks
using DelimitedFiles
using Dates

const year_seconds = 31556926
global const ctr = Ref{Int64}(1) 

include("ops_BP1-QD_structured.jl")
include("odefun_BP1-QD_structured.jl")
include("../utils_2D.jl")

function main()


  ### read input parameters from .dat file
  (pth, stride_space, stride_time, Lx, Lz, Hx, Hz, Nr, Ns, dx, 
  dz, el_r, el_s, sim_years, Vp, ρ, cs, σn, RSamin, RSamax, 
  RSb, RSDc, RSf0, RSV0, RSVinit, RSH1, RSH2, RSWf, SBPp) = read_params_BP1_CT(localARGS[1])

  
  try
    mkdir(pth)
  catch
    # folder already exists so make a new one.
    pth = pth*string(now())*"/"
    mkdir(pth)
  end

  year_seconds = 31556926
  #μ = cs^2 * ρ 
  μ = (x, z) -> (cs^2 * ρ .+ 0 .* x .+ 0 .* z)
  μshear = cs^2 * ρ
  η = μshear / (2 * cs)

  ################################## COORDINATE TRANSFORM ###################################
  rstar = Hx*(2/Nr)/dx - 1 # b/t -1 and 1; r_star = 1 is the limit of constant grid spacing
  sstar = Hz*(2/Ns)/dz -1 # b/t -1 and 1; s_star = 1 is the limit of constant grid spacing
  cx = log(Lx/Hx)/(1-rstar)
  cz = log(Lz/Hz)/(1-sstar)
  xt = (r,s) -> ((r .<= rstar) .* (Hx .* r ./ (rstar .+ 1) .+ Hx ./ (rstar .+ 1)) .+ (r .> rstar) .* (Hx .* exp.((r .- rstar) .* cx)) ,  (r .<= rstar) .* (Hx .* ones(size(r)) ./ (rstar .+ 1)) .+ (r .> rstar) .* (Hx .* cx .* exp.((r .- rstar) .* cx)), zeros(size(s)))
  zt = (r,s) -> ((s .<= sstar) .* (Hz .* s ./ (sstar .+ 1) .+ Hz ./ (sstar .+ 1)) .+ (s .> sstar) .* (Hz .* exp.((s .- sstar) .* cz)) , zeros(size(r)),  (s .<= sstar) .* (Hz .* ones(size(s)) ./ (sstar .+ 1)) .+ (s .> sstar) .* (Hz .* cz .* exp.((s .- sstar) .* cz)))


  metrics = create_metrics(Nr, Ns, μ, xt, zt) # create coordinate transform
    
     x = metrics.coord[1]
     z = metrics.coord[2]
     #p = scatter(x[2:end, 1], diff(x[:, 1]))
     #p = scatter(z[1, 2:end], diff(z[1, :]))
     #display(p)
     #print()

  ###################################################################### 
  # create finite difference operators on computational domain:
  (A, B, H̃, T, e, J, coord, facecoord, sJ, nx, ny) = get_operators(SBPp, Nr, Ns, metrics) 

  A = lu(A)  # matrix factorization

  # initialize time and vector b that stores boundary data (linear system will be Au = b, where b = B*g)
  t = 0
  b = zeros((Nr+1) * (Ns+1))

  # initial slip vector
  δ = zeros(Ns+1)
  

  # fill in initial boundary data into b
  bdry_vec!(b, B, δ, (t .* Vp./2)*ones(Ns+1), zeros(Nr+1), zeros(Nr+1), sJ)

  
  u = A \ b # solve linear system with a backsolve to obtain initial displacement u(x, z, 0)

  # Find z-index δNp corresponding to base of rate-and-state friction zone RSWf
  z = metrics.coord[2]

  (mm, δNp) = findmin(abs.(RSWf .- z[1,:]))
  @show z[1,δNp]
  #@assert z[1,δNp] ≈ RSWf

  # initialize change in shear stress due to quasi-static deformation
  Δτ = zeros(Ns+1)

  # Assemble fault variables/data
  RSa = zeros(δNp)
  for n = 1:δNp
      RSa[n] = RSamin - (RSamin - RSamax) *
          min(1, max(0, (RSH1 - z[1,n])/(RSH1 - RSH2)))
  end

  # Set pre-stress according to benchmark
  τ0 = σn * RSamax * asinh(RSVinit / (2 * RSV0) *
                                  exp.((RSf0 + RSb * log.(RSV0 / RSVinit)) /
                                      RSamax)) + η * RSVinit


  # Set initial state variable according to benchmark
  θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
      sinh.((τ0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)

  # Initialize psi version of state variable
  ψ = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)


  # Set initial condition for index 1 DAE - this is a stacked vector of psi, followed by slip
  ψδ = zeros(δNp + Ns + 1)  #because length(ψ) = δNp,  length(δ) = Nz+1
  ψδ[1:δNp] .= ψ
  ψδ[δNp+1:δNp + Ns + 1] .= δ

  # Set fault station locations (depths) specified in benchmark
  stations = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35] # km

  # Function that finds the depth-index corresponding to a station location
  function find_station_index(stations, grid_points)
      numstations = length(stations)
      station_ind = zeros(numstations)
      for i in range(1, stop=numstations)
        station_ind[i] = argmin(abs.(grid_points .- stations[i]))
        station_ind[i]
      end
      return Integer.(station_ind)
    end
    
  
  station_indices = find_station_index(stations, z[1,:])
  station_strings = ["000", "025", "005", "075", "100", "125", "150", "175", "200", "250", "300", "350"] # "125" corresponds to 12.5 km down dip; these are necessary for writing to files


  # Benchmark description also requests slip evolution output recorded at fault locations every ~stride_space spatial nodes:
  z0 = z[1,:] # TODO: this 
  flt_loc = z0[1:stride_space:δNp]  # physical stations (units of km)
  flt_loc_indices = find_station_index(flt_loc, z[1,:])

  # set up parameters sent to the right hand side of the DAE:
  odeparam = (
              Vp=Vp,
              A = A,
              sJ = sJ,
              u=u,
              Δτ = Δτ,
              τf = τ0*ones(Ns+1),
              b = b,
              μshear=μshear,
              RSa=RSa,
              RSb=RSb,
              σn=σn,
              η=η,
              RSV0=RSV0,
              τ0=τ0,
              RSDc=RSDc,
              RSf0=RSf0,
              δNp = δNp,
              Ns = Ns,
              Nr = Nr,
              N = Ns, # TODO: remove
              B = B,
              T = T,
              e = e,
              save_stride_fields = stride_time, # save every save_stride_fields time steps
              sim_years =  sim_years
              )

  # Set time span over which to solve:
  tspan = (0, sim_years * year_seconds)


  # Set up ODE problem corresponding to DAE
  prob = ODEProblem(odefun, ψδ, tspan, odeparam)

  # Set call-back function so that files are written to after successful time steps only.
  cb_fun = SavingCallback((ψδ, t, i) -> write_to_file(pth, ψδ, t, i, z, flt_loc, flt_loc_indices,station_strings, station_indices, odeparam, "BP1_", 0.1 * year_seconds), SavedValues(Float64, Float64))

  # Make text files to store on-fault time series and slip data,
  # Also initialize with initial data:
  create_text_files(pth, flt_loc, flt_loc_indices, stations, station_strings, station_indices, 0, RSVinit, δ, τ0, θ)

  # Solve DAE using Tsit5()
  sol = solve(prob, Tsit5(); dt=1e-5,
           abstol = 1e-6, reltol = 1e-6, save_everystep=true, gamma = 0.8,
           internalnorm=(x, _)->norm(x, Inf), callback=cb_fun)

  
  return (sol, z, δNp, pth)
end



(S, z, δNp, pth) = main();




# example of how to plot slip contours (uncomment if desired):
# plot_slip(pth*"slip.dat")

# examples of how ot plot times series of shear stress:
# plot_fault_time_series("slip", pth*"fltst_strk000.txt")
# plot_fault_time_series("slip_rate", pth*"fltst_strk000.txt")
# plot_fault_time_series("shear_stress", pth*"fltst_strk+10.txt")
# plot_fault_time_series("state", pth*"fltst_strk+25.txt")