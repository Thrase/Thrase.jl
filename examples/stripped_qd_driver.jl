# Stripped down version of Thrase code (for training purposes)
# Solves the SEAS benchmark problem BP1-QD:
# https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf

using Thrase
using LinearAlgebra
using DifferentialEquations
using DelimitedFiles



function main()


    ### input parameters
    (pth, stride_space, stride_time, xc, zc, Nx, Nz,
    sim_years, Vp, ρ, cs, σn, RSamin, RSamax, RSb, RSDc,
    RSf0, RSV0, RSVinit, RSH1,RSH2, RSWf, SBPp) = read_params(localARGS[1])

    try
        mkdir(pth)
    catch
        # folder already exists and data will be overwritten.
    end
    
    year_seconds = 31556926
    μ = cs^2 * ρ 
    μshear = cs^2 * ρ
    η = μshear / (2 * cs)

   

    x = Array(LinRange(xc[1], xc[2], Nx+1))  
    z = Array(LinRange(zc[1], zc[2], Nz+1))

    Lx = xc[2]
    Lz = zc[2]

    # create operators
    (M̃, F, τ, H̃, HfI_FT) = get_operators(SBPp, Nx, Nz, μ; xc = xc, zc = zc)

    # factor matrix with Cholesky
    A = cholesky(Symmetric(M̃))

    # initialize time and vector b that stores boundary data (linear system will be Au = b)
    t = 0
    b = zeros((Nx+1) * (Nz+1))

    # initial slip vector
    δ = zeros(Nz+1)
    
    # fill in initial boundary data into b
    bdry_vec_strip!(b, F, τ, x, z, δ ./ 2, (t .* Vp./2)*ones(size(z)), zeros(size(x)), Lx, Lz)
    
    u = A \ b # solve linear system with a backsolve to obtain initial displacement u(x, z, 0)
 
    # Find z-index δNp corresponding to base of rate-and-state friction zone RSWf
    (mm, δNp) = findmin(abs.(RSWf .- z))
    @assert z[δNp] ≈ RSWf
 
    # initialize change in shear stress due to quasi-static deformation
    Δτ = zeros(Nz+1)

    # Assemble fault variables/data
    RSa = zeros(δNp)
    for n = 1:δNp
        RSa[n] = RSamin - (RSamin - RSamax) *
            min(1, max(0, (RSH1 - z[n])/(RSH1 - RSH2)))
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
    ψδ = zeros(δNp + Nz + 1)  #because length(ψ) = δNp,  length(δ) = Nz+1
    ψδ[1:δNp] .= ψ
    ψδ[δNp+1:δNp + Nz + 1] .= δ

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
      
    station_indices = find_station_index(stations, z)
    station_strings = ["000", "025", "005", "075", "100", "125", "150", "175", "200", "250", "300", "350"] # "125" corresponds to 12.5 km down dip; these are necessary for writing to files
  
  
    # Benchmark description also requests slip evolution output recorded at fault locations every ~stride_space spatial nodes:
    flt_loc = z[1:stride_space:δNp]  # physical stations (units of km)
    flt_loc_indices = find_station_index(flt_loc, z)

    # set up parameters sent to the right hand side of the DAE:
    odeparam = (reject_step = [false],
                Vp=Vp,
                A = A,
                u=u,
                Δτ = Δτ,
                τf = τ0*ones(Nz+1),
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
                N = Nz,
                F = F,
                τ = τ,
                x = x ,
                z = z,
                HfI_FT = HfI_FT,
                Lx = Lx,
                Lz = Lz,
                save_stride_fields = stride_time # save every save_stride_fields time steps
                )

    # Set time span over which to solve:
    tspan = (0, sim_years * year_seconds)


    # Set up ODE problem corresponding to DAE
    prob = ODEProblem(odefun_stripped, ψδ, tspan, odeparam)
 
    # Set call-back function so that files are written to after successful time steps only.
    cb_fun = SavingCallback((ψδ, t, i) -> write_to_file(pth, ψδ, t, i, z, flt_loc, flt_loc_indices,station_strings, station_indices, odeparam, "BP1_", 0.1 * year_seconds), SavedValues(Float64, Float64))
  
    # Make text files to store on-fault time series and slip data,
    # Also initialize with initial data:
    create_text_files(pth, flt_loc, flt_loc_indices, stations, station_strings, station_indices, 0, RSVinit, δ, τ0, θ)

    # Solve DAE using Tsit5()
    sol = solve(prob, Tsit5(); dt=0.2,
             abstol = 1e-5, reltol = 1e-5, save_everystep=true, gamma = 0.2,
             internalnorm=(x, _)->norm(x, Inf), callback=cb_fun)

    
    return (sol, z, δNp, pth)
end



(S, z, δNp, pth) = main();




# example of how to plot slip contours (uncomment if desired):
# plot_slip(pth*"slip.dat")

# examples of how ot plot times series of shear stress:
# plot_fault_time_series("slip", pth*"fltst_strk000.txt")
# plot_fault_time_series("V", pth*"fltst_strk+10.txt")
# plot_fault_time_series("shear_stress", pth*"fltst_strk+10.txt")
# plot_fault_time_series("state", pth*"fltst_strk+25.txt")