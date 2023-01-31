# Stripped down version of Thrase code (for training purposes)
# Solves the SEAS benchmark problem BP1-QD:
# https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf

using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
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
    #(M̃, F, τ, H̃, HfI_FT) = get_operators(SBPp, Nr, Ns, μ)

    # factor with Cholesky
    M = cholesky(Symmetric(M̃))

    # initialize vector g that stores boundary data
    t = 0
    g = zeros((Nx+1) * (Nz+1))

    δ = zeros(Nz+1)
    
    bc_Dirichlet = (lf, x, z) -> (2-lf) .* (0 * x .+ 0 .* z) + (lf-1) .* (0 .* x .+ 0 .* z)
    bc_Neumann   = (lf, x, z, nx, ny) -> zeros(size(x))
    
    bdry_vec_mod!(g, F, τ, x, z, bc_Dirichlet, bc_Neumann, Lx, Lz)
    
    u = M \ g
 
    (mm, δNp) = findmin(abs.(RSWf .- z))
    @assert z[δNp] ≈ RSWf
 
    Δτ = zeros(Nz+1)

    # Assemble fault variables/data
    RSa = zeros(δNp)
    for n = 1:δNp
        RSa[n] = RSamin - (RSamin - RSamax) *
            min(1, max(0, (RSH1 - z[n])/(RSH1 - RSH2)))
    end


    τz0 = σn * RSamax * asinh(RSVinit / (2 * RSV0) *
                                    exp.((RSf0 + RSb * log.(RSV0 / RSVinit)) /
                                        RSamax)) + η * RSVinit


    θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
        sinh.((τz0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)


    ψ0 = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)


    ψδ = zeros(δNp + Nz + 1)  #because length(ψ) = δNp,  length(δ) = N+1
    ψδ[1:δNp] .= ψ0


    stations = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35] # km

    # Next part of code related to code output
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
    station_strings = ["000", "025", "005", "075", "100", "125", "150", "175", "200", "250", "300", "350"] # "125" corresponds to 12.5 km down dip. 
  
  
    # Fault locations to record slip evolution output:
    flt_loc = z[1:stride_space:δNp]  # physical stations (units of km)
    flt_loc_indices = find_station_index(flt_loc, z)

    # set up parameters sent to right hand side
    odeparam = (reject_step = [false],
                Vp=Vp,
                M = M,
                u=u,
                Δτ = Δτ,
                τf = τz0*ones(Nz+1),
                g = g,
                μshear=μshear,
                RSa=RSa,
                RSb=RSb,
                σn=σn,
                η=η,
                RSV0=RSV0,
                τz0=τz0,
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

    #dψV = zeros(δNp + Nz + 1)
    tspan = (0, sim_years * year_seconds)
    prob = ODEProblem(odefun, ψδ, tspan, odeparam)
 
    cb_fun = SavingCallback((ψδ, t, i) -> write_to_file(pth, ψδ, t, i, z, flt_loc, flt_loc_indices,station_strings, station_indices, odeparam, "BP1_", 0.1 * year_seconds), SavedValues(Float64, Float64))
  
    # Make text files to store on-fault time series and slip data,
    # Also initialize with initial data:
    create_text_files(pth, flt_loc, flt_loc_indices, stations, station_strings, station_indices, 0, RSVinit, δ, τz0, θ)

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