# On an unstructured mesh of quads:
# Solves a manufactured solution
# Geometry set up for BP3: dipping fault,  friction fault down dip to xd = 40km, then steady sliding
# this example is for a 30 degree thrust fault.

using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using DelimitedFiles
using Dates
using Plots

do_plotting = false # turn on if you want to plot the grid

const year_seconds = 31556926
global const ctr = Ref{Int64}(1) 

# parameters distinguishing RS fault from slowly creeping section
const RS_FAULT = 7
const VP_FAULT = 8
const EARTH_FREE = 3

include("ops_BP3-QD_unstructured.jl")
include("odefun_BP3fluid-QD_unstructured.jl")
include("../utils_2D.jl")
include("diffusions.jl")

# Before running this script, type
# localARGS = ["./examples/bp3-qd_unstructured.dat"]
# on command line 

function main()
    use_stage_limiter = false   # set false to skip pressure/permeability correction

    ### input parameters
    (pth, meshfile, Nx, Ny, sim_years, Vp, psi, ρ, ν, cs, σn_max, 
    RSamin, RSamax, RSbmin, RSbmax, RSDc,
    RSf0, RSV0, RSVinit, RSH1,RSH2, RSw, RSWf, SBPp, 
    kmin, kmax, T, L, σstar, q0, β, η_f, n, ρ_f, g) = read_params_planestrain_fluids(localARGS[1])

    
    try
        mkdir(pth)
    catch
        # folder already exists so make a new one.
        pth = pth*string(now())*"/"
        mkdir(pth)
    end

    # Define shear modulus μ and radiation damping η
    μ = cs^2 * ρ 
    μshear = cs^2 * ρ
    η = μ / (2 * cs)

    function exact_mu(x, y)
        return μ .+ 0 .* x
    end

    function exact_lambda(x, y)
        return 2 .* μ .* ν ./ (1 .- 2 .* ν) .+ 0 .* x
    end

    # Read in the unstructured mesh input file:
    (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/"*meshfile)


    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    @show (nelems, nfaces)

    # TODO: add assertions that lengthscales are sufficiently resolved.

    # Set grid points in each direction on each element in logical space - Nx and Ny read from input file.
    # currently we use the same number of nodes on each element.
    EToN0 = zeros(Int64, 2, nelems)
    EToN0[1, :] .= Nx
    EToN0[2, :] .= Ny 
    Nr =  EToN0[1, :] 
    Ns = EToN0[2, :]
    
    # Plot the grid
    if do_plotting
        p = Plots.scatter(verts[1,:], verts[2,:], markersize=2, markercolor=:blue, legend = false)
        for e = 1:nelems
            V = EToV[:, e]
            V[3], V[4] = V[4], V[3]
            push!(V, V[1])
            Plots.plot!(verts[1, V], verts[2, V],linecolor=:red, legend = false)
        end
        display(p)
    end

    # Secondary Grid Arrays:
    (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

    # Fault Arrays: 

    (fverts, fEToV, fEToF, fEToDomain, fFToB, Nf, fflip) = faultarrays(verts, EToV, EToF, FToB, Nr)
    (fFToE, fFToLF) = connectivityarrays_1D(fEToV, fEToF)
    
    F = (fverts, fEToV, fEToF, fEToDomain, fFToB, fFToE, fFToLF, Nf, fflip) #combine into array


    ######### create local operators on each block/element by applying coordinate transform
    # Create an empty dictionary to store the operators;
    # index via element number (an integer), return structure containing local ops.

    l = locoperator(2, 16, 16, exact_mu, exact_lambda)
    OPTYPE = typeof(l)
    lop = Dict{Int64, OPTYPE}()

    # create an empty dictionary to store each element's neighbor's penalty params:
    NTYPE = typeof([l[11][1], l[11][2], l[11][3], l[11][4]])
    neighborZ = Dict{Int64, NTYPE}()

    # Initialize array of indices corresponding to vertices as you move through blocks/elements
    vstarts = Array{Int64, 1}(undef, nelems + 1)
    vstarts[1] = 1   # start at 1. 
  
    # Loop over blocks/elements and create local operators:
    for e = 1:nelems

        Np = (Nr[e]+1)*(Ns[e]+1)  # total number of volume points on each element (these include boundary/interface nodes)
        vstarts[e+1] = vstarts[e] + Np # fill in vstarts

        # DOMAIN: ordered pairs that define the physical domain
        (x1,x2,x3,x4) = verts[1, EToV[:, e]] 
        (y1,y2,y3,y4) = verts[2, EToV[:, e]]

        # Details of transfinite interpolation: straight edges:
        # x-parameterization on each face
        ex = [(α) -> x1 * (1 .- α) / 2 + x3 * (1 .+ α) / 2,
        (α) -> x2 * (1 .- α) / 2 + x4 * (1 .+ α) / 2,
        (α) -> x1 * (1 .- α) / 2 + x2 * (1 .+ α) / 2,
        (α) -> x3 * (1 .- α) / 2 + x4 * (1 .+ α) / 2]

        # derivative of x-parameterization
        exα = [(α) -> -x1 / 2 + x3 / 2,
        (α) -> -x2 / 2 + x4 / 2,
        (α) -> -x1 / 2 + x2 / 2,
        (α) -> -x3 / 2 + x4 / 2]
        
        # y-parameterization on each face
        ey = [(α) -> y1 * (1 .- α) / 2 + y3 * (1 .+ α) / 2,
        (α) -> y2 * (1 .- α) / 2 + y4 * (1 .+ α) / 2,
        (α) -> y1 * (1 .- α) / 2 + y2 * (1 .+ α) / 2,
        (α) -> y3 * (1 .- α) / 2 + y4 * (1 .+ α) / 2]

        # derivative of y-parameterization
        eyα = [(α) -> -y1 / 2 + y3 / 2,
        (α) -> -y2 / 2 + y4 / 2,
        (α) -> -y1 / 2 + y2 / 2,
        (α) -> -y3 / 2 + y4 / 2]

        # Coordinate transformation with transfinite interpolation
        xt(r,s) = transfinite_blend(ex[1], ex[2], ex[3], ex[4], exα[1], exα[2], exα[3], exα[4], r, s)
        zt(r,s) = transfinite_blend(ey[1], ey[2], ey[3], ey[4], eyα[1], eyα[2], eyα[3], eyα[4], r, s)

     
        metrics = create_metrics(Nr[e], Ns[e], exact_mu, exact_lambda, xt, zt) # create coordinate transform


        ###################################################################### 
        # create local finite difference operators on computational domain:
        lop[e] = locoperator(SBPp, Nr[e], Ns[e], exact_mu, exact_lambda, metrics, FToB[EToF[:, e]]) 
    end

    # Calculate neighboring penalty parameters and add to lop 
    for e = 1:nelems
        lop[e].neighborZ  .= calculate_neighbors(lop, e, FToB, EToF, FToE, FToLF, EToO)
    end
    
    # Assemble the global volume operator and compute LU factorization:
    A = global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
    A = lu(A)

    # Get unique array indices for the faces corresponding to the fault/jump interface
    FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT), Nr, Ns)

    # Get unique array indices for the faces corresponding to Earth's free surface
    FToFREEstarts = bcstarts(FToB, FToE, FToLF, EARTH_FREE, Nr, Ns)

    ############ END COORDINATE TRANSFORM

    # Compute the total number of volume and jump (δ) points
    VNp = vstarts[nelems+1]-1
    δNp = FToδstarts[nfaces+1]-1
    FREENp = FToFREEstarts[nfaces+1]-1

    # Set boundary data functions:
    #{{{ Compute the boundary/interface functions:
    function creep(x,y,t, e, EToDomain)
        if EToDomain[e] == 1 # left hand side of fault
            return [(Vp/2) .* t * cosd(psi) .+ 0 .* x .+ 0 .* y, -(Vp/2) .* t * sind(psi) .+ 0 .* x .+ 0 .* y]
        elseif EToDomain[e] == 2
            return [-(Vp/2) .* t * cosd(psi) .+ 0 .* x .+ 0 .* y, (Vp/2) .* t * sind(psi) .+ 0 .* x .+ 0 .* y]
        else
            error("shouldn't get here")
        end
            
    end
    
    bc_Dirichlet = (lf, x, y, e, δ, t, EToDomain) -> creep(x,y,t,e,EToDomain)
    bc_Neumann   = (lf, x, y, nx, ny, e, δ, t, EToDomain) -> [zeros(size(x)), zeros(size(x))]

   
    # define the interface jump in displacement function
    in_jump      = (lf, x, y, e, δ, t, EToDomain) -> begin
      f = EToF[lf, e]  # Get global face number
      if EToS[lf, e] == 1  # check if face is on minus side
        if EToO[lf, e]     # check if on minus side, A&D define δ as minus side minus plus side, i do opposite.
          return -δ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else
          error("shouldn't get here")  # this is because "correct" orientation is always true of a face on the minus side. 
        end
      else                  # face on plus side, add minus sign
        if EToO[lf, e]      # check if orientation is correct
          return  δ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else                # if orientation is reversed, reverse the data
          return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f], :]
        end
      end
    end
    #}}}

    # initial time and rhs b of linear system (i.e. u = A\b)
    t = 0
    b = zeros(VNp, 2) # horizontal, vertical component contributions

    # initial slip vector
    δ = zeros(δNp, 2) # fault parallel (slip) followed by fault normal (opening)
  
  
      # fill in initial boundary data into b
    for e = 1:nelems
        loc_bdry_vec_v2!((@view b[vstarts[e]:vstarts[e+1]-1, :]), lop[e], FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, (e, δ, t, EToDomain))
    end
   
   
    U = A \ b[:] # solve linear system with a backsolve to obtain initial displacements u, w
    u = U[1:VNp]
    w = U[VNp+1:end]

    # initialize change in stress due to quasi-static deformation:
    Δτ = zeros(δNp)
    Δσ = zeros(δNp)

    # Assemble fault variables/data
    fault_nodes = zeros(δNp)
    RSa = zeros(δNp)
    RSb = zeros(δNp)
    ne = [0] # number of faces in 2D domain that are faults (i.e. number of 1D "elements")
    for f = 1:nfaces
    if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        ne[1] += 1
        (e1, _) = FToE[:, f]
        (lf1, _) = FToLF[:, f]
        xf = lop[e1].facecoord[1][lf1]
        yf = lop[e1].facecoord[2][lf1]
        xdf = yf ./ sind(psi)
        δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
        # a: piecewise linear with kink at RSH1
        #   RSamin at 0, RSamid at RSH1, RSamax at RSH1+RSH2
        # a-b: constant (RSamin-RSbmin) on [0,RSH1], then linear to (RSamax-RSbmax) at RSH1+RSH2
        # b = a - (a-b), continuous everywhere
        RSamid    = 0.035             # RSa at -RSH1 kink depth
        f_shallow = RSamin - RSbmin   # a-b on [0, RSH1]
        f_deep    = RSamax - RSbmax   # a-b at RSH1+RSH2
        for n = 1:length(δrng)
        x = abs(xdf[n])  # positive down-dip distance (km)
        if x <= RSH1
            RSa[δrng[n]] = RSamin + (RSamid - RSamin) / RSH1 * x
            f = f_shallow
        else
            RSa[δrng[n]] = RSamid + (RSamax - RSamid) / RSH2 * (x - RSH1)
            f = f_shallow + (f_deep - f_shallow) / RSH2 * (x - RSH1)
        end
        RSb[δrng[n]] = RSa[δrng[n]] - f
        fault_nodes[δrng[n]] = xdf[n]
        end
    end
    end

    #plt.scatter(fault_nodes, RSa.-RSb)
    #plt.scatter(fault_nodes, RSa)
   
    # initial surface displacements
    u1 = zeros(FREENp, 1) # horizontal
    u3 = zeros(FREENp, 1) # vertical
  

    # fill in surface displacements:
    compute_surface_displacements!(u1, u3, u, w, nfaces, FToB, FToE, FToLF, vstarts, lop, FToFREEstarts)

    surface_coords = zeros(FREENp)
    for f = 1:nfaces 
        if FToB[f] ∈ EARTH_SURFACE
            (e1, _) = FToE[:, f]
            (lf1, _) = FToLF[:, f]

            xs = lop[e1].facecoord[1][lf1]
            zs = lop[e1].facecoord[2][lf1]
            
           
                Frng = FToFREEstarts[f]:(FToFREEstarts[f+1]-1)
                for n = 1:length(Frng)
                      surface_coords[Frng[n]] = xs[n]

                end
        
        end
    end
                
   


    #plt.scatter(fault_nodes, RSa)


    # effective stress = normal_stress - P = intial normal stress + stress change - initial pressure - pressure change
    
    # σ0 if the initial effective normal stress
    # σ0: grows from σ0_surf at surface to σn_max at -RSH1 via tanh, then constant.
    # tanh(3x/RSH1)/tanh(3) is 0 at x=0, 1 at x=RSH1; min() caps it flat beyond RSH1.
    σn0_surf = 1                        # value at the surface (x = 0)
    σn0 = zeros(δNp)                      # THIS IS initial EFFECTIVE NORMAL STRESS
    σ0 = zeros(δNp)                       # initial normal strress
    for n = 1:δNp
        x = abs(fault_nodes[n])          # positive down-dip distance (km)
        σn0[n] = σn0_surf + (σn_max - σn0_surf) * min(tanh(3.0 * x / RSH1) / tanh(3.0), 1.0)
        σ0[n] = ρ * g * x
    end
    
    P = σ0 - σn0 # initial condition for pressure: total minus effective
    # plt.scatter(fault_nodes, P)
    # plt.scatter(fault_nodes, σn0)
    # plt.scatter(fault_nodes, σ0)
    
     hstar = μ .* RSDc .* RSb ./ (σn0 .* (RSb .- RSa) .^ 2)
    #plt.scatter(fault_nodes, hstar)
#@show hstar

    # Set pre-stress according to benchmark description
    # Recall that σn is the initial effective normal stress, so p will stand for pressure changes
    τ0 = RSamin .* σn0 .* asinh.(RSVinit ./ (2 .* RSV0) .*
             exp.((RSf0 .+ RSbmin .* log.(RSV0 ./ RSVinit)) ./ RSamin)) .+ η .* RSVinit

    


    kstar = kmax .* ones(δNp) # initial condition for permeability target; not sure what this should be
    k = kmin .+ (kstar .- kmin) .* exp.(-(σ0 .- P) ./ σstar) # initial condition for permeability

   
    # Set initial velocity (this cannot be set independently, it must be consistent with θ and ψ)
    V_0 = RSVinit .* ones(δNp)


    # Calculate total shear and total eff. normal stress at t = 0:
    τ = τ0 + Δτ
    σn = σn0 + Δσ
    σ = σ0

    #plt.scatter(fault_nodes, τ)

      # Set initial state variable according to benchmark           
    # θ = (RSDc ./ RSVinit) .*  ones(δNp)
    # # Initialize psi version of state variable    
    # ψ = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)
   
    # Set initial state variable        
    θ = (RSDc ./ RSV0) .* exp.((RSamin ./ RSbmin) .* log.((2 .* RSV0 ./ RSVinit) .*
        sinh.((τ0 .- η .* RSVinit) ./ (RSamin .* σn))) .- RSf0 ./ RSbmin)
    # Initialize psi version of state variable    
    ψ = RSf0 .+ RSbmin .* log.(RSV0 .* θ ./ RSDc)
    # Set initial condition for index 1 DAE - this is a stacked vector containing ψ, 
    # followed by slip (two components),
    # followed by pressure, followed by permeability 
    ψδ = zeros(6δNp)
    ψδ[1:δNp]       .= ψ        # state variable
    ψδ[δNp+1:2δNp]  .= δ[:, 1]  # horizontal jump
    ψδ[2δNp+1:3δNp] .= δ[:, 2]  # vertical jump
    ψδ[3δNp+1:4δNp] .= P        # pore pressure
    ψδ[4δNp+1:5δNp] .= k        # permeability
    ψδ[5δNp+1:6δNp] .= kstar    # permeability target

    # Set fault station locations every 2.5km down-dip) specified in benchmark
    stations_locations = [0 0
                        2.5*cosd(psi) -2.5*sind(psi)
                        5*cosd(psi) -5*sind(psi)
                        7.5*cosd(psi) -7.5*sind(psi)
                        10*cosd(psi) -10*sind(psi)
                        12.5*cosd(psi) -12.5*sind(psi)
                        15*cosd(psi) -15*sind(psi)
                        17.5*cosd(psi) -17.5*sind(psi)
                        20*cosd(psi) -20*sind(psi)
                        25*cosd(psi) -25*sind(psi)
                        30*cosd(psi) -30*sind(psi)
                        ]

     
    stations = setupfaultstations_withfluids(stations_locations, lop, FToB, FToE, FToLF,
                                (RS_FAULT, VP_FAULT))

    fault = setupfaultcoord(lop, FToB, FToE, FToLF,
                                (RS_FAULT, VP_FAULT), psi) 


    
    
    # Setup off-fault station locations as specified in benchmark (x_1, x3), x_2 = 0 for 2D problem
    offfault_stations_locations = [8.0 0.0 
                         16.0 0.0
                         32.0 0.0
                         -8.0 0.0
                         -16.0 0.0
                         -32.0 0.0
                         8.0 10.0
                         16.0 10.0
                         32.0 10.0
                         -8.0 10.0
                         -16.0 10.0
                         -32.0 10.0
                        ]

     
    offfault_stations = setup_off_faultstations(offfault_stations_locations, lop, FToB, FToE, FToLF,
                                EARTH_FREE)



    # set up parameters sent to the right hand side of the DAE:
    odeparam = (reject_step = [false],
    t_current = [0.0],
    SBPp = SBPp,
    P_old = P,
    P = P,
    k = k,
    u1 = u1,
    u3 = u3,
    Vp=Vp,
    A = A,
    lop=lop,
    neighborZ = neighborZ,
    EToF=EToF,
    EToS=EToS,
    FToE=FToE,
    FToLF=FToLF,
    EToO=EToO,
    FToB=FToB,
    FToδstarts=FToδstarts,
    EToDomain = EToDomain,
    b = b,
    u=u,
    w = w,
    τ=τ,
    σn=σn,
    σ = σ,
    Δτ = Δτ,
    Δσ = Δσ,
    τ0 = τ0,
    σ0 = σ0,
    σn0 = σn0,
    vstarts=vstarts,
    μshear=μshear,
    RSa=RSa,
    RSb=RSb,
    η=η,
    RSV0=RSV0, 
    RSDc=RSDc,
    RSf0=RSf0,
    fault_nodes = fault_nodes,
    psi = psi,
    kmin = kmin,
    kmax = kmax,
    T = T,
    L = L,
    σstar = σstar,
    q0 = q0,
    β = β,
    η_f = η_f,
    n = n,
    ρ_f = ρ_f,
    g = g,
    F = F
    )


    # Set time span over which to solve the DAE:
    tspan = (0, sim_years * year_seconds)

    # Set up ODE problem corresponding to DAE
    prob = ODEProblem(odefun, ψδ, tspan, odeparam)
 
    function stepcheck(_, p, _)
        if p.reject_step[1]
        p.reject_step[1] = false
        println("reject")
        return true
        end
        return false
    end
         

    # Set call-back function so that fields are written to text file after successful time step only.
    cb = SavingCallback((ψδ, t, i)->savedatafields_with_fluids(ψδ, t, i, stations, offfault_stations, fault, surface_coords, V_0, 
                                                    FToδstarts, odeparam,
                                                    pth*"BP3_N_$(Nr[1])_", pth*"off_fault_B39_N_$(Nr[1])_", pth*"Slip_BP3_N_$(Nr[1])_",  10year_seconds),
                        SavedValues(Float64, Float64))
 
    # Solve DAE using Tsit5(), an adaptive Runge-Kutta method
    limiter! = use_stage_limiter ? correct_permeability_and_pressure! : (ψδ, i, p, t) -> nothing
    sol = solve(prob, Tsit5(stage_limiter! = limiter!); isoutofdomain=stepcheck, dt=1,
                abstol = 1e-6, reltol = 1e-6, save_everystep=false, #gamma = 0.05,
                internalnorm=(x, _)->norm(x, Inf), callback=cb)
                
                
    return sol
end
        

S = main();

# example of how to plot slip contours (uncomment if desired):
# pth = "output/bp3-unstructured/2026-07-29T10:51:01.351/Slip_BP3_N_20_.dat"
# plot_slip(pth)

# examples of how ot plot times series of shear stress:
# pth = "output/bp3-unstructured/2026-07-29T10:51:01.351/BP3_N_20_0.0_0.0.dat"
# plot_fault_time_series("slip_rate", pth)
