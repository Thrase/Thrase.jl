# On an unstructured mesh of quads:
# Solves a manufactured solution
# Geometry set up for BP9: dipping fault,  friction fault down dip to xd = 40km, then steady sliding
# 30 degree thrust fault.

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
include("odefun_BP9-QD_unstructured.jl")
include("../utils_2D.jl")

# Before running this script, type
# localARGS = ["./examples/bp9-qd_unstructured.dat"]
# on command line 

function main()
    ### input parameters
    (pth, meshfile, Nx, Ny, sim_years, Vp, psi, ρ, ν, cs, σmin, zh, 
    RSamin, RSamax, RSb, RSDc,
    RSf0, RSV0, RSVinit, RSH ,RSHt, RSw, RSWf, SBPp) = read_params_planestrain(localARGS[1])

    ρw = 1
    g = 9.8

    try
        mkdir(pth)
    catch
        # folder already exists so make a new one.
        pth = pth*string(now())*"/"
        mkdir(pth)
    end

    # write input file to output folder
    source_pth = localARGS[1]
    idx = findlast('/', source_pth)
    filename = source_pth[idx+1:end]  
    destination_pth = joinpath(pth, filename)
    cp(source_pth, destination_pth, force=true)


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
#   verts = [-400 0 400 -400 0 400;
#               -800 -800 -800 0 0 0];

#                   EToV = [1 2; 2 3; 4 5; 5 6];
#     EToF = [1 2; 
#             2 3;
#             6 7;
#             4 5];


#     FToB = [1; 7; 1; 3; 3; 2; 2];
#    EToDomain = [1 2]

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
        p = Plots.scatter(verts[1,:], verts[2,:], markersize=2, markercolor=:blue, legend = false, size=(1600, 600))
        for e = 1:nelems
            V = EToV[:, e]
            V[3], V[4] = V[4], V[3]
            push!(V, V[1])
            Plots.plot!(verts[1, V], verts[2, V],linecolor=:red, legend = false, size=(1600, 600))
        end
        display(p)
    end


    # Secondary Grid Arrays:
    (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)


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

    # initialize backgrounds stresses and changes in stress due to quasi-static deformation:
    Δτ = zeros(δNp)
    Δσ = zeros(δNp)
    σ0 = zeros(δNp)
    τ0 = zeros(δNp)
    p0 = zeros(δNp)

    # Assemble fault variables/data
    fault_nodes = zeros(δNp)
    RSa = zeros(δNp)

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
                
   
    # compute stuff on fault:
    for f = 1:nfaces
        if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (e1, _) = FToE[:, f]
        (lf1, _) = FToLF[:, f]
        xf = lop[e1].facecoord[1][lf1]
        yf = lop[e1].facecoord[2][lf1]
        xdf = yf ./ sind(psi)
        δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
        for n = 1:length(δrng)
        
            if 0 >= xdf[n] >= -RSH
                 RSa[δrng[n]]  = RSamin
            elseif -RSH > xdf[n] >= -(RSH + RSHt)
                RSa[δrng[n]] = (-(RSamax - RSamin)/RSHt) * xdf[n] + RSamin - RSH *(RSamax - RSamin)/RSHt
            else             
                RSa[δrng[n]]  = RSamax
            end
            fault_nodes[δrng[n]] = xdf[n]

            x3 = -yf[n]

            p0[δrng[n]] = ρw*g*x3 .* (x3 <= zh) .+ (ρw*g*zh + ρ*g*(x3-zh)) .* (x3 > zh)
            # σ11_0[δrng[n]] = -(3.12)*(ρ - ρw)*g*x3 - p0[δrng[n]]
            # σ22_0[δrng[n]] = -ρ*g*x3
            # σ33_0[δrng[n]] = -ρ*g*x3

            # Set pre-stress according to benchmark description. σ0 is EFFECTIVE (substracts off pore pressure)
                        
            # σ0[δrng[n]]  = -sind(psi)*sind(psi) .* σ11_0[δrng[n]]  .- cosd(psi)*cosd(psi) .* σ33_0[δrng[n]]  .- p0[δrng[n]]
            
           σ0[δrng[n]]  = ρ*g*x3 - p0[δrng[n]]
            
           τ0[δrng[n]]  = 0.53*σ0[δrng[n]]

            
            # τ0[δrng[n]]  = -sind(psi)*cosd(psi) .* (σ11_0[δrng[n]]  .- σ33_0[δrng[n]] )

                            
            # σ0[δrng[n]]  = ρ*g*x3 .- p0[δrng[n]]
            # τ0[δrng[n]]  = 0.5 .* σ0[δrng[n]]
           
           
        end
    end
    end

   
    # DELETE Set pre-stress according to benchmark description
    # τ0 = fill(50 * RSamax * asinh(RSVinit / (2 * RSV0) *
    #                             exp.((RSf0 + RSb * log.(RSV0 / RSVinit)) /
    #                                     RSamax)) + η * RSVinit,
    #         δNp)
    # σ0 = 50 .* ones(δNp)

  
    # plt.scatter(sind(psi) .* fault_nodes, σ0)
    # plt.scatter(sind(psi) .* fault_nodes, τ0)

#    plt.scatter(sind(psi) .* fault_nodes, σ33_0)
   #plt.scatter(fault_nodes, σ33_0)


    
    
     # Calculate total shear/normal stress at t = 0:
    σ = σ0 + Δσ
    τ = τ0 + Δτ
    
    # Set initial state variable according to benchmark           
    θ = (RSDc / RSV0) .* ones(δNp)
    #θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
    #    sinh.((τ0 .- η .* RSVinit) ./ (RSa .* 50))) .- RSf0 ./ RSb)
 
    # Initialize psi version of state variable    
    ψ = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)
    
    # Set initial down-dip slip rate (this cannot be set independently, it must be consistent with θ and ψ)
    V_0 = Vp .* ones(δNp)
    
    # Set initial condition for index 1 DAE - this is a stacked vector containing ψ, followed by slip (two components)
    ψδ = zeros(3δNp)
    ψδ[1:δNp] .= ψ              # state variable
    ψδ[δNp+1:2δNp] .= δ[:, 1]   # fault parallel jump
    ψδ[2δNp+1:3δNp] .= δ[:, 2]  # fault normal jump

    # Set fault station locations every 2.5km down-dip) specified in benchmark (x_1, x_3); x_2 = 0 in 2D problem
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

     
    stations = setupfaultstations(stations_locations, lop, FToB, FToE, FToLF,
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




    # Length scales that need to be resolved by grid:
    Λ0 = μ .* RSDc ./ (RSb .* σ)
    hstar = (π/2) .* μ .* RSb .* RSDc ./ (σ .* (RSb .- RSa) .^ 2)
   
    # set up parameters sent to the right hand side of the DAE:
    odeparam = (reject_step = [false],
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
    FToFREEstarts = FToFREEstarts,
    EToDomain = EToDomain,
    b = b,
    u=u,
    w = w,
    τ=τ,
    σ=σ,
    Δτ = Δτ,
    Δσ = Δσ,
    vstarts=vstarts,
    μshear=μshear,
    RSa=RSa,
    RSb=RSb,
    σmin=σmin,
    zh=zh,
    η=η,
    RSV0=RSV0,
    τ0=τ0,
    σ0 = σ0, 
    RSDc=RSDc,
    RSf0=RSf0,
    fault_nodes = fault_nodes,
    psi = psi,
    u1 = u1,
    u3 = u3
    )

    # Set time span over which to solve the DAE:
    tspan = (0, sim_years * year_seconds)

    # Set up right-hand-side corresponding to DAE
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
    cb = SavingCallback((ψδ, t, i)->savedatafields(ψδ, t, i, stations, offfault_stations, fault, surface_coords, V_0, 
                                                    FToδstarts, odeparam,
                                                    pth*"BP9_N_$(Nr[1])_", pth*"off_fault_BP9_N_$(Nr[1])_", pth*"Slip_BP9_N_$(Nr[1])_",  10year_seconds),
                        SavedValues(Float64, Float64))
 
    # Solve DAE using Tsit5(), an adaptive Runge-Kutta method
    sol = solve(prob, Tsit5(); isoutofdomain=stepcheck, dt=1,
                abstol = 1e-7, reltol = 1e-7, save_everystep=false, #gamma = 0.05,
                internalnorm=(x, _)->norm(x, Inf), callback=cb)
                
                
    return sol
end
        

S = main();

# example of how to plot slip contours (uncomment if desired):
# pth = "output/bp9-unstructured/2026-07-29T10:51:01.351/Slip_BP9_N_20_.dat"
# plot_slip(pth)

# examples of how ot plot times series of shear stress:
# pth = "output/bp9-unstructured/2026-07-29T10:51:01.351/BP9_N_20_0.0_0.0.dat"
# plot_fault_time_series("slip_rate", pth)
