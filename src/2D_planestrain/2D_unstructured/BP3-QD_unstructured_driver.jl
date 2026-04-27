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

do_plotting = true # turn on if you want to plot the grid

const year_seconds = 31556926
global const ctr = Ref{Int64}(1) 

# parameters distinguishing RS fault from slowly creeping section
const RS_FAULT = 7
const VP_FAULT = 8

include("ops_BP3-QD_unstructured.jl")
include("odefun_BP3-QD_unstructured.jl")
include("../utils_2D.jl")

# Before running this script, type
# localARGS = ["./examples/bp3-qd_unstructured.dat"]
# on command line 

function main()
    ### input parameters
    (pth, meshfile, Nx, Ny, sim_years, Vp, psi, ρ, ν, cs, σn, 
    RSamin, RSamax, RSb, RSDc,
    RSf0, RSV0, RSVinit, RSH1,RSH2, RSWf, SBPp) = read_params_planestrain(localARGS[1])

    
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
   #(verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/"*meshfile)

   # 16 blocks total 
   Lx = 40
   Ly = 40
    x= [-Lx -.75*Lx -.5*Lx -.25*Lx 0 .25*Lx .5*Lx .75*Lx Lx]
    y = [-Ly -.75*Ly -.5*Ly -.25*Ly 0]
    
    old_verts = [x x x x x
            y[1]*ones(1,9) y[2]*ones(1,9) y[3]*ones(1,9) y[4]*ones(1,9) y[5]*ones(1,9)]
    
    verts = old_verts .+ 5 .* rand(2, 45)
    verts[:, 5] = old_verts[:, 5]
    verts[:, 14] = old_verts[:, 14]
    verts[:, 23] = old_verts[:, 23]
    verts[:, 32] = old_verts[:, 32]
    verts[:, 41] = old_verts[:, 41]
    # note the tranpose here:
    EToV = [1 2 10 11 
        2 3 11 12
        3 4 12 13
        4 5 13 14
        5 6 14 15
        6 7 15 16
        7 8 16 17
        8 9 17 18
        10 11 19 20
        11 12 20 21
        12 13 21 22 
        13 14 22 23
        14 15 23 24
        15 16 24 25
        16 17 25 26
        17 18 26 27
        19 20 28 29
        20 21 29 30
        21 22 30 31
        22 23 31 32
        23 24 32 33
        24 25 33 34
        25 26 34 35
        26 27 35 36
        28 29 37 38
        29 30 38 39
        30 31 39 40
        31 32 40 41
        32 33 41 42
        33 34 42 43
        34 35 43 44
        35 36 44 45]'
    
    # note transpose:
        EToF = [9 10 1 69
        10 11 2 70
        11 12 3 71
        12 13 4 72
        13 14 5 73
        14 15 6 74
        15 16 7 75
        16 17 8 76
        18 19 69 27
        19 20 70 28
        20 21 71 29
        21 22 72 30
        22 23 73 31
        23 24 74 32
        24 25 75 33
        25 26 76 34
        35 36 27 44
        36 37 28 45
        37 38 29 46
        38 39 30 47
        39 40 31 48
        40 41 32 49
        41 42 33 50
        42 43 34 51
        52 53 44 61
        53 54 45 62
        54 55 46 63
        55 56 47 64
        56 57 48 65
        57 58 49 66
        58 59 50 67
        59 60 51 68]'

        # first row ends on face 26, second row ends at 60
    FToB = [2 2 2 2 2 2 2 2 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0]'
    FToB[13] = 8
    FToB[22] = 8
    FToB[39] = 7
    FToB[56] = 7
    EToDomain = [1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2]
    #    verts =  4 .* [-20 -20 20 20 0 0;
    #             0 -20 -20 0 0 -20]

    # EToV = [2 6;
    #         6 3;
    #         1 5;
    #         5 4]            
    
    # EToF = 
    # [4 3;
    #  3 7;
    #  1 6;
    #  2 5]

    #  FToB = [2, 2, 7, 1, 2, 2, 1]

    #  EToDomain = [1, 2]

    #verts = [-20   0  -20 0  20 20 -20   0  20; 
          #   -20 -20    0 0 -20  0 -40 -40 -40];

#     xd = 20*cosd(psi)
#     yd = 20*sind(psi)
#     verts = [-20+xd  0+xd    -20  0  20+xd      20  -20+2*xd      0+2*xd     20+2*xd; 
#                 -yd  -yd       0    0    -yd     0  -2*yd  -2*yd  -2*yd];

#     FToB = [1; 7; 1; 2; 2; 0; 0; 1; 8; 1; 2; 2];
# #    #FToB = [1; 1; 1; 1; 1; 1; 1];


#     EToV = [1 7 2 8;
#     2 8 5 9;
#     3 1 4 2;
#     4 2 6 5]

#     EToF = [1 8 2 9;
#     2 9 3 10;
#     6 11 7 12;
#     4 6 5 7]


#    EToDomain = [1 1 2 2]

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

    ######### create local operators on each block/element by applying coordinate transform
    # Create an empty dictionary to store the operators;
    # index via element number (an integer), return structure containing local ops.

    l = locoperator(2, 16, 16, exact_mu, exact_lambda)
    OPTYPE = typeof(l)
    lop = Dict{Int64, OPTYPE}()

    # create an empty dictionary to store each element's neighbor's penalty params:
    NTYPE = typeof([l[12][1], l[12][2], l[12][3], l[12][4]])
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


    # Assemble the global volume operator and compute LU factorization:
    A = global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
    A = lu(A)

    # Get unique array indices for the faces corresponding to the fault/jump interface
    FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT), Nr, Ns)


    # Calculate interface jump displacement penalty parameters:
    for e = 1:nelems
        LFToB = FToB[EToF[:,e]]

        nz = [similar(lop[e].IsJZ[1]), similar(lop[e].IsJZ[2]), similar(lop[e].IsJZ[3]), similar(lop[e].IsJZ[4])]
            
        # loop over the faces to get the neighboring penalty parameters:
        for lf = 1:4   
            if LFToB[lf] == BC_LOCKED_INTERFACE || LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
                f = EToF[lf, e]  # get global face number
                (em, ep) = FToE[:, f]  # find the two elements that share global face f.
                (fm, fp) = FToLF[:, f]

                if em == e
                    eo = ep 
                    nf = fp
                else
                    eo = em 
                    nf = fm
                end
                
                nz[nf] .= lop[eo].IsJZ[nf] # need to store the Z from the other element's local face
                
            else  
            end
        
        end
        
        neighborZ[e] = nz
    end
    ############ END COORDINATE TRANSFORM

    # Compute the total number of volume and jump (δ) points
    VNp = vstarts[nelems+1]-1
    δNp = FToδstarts[nfaces+1]-1

    # Set boundary data functions:
    #{{{ Compute the boundary/interface functions:
    function creep(x,y,t, e, EToDomain)
        if EToDomain[e] == 1 # left hand side of fault
            return [(Vp/2) .* t * cosd(psi) .+ 0 .* x .+ 0 .* y, (Vp/2) .* t * sind(psi) .+ 0 .* x .+ 0 .* y]
        elseif EToDomain[e] == 2
            return [-(Vp/2) .* t * cosd(psi) .+ 0 .* x .+ 0 .* y, -(Vp/2) .* t * sind(psi) .+ 0 .* x .+ 0 .* y]
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
        if EToO[lf, e]     # check if on minus side, A&D define δ as minus side minus plus side
          return δ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else
          error("shouldn't get here")  # this is because "correct" orientation is always true of a face on the minus side. 
        end
      else                  # face on plus side, add minus sign
        if EToO[lf, e]      # check if orientation is correct
          return  -δ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else                # if orientation is reversed, reverse the data
          return  -δ[(FToδstarts[f+1]-1):-1:FToδstarts[f], :]
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
        loc_bdry_vec_v2!((@view b[vstarts[e]:vstarts[e+1]-1, :]), lop[e], neighborZ[e], FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, (e, δ, t, EToDomain))
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
    for f = 1:nfaces
    if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (e1, _) = FToE[:, f]
        (lf1, _) = FToLF[:, f]
        xf = lop[e1].facecoord[1][lf1]
        yf = lop[e1].facecoord[2][lf1]
        xdf = yf ./ sind(psi)
        δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
        for n = 1:length(δrng)
        RSa[δrng[n]] = RSamin - (RSamin - RSamax) *
            min(1, max(0, (RSH1 + xdf[n])/(RSH1 - RSH2)))
        
        fault_nodes[δrng[n]] = xdf[n]
        end
    end
    end
   # plt.scatter(fault_nodes, RSa)
   
    # Set pre-stress according to benchmark description
    τ0 = fill(σn * RSamax * asinh(RSVinit / (2 * RSV0) *
                                exp.((RSf0 + RSb * log.(RSV0 / RSVinit)) /
                                        RSamax)) + η * RSVinit,
            δNp)
    σ0 = σn .* ones(δNp)

    # for f = 1:nfaces
    #     if FToB[f] == RS_FAULT
    #         (e1, e2) = FToE[:, f]
    #         (lf1, lf2) = FToLF[:, f]
    #         nx = lop[e1].nx
    #         δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
    #         for n = 1:length(δrng)
    #             τ0[δrng[n]] = sign(nx[lf1][n])*abs(τ0[δrng[n]]) # TODO WTF is this for.
    #         end
    #     end
    # end
   
    # Set initial state variable according to benchmark           
    θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
        sinh.((τ0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)
    # Initialize psi version of state variable    
    ψ = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)
    # Set initial velocity (this cannot be set independently, it must be consistent with θ and ψ)
    V_0 = Vp .* ones(δNp)


    # Calculate total shear/normal stress at t = 0:
    τ = τ0 + Δτ
    σ = σ0 + Δσ

    # Set initial condition for index 1 DAE - this is a stacked vector containing ψ, followed by slip (two components)
    ψδ = zeros(3δNp)
    ψδ[1:δNp] .= ψ              # state variable
    ψδ[δNp+1:2δNp] .= δ[:, 1]   # horizontal jump
    ψδ[2δNp+1:3δNp] .= δ[:, 2]  # vertical jump

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

     
    stations = setupfaultstations(stations_locations, lop, FToB, FToE, FToLF,
                                (RS_FAULT, VP_FAULT))

    fault = setupfaultcoord(lop, FToB, FToE, FToLF,
                                (RS_FAULT, VP_FAULT)) 




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
    σn=σn,
    η=η,
    RSV0=RSV0,
    τ0=τ0,
    σ0 = σ0, 
    RSDc=RSDc,
    RSf0=RSf0,
    fault_nodes = fault_nodes,
    psi = psi
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
         
    @show Nr[1]
    # Set call-back function so that fields are written to text file after successful time step only.
    cb = SavingCallback((ψδ, t, i)->savedatafields(ψδ, t, i, stations, fault, V_0, 
                                                    FToδstarts, odeparam,
                                                    pth*"BP3_N_$(Nr[1])_", pth*"Slip_BP3_N_$(Nr[1])_",  10year_seconds),
                        SavedValues(Float64, Float64))
 
    # Solve DAE using Tsit5(), an adaptive Runge-Kutta method
    sol = solve(prob, Tsit5(); isoutofdomain=stepcheck, dt=1,
                abstol = 1e-6, reltol = 1e-6, save_everystep=false, #gamma = 0.05,
                internalnorm=(x, _)->norm(x, Inf), callback=cb)
                
                
    return sol
end
        

S = main();

# example of how to plot slip contours (uncomment if desired):
# plot_slip(pth*"slip.dat")

# examples of how ot plot times series of shear stress:
#pth = "BP3_N_40_0.0_-7.5.dat"
#plot_fault_time_series("shear_stress", pth)
# plot_fault_time_series("slip_rate", pth*"fltst_strk000.txt")
# plot_fault_time_series("shear_stress", pth*"fltst_strk+10.txt")
# plot_fault_time_series("state", pth*"fltst_strk+25.txt")