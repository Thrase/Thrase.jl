# On an unstructured mesh of quads:
# Solves a manufactured solution
# Geometry set up to resemble BP1 - so a vertical friction fault down to z = 40km, then steady sliding

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

include("ops_BP1-QD_unstructured.jl")
include("odefun_BP1-QD_unstructured.jl")
include("../utils_2D.jl")

# Before running this script, type
# localARGS = ["./examples/bp1-qd_unstructured.dat"]
# on command line 

function main()
    ### input parameters
    (pth, meshfile, Nx, Ny, sim_years, Vp, ρ, cs, σn, 
    RSamin, RSamax, RSb, RSDc,
    RSf0, RSV0, RSVinit, RSH1,RSH2, RSWf, SBPp) = read_params(localARGS[1])

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
        return μ 
    end

    # Read in the unstructured mesh input file:
    (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/"*meshfile)
    
    # Domain size in x-direction (used to specify loading conditions):
    Lx = maximum(verts[1,:])
    
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
    l = locoperator(2, 16, 16, exact_mu)
    OPTYPE = typeof(l)
    lop = Dict{Int64, OPTYPE}()

    # create an empty dictionary to store each element's neighbor's penalty params:
    NTYPE = typeof([l[13][1], l[13][2], l[13][3], l[13][4]])
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

     
        metrics = create_metrics(Nr[e], Ns[e], exact_mu, xt, zt) # create coordinate transform


        ###################################################################### 
        # create local finite difference operators on computational domain:
        lop[e] = locoperator(SBPp, Nr[e], Ns[e], exact_mu, metrics, FToB[EToF[:, e]]) 
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
            if LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
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
    creep(x, y, t) = (x ./ Lx) .* (Vp/2) .* t
    bc_Dirichlet = (lf, x, y, e, δ, t) -> creep(x,y,t)
    bc_Neumann   = (lf, x, y, nx, ny, e, δ, t) -> zeros(size(x))

    in_jump      = (lf, x, y, e, δ, t) -> begin
    f = EToF[lf, e]
    if EToS[lf, e] == 1
      if EToO[lf, e]
        return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else
        return -δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
      end
    else
      if EToO[lf, e]
        return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else
        return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
      end
    end
    end
    #}}}

    # initial time and rhs b of linear system (i.e. u = A\b)
    t = 0
    b = zeros(VNp)

    # initial slip vector
    δ = zeros(δNp)
  
    # fill in initial boundary data into b
    for e = 1:nelems
        loc_bdry_vec_v2!((@view b[vstarts[e]:vstarts[e+1]-1]), lop[e], neighborZ[e], FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, (e, δ, t))
    end
   
    u = A \ b # solve linear system with a backsolve to obtain initial displacement u(x, z, 0)

    # initialize change in shear stress due to quasi-static deformation:
    Δτ = zeros(δNp)

    # Assemble fault variables/data
    fault_nodes = zeros(δNp)
    RSa = zeros(δNp)
    for f = 1:nfaces
    if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (e1, _) = FToE[:, f]
        (lf1, _) = FToLF[:, f]
        xf = lop[e1].facecoord[1][lf1]
        yf = lop[e1].facecoord[2][lf1]
        δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
        for n = 1:length(δrng)
        RSa[δrng[n]] = RSamin - (RSamin - RSamax) *
            min(1, max(0, (RSH1 + yf[n])/(RSH1 - RSH2)))
        
        fault_nodes[δrng[n]] = yf[n]
        end
    end
    end

    # Set pre-stress according to benchmark description
    τ0 = fill(σn * RSamax * asinh(RSVinit / (2 * RSV0) *
                                exp.((RSf0 + RSb * log.(RSV0 / RSVinit)) /
                                        RSamax)) + η * RSVinit,
            δNp)
    for f = 1:nfaces
        if FToB[f] == RS_FAULT
            (e1, e2) = FToE[:, f]
            (lf1, lf2) = FToLF[:, f]
            nx = lop[e1].nx
            δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
            for n = 1:length(δrng)
            τ0[δrng[n]] = sign(nx[lf1][n])*abs(τ0[δrng[n]])
            end
        end
    end
     
    # Set initial state variable according to benchmark           
    θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
        sinh.((τ0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)
    # Initialize psi version of state variable    
    ψ = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)

    # Calculate total shear stress at t = 0:
    τ = τ0 + Δτ

    # Set initial condition for index 1 DAE - this is a stacked vector containing ψ, followed by slip
    ψδ = zeros(2δNp)
    ψδ[1:δNp] .= ψ
    ψδ[δNp+1:2δNp] .= δ

    # Set fault station locations (depths) specified in benchmark
    stations_locations = [0 0
                        0 -0.05
                        0 -2.5
                        0 -5
                        0 -7.5
                        0 -10
                        0 -12.5
                        0 -15
                        0 -17.5
                        0 -20
                        0 -25
                        0 -30
                        ]
    stations = setupfaultstations(stations_locations, lop, FToB, FToE, FToLF,
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
    b = b,
    u=u,
    τ=τ,
    Δτ = Δτ,
    vstarts=vstarts,
    μshear=μshear,
    RSa=RSa,
    RSb=RSb,
    σn=σn,
    η=η,
    RSV0=RSV0,
    τ0=τ0,
    RSDc=RSDc,
    RSf0=RSf0,
    Lx = Lx,
    fault_nodes
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
    cb = SavingCallback((ψδ, t, i)->savefaultstation(ψδ, t, i, stations,
                                                    FToδstarts, odeparam,
                                                    pth*"BP1_N_$(Nr[1])_$(Lx)", 10year_seconds),
                        SavedValues(Float64, Float64))

   
    # Solve DAE using Tsit5(), an adaptive Runge-Kutta method
    sol = solve(prob, Tsit5(); isoutofdomain=stepcheck, dt=0.1*year_seconds,
                abstol = 1e-9, reltol = 1e-9, save_everystep=false, #gamma = 0.05,
                internalnorm=(x, _)->norm(x, Inf), callback=cb)
                
                
    return sol
end
        

S = main();

# example of how to plot slip contours (uncomment if desired):
# plot_slip(pth*"slip.dat")

# examples of how ot plot times series of shear stress:
pth = "BP1_N_40_0.0_-7.5.dat"
plot_fault_time_series("shear_stress", pth)
# plot_fault_time_series("slip_rate", pth*"fltst_strk000.txt")
# plot_fault_time_series("shear_stress", pth*"fltst_strk+10.txt")
# plot_fault_time_series("state", pth*"fltst_strk+25.txt")