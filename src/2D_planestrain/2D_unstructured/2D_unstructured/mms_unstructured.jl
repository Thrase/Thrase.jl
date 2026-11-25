# On an unstructured mesh of quads:
# Solves a manufactured solution
# Geometry set up to resemble BP1 - so a vertical friction fault down to z = 40km, then steady sliding

using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using DelimitedFiles
using PyPlot

const year_seconds = 31556926
global const ctr = Ref{Int64}(1) 

const RS_FAULT = 7
const VP_FAULT = 8

include("ops_BP1-QD_unstructured.jl")
include("../utils_2D.jl")

const ϵ = 1
const γ = 0.01
const α = 0.02
const β = 0.03
const ν = 0
const c = 1
const L = 160

# manufactured solution
function exact_mu(x, z)
    # Define shear modulus μ
    μ = 2 .+ c .* sin.(x ./ L .+ z ./ L)
    return μ
end

function exact_mu_x(x, z)
    # Define shear modulus derivative μ_x
    μ_x = c ./ L .* cos.(x ./ L .+ z ./ L)
    return μ_x
end

function exact_mu_z(x, z)
    # Define shear modulus derivative μ_z
    μ_z = c ./ L .* cos.(x ./ L .+ z ./ L)
    return μ_z
end

function exact(x, z, e, EToDomain)
    if EToDomain[e] == 1
        uexact = 1 .+ β .* x .+ γ .* z .+ ϵ .* sin.(pi * x ./ L) .* cos.(2 .* pi * z ./ L)
    else
        uexact = 2 .+ ν .* x .+ α .* z .+ ϵ .* cos.(pi * x ./ L) .* cos.(2 .*pi * z ./ L) # never getting here.
    end
    return uexact[:]
end

function exact_x(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return β .+ (pi ./ L) * ϵ .* cos.(pi * x ./ L) .* cos.(2 .* pi * z ./ L)
    else
        return ν .+ (-pi ./ L) * ϵ .* sin.(pi * x ./ L) .* cos.(2 .* pi * z ./ L)
    end
end
  
function exact_xx(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ -pi^2 ./ L^2 * ϵ .* sin.(pi * x ./ L) .* cos.(2 .* pi * z ./ L)
    else
        return 0 .- pi^2 ./ L^2 * ϵ .* cos.(pi * x ./ L) .* cos.(2 .* pi * z ./ L)
    end
end
  
function exact_z(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ γ  .+ -2 .* pi ./ L * ϵ .* sin.(pi * x ./ L) .* sin.(2 .* pi * z ./ L)
    else
        return 0 .+ α  .+ -2 .* pi ./ L * ϵ .* cos.(pi * x ./ L) .* sin.(2 .* pi * z ./ L)
    end
end
  
function exact_zz(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0  .+ -4 .* pi^2 ./ L^2 * ϵ .* sin.(pi * x ./ L) .* cos.(2 .* pi * z ./ L)
    else
        return 0 .+ -4 .* pi^2 ./ L^2 * ϵ .* cos.(pi * x ./ L) .* cos.(2 .* pi * z ./ L)
    end       
end



function source(x, z, e, EToDomain)
    f = -exact_mu_x(x, z) .* exact_x(x, z, e, EToDomain) .- exact_mu_z(x, z) .* exact_z(x, z, e, EToDomain) .- exact_mu(x, z) .* (exact_xx(x, z, e, EToDomain) .+ exact_zz(x, z, e, EToDomain))  # source data
    return f[:]
end




function antiplane_solve(Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)

    nelems = length(Nr)
    
    ################################## COORDINATE TRANSFORM ###################################
    #
    # Build the local volume operators
    #

    # Create an empty dictionary to store the operators;
    # index via element number, return structure containing local ops.
    OPTYPE = typeof(locoperator(2, 16, 16, exact_mu))
    lop = Dict{Int64, OPTYPE}()
    
    # Indices corresponding to vertices as you move through blocks/elements
    vstarts = Array{Int64, 1}(undef, nelems + 1)
    vstarts[1] = 1   # start at 1. 

    # Loop over blocks/elements and create local operators:
    for e = 1:nelems

        Np = (Nr[e]+1)*(Ns[e]+1)  # total number of volume points on each element
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

        # Uncomment if coordinate transform user-specified, for example:
        #xt(r,s) = (2 * r, 2 * ones(size(r)), 0 * ones(size(s)))
        #zt(r,s) = (2 * s, 0 * ones(size(r)), 2 * ones(size(s)))

        metrics = create_metrics(Nr[e], Ns[e], exact_mu, xt, zt) # create coordinate transform


        ###################################################################### 
        # create local finite difference operators on computational domain:
        lop[e] = locoperator(SBPp, Nr[e], Ns[e], exact_mu, metrics, FToB[EToF[:, e]]) 
    end

  
    # Assemble the global volume operator and factor:
    M = global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
    M = lu(M)

    # Get a unique array indices for the faces corresponding to the fault/jump interface
    FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT), Nr, Ns)
   
    # Compute the total number of volume and fault/jump (δ) points
    VNp = vstarts[nelems+1]-1
    δNp = FToδstarts[nfaces+1]-1
    
    # define the data for the jump in displacement on fault interface
    δ = zeros(δNp)
    for f = 1:nfaces
      if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (xf, yf) = lop[e1].facecoord
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] =
        exact(xf[lf1], yf[lf1], e2, EToDomain) - exact(xf[lf1], yf[lf1], e1, EToDomain)
        
      end
    end
  
    # define the data for the jump in traction data on fault interface
    θ = zeros(δNp)
    for f = 1:nfaces
        if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (xf, yf) = lop[e1].facecoord
        
        tf1 = lop[e1].nx[lf1] .* exact_mu(xf[lf1], yf[lf1]) .* exact_x(xf[lf1], yf[lf1], e1, EToDomain) + lop[e1].ny[lf1] .* exact_mu(xf[lf1], yf[lf1]) .* exact_z(xf[lf1], yf[lf1], e1, EToDomain)
        tf2 = lop[e2].nx[lf2] .* exact_mu(xf[lf1], yf[lf1]) .* exact_x(xf[lf1], yf[lf1], e2, EToDomain) + lop[e2].ny[lf2] .* exact_mu(xf[lf1], yf[lf1]) .* exact_z(xf[lf1], yf[lf1], e2, EToDomain)
        
        @views θ[FToδstarts[f]:(FToδstarts[f+1]-1)] = tf1 + tf2
        end
    end

    # define the interface jump in displacement function
    in_jump      = (lf, x, y, e, δ, θ) -> begin
      f = EToF[lf, e]  # Get global face number
      if EToS[lf, e] == 1  # check if face is on minus side
        if EToO[lf, e]     # check is correct orientation
          return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else
          error("shouldn't get here")  # this is because "correct" orientation is always true of a face on the minus side. 
        end
      else                  # face on plus side
        if EToO[lf, e]      # check if orientation is correct
          return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else                # if orientation is reversed, reverse the data
          return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
        end
      end
    end

    # define the interface jump in traction fucntion
    in_tractionjump      = (lf, x, y, e, δ, θ) -> begin
      f = EToF[lf, e]  # Get global face number
      if EToS[lf, e] == 1  # check if face is on minus side
        if EToO[lf, e]     # check is correct orientation
          return θ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else
          error("shouldn't get here")  # this is bcasue "correct" orientation is always true of a face on teh minus side. 
        end
      else                  # face on plus side
        if EToO[lf, e]      # check if orientation is correct
          return  θ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else                # if orientation is reversed, reverse the data
          return  θ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
        end
      end
    end

    # Set boundary data:
    bc_Dirichlet = (lf, x, y, e, δ, θ) -> exact(x, y, e, EToDomain)
    bc_Neumann   = (lf, x, y, nx, ny, e, δ, θ) -> (nx .* exact_mu(x, y) .* exact_x(x, y, e, EToDomain)
                                                + ny .* exact_mu(x,y) .* exact_z(x, y, e, EToDomain))

    # initialize RHS array that stores boundary data (linear system will be Au = b, where b = B*g)
    b = zeros(VNp)
    uexact = zeros(VNp)
    Jf = zeros(VNp)

    
    # # modify b to incorporate BC.
    for e = 1:nelems
        LFToB = FToB[EToF[:,e]]

        neighborZ = [similar(lop[e].IsJZ[1]), similar(lop[e].IsJZ[2]), similar(lop[e].IsJZ[3]), similar(lop[e].IsJZ[4])]
            
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

                neighborZ[nf] .= lop[eo].IsJZ[nf] # need to store the Z from the other element's local face
                
            else  
            end
        
        end
    
   
    loc_bdry_vec!((@view b[vstarts[e]:vstarts[e+1]-1]), lop[e], neighborZ, FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, in_tractionjump, (e, δ, θ))
    uexact[vstarts[e]:vstarts[e+1]-1] = exact(lop[e].coord[1], lop[e].coord[2], e, EToDomain)
    locsourcearray!((@view Jf[vstarts[e]:vstarts[e+1]-1]), source, lop[e].coord, lop[e].J, (e, EToDomain))

    end
    
    # solve linear system with a backsolve to obtain displacement.

    u = M \ (b-Jf) 
    
    

    return (u, uexact, lop, vstarts)

end


let 
    # SBP interior order of accuracy
    SBPp = 2 

    # This actually has to go with .inp file. Mesh file side set type to actually boundary condition type. Intially, every face is locked interface unless specified otherwise in sideset info.
    #REbc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN, BC_NEUMANN, BC_NEUMANN, BC_JUMP_INTERFACE, BC_JUMP_INTERFACE]
    # read in mesh from an .inp file and put in boundary condition type
    
   (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/BP1_v1.inp")


   #uncomment for two blocks
#  verts = [-100 0 100 -100 0 100;
#              0 0 0 100  100 100];

#    FToB = [1; 7; 1; 1; 1; 1; 1];

#   #  same orientation
#     EToV = [1 2; 2 3; 4 5; 5 6];
#     EToF = [1 2; 
#             2 5;
#             4 3;
#             6 7];

    
#    EToDomain = [1 2]
    #EToV defines the element by its vertices
    #EToF defines element by its four faces, in global face number
    #FToB defines whether face is Dirichlet (1), Neumann (2), interior jump (7 for RS fault, 8 for steady sliding)
     #    or just an interior interface (0)
    #EToDomain specific to domain, e.g. 1 if inside circle, 2 otherwise

    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    #@show (nelems, nfaces)

    # Determine secondary arrays
    # FToE : Unique Global Face to Element Number
    #        (the i'th column of this stores the element numbers that share the
    #        global face number i)
    # FToLF: Unique Global Face to Element local face number
    #        (the i'th column of this stores the element local face numbers that
    #        shares the global face number i)
    # EToO : Element to Unique Global Faces Orientation
    #        (the i'th column of this stores the whether the element and global
    #        face are oriented in the same way in physical memory or need to be
    #        rotated)
    # EToS : Element to Unique Global Face Side
    #        (the i'th column of this stores whether an element face is on the
    #        plus side or minus side of the global face)
    (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

    # This is the base mesh size in each dimension
    N1 = N0 = 17

    # EToN0 is the base mesh size (e.g., before refinement)
    EToN0 = zeros(Int64, 2, nelems)
    EToN0[1, :] .= N0
    EToN0[2, :] .= N1


    # create arrays for storing errors and rates of convergence
    ϵ = zeros(2)
    rates = zeros(length(ϵ)-1)

    # successively refine 
    for lvl = 1:length(ϵ)

       
        # Set up the local grid dimensions (i.e. number of nodes in each direction on each element)
        Nr = EToN0[1, :] * (2^(lvl-1))
        Ns = EToN0[2, :] * (2^(lvl-1))

        # solve for displacements u and compute exact solution uexact. lop are local ops and vstars are indices for each element.
        (u, uexact, lop, vstarts) = antiplane_solve(Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)
                                                   
        if lvl == 1 # plot approximation, exact solution, or error on coarsest grid
            
            #better_plot_solution(uexact, nelems, vstarts, Nr, Ns, lop) 
            
            better_plot_solution(u, nelems, vstarts, Nr, Ns, lop)
            #better_plot_solution(uexact-u, nelems, vstarts, Nr, Ns, lop)
            #poo
          
        end

        H̃ = blockdiag(lop[1].H̃)
        for e = 2:length(lop)
            H̃ = blockdiag(H̃, lop[e].H̃)
        end
        diff = u .- uexact

        ϵ[lvl] = sqrt(diff' * H̃ * diff) ./ sqrt(uexact' * H̃ * uexact)   # relative error
        
    end

    for lvl = 2:length(ϵ)
        rates[lvl-1] = log2(ϵ[lvl-1]/ϵ[lvl])
    end

    @show ϵ
    @show rates
end