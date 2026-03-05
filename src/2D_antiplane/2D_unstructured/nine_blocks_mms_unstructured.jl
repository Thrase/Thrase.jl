# On an unstructured mesh of quads:
# Solves a manufactured solution

using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using DelimitedFiles
using MeshGrid


const year_seconds = 31556926
global const ctr = Ref{Int64}(1) 

include("ops_BP1-QD_unstructured.jl")
include("../utils_2D.jl")

const ϵ = 1
const γ = 1
const α = 1
const β = 1

# manufactured solution
function exact_mu(x, z)
    # Define shear modulus μ
    μ = 2 .+ 1 .* sin.(x .+ z)
    return μ
end

function exact_mu_x(x, z)
    # Define shear modulus derivative μ_x
    μ_x = 1 .* cos.(x .+ z)
    return μ_x
end

function exact_mu_z(x, z)
    # Define shear modulus derivative μ_z
    μ_z = 1 .* cos.(x .+ z)
    return μ_z
end

function exact(x, z, e, EToDomain)
    if EToDomain[e] == 1
        uexact = 1 .+ β .* x .+ γ .* z .+ ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    else
        uexact = 1 .+ β .* x .+ α .* z .+ ϵ .* sin.(pi * x) .* cos.(2 .*pi * z)
    end
    return uexact[:]
end

function exact_x(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return β .+ pi * ϵ .* cos.(pi * x) .* cos.(2 .* pi * z)
    else
        return β .+ pi * ϵ .* cos.(pi * x) .* cos.(2 .* pi * z)
    end
end
  
function exact_xx(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ -pi^2 * ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    else
        return 0 .+ -pi^2 * ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    end
end
  
function exact_z(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ γ  .+ -2 .* pi * ϵ .* sin.(pi * x) .* sin.(2 .* pi * z)
    else
        return 0 .+ α  .+ -2 .* pi * ϵ .* sin.(pi * x) .* sin.(2 .* pi * z)
    end
end
  
function exact_zz(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0  .+ -4 .* pi^2 * ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    else
        return 0 .+ -4 .* pi^2 * ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
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
    
    # I think these are the indices corresponding to vertices as you move through elements
    vstarts = Array{Int64, 1}(undef, nelems + 1)
    vstarts[1] = 1

    # Loop over blocks and create local operators
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
    
    #M = lu(M)

    # Get a unique array indexes for the face to jumps map
    FToδstarts = bcstarts(FToB, FToE, FToLF, BC_JUMP_INTERFACE, Nr, Ns)
   
    # Compute the number of volume and jump (δ) points
    VNp = vstarts[nelems+1]-1
    δNp = FToδstarts[nfaces+1]-1
    
    # define the jump data
    δ = zeros(δNp)
    for f = 1:nfaces
      if FToB[f] == BC_JUMP_INTERFACE
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (xf, yf) = lop[e1].facecoord
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] =
        exact(xf[lf1], yf[lf1], e2, EToDomain) - exact(xf[lf1], yf[lf1], e1, EToDomain)
        
      end
    end

   
    # define the interface function
    in_jump      = (lf, x, y, e, δ) -> begin
      f = EToF[lf, e]  # Get global face number
      if EToS[lf, e] == 1  # check if face is on minus side
        if EToO[lf, e]     # check is correct orientation
          return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else
          error("shouldn't get here")  # this is bcasue "correct" orientation is always true of a face on teh minus side. 
        end
      else                  # face on plus side
        if EToO[lf, e]      # check if orientation is correct
          return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else                # if orientation is reversed, reverse the data
          return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
        end
      end
    end

    # Set boundary data:
    bc_Dirichlet = (lf, x, y, e, δ) -> exact(x, y, e, EToDomain)
    bc_Neumann   = (lf, x, y, nx, ny, e, δ) -> (nx .* exact_mu(x, y) .* exact_x(x, y, e, EToDomain)
                                                + ny .* exact_mu(x,y) .* exact_z(x, y, e, EToDomain))

    # initialize RHS array that stores boundary data (linear system will be Au = b, where b = B*g)
    b = zeros(VNp)
    uexact = zeros(VNp)
    Jf = zeros(VNp)



    # # modify b to incorporate BC - check this. Do we have to reverse Dirichlet/Neumann boundary data if reversed?
    for e = 1:nelems
        loc_bdry_vec!((@view b[vstarts[e]:vstarts[e+1]-1]), lop[e], FToB[EToF[:,e]], bc_Dirichlet, bc_Neumann, in_jump, (e, δ))
        uexact[vstarts[e]:vstarts[e+1]-1] = exact(lop[e].coord[1], lop[e].coord[2], e, EToDomain)
        locsourcearray!((@view Jf[vstarts[e]:vstarts[e+1]-1]), source, lop[e].coord, lop[e].J, (e, EToDomain))
    end

    
    # solve linear system with a backsolve to obtain displacement.

    u = M \ (b-Jf) 
    

    return (u, uexact, lop, vstarts)
end


let

    ###### PHYSICAL DOMAIN #######


    # mesh file side set type to actually boundary condition type
    bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN, BC_JUMP_INTERFACE]
   
    (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/plate_test.inp")

    # modify a bit to change orientation
    EToV[:, 7] = [10; 5; 15; 7]
    EToF[:, 7] = [19; 18; 20; 13]
    #Can be misoriented, but faces and vertices need to follow! 

   

    # modify a bit to put fault in - this info should really be in .inp file tho
    EToDomain[1] = 2
    EToDomain[4] = 2
    EToDomain[7] = 2
    
    FToB = [1, 7, 1, 0, 0, 1, 0, 1, 1, 0, 1, 7, 0, 0, 0, 1, 0, 1, 7, 1, 0, 1, 1, 1]
   
 
    (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)
   

    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    @show (nelems, nfaces)
    
    # Plot the connectivity
    #plot_connectivity(verts, EToV)#- this needs some work.



    # Convergence tests:

    # Define spatial order of accuracy
    SBPp = 2

    # This is the base mesh size in each dimension
    N1 = N0 = 17

    # EToN0 is the base mesh size (e.g., before refinement)
    EToN0 = zeros(Int64, 2, nelems)
    EToN0[1, :] .= N0
    EToN0[2, :] .= N1


    ϵ = zeros(4)
    rates = zeros(length(ϵ)-1)


    for lvl = 1:length(ϵ)
        # Set up the local grid dimensions
        Nr = EToN0[1, :] * (2^(lvl-1))
        Ns = EToN0[2, :] * (2^(lvl-1))

        (u, uexact, lop, vstarts) = antiplane_solve(Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)
        
        H̃ = blockdiag(lop[1].H̃, lop[2].H̃, lop[3].H̃, lop[4].H̃, lop[5].H̃, lop[6].H̃, lop[7].H̃, lop[8].H̃, lop[9].H̃)
        diff = u .- uexact

        ϵ[lvl] = sqrt(diff' * H̃ * diff) ./ sqrt(uexact' * H̃ * uexact)   # relative error
        
    end

    for lvl = 2:length(ϵ)
        rates[lvl-1] = log2(ϵ[lvl-1]/ϵ[lvl])
    end

    @show ϵ
    @show rates

end

