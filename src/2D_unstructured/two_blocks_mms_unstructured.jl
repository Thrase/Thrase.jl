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

function exact(x, z, e)
    if e == 1
        uexact = 1 .+ β .* x .+ γ .* z .* (z .- 1) .+ ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    else
        uexact = 1 .+ β .* x .+ α .* z .* (z .- 1) .+ ϵ .* sin.(pi * x) .* cos.(2 .*pi * z)
    end
    return uexact[:]
end

function exact_x(x, z, e)
    if e == 1
        return β .+ pi * ϵ .* cos.(pi * x) .* cos.(2 .* pi * z)
    else
        return β .+ pi * ϵ .* cos.(pi * x) .* cos.(2 .* pi * z)
    end
end
  
function exact_xx(x, z, e)
    if e == 1
        return 0 .+ -pi^2 * ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    else
        return 0 .+ -pi^2 * ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    end
end
  
function exact_z(x, z, e)
    if e == 1
        return 0 .+ γ .* (2 .* z .- 1) .+ -2 .* pi * ϵ .* sin.(pi * x) .* sin.(2 .* pi * z)
    else
        return 0 .+ α .* (2 .* z .- 1).+ -2 .* pi * ϵ .* sin.(pi * x) .* sin.(2 .* pi * z)
    end
end
  
function exact_zz(x, z, e)
    if e == 1
        return γ .* 2  .+ -4 .* pi^2 * ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    else
        return α .* 2 .+ -4 .* pi^2 * ϵ .* sin.(pi * x) .* cos.(2 .* pi * z)
    end       
end



function source(x, z, e)
    f = -exact_mu_x(x, z) .* exact_x(x, z, e) .- exact_mu_z(x, z) .* exact_z(x, z, e) .- exact_mu(x, z) .* (exact_xx(x, z, e) .+ exact_zz(x, z, e))  # source data
    return f[:]
end




function antiplane_solve(Nr, Ns, SBPp)

    ###### PHYSICAL DOMAIN #######

    
    #verts = [-1.1 1 -1.2 0.9;
    #        -.29 -1 .891 1.12];

    #verts = [-1 1 -1 1;
    #        -1 -1 1 1];

    # mesh file side set type to actually boundary condition type
    bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN, BC_JUMP_INTERFACE]
    #(verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/square_circle.inp";
    #                                        bc_map = bc_map)
    # EToV defines the element by its vertices
    # EToF defines element by its four faces, in global face number
    # FToB defines whether face is Dirichlet (1), Neumann (2), interior jump (7)
    #      or just an interior interface (0)

    # TWO ELEMENTS
    # verts = [-1 0 1 -1 0 1;
    #           0 0 0 1  1 1];

    # FToB = [1; 7; 1; 1; 1; 1; 1];

    # same orientation
    # EToV = [1 2; 2 3; 4 5; 5 6];
    # EToF = [1 2; 
    #         2 5;
    #         4 3;
    #         6 7];
    # opposite orientation
    # EToV = [1 6; 2 5; 4 3; 5 2];
    # EToF = [1 3; 
    #         2 2;
    #         4 7;
    #         6 5];


  # TWO ELEMENTS
    # verts = [-1 0 1 -1 0 1;
    #           0 0 0 1  1 1];

    # FToB = [1; 7; 1; 1; 1; 1; 1];
    
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
   
    @show EToO

    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    @show (nelems, nfaces)
    
    # Plot the connectivity
    #plot_connectivity(verts, EToV)#- this needs some work.



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

        Np = (Nr+1)*(Ns+1)  # TODO: modify to allow different numbers on different elements
        vstarts[e+1] = vstarts[e] + Np

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

        metrics = create_metrics(Nr, Ns, exact_mu, xt, zt) # create coordinate transform


        ###################################################################### 
        # create local finite difference operators on computational domain:
        lop[e] = locoperator(SBPp, Nr, Ns, exact_mu, metrics, FToB[EToF[:, e]]) 
    end

  
    # Assemble the global volume operator and factor:
    M = global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
    #H̃ = blockdiag(lop[1].H̃, lop[2].H̃)
    #mMH = -H̃*M 
    #@show extrema(mMH-mMH')
    #@show UnicodePlots.spy(mMH-mMH')
    
    M = lu(M)

    # Get a unique array indexes for the face to jumps map
    FToδstarts = bcstarts(FToB, FToE, FToLF, BC_JUMP_INTERFACE, Nr, Ns)
    # Compute the number of volume and jump (δ) points
    VNp = vstarts[nelems+1]-1
    δNp = FToδstarts[nfaces+1]-1
    
    # define the jump data
    δ = zeros(δNp)
    for f = 1:nfaces
      if FToB[f] == BC_JUMP_INTERFACE
        @show (e1, e2) = FToE[:, f]
        @show (lf1, lf2) = FToLF[:, f]
        
        (xf, yf) = lop[e1].facecoord
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] =
        exact(xf[lf1], yf[lf1], e2) - exact(xf[lf1], yf[lf1], e1)
        
        
        #@show δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      end
    end

    #@show δ
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
    bc_Dirichlet = (lf, x, y, e, δ) -> exact(x, y, e)
    bc_Neumann   = (lf, x, y, nx, ny, e, δ) -> (nx .* exact_mu(x, y) .* exact_x(x, y, e)
                                                + ny .* exact_mu(x,y) .* exact_z(x, y, e))

    # initialize RHS array that stores boundary data (linear system will be Au = b, where b = B*g)
    b = zeros(2*(Nr+1) *(Ns+1))
    #vstarts = [1, (Nr+1)*(Ns+1)+1, 2*(Nr+1) * (Ns+1)+1]
    uexact = zeros(2*(Nr+1) *(Ns+1))
    Jf = zeros(2*(Nr+1) *(Ns+1))
    # modify b to incorporate BC
    for e = 1:nelems
        loc_bdry_vec!((@view b[vstarts[e]:vstarts[e+1]-1]), lop[e], FToB[EToF[:,e]], bc_Dirichlet, bc_Neumann, in_jump, (e, δ))
        uexact[vstarts[e]:vstarts[e+1]-1] = exact(lop[e].coord[1], lop[e].coord[2], e)
        locsourcearray!((@view Jf[vstarts[e]:vstarts[e+1]-1]), source, lop[e].coord, lop[e].J, (e))
    end

    
    # solve linear system with a backsolve to obtain displacement.

    u = M \ (b-Jf) 
    

    return (u, uexact, lop[1].H̃, lop[2].H̃, lop[1], lop[2], vstarts)
end

# Convergence tests:

# Define spatial order of accuracy
SBPp = 2

# Number of grid points on each element
Nr = [2^3, 2^4, 2^5, 2^6, 2^7, 2^8]
Ns = Nr

errors = zeros(length(Nr))
rates = zeros(length(Nr)-1)

for n = 1:length(Nr)
#n = 3
    (u, uexact, H̃1, H̃2, lop1, lop2, vstarts) = antiplane_solve(Nr[n], Ns[n], SBPp)
    #pp = lop1.eRST[2]*lop1.T[2]*u[vstarts[1]:vstarts[2]-1]
    #mm = lop2.eRST[2]*lop2.T[2]*u[vstarts[2]:vstarts[3]-1]
    # @show pp
    # @show mm
    # @show lop1.eRST[2]*u[vstarts[1]:vstarts[2]-1]
    # @show lop2.eRST[2]*u[vstarts[2]:vstarts[3]-1]
    # @show pp
    # poo
    #   u1 = reshape(u[vstarts[1]:vstarts[2]-1], Nr[n]+1, Ns[n]+1)
    #   u2 = reshape(u[vstarts[2]:vstarts[3]-1], Nr[n]+1, Ns[n]+1)

    #   u1exact = reshape(uexact[vstarts[1]:vstarts[2]-1], Nr[n]+1, Ns[n]+1)
    #   u2exact = reshape(uexact[vstarts[2]:vstarts[3]-1], Nr[n]+1, Ns[n]+1)
   
     
     #surface(lop1.coord[1], lop1.coord[2], u1, title="My Surface Plot (with MeshGrid.jl)", camera = (15, 35))
     #surface!(lop2.coord[1], lop2.coord[2], u2, title="My Surface Plot (with MeshGrid.jl)", camera = (15, 35))

    #  surface(lop1.coord[1], lop1.coord[2], u1exact-u1, title="My Surface Plot (with MeshGrid.jl)", camera = (15, 35))
    #  surface!(lop2.coord[1], lop2.coord[2], u2exact-u2, title="My Surface Plot (with MeshGrid.jl)", camera = (15, 35))

    #poo
    H̃ = blockdiag(H̃1, H̃2)
    diff = u .- uexact

    errors[n] = sqrt(diff' * H̃ * diff) ./ sqrt(uexact' * H̃ * uexact)   # relative error
    
end

 for n = 2:length(Nr)
     rates[n-1] = log2(errors[n-1]/errors[n])
 end

 @show errors
 @show rates

# p = (2, 4, 6) should yield rates of (2, 3, 4), respectively.

# Plot error on coarsest grid
#using Plots
#pyplot();
#(u, uexact, H̃, coord) = antiplane_solve(Nr[1], Ns[1], SBPp)
#u = reshape(u, Nr[1]+1, Ns[1]+1)
#uexact = reshape(uexact, Nr[1]+1, Ns[1]+1)
#p1 = surface(coord[1], coord[2], uexact - u, xlabel = "x", ylabel = "z", showaxis = true,camera=(30,30))
#display(p1)


