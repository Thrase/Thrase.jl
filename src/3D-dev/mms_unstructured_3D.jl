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

include("ops_BP9-QD_unstructured.jl")
include("utils_3D.jl")

const ϵ = 1
const γ = 0.01
const α = 0.02
const β = 0.03
const ν = 0
const c = 1
const L = 160

# manufactured solution
# BASIN STRUCTURE
B_p = (λ = 1.0, μ_in = 2, μ_out = 5, c = 1, r̄ = 144, r_w = 20, on = true)

function μ(x, y, z, B_p)

    c = B_p.c
    μ_in = B_p.μ_in
    μ_out = B_p.μ_out
    r̄ = B_p.r̄
    r_w = B_p.r_w
    on = B_p.on

    if on == false
        return repeat([μ_out], outer=size(x))
    else
        return (μ_out - μ_in)/2 *
            (tanh.((x .^ 2 .+ c^2 * y .^ 2 .- r̄) ./ r_w) .+ 1) .+ μ_in
    end
end

function μ_x(x, y, z, B_p)
    
    c = B_p.c
    μ_in = B_p.μ_in
    μ_out = B_p.μ_out
    r̄ = B_p.r̄
    r_w = B_p.r_w
    on = B_p.on

    if on == false
        return zeros(size(x))  
    else
        return ((μ_out - μ_in) .* x .*
                sech.((x .^ 2 .+ c^2 * y .^ 2 .- r̄) ./ r_w) .^ 2) ./ r_w
    end
end

function μ_y(x, y, z, B_p)

    c = B_p.c
    μ_in = B_p.μ_in
    μ_out = B_p.μ_out
    r̄ = B_p.r̄
    r_w = B_p.r_w
    on = B_p.on

    if on == false
        
        return zeros(size(x))

    else    
        return ((μ_out - μ_in) .* (c^2 * y) .*
            sech.((x .^ 2 + c^2 * y .^ 2 .- r̄) ./ r_w) .^ 2) ./ r_w
    end
end

function μ_z(x, y, z, B_p)

    c = B_p.c
    μ_in = B_p.μ_in
    μ_out = B_p.μ_out
    r̄ = B_p.r̄
    r_w = B_p.r_w
    on = B_p.on

    if on == false
        
        return zeros(size(x))

    else    
        return 0 .* ((μ_out - μ_in) .* (c^2 * y) .*
            sech.((x .^ 2 + c^2 * y .^ 2 .- r̄) ./ r_w) .^ 2) ./ r_w   #zeros!
    end
end

function λ(x, y, z, B_p)
    return B_p.λ .+ zeros(size(x))
end

function λ_x(x, y, z, B_p)
    return zeros(size(x))
end

function λ_y(x, y, z, B_p)
    return zeros(size(x))
end

function λ_z(x, y, z, B_p)
    return zeros(size(x))
end

function K(x, y, z, B_p)
    return λ(x, y, z, B_p) .+ 2 .* μ(x, y, z, B_p) ./ 3
end

function K_x(x, y, z, B_p)
    return λ_x(x, y, z, B_p) .+ 2 .* μ_x(x, y, z, B_p) ./ 3
end

function K_y(x, y, z, B_p)
    return λ_y(x, y, z, B_p) .+ 2 .* μ_y(x, y, z, B_p) ./ 3
end

function K_z(x, y, z, B_p)
    return λ_z(x, y, z, B_p) .+ 2 .* μ_z(x, y, z, B_p) ./ 3
end

function exact_u1(x, y, z)
    a = 1
    b = 1
    u1 =    a .+ b .* sin.(x .+ y .+  z)
    u1_x =   b * cos.(x .+ y .+  z)
    u1_y =   b * cos.(x .+ y .+  z)
    u1_z =   b * cos.(x .+ y .+  z)
    u1_xx = -b * sin.(x .+ y .+  z)
    u1_xy = -b * sin.(x .+ y .+  z)
    u1_xz = -b * sin.(x .+ y .+  z)
    u1_yy = -b * sin.(x .+ y .+  z)
    u1_yz = -b * sin.(x .+ y .+  z)
    u1_zz = -b * sin.(x .+ y .+  z)

    return (u1, u1_x, u1_y, u1_z, u1_xx, u1_xy, u1_xz, u1_yy, u1_yz, u1_zz)
end

function exact_u2(x, y, z)
    
    a = 1
    b = 1
    u2 =    a .+ b .* cos.(x .+ y .+  z)
    u2_x =   -b * sin.(x .+ y .+  z)
    u2_y =   -b * sin.(x .+ y .+  z)
    u2_z =   -b * sin.(x .+ y .+  z)
    u2_xx = -b * cos.(x .+ y .+  z)
    u2_xy = -b * cos.(x .+ y .+  z)
    u2_xz = -b * cos.(x .+ y .+  z)
    u2_yy = -b * cos.(x .+ y .+  z)
    u2_yz = -b * cos.(x .+ y .+  z)
    u2_zz = -b * cos.(x .+ y .+  z)


    
    return (u2, u2_x, u2_y, u2_z, u2_xx, u2_xy, u2_xz, u2_yy, u2_yz, u2_zz)
end

function exact_u3(x, y, z)
    
    a = 1
    b = 1
    u3 =    a .+ b .* sin.(0.2x .+ 0.3y .+  0.4z)
    u3_x =   0.2b  * cos.(0.2x .+ 0.3y .+  0.4z)
    u3_y =   0.3b  * cos.(0.2x .+ 0.3y .+  0.4z)
    u3_z =   0.4b  * cos.(0.2x .+ 0.3y .+  0.4z)
    u3_xx = -.04b  * sin.(0.2x .+ 0.3y .+  0.4z)
    u3_xy = -.06b  * sin.(0.2x .+ 0.3y .+  0.4z)
    u3_xz = -.08b * sin.(0.2x .+ 0.3y .+  0.4z)
    u3_yy = -.09b  * sin.(0.2x .+ 0.3y .+  0.4z)
    u3_yz = -.12b * sin.(0.2x .+ 0.3y .+  0.4z)
    u3_zz = -.16b * sin.(0.2x .+ 0.3y .+  0.4z)

   
    return (u3, u3_x, u3_y, u3_z, u3_xx, u3_xy, u3_xz, u3_yy, u3_yz, u3_zz)

end

function source(x, y, z)
    (u1, u1_x, u1_y, u1_z, u1_xx, u1_xy, u1_xz, u1_yy, u1_yz, u1_zz) = exact_u1(x, y, z)
    (u2, u2_x, u2_y, u2_z, u2_xx, u2_xy, u2_xz, u2_yy, u2_yz, u2_zz) = exact_u2(x, y, z)
    (u3, u3_x, u3_y, u3_z, u3_xx, u3_xy, u3_xz, u3_yy, u3_yz, u3_zz) = exact_u3(x, y, z)

   
    σxx_x = λ_x(x, y, z, B_p) .* (u1_x + u2_y .+ u3_z) + 2μ_x(x, y, z, B_p) .* u1_x +  λ(x, y, z, B_p) .* (u1_xx .+ u2_xy .+ u3_xz) .+ 2μ(x, y, z, B_p) .* u1_xx
    σxy_y = μ_y(x, y, z, B_p) .* (u1_y + u2_x) .+ μ(x, y, z, B_p) .* (u1_yy .+ u2_xy)
    σxz_z = μ_z(x, y, z, B_p) .* (u1_z + u3_x) .+ μ(x, y, z, B_p) .* (u1_zz .+ u3_xz)

    σxy_x = μ_x(x, y, z, B_p) .* (u1_y + u2_x) .+ μ(x, y, z, B_p) .* (u1_xy .+ u2_xx)
    σyy_y = λ_y(x, y, z, B_p) .* (u1_x + u2_y .+ u3_z) + 2μ_y(x, y, z, B_p) .* u2_y + λ(x, y, z, B_p) .* (u1_xy .+ u2_yy .+ u3_yz) + 2μ(x, y, z, B_p) .* u2_yy
    σyz_z = μ_z(x, y, z, B_p) .* (u2_z + u3_y) .+ μ(x, y, z, B_p) .* (u2_zz .+ u3_yz)

    σxz_x = μ_x(x, y, z, B_p) .* (u1_z + u3_x) .+ μ(x, y, z, B_p) .* (u1_xz .+ u3_xx)
    σyz_y = μ_y(x, y, z, B_p) .* (u2_z + u3_y) .+ μ(x, y, z, B_p) .* (u2_yz .+ u3_yy)
    σzz_z = λ_z(x, y, z, B_p) .* (u1_x + u2_y + u3_z) .+ 2μ_z(x, y, z, B_p) .* u3_z + λ(x, y, z, B_p) .* (u1_xz .+ u2_yz + u3_zz) + 2μ(x, y, z, B_p) .* u3_zz

    
    s1 = σxx_x .+ σxy_y .+ σxz_z
    s2 = σxy_x .+ σyy_y .+ σyz_z
    s3 = σxz_x .+ σyz_y .+ σzz_z
    
    
    return [s1[:];s2[:];s3[:]]
end






function elasticity_solve(Nq, Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)

    nelems = length(Nr)
    
    ################################## COORDINATE TRANSFORM ###################################
    #
    # Build the local volume operators
    #

    # Create an empty dictionary to store the operators;
    # index via element number, return structure containing local ops.
    OPTYPE = typeof(locoperator(2, 16, 16, 16, exact_mu))
    lop = Dict{Int64, OPTYPE}()
    
    # Indices corresponding to vertices as you move through blocks/elements
    vstarts = Array{Int64, 1}(undef, nelems + 1)
    vstarts[1] = 1   # start at 1. 

    # Loop over blocks/elements and create local operators:
    for e = 1:nelems

        Np = (Nq[e]+1)*(Nr[e]+1)*(Ns[e]+1)  # total number of volume points on each element
        vstarts[e+1] = vstarts[e] + Np # fill in vstarts

        # DOMAIN: ordered pairs that define the physical domain
        xcorners = verts[1, EToV[:, e]] 
        ycorners = verts[2, EToV[:, e]]
        zcornder = verts[3, EToV[:, e]]

        xt, yt, zt = transfinite_blend_3D(xcorners, ycorners, zcorners)


        metrics = create_metrics_3D(Nq[e], Nr[e], Ns[e], exact_mu, xt, yt, zt) # create coordinate transform


        ###################################################################### 
        # create local finite difference operators on computational domain:
        lop[e] = locoperator_3D(SBPp, Nq[e], Nr[e], Ns[e], exact_mu, metrics, FToB[EToF[:, e]]) 
    end

  
    # Assemble the global volume operator and factor:
    M = global_operator_3D(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
    M = lu(M)

    # Get a unique array indices for the faces corresponding to the fault/jump interface
    FToδstarts = bcstarts_3D(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT), Nq, Nr, Ns)
   
    # Compute the total number of volume and fault/jump (δ) points
    VNp = vstarts[nelems+1]-1
    δNp = FToδstarts[nfaces+1]-1
    
    # define the data for the jump in displacement on fault interface
    δ = zeros(δNp)
    for f = 1:nfaces
      if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (xf, yf, zf) = lop[e1].facecoord
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] =
        exact(xf[lf1], yf[lf1], zf[lf1], e2, EToDomain) - exact(xf[lf1], yf[lf1], zf[lf1], e1, EToDomain)
        
      end
    end
  
    # define the data for the jump in traction data on fault interface
    θ = zeros(δNp)
    for f = 1:nfaces
        if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (xf, yf, zf) = lop[e1].facecoord
        
        tf1 = lop[e1].nx[lf1] .* exact_mu(xf[lf1], yf[lf1], zf[lf1]) .* exact_x(xf[lf1], yf[lf1], zf[lf1], e1, EToDomain) + lop[e1].ny[lf1] .* exact_mu(xf[lf1], yf[lf1], zf[lf1]) .* exact_y(xf[lf1], yf[lf1], zf[lf1], e1, EToDomain) + lop[e1].nz[lf1] .* exact_mu(xf[lf1], yf[lf1], zf[lf1]) .* exact_z(xf[lf1], yf[lf1], zf[lf1], e1, EToDomain)
        tf2 = lop[e2].nx[lf2] .* exact_mu(xf[lf1], yf[lf1], zf[lf1]) .* exact_x(xf[lf1], yf[lf1], zf[lf1], e2, EToDomain) + lop[e2].ny[lf2] .* exact_mu(xf[lf1], yf[lf1], zf[lf1]) .* exact_y(xf[lf1], yf[lf1], zf[lf1], e2, EToDomain) + lop[e2].ny[lf2] .* exact_mu(xf[lf1], yf[lf1], zf[lf1]) .* exact_z(xf[lf1], yf[lf1], zf[lf1], e2, EToDomain)
        
        @views θ[FToδstarts[f]:(FToδstarts[f+1]-1)] = tf1 + tf2
        end
    end

    # define the interface jump in displacement function
    in_jump      = (lf, x, y, z, e, δ, θ) -> begin
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
    in_tractionjump      = (lf, x, y, z, e, δ, θ) -> begin
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
    bc_Dirichlet = (lf, x, y, z, e, δ, θ) -> exact(x, y, z, e, EToDomain)
    bc_Neumann   = (lf, x, y, z, nx, ny, e, δ, θ) -> (nx .* exact_mu(x, y, z) .* exact_x(x, y, z, e, EToDomain)
                                                + ny .* exact_mu(x,y,z) .* exact_y(x, y, z, e, EToDomain)
                                                + nz .* exact_mu(x,y,z) .* exact_z(x, y, z, e, EToDomain))

    # initialize RHS array that stores boundary data (linear system will be Au = b, where b = B*g)
    b = zeros(VNp)
    uexact = zeros(VNp)
    Jf = zeros(VNp)

    
    # # modify b to incorporate BC.
    for e = 1:nelems
        LFToB = FToB[EToF[:,e]]

        neighborZ = [similar(lop[e].IsJZ[1]), similar(lop[e].IsJZ[2]), similar(lop[e].IsJZ[3]), similar(lop[e].IsJZ[4]), similar(lop[e].IsJZ[5]), similar(lop[e].IsJZ[6])]
            
        # loop over the faces to get the neighboring penalty parameters:
        for lf = 1:6   
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
    
   
    loc_bdry_vec_3D!((@view b[vstarts[e]:vstarts[e+1]-1]), lop[e], neighborZ, FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, in_tractionjump, (e, δ, θ))
    uexact[vstarts[e]:vstarts[e+1]-1] = exact(lop[e].coord[1], lop[e].coord[2], lop[e].coord[3], e, EToDomain)
    locsourcearray_3D!((@view Jf[vstarts[e]:vstarts[e+1]-1]), source, lop[e].coord, lop[e].J, (e, EToDomain))

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
    
   #(verts, EToV, EToF, FToB, EToDomain) = read_inp_3d("meshes/BP9.inp")


    # two blocks:

    EToV = [8 4;
        7 6;
        10 3;
        6 2;
        12 5;
        5 7;
        1 9;
        4 11]


    EToF = [1 7;
            2 8;
            3 9;
            4 10;
            5 2;
            6 11]

    (FToE, FToLF, EToO, EToS) = connectivityarrays_3d(EToV, EToF)


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
   # (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

    # This is the base mesh size in each dimension
    N2 = N1 = N0 = 17

    # EToN0 is the base mesh size (e.g., before refinement)
    EToN0 = zeros(Int64, 3, nelems)
    EToN0[1, :] .= N0
    EToN0[2, :] .= N1
    EToN0[3, :] .= N2


    # create arrays for storing errors and rates of convergence
    ϵ = zeros(2)
    rates = zeros(length(ϵ)-1)

    # successively refine 
    for lvl = 1:length(ϵ)

       
        # Set up the local grid dimensions (i.e. number of nodes in each direction on each element)
        Nq = EToN0[1, :] * (2^(lvl-1))
        Nr = EToN0[2, :] * (2^(lvl-1))
        Ns = EToN0[3, :] * (2^(lvl-1))

        # solve for displacements u and compute exact solution uexact. lop are local ops and vstars are indices for each element.
        (u, uexact, lop, vstarts) = elasticity_solve(Nq, Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)
                                                   
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