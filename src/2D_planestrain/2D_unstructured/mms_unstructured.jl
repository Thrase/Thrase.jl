# On an unstructured mesh of quads:
# Solves a manufactured solution
# Geometry set up to resemble BP3 (or BP9-2D)

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

include("ops_BP3-QD_unstructured.jl")
include("../utils_2D.jl")

const ϵ = 1
const γ = 0.1
const α = 0.2
const β = 0.3
const ν = 0.4
const c = 0
const c1 = 0
const b1 = 1
const b2 = 1
const L = 40

# manufactured solution - need to test for discontinouous coefficients

function exact_mu(x, z)
    # Define shear modulus μ
    μ = 2 .+ c .* sin.(b2 .* x ./ L .+ b1 .* z ./ L)
    return μ
end

function exact_mu_x(x, z)
    # Define shear modulus derivative μ_x
    μ_x = b2 .* c ./ L .* cos.(b2 .* x ./ L .+ b1 .* z ./ L)
    return μ_x
end

function exact_mu_z(x, z)
    # Define shear modulus derivative μ_z
    μ_z = b1 .* c ./ L .* cos.(b2 .* x ./ L .+ b1 .* z ./ L)
    return μ_z
end



function exact_lambda(x, z)
    lam = 1 .+ c1 .* cos.(b2 .* x ./ L .+ b1 .* z ./ L)
    return lam
end

function exact_lambda_x(x, z)
    lam_x = -c1 .* (b2 ./ L) .* sin.(b2 .* x ./ L .+ b1 .* z ./ L)
    return lam_x
end

function exact_lambda_z(x, z)
    
    lam_z = -c1 .* (b1 ./ L) .* sin.(b2 .* x ./ L .+ b1 .* z ./ L)

    return lam_z
end

function exact_nu(x, z) # nu should range strictly between 0 and 1/2 
    # Define Poisson ratio
    nu = exact_lambda ./ (2 .* exact_lambda .+ 2 .* exact_mu)
    return nu
end


function exactu(x, z, e, EToDomain)
    if EToDomain[e] == 1
        uexact = 0 .+ β .* x .+ γ .* z .+ ϵ .* sin.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    else
        uexact = 0 .+ ν .* x .+ α .* z .+ ϵ .* cos.(pi * x ./ L) .* cos.(b1 .* 2 .*pi * z ./ L) # never getting here.
    end
    return uexact[:]
end

function exactu_x(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return β .+ (pi ./ L) * ϵ .* cos.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    else
        return ν .- (pi ./ L) * ϵ .* sin.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    end
end
  
function exactu_xx(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ -pi^2 ./ L^2 * ϵ .* sin.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    else
        return 0 .- pi^2 ./ L^2 * ϵ .* cos.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    end
end

function exactu_z(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ γ  .+ 2 .* pi ./ L * ϵ .* b1 .* sin.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    else
        return 0 .+ α  .+ -2 .* pi ./ L * ϵ .* b1 .* cos.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    end
end
  
function exactu_zz(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0  .+ -4 .* pi^2 ./ L^2 * ϵ .* b1 .^ 2 .* sin.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    else
        return 0 .+ -4 .* pi^2 ./ L^2 * ϵ .* b1 .^ 2 .* cos.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    end       
end

function exactu_xz(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ (2*pi/L) .* (pi ./ L) * ϵ .* b1 .* cos.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    else
        return 0 .+ (2*pi/L) * (pi ./ L) * ϵ .* b1 .* sin.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    end
end


function exactw(x, z, e, EToDomain)
    if EToDomain[e] == 1
        wexact = 0 .+ β .* x .+ γ .* z .+ ϵ .* cos.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    else
        wexact = 0 .+ ν .* x .+ α .* z .+ ϵ .* sin.(pi * x ./ L) .* cos.(b1 .* 2 .*pi * z ./ L) # never getting here.
    end
    return wexact[:]
end

function exactw_x(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return β .+ (-pi ./ L) * ϵ .* sin.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    else
        return ν .+ (pi ./ L) * ϵ .* cos.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    end
end
  
function exactw_xx(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ -pi^2 ./ L^2 * ϵ .* cos.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    else
        return 0 .- pi^2 ./ L^2 * ϵ .* sin.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    end
end
  
function exactw_z(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ γ  .+ -2 .* pi ./ L * ϵ .* b1 .* cos.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    else
        return 0 .+ α  .+ -2 .* pi ./ L * ϵ .* b1 .* sin.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    end
end
  
function exactw_zz(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0  .+ -4 .* pi^2 ./ L^2 * ϵ .* b1 .^ 2 .* cos.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    else
        return 0 .+ -4 .* pi^2 ./ L^2 * ϵ .* b1 .^ 2 .* sin.(pi * x ./ L) .* cos.(b1 .* 2 .* pi * z ./ L)
    end       
end

function exactw_xz(x, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .+ (2*pi/L) * (pi ./ L) * ϵ .* b1 .* sin.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    else
        return 0 .- (2*pi/L) * (pi ./ L) * ϵ .* b1 .* cos.(pi * x ./ L) .* sin.(b1 .* 2 .* pi * z ./ L)
    end
end

function source(x, z, e, EToDomain)
    c = exact_lambda(x, z) .+ 2 .* exact_mu(x, z)
    cx = exact_lambda_x(x, z) .+ 2 .* exact_mu_x(x, z)
    cz = exact_lambda_z(x, z) .+ 2 .* exact_mu_z(x, z)

    lam = exact_lambda(x, z)
    lamx = exact_lambda_x(x, z)
    lamz = exact_lambda_z(x, z)

    mu = exact_mu(x, z)
    mux = exact_mu_x(x, z)
    muz = exact_mu_z(x, z)
   
    ux = exactu_x(x, z, e, EToDomain)
    uz = exactu_z(x, z, e, EToDomain)
    uxx = exactu_xx(x, z, e, EToDomain)
    uxz = exactu_xz(x, z, e, EToDomain)
    uzz = exactu_zz(x, z, e, EToDomain)

    wx = exactw_x(x, z, e, EToDomain)
    wz = exactw_z(x, z, e, EToDomain)
    wxx = exactw_xx(x, z, e, EToDomain)
    wxz = exactw_xz(x, z, e, EToDomain)
    wzz = exactw_zz(x, z, e, EToDomain)

    f1 = -(cx .* ux .+ c .* uxx .+ lamx .* wz .+ lam .* wxz .+ muz .* uz .+ mu .* uzz .+ muz .* wx .+ mu .* wxz)
    f2 = -(mux .* uz .+ mu .* uxz .+ mux .* wx .+ mu .* wxx .+ lamz .* ux .+ lam .* uxz .+ cz .* wz .+ c .* wzz)
    
    #.- (exact_mu(x, z) + exact_lambda(x, z)) .* exactw_xz(x, z, e, EToDomain) .- exact_mu(x, z) .* exactu_zz(x, z, e, EToDomain)  # source data
    #f2 = -exact_mu(x, z) .* exactw_xx(x, z, e, EToDomain) .- (exact_mu(x, z) + exact_lambda(x, z)) .* exactu_xz(x, z, e, EToDomain)  .- (exact_lambda(x, z) + 2*exact_mu(x, z)) .* exactw_zz(x, z, e, EToDomain)  # source data

    f = [f1[:] f2[:]]
    return f
end

function exact(x, z, e, EToDomain)
    e1 = exactu(x, z, e, EToDomain)
    e2 = exactw(x, z, e, EToDomain)

    return [e1 e2]
end

function exact_traction(x, y, nx, ny, e, EToDomain)


    c = exact_lambda(x, y) + 2 * exact_mu(x, y)
    mu = exact_mu(x, y)
    lam = exact_lambda(x, y)

    sig_xx = c .*  exactu_x(x, y, e, EToDomain) .+ lam .*  exactw_z(x, y, e, EToDomain) 
    sig_xz = mu .* exactu_z(x, y,  e, EToDomain) .+ mu .*  exactw_x(x, y, e, EToDomain) 
    sig_zz = lam .*  exactu_x(x, y,  e, EToDomain) .+ c .* exactw_z(x, y,  e, EToDomain) 
        

    T1 = nx .* sig_xx + ny .* sig_xz
    T2 = nx .* sig_xz + ny .* sig_zz

    T = [T1[:] T2[:]]

    return T 

end


function planestrain_solve(Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)

    nelems = length(Nr)
    
    ################################## COORDINATE TRANSFORM ###################################
    #
    # Build the local volume operators
    #

    # Create an empty dictionary to store the operators;
    # index via element number, return structure containing local ops.
    OPTYPE = typeof(locoperator(2, 16, 16, exact_mu, exact_lambda))
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

        metrics = create_metrics(Nr[e], Ns[e], exact_mu, exact_lambda, xt, zt) # create coordinate transform


        ###################################################################### 
        # create local finite difference operators on computational domain:
        lop[e] = locoperator(SBPp, Nr[e], Ns[e], exact_mu, exact_lambda, metrics, FToB[EToF[:, e]]) 
    end

  
    # Assemble the global volume operator and factor:
    # HERE:
    Mo = global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
    M = lu(Mo)

    # Get a unique array indices for the faces corresponding to the fault/jump interface
    FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT), Nr, Ns)
   
    # Compute the total number of volume and fault/jump (δ) points
    VNp = vstarts[nelems+1]-1
    δNp = FToδstarts[nfaces+1]-1
    
    # define the data for the jump in displacement on fault interface
    δ = zeros(δNp, 2)
    for f = 1:nfaces
      if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (em, ep) = FToE[:, f]
        (lfm, lfp) = FToLF[:, f]
        (xf, yf) = lop[em].facecoord

        exactm = exact(xf[lfm], yf[lfm], em, EToDomain)
        exactp = exact(xf[lfm], yf[lfm], ep, EToDomain)
        um = exactm[:, 1]
        wm = exactm[:, 2]
        up = exactp[:, 1]
        wp = exactp[:, 2]


        # always define δ as the displacement on minus side minus displacement on plus, to match A&D
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1), 1] = um - up
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1), 2] = wm - wp
            
        
      end
    end
  
    # define the data for the jump in traction data on fault interface
    θ = zeros(δNp, 2)
    for f = 1:nfaces
        if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (em, ep) = FToE[:, f]
        (lfm, lfp) = FToLF[:, f]
        (xf, yf) = lop[em].facecoord

        cm = exact_lambda(xf[lfm], yf[lfm]) .+ 2 .* exact_mu(xf[lfm], yf[lfm])
        mum = exact_mu(xf[lfm], yf[lfm])
        lamm = exact_lambda(xf[lfm], yf[lfm])

        cp = exact_lambda(xf[lfm], yf[lfm]) .+ 2 .* exact_mu(xf[lfm], yf[lfm])
        mup = exact_mu(xf[lfm], yf[lfm])
        lamp = exact_lambda(xf[lfm], yf[lfm])

        sig_xx_minus = cm .*  exactu_x(xf[lfm], yf[lfm], em, EToDomain) .+ lamm .*  exactw_z(xf[lfm], yf[lfm], em, EToDomain) 
        sig_xx_plus = cp .*  exactu_x(xf[lfm], yf[lfm], ep, EToDomain) .+ lamp .*  exactw_z(xf[lfm], yf[lfm], ep, EToDomain) 

        sig_xz_minus = mum .*  exactu_z(xf[lfm], yf[lfm], em, EToDomain) .+ mum .*  exactw_x(xf[lfm], yf[lfm], em, EToDomain) 
        sig_xz_plus =  mup .*  exactu_z(xf[lfm], yf[lfm], ep, EToDomain) .+ mup .*  exactw_x(xf[lfm], yf[lfm], ep, EToDomain) 


        sig_zz_minus = lamm .*  exactu_x(xf[lfm], yf[lfm], em, EToDomain) .+ cm .*  exactw_z(xf[lfm], yf[lfm], em, EToDomain) 
        sig_zz_plus =  lamp .*  exactu_x(xf[lfm], yf[lfm], ep, EToDomain) .+ cp .*  exactw_z(xf[lfm], yf[lfm], ep, EToDomain) 
        
    

        T1m = lop[em].nx[lfm] .* sig_xx_minus + lop[em].ny[lfm] .* sig_xz_minus
        T2m = lop[em].nx[lfm] .* sig_xz_minus + lop[em].ny[lfm] .* sig_zz_minus

        T1p = lop[ep].nx[lfp] .* sig_xx_plus + lop[ep].ny[lfp] .* sig_xz_plus
        T2p = lop[ep].nx[lfp] .* sig_xz_plus + lop[ep].ny[lfp] .* sig_zz_plus
        
        @views θ[FToδstarts[f]:(FToδstarts[f+1]-1), 1] = T1m + T1p
                
        @views θ[FToδstarts[f]:(FToδstarts[f+1]-1), 2] = T2m + T2p

        end
    end

    # define the interface jump in displacement function
    in_jump      = (lf, x, y, e, δ, θ) -> begin
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

    # define the interface jump in traction fucntion
    in_tractionjump      = (lf, x, y, e, δ, θ) -> begin
      f = EToF[lf, e]  # Get global face number
      if EToS[lf, e] == 1  # check if face is on minus side
        if EToO[lf, e]     # check is correct orientation
          return θ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else
          error("shouldn't get here")  # this is bcasue "correct" orientation is always true of a face on teh minus side. 
        end
      else                  # face on plus side
        if EToO[lf, e]      # check if orientation is correct
          return  θ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else                # if orientation is reversed, reverse the data
          return  θ[(FToδstarts[f+1]-1):-1:FToδstarts[f], :]
        end
      end
    end

    # Set boundary data:
    bc_Dirichlet = (lf, x, y, e, δ, θ) -> exact(x, y, e, EToDomain)
    bc_Neumann   = (lf, x, y, nx, ny, e, δ, θ) -> exact_traction(x, y, nx, ny, e, EToDomain)
    
                                              
    # initialize RHS array that stores boundary data (linear system will be Au = b, where b = B*g)
    b = zeros(VNp, 2)
    uexact = zeros(VNp, 2)
    Jf = zeros(VNp, 2)

    
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
    
   
    
    loc_bdry_vec!((@view b[vstarts[e]:vstarts[e+1]-1, :]), lop[e], neighborZ, FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, in_tractionjump, (e, δ, θ))
    uexact[vstarts[e]:vstarts[e+1]-1, :] = exact(lop[e].coord[1], lop[e].coord[2], e, EToDomain)
    locsourcearray!((@view Jf[vstarts[e]:vstarts[e+1]-1, :]), source, lop[e].coord, lop[e].J, (e, EToDomain))

    end
    
    # solve linear system with a backsolve to obtain displacement.
 
    u = M \ (b[:]-Jf[:]) 
    
    

    return (u, uexact[:], lop, vstarts, Mo)

end


let 
    # SBP interior order of accuracy
    SBPp = 2

    # This actually has to go with .inp file. Mesh file side set type to actually boundary condition type. Intially, every face is locked interface unless specified otherwise in sideset info.
    #REbc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN, BC_NEUMANN, BC_NEUMANN, BC_JUMP_INTERFACE, BC_JUMP_INTERFACE]
    # read in mesh from an .inp file and put in boundary condition type
    
   #(verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/BP1_v1.inp")

   #uncomment for one block
    # verts = [-10.5 9.9 -10.1 10.5;
    #         -10.2 -9.8 9.95 10];
    # #verts = 1 .* [-1 1.5 -1.5 1;
    #  #   -1 -2 0.5 1];
    # FToB = [1; 1; 1; 1];
    # #  same orientation
    # EToV = [1;2;3;4];
    # EToF = [1; 
    #         2;
    #         3;
    #         4];
    # EToDomain = [1]

   #uncomment for two blocks with same orientation
#  verts = [-120 1 100 -111 1.4 100;
#              1 -1 2 100 102 97];

#    #          verts = [-10 0 10 -10 0 10;
#    #      0 0 0 10 10 10];

#    FToB = [1; 7; 1; 2; 2; 2; 2];
#    #FToB = [1; 1; 1; 1; 1; 1; 1];


#    #  same orientation
#     EToV = [1 2; 2 3; 4 5; 5 6];
#     EToF = [1 2; 
#             2 5;
#             4 3;
#             6 7];

#   #  different orientation
#     # EToV = [1 6; 2 5; 4 3; 5 2];
#     # EToF = [1 5; 
#     #         2 2;
#     #         4 7;
#     #         6 3];
    
#   EToDomain = [1 2]

   #uncomment for four blocks with same orientation
#   verts = [-40   0  -40 0  40 40 -40   0  40; 
#            -40 -40    0 0 -40  0 -80 -80 -80];

#     FToB = [1; 7; 1; 2; 2; 0; 0; 1; 8; 1; 2; 2];
# #    #FToB = [1; 1; 1; 1; 1; 1; 1];


#     EToV = [1 2 7 8;
#     2 5 8 9;
#     3 4 1 2;
#     4 6 2 5]

#     EToF = [1 2 8 9;
#     2 3 9 10;
#     6 7 11 12;
#     4 5 6 7]


#    EToDomain = [1 2 1 2]


# 9 blocks:
# (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("../../../meshes/plate_test.inp")

#     # modify a bit to change orientation
#     EToV[:, 7] = [10; 5; 15; 7]
#     EToF[:, 7] = [19; 18; 20; 13]
#     #Can be misoriented, but faces and vertices need to follow! 

   
#     # modify a bit to put fault in - this info should really be in .inp file tho
#     EToDomain[1] = 2
#     EToDomain[4] = 2
#     EToDomain[7] = 2
    
#     #FToB = [1, 7, 2, 0, 0, 2, 0, 1, 2, 0, 1, 7, 0, 0, 0, 1, 0, 1, 7, 2, 0, 2, 1, 2]
#      FToB = [1, 7, 2, 0, 0, 2, 0, 1, 2, 0, 1, 7, 7, 7, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2]
   
 
#     (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)
    
    psi = 30
   
    xd = 20*cosd(psi)
    yd = 20*sind(psi)
    verts = [-20+xd  0+xd    -20  0  20+xd      20  -20+2*xd      0+2*xd     20+2*xd; 
                -yd  -yd       0    0    -yd     0  -2*yd  -2*yd  -2*yd];

    FToB = [1; 7; 1; 2; 2; 0; 0; 1; 8; 1; 2; 2];
#    #FToB = [1; 1; 1; 1; 1; 1; 1];


    EToV = [1 7 2 8;
    2 8 5 9;
    3 1 4 2;
    4 2 6 5]

    EToF = [1 8 2 9;
    2 9 3 10;
    6 11 7 12;
    4 6 5 7]


   EToDomain = [1 1 2 2]



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
    N1 = N0 = 11

    # EToN0 is the base mesh size (e.g., before refinement)
    EToN0 = zeros(Int64, 2, nelems)
    EToN0[1, :] .= N0
    EToN0[2, :] .= N1


    # create arrays for storing errors and rates of convergence
    ϵ = zeros(6)
    rates = zeros(length(ϵ)-1)

    # successively refine 
    for lvl = 1:length(ϵ)

       
        # Set up the local grid dimensions (i.e. number of nodes in each direction on each element)
        Nr = EToN0[1, :] * (2^(lvl-1))
        Ns = EToN0[2, :] * (2^(lvl-1))

        # solve for displacements u and compute exact solution uexact. lop are local ops and vstars are indices for each element.
        (u, uexact, lop, vstarts, M) = planestrain_solve(Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)
                  
        N = Integer(length(u)/2)

        u1_exact = uexact[1:N]
        u2_exact = uexact[N+1:2*N]

        u1 = u[1:N]
        u2 = u[N+1:2*N]
        if lvl == 1 # plot approximation, exact solution, or error on coarsest grid
            
            better_plot_solution(u2, nelems, vstarts, Nr, Ns, lop) 
            
            better_plot_solution(u2_exact, nelems, vstarts, Nr, Ns, lop)
            better_plot_solution(u2_exact-u2, nelems, vstarts, Nr, Ns, lop)
            #poo
          
        end

        H̃ = blockdiag(lop[1].H̃, lop[1].H̃)
        for e = 2:length(lop)
            H̃ = blockdiag(H̃, lop[e].H̃, lop[e].H̃)
        end
        # A = -H̃ * M
        # @show extrema(real(eigvals(Matrix(A))))
        # @show extrema(imag(eigvals(Matrix(A))))
        # @show extrema(A - A')
        # spy(A)
        # spy(A-A')
        # poo
        diff = u .- uexact

        ϵ[lvl] = sqrt(diff' * H̃ * diff) ./ sqrt(uexact' * H̃ * uexact)   # relative error
        
    end

    for lvl = 2:length(ϵ)
        rates[lvl-1] = log2(ϵ[lvl-1]/ϵ[lvl])
    end

    @show ϵ
    @show rates
end