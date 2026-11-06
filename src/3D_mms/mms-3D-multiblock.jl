# On an single block:
# Solves a manufactured solution
# Geometry set up to resemble BP3-3D 

using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using DelimitedFiles
using PyPlot
using IterativeSolvers


const year_seconds = 31556926
global const ctr = Ref{Int64}(1) 

const RS_FAULT = 7
const VP_FAULT = 8

include("ops_BP9-QD-3D_unstructured.jl")
include("utils_3D.jl")
include("transfinite.jl")
include("simple_exact_sol.jl")



function solve_3D(Nq, Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)

        
        nelems = length(Nq)

        ################################## COORDINATE TRANSFORM ###################################
        #
        # Build the local volume operators
        #

        # Create an empty dictionary to store the operators;
        # index via element number, return structure containing local ops.
        OPTYPE = typeof(locoperator(2, 16, 16, 16, exact_mu, exact_lambda))
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
            zcorners = verts[3, EToV[:, e]]

           (xt, yt, zt) = transfinite_blend_3D(xcorners, ycorners, zcorners)
            metrics = create_metrics(SBPp, Nq[e], Nr[e], Ns[e], exact_lambda, exact_mu, xt, yt, zt) # create coordinate transform

            ###################################################################### 
            # create local finite difference operators on computational domain:
            lop[e] = locoperator(SBPp, Nq[e], Nr[e], Ns[e], exact_mu, exact_lambda, metrics, FToB[EToF[:, e]]) 

        end


    #Calculate neighboring penalty parameters and add to lop 
    for e = 1:nelems
        lop[e].neighborZ  .= calculate_neighbors(lop, e, FToB, EToF, FToE, FToLF, EToO)
    end

    # Assemble the global volume operator:
    Mo = global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nq, Nr, Ns)
   

    # Get a unique array indices for the faces corresponding to the fault/jump interface
    FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT), Nq, Nr, Ns)

    # Compute the total number of volume and fault/jump (δ) points
    VNp = vstarts[nelems+1]-1
    δNp = FToδstarts[nfaces+1]-1

    # Set boundary/interface data:
    δ = zeros(δNp, 3)
    θ = zeros(δNp, 3)

    # Fill δ (displacement jump) and θ (traction jump) from the exact solution.
    # Both arrays are aligned with the minus side.
    for f = 1:nfaces
      if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (em, ep) = FToE[:, f]
        (lfm, lfp) = FToLF[:, f]
        (xfm, yfm, zfm) = lop[em].facecoord
        (xfp, yfp, zfp) = lop[ep].facecoord

        exactm = exact(xfm[lfm], yfm[lfm], zfm[lfm], em, EToDomain)
        exactp = exact(xfp[lfp], yfp[lfp], zfp[lfp], ep, EToDomain)

        um = exactm[:, 1]; vm = exactm[:, 2]; wm = exactm[:, 3]
        up = exactp[:, 1]; vp = exactp[:, 2]; wp = exactp[:, 3]

        ori = EToO[lfp, ep]
        if ori == (true, true, true)
          @views δ[FToδstarts[f]:(FToδstarts[f+1]-1), 1] = up - um
          @views δ[FToδstarts[f]:(FToδstarts[f+1]-1), 2] = vp - vm
          @views δ[FToδstarts[f]:(FToδstarts[f+1]-1), 3] = wp - wm
        else
          @views δ[FToδstarts[f]:(FToδstarts[f+1]-1), 1] = up[end:-1:1] - um
          @views δ[FToδstarts[f]:(FToδstarts[f+1]-1), 2] = vp[end:-1:1] - vm
          @views δ[FToδstarts[f]:(FToδstarts[f+1]-1), 3] = wp[end:-1:1] - wm
        end
      end
    end

    for f = 1:nfaces
      if FToB[f] ∈ (RS_FAULT, VP_FAULT)
        (em, ep) = FToE[:, f]
        (lfm, lfp) = FToLF[:, f]
        (xfm, yfm, zfm) = lop[em].facecoord
        (xfp, yfp, zfp) = lop[ep].facecoord

        nxm = lop[em].nx[lfm]; nym = lop[em].ny[lfm]; nzm = lop[em].nz[lfm]
        nxp = lop[ep].nx[lfp]; nyp = lop[ep].ny[lfp]; nzp = lop[ep].nz[lfp]

        Tm = exact_traction(xfm[lfm], yfm[lfm], zfm[lfm], nxm, nym, nzm, em, EToDomain)
        Tp = exact_traction(xfp[lfp], yfp[lfp], zfp[lfp], nxp, nyp, nzp, ep, EToDomain)

        T1m = Tm[:, 1]; T2m = Tm[:, 2]; T3m = Tm[:, 3]
        T1p = Tp[:, 1]; T2p = Tp[:, 2]; T3p = Tp[:, 3]

        ori = EToO[lfp, ep]
        if ori == (true, true, true)
          @views θ[FToδstarts[f]:(FToδstarts[f+1]-1), 1] = T1m + T1p
          @views θ[FToδstarts[f]:(FToδstarts[f+1]-1), 2] = T2m + T2p
          @views θ[FToδstarts[f]:(FToδstarts[f+1]-1), 3] = T3m + T3p
        else
          @views θ[FToδstarts[f]:(FToδstarts[f+1]-1), 1] = T1m + T1p[end:-1:1]
          @views θ[FToδstarts[f]:(FToδstarts[f+1]-1), 2] = T2m + T2p[end:-1:1]
          @views θ[FToδstarts[f]:(FToδstarts[f+1]-1), 3] = T3m + T3p[end:-1:1]
        end
      end
    end

    # δ and θ are the global data, and they are always aligned with minues side
    # in_jump tranforms the data locally to match local orientation
    # in_jump and in_traction jump just take the global data and negate/reverse it if need be.
    # define the interface jump in displacement function
    in_jump      = (lf, x, y, z, e, δ, θ) -> begin
      f = EToF[lf, e]  # Get global face number
      if EToS[lf, e] == 1  # minus side — always canonical orientation
        return -δ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
      else                  # plus side
        ori = EToO[lf, e]   # Tuple{Bool,Bool,Bool}
        @show ori
        if ori == (true, true, true)
          return  δ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else                # reversed orientation — full 2D face permutation
          em = FToE[1, f]
          fm = FToLF[1, f]  # local face of minus element — determines face grid size
          
          # face dimensions in minus-side ordering
          if fm <= 2
            n1, n2 = Nr[em]+1, Ns[em]+1
          elseif fm <= 4
            n1, n2 = Nq[em]+1, Ns[em]+1
          else
            n1, n2 = Nq[em]+1, Nr[em]+1
          end
          
          perm = zeros(Int, n1*n2)
          for jj in 1:n2, ii in 1:n1
            k = ii + n1*(jj - 1)
            if ori[1]
              ni = ori[2] ? ii : n1 - ii + 1
              nj = ori[3] ? jj : n2 - jj + 1
              perm[k] = ni + n1*(nj - 1)
            else
              ni = ori[2] ? jj : n2 - jj + 1
              nj = ori[3] ? ii : n1 - ii + 1
              perm[k] = ni + n2*(nj - 1)
            end
          end
          
          return δ[FToδstarts[f] .+ invperm(perm) .- 1, :]
        end
      end
    end

    # define the interface jump in traction function
    in_tractionjump      = (lf, x, y, z, e, δ, θ) -> begin
      f = EToF[lf, e]  # Get global face number
      if EToS[lf, e] == 1  # minus side — always canonical orientation
        return θ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
      else                  # plus side
        ori = EToO[lf, e]   # Tuple{Bool,Bool,Bool}
        if ori == (true, true, true)
          return  θ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else                # reversed orientation — full 2D face permutation
          em = FToE[1, f]
          fm = FToLF[1, f]  # local face of minus element — determines face grid size
          
          # face dimensions in minus-side ordering
          if fm <= 2
            n1, n2 = Nr[em]+1, Ns[em]+1
          elseif fm <= 4
            n1, n2 = Nq[em]+1, Ns[em]+1
          else
            n1, n2 = Nq[em]+1, Nr[em]+1
          end
          
          perm = zeros(Int, n1*n2)
          for jj in 1:n2, ii in 1:n1
            k = ii + n1*(jj - 1)
            if ori[1]
              ni = ori[2] ? ii : n1 - ii + 1
              nj = ori[3] ? jj : n2 - jj + 1
              perm[k] = ni + n1*(nj - 1)
            else
              ni = ori[2] ? jj : n2 - jj + 1
              nj = ori[3] ? ii : n1 - ii + 1
              perm[k] = ni + n2*(nj - 1)
            end
          end
          
          return θ[FToδstarts[f] .+ invperm(perm) .- 1, :]
        end
      end
    end
       
    bc_Dirichlet = (lf, x, y, z, e, δ, θ) -> exact(x, y, z, e, EToDomain)
    bc_Neumann   = (lf, x, y, z, nx, ny, nz, e, δ, θ) -> exact_traction(x, y, z, nx, ny, nz, e, EToDomain)
    

    b = zeros(VNp, 3)
    uexact = zeros(VNp, 3)
    Jf = zeros(VNp, 3)
    
      for e = 1:nelems
        loc_bdry_vec!((@view b[vstarts[e]:vstarts[e+1]-1, :]), lop[e], FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, in_tractionjump, (e, δ, θ))
        uexact[vstarts[e]:vstarts[e+1]-1, :] = exact(lop[e].coord[1], lop[e].coord[2], lop[e].coord[3], e, EToDomain)
        locsourcearray!((@view Jf[vstarts[e]:vstarts[e+1]-1, :]), source, lop[e].coord, lop[e].J, (e, EToDomain))

    end
    
    # Direct Solve:
    # M = lu(Mo)
    # u = M \ (b[:]-Jf[:]) 

   

    # Iterative Solve:
    Hblock = blockdiag([lop[e].H̃ for e in 1:length(lop)]...)
    H̃ = blockdiag(Hblock, Hblock, Hblock)
    A = -H̃*Mo 
    # To check for SPD:
    # @show extrema(A-A')
    # F = eigen(Matrix(A))
    # @show F.values
    rhs = -H̃*(b[:]-Jf[:])
    u, stats = IterativeSolvers.cg(A, rhs; reltol=1e-12, maxiter=10000, log=true)
  
    @show stats
    return (u, uexact[:], lop)

  end


let 

    # SBP interior order of accuracy
    SBPp = 2

    (verts, EToV, EToF, FToB, EToDomain) = read_inp_3d("src/3D/mesh3D_32blocks_30degreeDip.inp")
   
    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    @show (nelems, nfaces)

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


    # N_base is the base mesh size in each dimension, for all the elements
    N_base = 16
    Nq = Nr = Ns =  Vector{Int64}(N_base*ones(nelems))

    (u, uexact, lop) = solve_3D(Nq, Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)

    N = Integer(length(u)/3)

    u1_exact = uexact[1:N]
    u2_exact = uexact[N+1:2*N]
    u3_exact = uexact[2N+1:3*N]

    u1 = u[1:N]
    u2 = u[N+1:2*N]
    u3 = u[2N+1:3*N]


    Hblock = blockdiag([lop[e].H̃ for e in 1:length(lop)]...)
    H̃ = blockdiag(Hblock, Hblock, Hblock)


    diff = u .- uexact
     
    @show sqrt(diff' * H̃ * diff) ./ sqrt(uexact' * H̃ * uexact)   # relative error
        

    
end