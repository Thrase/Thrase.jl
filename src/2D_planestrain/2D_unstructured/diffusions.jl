using Plots
using SparseArrays
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks

const BC_DIRICHLET       = 1
const BC_LOCKED_INTERFACE = 0
const BC_NEUMANN         = 2
const EARTH_SURFACE      = 3


include("../../diagonal_sbp.jl")




function calculate_neighbors_diffusion(lop, e, FToB, EToF, FToE, FToLF)

  neighborZ = [similar(lop[e].sJZ[1]), similar(lop[e].sJZ[2])]
  
  LFToB = FToB[EToF[:,e]]            
  # loop over the faces to get the neighboring penalty parameters:
  for lf = 1:2
      if LFToB[lf] == BC_LOCKED_INTERFACE 
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


        neighborZ[lf] = lop[eo].sJZ[nf]
        
        end
    end
  return neighborZ
end



function create_diffusion_metrics(N, exact_c, xf=(r)->(r, ones(size(r))))

  Np = N + 1
  r = range(-1, stop=1, length=Np)

  (x, xr) = xf(r)

  c = exact_c(x)
  J = xr

  @assert minimum(J) > 0

  JI = 1 ./ J
  rx = JI

  # variable coefficient matrix components
  crr = J .* (rx .* c .* rx)

  # surface matrices
  xf1 = x[1]       # FIX: view(x[1]) → x[1]  (view needs index args; scalar is fine here)
  nx1 = -xr[1]
  sJ1 = abs.(nx1)
  nx1 = nx1 ./ sJ1

  xf2 = x[Np]      # FIX: (xf2, yf2) = view(x[Np]) → xf2 = x[Np]  (can't destructure scalar; no yf2 in 1D)
  nx2 = xr[end]
  sJ2 = abs.(nx2)
  nx2 = nx2 ./ sJ2

  (coord = x,
  facecoord = (xf1, xf2),
  crr = crr, J = J, JI = JI,
  sJ = (sJ1, sJ2),
  nx = (nx1, nx2),   # FIX: removed ny = (ny1, ny2) — ny1/ny2 undefined and meaningless in 1D
  rx = rx)
end

function loc_bdry_vec_diffusion!(ge, t, lop, LFToB, bc_Dirichlet, bc_Neumann, bcargs = ())

  dB = lop.dB
  nB = lop.nB

  xf = lop.facecoord
  sJ = lop.sJ
  nx = lop.nx
  ge[:] .= 0

  for lf = 1:2
    if LFToB[lf] == BC_DIRICHLET
      gD = bc_Dirichlet(lf, t, xf[lf], bcargs...)
      ge[:] .+= vec(dB[lf]*gD)

    elseif LFToB[lf] == BC_NEUMANN
     
      gN = bc_Neumann(lf, t, xf[lf], nx[lf], bcargs...)
      ge[:] .+= vec(nB[lf] * (sJ[lf] .* gN))
     

    elseif LFToB[lf] == BC_LOCKED_INTERFACE
        # no data since locked
    else
      error("invalid bc")
    end
    
  end

end

function loc_diffusion_ops(p, N, c, diffusion_metrics = create_diffusion_metrics(N, c),
             LFToB = (BC_DIRICHLET, BC_DIRICHLET); C = diffusion_metrics.crr)

    Np = N+1
    J = spdiagm(0 => reshape(diffusion_metrics.J, Np))  # FIX: Np*Np → Np  (1D: J has Np elements)
    JI = spdiagm(0 => reshape(diffusion_metrics.JI, Np)) 

    (D, HI, H, r) = diagonal_sbp_D1(p, N)
    D = JI * D
    
    (A, CS0, CSN, HI, H, _, _) = variable_diagonal_sbp_D2(p, N, C)
    A = JI * A

    Id = sparse(I, Np, Np)

    e1 = sparse([1  ], [1], [1], Np, 1)
    e2 = sparse([Np], [1], [1], Np, 1)
    e  = (e1, e2)
    eT = (e1', e2')

    Cmat = spdiagm(0 => C[:])  

    d  = 1
    h1 = H[1, 1]
    h2 = H[end,end]

    β  = 10
    sJZ_1 = β * (d ./ h1) * Cmat
    sJZ_2 = β * (d ./ h2) * Cmat
    sJZ = (sJZ_1, sJZ_2)

    T1 = -CS0
    T2 = CSN
    T  = (T1, T2)
    TT = (T1', T2')

    # Dirichlet SAT operators (added to A, and stored for RHS boundary vector)
    dSAT_1 = JI * HI * (T1 .- sJZ_1)' * e1 * e1'
    dSAT_2 = JI * HI * (T2 .- sJZ_2)' * e2 * e2'
    dSAT = (dSAT_1, dSAT_2)

    # Neumann SAT operators
    nSAT_1 = -JI * HI * e1 * e1' * T1
    nSAT_2 = -JI * HI * e2 * e2' * T2
    nSAT = (nSAT_1, nSAT_2)

    # RHS boundary contribution matrices
    dB1 = JI * HI * (T1 .- sJZ_1)' * e1
    dB2 = JI * HI * (T2 .- sJZ_2)' * e2
    dB = (dB1, dB2)

    nB1 = -JI * HI * e1
    nB2 = -JI * HI * e2
    nB = (nB1, nB2)

    neighborZ = [similar(sJZ_1), similar(sJZ_2)]

    for lf = 1:2
        if LFToB[lf] == BC_DIRICHLET
            A += dSAT[lf]

        elseif (LFToB[lf] == BC_NEUMANN || LFToB[lf] == EARTH_SURFACE)
            A += nSAT[lf]
            

        elseif LFToB[lf] == BC_LOCKED_INTERFACE   
            # do in global assembly
        else
            print("shouldn't get here")
        end
    end

    return (A = A, D = D, dB = dB, nB = nB, e = e, eT = eT, H = H, HI = HI, T = T, TT = TT, J = J, JI = JI,
            coord      = diffusion_metrics.coord,
            facecoord  = diffusion_metrics.facecoord,
            sJ         = diffusion_metrics.sJ,
            nx         = diffusion_metrics.nx,
            sJZ        = sJZ,
            neighborZ  = neighborZ)
end

function compute_glob_diff_op(lop, fstarts, FToB, FToE, FToLF)

    n_elems = length(lop)
    n_faces = length(FToB)

    Nxd = fstarts[n_elems+1] - 1
    M = spzeros(Nxd, Nxd)

    D = spzeros(Nxd, Nxd)



    for e = 1:n_elems
        idx = fstarts[e]:fstarts[e+1]-1   
 
        M[idx, idx] = lop[e].A
        D[idx, idx] = lop[e].D
    end

    for f = 1:n_faces

        if FToB[f] == BC_LOCKED_INTERFACE   
            (em, ep) = FToE[:, f]
            (fm, fp) = FToLF[:, f]

            sJZm = 0.5 .* ((lop[em].sJZ[fm])' + (lop[em].neighborZ[fm])')
            sJZp = 0.5 .* ((lop[ep].sJZ[fp])' + (lop[ep].neighborZ[fp])')

            Tm  = lop[em].T[fm]
            Tp  = lop[ep].T[fp]
            TTm = lop[em].TT[fm]
            TTp = lop[ep].TT[fp]

            idxm = fstarts[em]:fstarts[em+1]-1
            idxp = fstarts[ep]:fstarts[ep+1]-1

            JI_em = lop[em].JI
            JI_ep = lop[ep].JI

            M[idxm, idxm] += JI_em * (-lop[em].HI * (sJZm - 0.5 * TTm) * lop[em].e[fm] * lop[em].eT[fm] +
                               -0.5 * lop[em].HI * lop[em].e[fm] * lop[em].eT[fm] * Tm)
            M[idxm, idxp] += JI_em *  (lop[em].HI * (sJZm - 0.5 * TTm) * lop[em].e[fm] * lop[ep].eT[fp] +
                               -0.5 * lop[em].HI * lop[em].e[fm] * lop[ep].eT[fp] * Tp)
            M[idxp, idxm] += JI_ep *  (lop[ep].HI * (sJZp - 0.5 * TTp) * lop[ep].e[fp] * lop[em].eT[fm] +
                               -0.5 * lop[ep].HI * lop[ep].e[fp] * lop[em].eT[fm] * Tm)
            M[idxp, idxp] += JI_ep * (-lop[ep].HI * (sJZp - 0.5 * TTp) * lop[ep].e[fp] * lop[ep].eT[fp] +
                               -0.5 * lop[ep].HI * lop[ep].e[fp]  * lop[ep].eT[fp] * Tp)
        end
    end
    return D, M
end


function diffusion_rhs(t, N, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, nfaces, EToDomain, p, c)

    nelems = length(N)

    OPTYPE = typeof(loc_diffusion_ops(2, 16, x -> c(x, 1, EToDomain)))
    lop = Dict{Int64, OPTYPE}()

    fstarts = Array{Int64, 1}(undef, nelems + 1)
    fstarts[1] = 1

    for e = 1:nelems

        Np = N[e]+1
        fstarts[e+1] = fstarts[e] + Np

        (x1, x2) = -verts[2, EToV[:, e]]  # depth = -y, positive downward

        m = (x2 - x1) / 2
        b = (x1 + x2) / 2   # FIX: (x1 + x)/2 → (x1 + x2)/2  (x undefined)
        xt(r) = (m*r .+ b, m * ones(size(r)))


        metrics = create_diffusion_metrics(N[e], x -> c(x, e,  EToDomain), xt)

        lop[e] = loc_diffusion_ops(SBPp, N[e], x -> c(x, e,  EToDomain), metrics, FToB[EToF[:, e]])
    end

    for e = 1:nelems
        lop[e].neighborZ .= calculate_neighbors_diffusion(lop, e, FToB, EToF, FToE, FToLF)        # TODO: calculate_neighbors must be defined (not present in this file)
    end

    M = compute_glob_diff_op(lop, fstarts, FToB, FToE, FToLF)

    fNp = fstarts[nelems+1] - 1

    b = zeros(fNp)
  

    for e = 1:nelems
      

        bc_Dirichlet = (lf, t, x, e, EToDomain) -> 0
        # c = kprime / eta_f
        # Divide by sJ because loc_bdry_vec_diffusion! multiplies gN by sJ[lf].
        # For the 1D fault mesh sJ = J (the volume Jacobian at the endpoint, not
        # a surface-area ratio), so the factor must be cancelled here.
        bc_Neumann   = (lf, t, x, nx, e, EToDomain) ->
            (c(x, e, EToDomain) * p.ρ_f * p.g + p.q0/(p.n*p.β)) / lop[e].sJ[lf]
        bcargs = (e, EToDomain)

        LFToB = FToB[EToF[:, e]]

        loc_bdry_vec_diffusion!((@view b[fstarts[e]:fstarts[e+1]-1]), t, lop[e],
                                 LFToB, bc_Dirichlet, bc_Neumann, bcargs)
   

    end

    return M, b, lop, fstarts

end


#{{{ connectivityarrays
function connectivityarrays_1D(EToV, EToF)
  nelems = size(EToV, 2)
  nfaces = maximum(maximum(EToF))

  FToE  = zeros(Int64, 2, nfaces)
  FToLF = zeros(Int64, 2, nfaces)

  LFToLV = [1, 2]   
  for e = 1:nelems
    for lf = 1:2
      gf = EToF[lf, e]
      if FToE[1, gf] == 0
        @assert FToLF[1, gf] == 0
        FToE[1, gf] = e
        FToLF[1, gf] = lf
      else
        @assert FToE[2, gf] == 0
        @assert FToLF[2, gf] == 0
        FToE[2, gf] = e
        FToLF[2, gf] = lf
      end
    end
  end
  (FToE, FToLF)
end
#}}}


export diffusion_rhs