using SparseArrays
using LinearAlgebra


⊗(A,B) = kron(A, B)

#{{{ Transfinite Blend
function transfinite_blend(α1, α2, α3, α4, α1s, α2s, α3r, α4r, r, s)
  # +---4---+
  # |       |
  # 1       2
  # |       |
  # +---3---+
  @assert [α1(-1) α2(-1) α1( 1) α2( 1)] ≈ [α3(-1) α3( 1) α4(-1) α4( 1)]


  x = (1 .+ r) .* α2(s)/2 + (1 .- r) .* α1(s)/2 +
      (1 .+ s) .* α4(r)/2 + (1 .- s) .* α3(r)/2 -
     ((1 .+ r) .* (1 .+ s) .* α2( 1) +
      (1 .- r) .* (1 .+ s) .* α1( 1) +
      (1 .+ r) .* (1 .- s) .* α2(-1) +
      (1 .- r) .* (1 .- s) .* α1(-1)) / 4

  xr =  α2(s)/2 - α1(s)/2 +
        (1 .+ s) .* α4r(r)/2 + (1 .- s) .* α3r(r)/2 -
      (+(1 .+ s) .* α2( 1) +
       -(1 .+ s) .* α1( 1) +
       +(1 .- s) .* α2(-1) +
       -(1 .- s) .* α1(-1)) / 4


  xs = (1 .+ r) .* α2s(s)/2 + (1 .- r) .* α1s(s)/2 +
       α4(r)/2 - α3(r)/2 -
      (+(1 .+ r) .* α2( 1) +
       +(1 .- r) .* α1( 1) +
       -(1 .+ r) .* α2(-1) +
       -(1 .- r) .* α1(-1)) / 4

  return (x, xr, xs)
end

function transfinite_blend(α1, α2, α3, α4, r, s, p)
  (Nrp, Nsp) = size(r)
  (Dr, _, _, _) = diagonal_sbp_D1(p, Nrp-1; xc = (-1,1))
  (Ds, _, _, _) = diagonal_sbp_D1(p, Nsp-1; xc = (-1,1))

  α2s(s) = α2(s) * Ds'
  α1s(s) = α1(s) * Ds'
  α4r(s) = Dr * α4(r)
  α3r(s) = Dr * α3(r)

  transfinite_blend(α1, α2, α3, α4, α1s, α2s, α3r, α4r, r, s)
end

function transfinite_blend(v1::T1, v2::T2, v3::T3, v4::T4, r, s
                          ) where {T1 <: Number, T2 <: Number,
                                   T3 <: Number, T4 <: Number}
  e1(α) = v1 * (1 .- α) / 2 + v3 * (1 .+ α) / 2
  e2(α) = v2 * (1 .- α) / 2 + v4 * (1 .+ α) / 2
  e3(α) = v1 * (1 .- α) / 2 + v2 * (1 .+ α) / 2
  e4(α) = v3 * (1 .- α) / 2 + v4 * (1 .+ α) / 2
  e1α(α) = -v1 / 2 + v3 / 2
  e2α(α) = -v2 / 2 + v4 / 2
  e3α(α) = -v1 / 2 + v2 / 2
  e4α(α) = -v3 / 2 + v4 / 2
  transfinite_blend(e1, e2, e3, e4, e1α, e2α, e3α, e4α, r, s)
end
#}}}



function create_metrics(Nr, Ns, μ,
  xf=(r,s)->(r, ones(size(r)), zeros(size(r))),
  yf=(r,s)->(s, zeros(size(s)), ones(size(s))))


  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  r = range(-1, stop=1, length=Nrp)
  s = range(-1, stop=1, length=Nsp)

  # Create the mesh
  r = ones(1, Nsp) ⊗ r
  s = s' ⊗ ones(Nrp)

  (x, xr, xs) = xf(r, s)
  (y, yr, ys) = yf(r, s)

  J = xr .* ys - xs .* yr

  @assert minimum(J) > 0

  JI = 1 ./ J

  rx =  ys ./ J
  sx = -yr ./ J
  ry = -xs ./ J
  sy =  xr ./ J

  # variable coefficient matrix components
  crr = J .* (rx .* μ .* rx + ry .* μ .* ry)
  crs = J .* (sx .* μ .* rx + sy .* μ .* ry)
  css = J .* (sx .* μ .* sx + sy .* μ .* sy)

  # surface matrices
  (xf1, yf1) = (view(x, 1, :), view(y, 1, :))
  nx1 = -ys[1, :]
  ny1 =  xs[1, :]
  sJ1 = hypot.(nx1, ny1)
  nx1 = nx1 ./ sJ1
  ny1 = ny1 ./ sJ1

  (xf2, yf2) = (view(x, Nrp, :), view(y, Nrp, :))
  nx2 =  ys[end, :]
  ny2 = -xs[end, :]
  sJ2 = hypot.(nx2, ny2)
  nx2 = nx2 ./ sJ2
  ny2 = ny2 ./ sJ2

  (xf3, yf3) = (view(x, :, 1), view(y, :, 1))
  nx3 =  yr[:, 1]
  ny3 = -xr[:, 1]
  sJ3 = hypot.(nx3, ny3)
  nx3 = nx3 ./ sJ3
  ny3 = ny3 ./ sJ3

  (xf4, yf4) = (view(x, :, Nsp), view(y, :, Nsp))
  nx4 = -yr[:, end]
  ny4 =  xr[:, end]
  sJ4 = hypot.(nx4, ny4)
  nx4 = nx4 ./ sJ4
  ny4 = ny4 ./ sJ4

  (coord = (x,y),
  facecoord = ((xf1, xf2, xf3, xf4), (yf1, yf2, yf3, yf4)),
  crr = crr, css = css, crs = crs,
  J=J,
  JI = JI,
  sJ = (sJ1, sJ2, sJ3, sJ4),
  nx = (nx1, nx2, nx3, nx4),
  ny = (ny1, ny2, ny3, ny4),
  rx = rx, ry = ry, sx = sx, sy = sy)
end

function get_operators(p, Nr, Ns, metrics=create_metrics(p,Nr,Ns),
                     crr = metrics.crr,
                     css = metrics.css,
                     crs = metrics.crs)
  csr = crs
  J = metrics.J

  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  # Derivative operators in logical space
  (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
  (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))

  # Identity matrices for the computation
  Ir = sparse(I, Nrp, Nrp)
  Is = sparse(I, Nsp, Nsp)

  #{{{ Set up the rr derivative matrix
  ISr0 = Array{Int64,1}(undef,0)
  JSr0 = Array{Int64,1}(undef,0)
  VSr0 = Array{Float64,1}(undef,0)
  ISrN = Array{Int64,1}(undef,0)
  JSrN = Array{Int64,1}(undef,0)
  VSrN = Array{Float64,1}(undef,0)

  (D2e, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Nr, rand(Nrp))
  IArr = Array{Int64,1}(undef,Nsp * length(D2e.nzval))
  JArr = Array{Int64,1}(undef,Nsp * length(D2e.nzval))
  VArr = Array{Float64,1}(undef,Nsp * length(D2e.nzval))
  stArr = 0

  ISr0 = Array{Int64,1}(undef,Nsp * length(S0e.nzval))
  JSr0 = Array{Int64,1}(undef,Nsp * length(S0e.nzval))
  VSr0 = Array{Float64,1}(undef,Nsp * length(S0e.nzval))
  stSr0 = 0

  ISrN = Array{Int64,1}(undef,Nsp * length(SNe.nzval))
  JSrN = Array{Int64,1}(undef,Nsp * length(SNe.nzval))
  VSrN = Array{Float64,1}(undef,Nsp * length(SNe.nzval))
  stSrN = 0
  for j = 1:Nsp
    rng = (j-1) * Nrp .+ (1:Nrp)
    (D2e, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Nr, crr[rng])
    (Ie, Je, Ve) = findnz(D2e)
    IArr[stArr .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JArr[stArr .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VArr[stArr .+ (1:length(Ve))] = Ve
    stArr += length(Ve)

    (Ie, Je, Ve) = findnz(S0e)
    ISr0[stSr0 .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JSr0[stSr0 .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VSr0[stSr0 .+ (1:length(Ve))] = Ve
    stSr0 += length(Ve)

    (Ie, Je, Ve) = findnz(SNe)
    ISrN[stSrN .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JSrN[stSrN .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VSrN[stSrN .+ (1:length(Ve))] = Ve
    stSrN += length(Ve)
  end
  # TODO: is all this for Sr0 etc necessayr? can't we just do kronecker with 1D op?
  Drr = sparse(IArr[1:stArr], JArr[1:stArr], VArr[1:stArr], Np, Np)
  Sr0 = sparse(ISr0[1:stSr0], JSr0[1:stSr0], VSr0[1:stSr0], Np, Np)  
  SrN = sparse(ISrN[1:stSrN], JSrN[1:stSrN], VSrN[1:stSrN], Np, Np)
  Sr0T = sparse(JSr0[1:stSr0], ISr0[1:stSr0], VSr0[1:stSr0], Np, Np)
  SrNT = sparse(JSrN[1:stSrN], ISrN[1:stSrN], VSrN[1:stSrN], Np, Np)
  
  #= affine mesh test
  (D2, S0, SN, _, _, _) = diagonal_sbp_D2(p, Nr)
  @assert Drr ≈ Is ⊗ D2

  @assert Sr0 ≈ ((sparse(Diagonal(crr[1   .+ Nrp*(0:Ns)])) * Hs) ⊗ S0)
  @assert SrN ≈ ((sparse(Diagonal(crr[Nrp .+ Nrp*(0:Ns)])) * Hs) ⊗ SN)
 =#

  #{{{ Set up the ss derivative matrix
  (D2e, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Ns, rand(Nsp))
  IAss = Array{Int64,1}(undef,Nrp * length(D2e.nzval))
  JAss = Array{Int64,1}(undef,Nrp * length(D2e.nzval))
  VAss = Array{Float64,1}(undef,Nrp * length(D2e.nzval))
  stAss = 0

  ISs0 = Array{Int64,1}(undef,Nrp * length(S0e.nzval))
  JSs0 = Array{Int64,1}(undef,Nrp * length(S0e.nzval))
  VSs0 = Array{Float64,1}(undef,Nrp * length(S0e.nzval))
  stSs0 = 0

  ISsN = Array{Int64,1}(undef,Nrp * length(SNe.nzval))
  JSsN = Array{Int64,1}(undef,Nrp * length(SNe.nzval))
  VSsN = Array{Float64,1}(undef,Nrp * length(SNe.nzval))
  stSsN = 0
  for i = 1:Nrp
    rng = i .+ Nrp * (0:Ns)
    (D2e, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Ns, css[rng])

    (Ie, Je, Ve) = findnz(D2e)
    IAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VAss[stAss .+ (1:length(Ve))] = Ve
    stAss += length(Ve)

    (Ie, Je, Ve) = findnz(S0e)
    ISs0[stSs0 .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JSs0[stSs0 .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VSs0[stSs0 .+ (1:length(Ve))] = Ve
    stSs0 += length(Ve)

    (Ie, Je, Ve) = findnz(SNe)
    ISsN[stSsN .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JSsN[stSsN .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VSsN[stSsN .+ (1:length(Ve))] = Ve
    stSsN += length(Ve)
  end
  # TODO: is all this for Ss0 etc necessayr? can't we just do kronecker with 1D op?
  Dss = sparse(IAss[1:stAss], JAss[1:stAss], VAss[1:stAss], Np, Np)
  Ss0 = sparse(ISs0[1:stSs0], JSs0[1:stSs0], VSs0[1:stSs0], Np, Np)
  SsN = sparse(ISsN[1:stSsN], JSsN[1:stSsN], VSsN[1:stSsN], Np, Np)
  Ss0T = sparse(JSs0[1:stSs0], ISs0[1:stSs0], VSs0[1:stSs0], Np, Np)
  SsNT = sparse(JSsN[1:stSsN], ISsN[1:stSsN], VSsN[1:stSsN], Np, Np)

  #= affine mesh test
  (D2, S0, SN, _, _, _) = diagonal_sbp_D2(p, Ns)
  @assert Dss ≈ D2 ⊗ Ir
  @assert Ss0 ≈ (S0 ⊗ (Hr * sparse(Diagonal(css[1:Nrp]))))
  @assert SsN ≈ (SN ⊗ (Hr * sparse(Diagonal(css[Nrp*Ns .+ (1:Nrp)]))))
  =#

  #{{{ Set up the sr and rs derivative matrices
  Dsr = (Ds ⊗ Ir) * sparse(1:length(crs), 1:length(crs), view(crs, :)) * (Is ⊗ Dr)
  Drs = (Is ⊗ Dr) * sparse(1:length(csr), 1:length(csr), view(csr, :)) * (Ds ⊗ Ir)
  #}}}

  #
  # Boundary point matrices

  er0 = sparse([1  ], [1], [1], Nrp, 1)
  erN = sparse([Nrp], [1], [1], Nrp, 1)
  es0 = sparse([1  ], [1], [1], Nsp, 1)
  esN = sparse([Nsp], [1], [1], Nsp, 1)

  er0T = sparse([1], [1  ], [1], 1, Nrp)
  erNT = sparse([1], [Nrp], [1], 1, Nrp)
  es0T = sparse([1], [1  ], [1], 1, Nsp)
  esNT = sparse([1], [Nsp], [1], 1, Nsp)

  
  #cmax = maximum([maximum(crr), maximum(crs), maximum(css)])

  #
  # Store boundary coefficients as matrices - waht are these for?
  #
  #crs0 = sparse(Diagonal(crs[1:Nrp]))
  #crsN = sparse(Diagonal(crs[Nrp*Ns .+ (1:Nrp)]))
  #csr0 = sparse(Diagonal(csr[1   .+ Nrp*(0:Ns)]))
  #csrN = sparse(Diagonal(csr[Nrp .+ Nrp*(0:Ns)]))

  #
  # Surface mass matrices
  #
  H1 = Hs
  H1I = HsI

  H2 = Hs
  H2I = HsI

  H3 = Hr
  H3I = HrI

  H4 = Hr
  H4I = HrI


  # Surface quadtrature matrices
  H1 = H2 = Hs 
  H3 = H4 = Hr

  H = (Hs, Hs, Hr, Hr)
  HI = (HsI, HsI, HrI, HrI)

  J = spdiagm(0 => reshape(metrics.J, Nrp*Nsp))

  # Face lift operators (their transpose gives the face extraction operators)
  e = (convert(SparseMatrixCSC{Float64, Int64}, kron(Is, er0)),
       convert(SparseMatrixCSC{Float64, Int64}, kron(Is, erN)),
       convert(SparseMatrixCSC{Float64, Int64}, kron(es0, Ir)),
       convert(SparseMatrixCSC{Float64, Int64}, kron(esN, Ir)))
  
  # coefficent matrices on each of the four faces
  Crr1 = spdiagm(0 => crr[1, :])
  Crs1 = spdiagm(0 => crs[1, :])
  Csr1 = spdiagm(0 => crs[1, :])
  Css1 = spdiagm(0 => css[1, :])
  
  Crr2 = spdiagm(0 => crr[Nrp, :])
  Crs2 = spdiagm(0 => crs[Nrp, :])
  Csr2 = spdiagm(0 => crs[Nrp, :])
  Css2 = spdiagm(0 => css[Nrp, :])
  
  Css3 = spdiagm(0 => css[:, 1])
  Crs3 = spdiagm(0 => crs[:, 1])
  Csr3 = spdiagm(0 => crs[:, 1])
  Crr3 = spdiagm(0 => crr[:, 1])
  
  Css4 = spdiagm(0 => css[:, Nsp])
  Crs4 = spdiagm(0 => crs[:, Nsp])
  Csr4 = spdiagm(0 => crs[:, Nsp])
  Crr4 = spdiagm(0 => crr[:, Nsp])
 
  # penalty matrices:
  β = 1
  d = 1
  h1 = Hr[1, 1] #TODO: ???
  Z_1 = β * (d/h1) * Crr1 # TODO: are these correct?
  Z_2 = β * (d/h1) * Crr2
  
  # discrete boundary traction operators:
  T_1 = nCrr1 * (Iz ⊗ S0x) - Crs1 * (Dz ⊗ Ix)
  T_2 = Crr2 * (Iz ⊗ SNx) + Crs2 * (Dz ⊗ Ix)
  T_3 = -Css3 * (S0z ⊗ Ix) - Csr3 * (Iz ⊗ Dx)
  T_4 = -Css4 * (SNz ⊗ Ix) - Csr4 * (Iz ⊗ Dx)  # TODO: I'm sure these are wrong?   Do they not build in normal derivative? maybe i did. 
  T = (T_1, T_2, T_3, T_4)
  
        
  # Dirichlet on faces 1 and 2:
  SAT_1 = Hinv * (T_1 .- Z_1)' * e[1] * H1 * e[1]'   #TODO: fix these with A&D paper
  SAT_2 = Hinv * (T_2 .- Z_2)' * e[2] * H2 * e[2]'   #TODO: fix these with A&D paper

  # Neumann on faces 3 and 4:
  SAT_3 = -Hinv * e[3] * H3 * e[3]' * T_3        #TODO: fix these with A&D paper
  SAT_4 = -Hinv * e[4] * H4 * e[4]' * T_4        #TODO: fix these with A&D paper

  # boundary data operators for Dirichlet faces 1 and 2: 
  B1 = Hinv * (T_1 .- Z_1)' * e[1] * H1   #TODO: fix this with A&D paper
  B2 = Hinv * (T_2 .- Z_2)' * e[2] * H2   #TODO: fix this with A&D paper

  # boundary data operators for Neumann faces 3 and 4: 
  B3 = -Hinv * e[3] * H3  #TODO: fix this with A&D paper
  B4 = -Hinv * e[4] * H4  #TODO: fix this with A&D paper

  D2 = Drr + Dss + Drs + Dsr
  A = D2 + SAT_1 + SAT_2 + SAT_3 + SAT_4  # Creation of LHS matrix

  B = (B1, B2, B3, B4)  # These will multiply the boundary data when forming the RHS.

  H̃ = Hz ⊗ Hx            # Need this for convergence tests.


  return (A, B, H̃, T, e, J, coord = metrics.coord, 
  facecoord = metrics.facecoord, 
  sJ = metrics.sJ,
  nx = metrics.nx,
  ny = metrics.ny)


end



function bdry_vec_strip!(g, B, x, z, slip_data, remote_data, free_surface_data, Lx, Lz)

    
  g[:] .= 0

  # fault (Dirichlet):
  vf = slip_data
  g[:] += B[1] * vf

  # FACE 2 (Dirichlet):
  vf = remote_data
  g[:] += B[2] * vf

  # FACE 3 (Neumann):
  gN = free_surface_data
  vf = gN
  g[:] += B[3] * sJ[3] * vf  #TODO: prob error

  # FACE 4 (Neumann):
  gN = free_surface_data
  vf = gN
  g[:] += B[4] * sJ[4] * vf #TODO: prob error
  

end


function computetraction_stripped(T, u, e)
  e1 = e[1]
  tau = (-e1' * T[1]* u) ./ sJ #TODO: check for correctness
  return tau
end
  

function rateandstate(V, psi, σn, ϕ, η, a, V0)
  Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
  f = a .* asinh.(V .* Y)
  dfdV  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* Y
  
  g    = σn .* f    + η .* V - ϕ
  dgdV = σn .* dfdV + η
  (g, dgdV)
end
  
  function newtbndv(func, xL, xR, x; ftol = 1e-6, maxiter = 500, minchange=0,
                    atolx = 1e-4, rtolx = 1e-4)
    (fL, _) = func(xL)
    (fR, _) = func(xR)
    if fL .* fR > 0
      return (typeof(x)(NaN), typeof(x)(NaN), -maxiter)
    end
  
    (f, df) = func(x)
    dxlr = xR - xL
  
    for iter = 1:maxiter
      dx = -f / df
      x  = x + dx
  
      if x < xL || x > xR || abs(dx) / dxlr < minchange
        x = (xR + xL) / 2
        dx = (xR - xL) / 2
      end
  
      (f, df) = func(x)
  
      if f * fL > 0
        (fL, xL) = (f, x)
      else
        (fR, xR) = (f, x)
      end
      dxlr = xR - xL
  
      if abs(f) < ftol && abs(dx) < atolx + rtolx * (abs(dx) + abs(x))
        return (x, f, iter)
      end
    end
    return (x, f, -maxiter)
  end



export transfinite_blend, create_metrics, get_operators, computetraction
export rateandstate, newtbndv
