using SparseArrays
using LinearAlgebra
#using UnicodePlots

⊗(A,B) = kron(A, B)

const BC_DIRICHLET        = 1
const BC_NEUMANN          = 2
const BC_LOCKED_INTERFACE = 0
const BC_JUMP_INTERFACE   = 7

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


function create_metrics(Nr, Ns, exact_mu,
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

 

  μ = exact_mu(x, y)
  J = xr .* ys - xs .* yr
  # The Jacobian determinant of mapping from logical to physical.  (large if mapping from domain larger than 1. )
 
  @assert minimum(J) > 0

  JI = 1 ./ J

  rx =  ys ./ J
  sx = -yr ./ J
  ry = -xs ./ J
  sy =  xr ./ J

  # variable coefficient matrix components
  crr = J .* (rx .* μ .* rx .+ ry .* μ .* ry)
  crs = J .* (sx .* μ .* rx .+ sy .* μ .* ry)
  css = J .* (sx .* μ .* sx .+ sy .* μ .* sy)


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

function neighborZ(e, lf, lop)

  # takes in an element e, which is assumed to have an INTERFACE on its local face lf
  # returns the sJZ value on opposing face. 

  
      f = EToF[lf, e] # get global face number corresponding to this elements local face lf. CHECK THIS.

      (em, ep) = FToE[:, f]  # find the two elements that share global face f.
      @assert em == e
      (fe, fp) = FToLF[:, f]  # find corresponding local faces corresponding to global face f.
      @assert fe == lf

      return lop[ep].sJZ[fp]
end



function locoperator(p, Nr, Ns, μ, metrics=create_metrics(Nr,Ns,μ),  
            LFToB = (BC_DIRICHLET, BC_DIRICHLET,
                     BC_DIRICHLET, BC_DIRICHLET);
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

  # Identity matrices
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

  #{{{ Set up the rs and sr derivative matrices
  Drs = (Is ⊗ Dr) * sparse(1:length(crs), 1:length(crs), view(crs, :)) * (Ds ⊗ Ir)
  Dsr = (Ds ⊗ Ir) * sparse(1:length(csr), 1:length(csr), view(csr, :)) * (Is ⊗ Dr)
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


  # Surface mass matrices
  H1 = Hs
  H1I = HsI

  H2 = Hs
  H2I = HsI

  H3 = Hr
  H3I = HrI

  H4 = Hr
  H4I = HrI

  # Surface quadrature matrices
  H1 = H2 = Hs 
  H3 = H4 = Hr

  H = (Hs, Hs, Hr, Hr)
  HI = (HsI, HsI, HrI, HrI)

  sJ1 = spdiagm(0 => metrics.sJ[1])
  sJ2 = spdiagm(0 => metrics.sJ[2])
  sJ3 = spdiagm(0 => metrics.sJ[3])
  sJ4 = spdiagm(0 => metrics.sJ[4])

  sJ1_J1 = spdiagm(0 => metrics.sJ[1] ./ metrics.J[1, :])
  sJ2_J2 = spdiagm(0 => metrics.sJ[2] ./ metrics.J[Nrp, :])
  sJ3_J3 = spdiagm(0 => metrics.sJ[3] ./ metrics.J[:, 1])
  sJ4_J4 = spdiagm(0 => metrics.sJ[4] ./ metrics.J[:, Nsp])


  J = spdiagm(0 => reshape(metrics.J, Nrp*Nsp))

  # Face lift operators 
  eRS = (convert(SparseMatrixCSC{Float64, Int64}, kron(Is, er0)),
       convert(SparseMatrixCSC{Float64, Int64}, kron(Is, erN)),
       convert(SparseMatrixCSC{Float64, Int64}, kron(es0, Ir)),
       convert(SparseMatrixCSC{Float64, Int64}, kron(esN, Ir)))
  # face extraction ops
  eRST = (convert(SparseMatrixCSC{Float64, Int64}, kron(Is, er0')),
       convert(SparseMatrixCSC{Float64, Int64}, kron(Is, erN')),
       convert(SparseMatrixCSC{Float64, Int64}, kron(es0', Ir)),
       convert(SparseMatrixCSC{Float64, Int64}, kron(esN', Ir)))
  
 
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
 

  Hinv = HsI ⊗ HrI
  # penalty matrices:
  β = 10
  d = 2 # dimension

  # h1 and h2 in A&D are in logical space
  
  h1 = Hr[1, 1] 
  h2 = Hr[Nrp, Nrp] 
  h3 = Hs[1, 1] 
  h4 = Hs[Nsp, Nsp] 

  
  #  Z_1 = β * (d/h1) * (Is ⊗ (sJ1_J1)) * (Is ⊗ Crr1) # TODO: check; this is al messed up
  #  Z_2 = β * (d/h1) * (Is ⊗ (sJ2_J2)) * (Is ⊗ Crr2)
  #  Z_3 = β * (d/h2) * ((sJ3_J3) ⊗ Ir) * (Css3 ⊗ Ir) # TODO: check
  #  Z_4 = β * (d/h2) * ((sJ4_J4) ⊗ Ir) * (Css4 ⊗ Ir)

  # Z' penalty parameters (for Dirichlet BC) in logical space:
  sJZ_1 = β * (d ./ h1) * (Is ⊗ Crr1) 
  sJZ_2 = β * (d ./ h2) * (Is ⊗ Crr2)
  sJZ_3 = β * (d ./ h3) * (Css3 ⊗ Ir) 
  sJZ_4 = β * (d ./ h4) * (Css4 ⊗ Ir)
  sJZ = (sJZ_1, sJZ_2, sJZ_3, sJZ_4)

  # discrete boundary traction operators in logical space:
  T_1 = -Sr0 - ((Crs1*Ds) ⊗ Ir)
  T_2 = SrN + ((Crs2*Ds) ⊗ Ir)
  T_3 = -(Is ⊗ (Csr3*Dr)) - Ss0
  T_4 = (Is ⊗ (Csr4*Dr)) + SsN
  T = (T_1, T_2, T_3, T_4)  
  TT = (T_1', T_2', T_3', T_4')
    
  # Dirichlet faces:
  dSAT_1 = Hinv * (T_1 .- sJZ_1)' * eRS[1] * H1 * eRS[1]'   
  dSAT_2 = Hinv * (T_2 .- sJZ_2)' * eRS[2] * H2 * eRS[2]'  
  dSAT_3 = Hinv * (T_3 .- sJZ_3)' * eRS[3] * H3 * eRS[3]'   
  dSAT_4 = Hinv * (T_4 .- sJZ_4)' * eRS[4] * H4 * eRS[4]'  
  dSAT = (dSAT_1, dSAT_2, dSAT_3, dSAT_4)

  # Neumann faces:
  nSAT_1 = -Hinv * eRS[1] * H1 * eRS[1]' * T_1      
  nSAT_2 = -Hinv * eRS[2] * H2 * eRS[2]' * T_2 
  nSAT_3 = -Hinv * eRS[3] * H3 * eRS[3]' * T_3      
  nSAT_4 = -Hinv * eRS[4] * H4 * eRS[4]' * T_4  
  nSAT = (nSAT_1, nSAT_2, nSAT_3, nSAT_4)     


   # Interfaces (these are missing the face extraction and latter traction ops on purpose): 
  #  iSAT_1 = Hinv * (0.5*T_1 .- Z_1)' * eRS[1] * H1 * eRS[1]' - 0.5*Hinv * eRS[1] * H1 * eRS[1]' * T_1      
  #  iSAT_2 = Hinv * (0.5*T_2 .- Z_2)' * eRS[2] * H2 * eRS[2]' - 0.5*Hinv * eRS[2] * H2 * eRS[2]' * T_2 
  #  iSAT_3 = Hinv * (0.5*T_3 .- Z_3)' * eRS[3] * H3 * eRS[3]' - 0.5*Hinv * eRS[3] * H3 * eRS[3]' * T_3      
  #  iSAT_4 = Hinv * (0.5*T_4 .- Z_4)' * eRS[4] * H4 * eRS[4]' - 0.5*Hinv * eRS[4] * H4 * eRS[4]' * T_4  
  
  

  # Z' penalty parameters (for interface Dirichlet conditions) in logical space:
  IsJZ_1 = β * (d ./ (4*h1)) * (Is ⊗ Crr1) 
  IsJZ_2 = β * (d ./ (4*h2)) * (Is ⊗ Crr2)
  IsJZ_3 = β * (d ./ (4*h3)) * (Css3 ⊗ Ir) 
  IsJZ_4 = β * (d ./ (4*h4)) * (Css4 ⊗ Ir)
  IsJZ = (IsJZ_1, IsJZ_2, IsJZ_3, IsJZ_4)


  # iSAT_1 = Hinv * (0.5*T_1 .- sJZ_1)' * eRS[1] * H1    
  # iSAT_2 = Hinv * (0.5*T_2 .- sJZ_2)' * eRS[2] * H2
  # iSAT_3 = Hinv * (0.5*T_3 .- sJZ_3)' * eRS[3] * H3 
  # iSAT_4 = Hinv * (0.5*T_4 .- sJZ_4)' * eRS[4] * H4
  # iSAT = (iSAT_1, iSAT_2, iSAT_3, iSAT_4) 
 
  # these are for the jump in traction
  tSAT_1 = -0.5 * Hinv * eRS[1] * H1
  tSAT_2 = -0.5 * Hinv * eRS[2] * H2
  tSAT_3 = -0.5 * Hinv * eRS[3] * H3
  tSAT_4 = -0.5 * Hinv * eRS[4] * H4
  tSAT = (tSAT_1, tSAT_2, tSAT_3, tSAT_4) 
  
  
  ########  BELOW ARE FOR RHS VECTOR ########
  # boundary data operators for Dirichlet faces
  dB1 = Hinv * (T_1 .- sJZ_1)' * eRS[1] * H1   
  dB2 = Hinv * (T_2 .- sJZ_2)' * eRS[2] * H2  
  dB3 = Hinv * (T_3 .- sJZ_3)' * eRS[3] * H3
  dB4 = Hinv * (T_4 .- sJZ_4)' * eRS[4] * H4  
  dB = (dB1, dB2, dB3, dB4)

  # boundary data operators for Neumann faces
  nB1 = -Hinv * eRS[1] * H1
  nB2 = -Hinv * eRS[2] * H2 
  nB3 = -Hinv * eRS[3] * H3 
  nB4 = -Hinv * eRS[4] * H4 
  nB = (nB1, nB2, nB3, nB4)

  
   # boundary data operators for INTERFACES 
   # code set up to impose continuity of traction only, which is why these ops don't include those for enforcing non-zero traciton jumps.
   iB1 = Hinv * (0.5*T_1 .- sJZ_1)' * eRS[1] * H1  
   iB2 = Hinv * (0.5*T_2 .- sJZ_2)' * eRS[2] * H2  
   iB3 = Hinv * (0.5*T_3 .- sJZ_3)' * eRS[3] * H3
   iB4 = Hinv * (0.5*T_4 .- sJZ_4)' * eRS[4] * H4  
   iB = (iB1, iB2, iB3, iB4)

  # Build operator for element local
  A = Drr + Dss + Drs + Dsr

  # Build operators to connect interface
  # B = spzeros(size(A))

  # Modify operator to handle the boundary conditions
  bctype=(BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE)
  
  for lf = 1:4
    if LFToB[lf] == BC_NEUMANN
      A += nSAT[lf]
    elseif LFToB[lf] == BC_DIRICHLET
      A += dSAT[lf]
    elseif (LFToB[lf] == BC_LOCKED_INTERFACE) || (LFToB[lf] == BC_JUMP_INTERFACE) || (LFToB[lf] == RS_FAULT) || (LFToB[lf] == VP_FAULT)
      #A += iSAT[lf][1]*eRST[lf] + iSAT[lf][2]*tRS[lf]
    else  
      print("invalid bc")
     end
  end
  bctype=(LFToB[1], LFToB[2], LFToB[3], LFToB[4])


  H̃ = Hs ⊗ Hr            # Need this for convergence tests.


  return (A, dB, nB, iB, eRST, tSAT, H̃, T, TT, eRS, J, coord = metrics.coord, IsJZ, Hinv, H,
  facecoord = metrics.facecoord, 
  sJ = metrics.sJ,
  nx = metrics.nx,
  ny = metrics.ny,
  sJZ)


end

#{{{ volsourcearray()
function locsourcearray!(ge, source, coord, J, volargs = ())

  (xloc, yloc) = coord
  Jloc = J
  ge[:] += Jloc * source(xloc[:], yloc[:], volargs...)

end

#}}}

function loc_bdry_vec!(ge, lop, neighborZ, LFToB, EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, in_tractionjump, bcargs = ())

  e = bcargs[1]
  dB = lop.dB 
  nB = lop.nB 
  iB = lop.iB

  (xf, yf) = lop.facecoord
  sJ = lop.sJ 
  nx = lop.nx
  ny = lop.ny
  ge[:] .= 0
  for lf = 1:4
    if LFToB[lf] == BC_DIRICHLET
      # fault (Dirichlet):
      gD = bc_Dirichlet(lf, xf[lf], yf[lf], bcargs...) # this should return B[1]*g
      vf = dB[lf]*gD
  
 
    elseif LFToB[lf] == BC_NEUMANN
      gN = bc_Neumann(lf, xf[lf], yf[lf], nx[lf], ny[lf], bcargs...) #this should return B[3]*g3
      vf = nB[lf] * (sJ[lf] .* gN)
      

    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here

    elseif LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
      
      f = EToF[lf, e]  # get global face number

      (em, ep) = FToE[:, f]  # find the two elements that share global face f; ep is on plus side
      (fm, fp) = FToLF[:, f]

      if em == e
        #eo = ep 
        nf = fp
      else
        #eo = em 
        nf = fm
      end

      sJZm = 1*(lop.IsJZ[lf])' + (neighborZ[nf])'
  
      # Jump in displacement
      gJ = in_jump(lf, xf[lf], yf[lf], bcargs...) 
      A = -lop.Hinv * sJZm * lop.eRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf] * lop.eRS[lf] * lop.H[lf]
      
      # Jump in traction
      gT = in_tractionjump(lf, xf[lf], yf[lf], bcargs...) 
    
      B = -0.5 * lop.Hinv * lop.eRS[lf] * lop.H[lf]
    
      vf = A*gJ + B*(sJ[lf] .* gT)

    else
      error("invalid bc")
    end
    ge[:] += vf
  end

end

function loc_bdry_vec_v2!(ge, lop, neighborZ, LFToB, EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, bcargs = ())

  e = bcargs[1]
  dB = lop.dB 
  nB = lop.nB 
  iB = lop.iB

  (xf, yf) = lop.facecoord
  sJ = lop.sJ 
  nx = lop.nx
  ny = lop.ny
  ge[:] .= 0
  for lf = 1:4
    if LFToB[lf] == BC_DIRICHLET
      # fault (Dirichlet):
      gD = bc_Dirichlet(lf, xf[lf], yf[lf], bcargs...) # this should return B[1]*g
      vf = dB[lf]*gD
  
 
    elseif LFToB[lf] == BC_NEUMANN
      gN = bc_Neumann(lf, xf[lf], yf[lf], nx[lf], ny[lf], bcargs...) #this should return B[3]*g3
      vf = nB[lf] * (sJ[lf] .* gN)
      

    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here

    elseif LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
      
      f = EToF[lf, e]  # get global face number

      (em, ep) = FToE[:, f]  # find the two elements that share global face f; ep is on plus side
      (fm, fp) = FToLF[:, f]

      if em == e
        #eo = ep 
        nf = fp
      else
        #eo = em 
        nf = fm
      end

      sJZm = 1*(lop.IsJZ[lf])' + (neighborZ[nf])'
  
      # Jump in displacement
      gJ = in_jump(lf, xf[lf], yf[lf], bcargs...) 
      A = -lop.Hinv * sJZm * lop.eRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf] * lop.eRS[lf] * lop.H[lf]
      
      B = -0.5 * lop.Hinv * lop.eRS[lf] * lop.H[lf]
    
      vf = A*gJ

    else
      error("invalid bc")
    end
    ge[:] += vf
  end

end

# function computetraction(T, u, e, sJ)
#   e1 = e[1]
#   tau = (-e1' * T[1]* u) ./ sJ[1] #TODO: check for correctness and generality
#   return tau
# end
  
# function computetraction(lop, lf, u, λ, δ)
#   HfI_FT = lop.HfI_FT[lf]
#   τf = lop.τ[lf]
#   sJ = lop.sJ[lf]

#   return (HfI_FT * u + τf * (λ .- δ / 2)) ./ sJ
# end

function computetraction(lop, lf, u)
  T = lop.T[lf]
  sJ = lop.sJ[lf]
  eRST = lop.eRST[lf]
  return (eRST * T * u) ./ sJ 
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



  function plot_connectivity(verts, EToV)
    Lx = extrema(verts[1,:])
    Lx = (floor(Int, Lx[1]), ceil(Int, Lx[2]))
    Ly = extrema(verts[2,:])
    Ly = (floor(Int, Ly[1]), ceil(Int, Ly[2]))
    width = Lx[2] - Lx[1]
    height = Ly[2] - Ly[1]
    plt = Plot(BrailleCanvas(80, ceil(Int, 40 * height / width),
                             origin_x = Lx[1], origin_y = Ly[1],
                             width = width, height = height))
  
  
    #annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
    #annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
    #annotate!(plt, :bl, string(Lx[1]), color = :light_black)
    #annotate!(plt, :br, string(Lx[2]), color = :light_black)
    for e = 1:size(EToV, 2)
      (v1, v2, v3, v4) = EToV[1:4, e]
      x = verts[1, [v1 v2 v4 v3 v1]][:]
      y = verts[2, [v1 v2 v4 v3 v1]][:]
      lineplot!(plt, x, y)
    end
    #title!(plt, "connectivity")
    display(plt)
  end
  
  function plot_blocks(lop)
    Lx = (floatmax(), -floatmax())
    Ly = (floatmax(), -floatmax())
    for e = 1:length(lop)
      (x, y) = lop[e].coord
      Lxe = extrema(x)
      Lye = extrema(y)
      Lx = (min(Lx[1], Lxe[1]), max(Lx[2], Lxe[2]))
      Ly = (min(Ly[1], Lye[1]), max(Ly[2], Lye[2]))
    end
  
    Lx = (floor(Int, Lx[1]), ceil(Int, Lx[2]))
    Ly = (floor(Int, Ly[1]), ceil(Int, Ly[2]))
  
    width = Lx[2] - Lx[1]
    height = Ly[2] - Ly[1]
    plt = Plot(BrailleCanvas(80, ceil(Int, 40 * height / width),
                             origin_x = Lx[1], origin_y = Ly[1],
                             width = width, height = height))
  
  
    #annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
    #annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
    #annotate!(plt, :bl, string(Lx[1]), color = :light_black)
    #annotate!(plt, :br, string(Lx[2]), color = :light_black)
  
    for e = 1:length(lop)
      (xf, yf) = lop[e].facecoord
      bctype = lop[e].bctype
      for lf = 1:length(xf)
        if bctype[lf] == BC_LOCKED_INTERFACE
          lineplot!(plt, xf[lf], yf[lf], color=:blue)
        elseif bctype[lf] == BC_DIRICHLET
          lineplot!(plt, xf[lf], yf[lf], color=:green)
        elseif bctype[lf] == BC_DIRICHLET
          lineplot!(plt, xf[lf], yf[lf], color=:yellow)
        else
          lineplot!(plt, xf[lf], yf[lf], color=:red)
        end
      end
    end
    #title!(plt, "mesh")
    display(plt)
  end



#{{{ global_operator: Build the global operator
function global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
  
  nelems = length(lop)
  nfaces = length(FToB)

  # Get a unique array indexes for the face to jumps map
  #FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT, BC_JUMP_INTERFACE), Nr, Ns) # might want to send this in?
  
  # Compute the total number of volume and jump (δ) points
  VNp = vstarts[nelems+1]-1
  #δNp = FToδstarts[nfaces+1]-1

  M = spzeros(VNp, VNp)
  for e = 1:nelems
    M[vstarts[e]:vstarts[e+1]-1, vstarts[e]:vstarts[e+1]-1] = lop[e].A  # put local A along block diagonal
  end

  for f = 1:nfaces  # Loop over the global faces to put in off-diagonal interface ops.
    
    if (FToB[f] == BC_JUMP_INTERFACE) || (FToB[f] == BC_LOCKED_INTERFACE) || (FToB[f] == RS_FAULT) || (FToB[f] == VP_FAULT)  # check if face is an interface
   
      (em, ep) = FToE[:, f]  # find the two elements that share global face f.
      (fm, fp) = FToLF[:, f]  # find corresponding local faces corresponding to global face f.
     
      #M[vstarts[em]:vstarts[em+1]-1, vstarts[em]:vstarts[em+1]-1] += lop[em].iSAT[fm]*lop[em].eRST[fm] + lop[em].tSAT[fm]*lop[em].eRST[fm]*lop[em].T[fm]
      #M[vstarts[ep]:vstarts[ep+1]-1, vstarts[ep]:vstarts[ep+1]-1] += lop[ep].iSAT[fp]*lop[ep].eRST[fp] + lop[ep].tSAT[fp]*lop[ep].eRST[fp]*lop[ep].T[fp]
      sJZm = (lop[em].IsJZ[fm])' + (lop[ep].IsJZ[fp])'

      c = 1
      M[vstarts[em]:vstarts[em+1]-1, vstarts[em]:vstarts[em+1]-1] += (-lop[em].Hinv * sJZm * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] + 
                                                                      + 0.5 * lop[em].Hinv * lop[em].TT[fm] * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] + 
                                                                      -  c * 0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] * lop[em].T[fm])
      M[vstarts[ep]:vstarts[ep+1]-1, vstarts[ep]:vstarts[ep+1]-1] += (-lop[ep].Hinv * sJZm * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] + 
                                                                      + 0.5 * lop[ep].Hinv * lop[ep].TT[fp] * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] + 
                                                                      -  c * 0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] * lop[ep].T[fp])

      # If face orientation on plus side does not match minus side, then flip:
      Np = (fm <= 2 ? Ns[em]+1 : Nr[em]+1) # Set Np to be either Nr or Ns, depending on face
      if EToO[fp, ep]                     # orientation matches
        R = sparse(I, Np, Np)         
      else  
        R = rotr90(sparse(I, Np, Np)) # orientation doesn't match
      end

      # M[vstarts[em]:vstarts[em+1]-1, vstarts[ep]:vstarts[ep+1]-1] += -lop[em].iSAT[fm]*R*lop[ep].eRST[fp] + lop[em].tSAT[fm]*R*lop[ep].eRST[fp]*lop[ep].T[fp]
      # M[vstarts[ep]:vstarts[ep+1]-1, vstarts[em]:vstarts[em+1]-1] += -lop[ep].iSAT[fp]*R*lop[em].eRST[fm] + lop[ep].tSAT[fp]*R*lop[em].eRST[fm]*lop[em].T[fm]
   
      M[vstarts[em]:vstarts[em+1]-1, vstarts[ep]:vstarts[ep+1]-1] += (lop[em].Hinv * sJZm * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] + 
                                                                    - 0.5 * lop[em].Hinv * lop[em].TT[fm] * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] + 
                                                                    - c * 0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] * lop[ep].T[fp])

      M[vstarts[ep]:vstarts[ep+1]-1, vstarts[em]:vstarts[em+1]-1] += (lop[ep].Hinv * sJZm * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] + 
                                                                    - 0.5 * lop[ep].Hinv * lop[ep].TT[fp] * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] + 
                                                                    - c * 0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] * lop[em].T[fm])
    end
  end
  return M
end
#}}}


function bcstarts(FToB, FToE, FToLF, bctype, Nr, Ns)
    # compute an array of indices to loop over vertices on element faces
    nfaces = length(FToB)   # get total number of faces
    bcstarts = Array{Int64, 1}(undef, nfaces + 1) # initialize array
    bcstarts[1] = 1                               # array starts at 1
    for f = 1:nfaces                              # loop over the global face numbers
      if FToB[f] ∈ bctype                         # determine if face is of type bctype, e.g. BC_JUMP_INTERFACE
        e  = FToE[1,f]                            # if so, set e = first element shared by face (plus side)
        lf = FToLF[1,f]                           # get the local face number for first element
        bcstarts[f+1] = bcstarts[f] + (lf ∈ (1,2) ? Ns[e] : Nr[e]) + 1 # fill in next indices to be previous index plus either Ns (if local face is 1 or 2) or Nr (local face 3 or 4)
      else
        bcstarts[f+1] = bcstarts[f]               # if face is not of bctype, just insert a 1. This will imply that we skipped over a face that wasn't of bctype.
      end
    end
    bcstarts                                      # return array of indices
end





#{{{ connectivityarrays
function connectivityarrays(EToV, EToF)
  # number of elements
  nelems = size(EToV, 2)
  nfaces = maximum(maximum(EToF))

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side

  FToE  = zeros(Int64, 2, nfaces)
  FToLF = zeros(Int64, 2, nfaces)
  EToO  = Array{Bool,2}(undef, 4, nelems)
  EToS  = zeros(Int64, 4, nelems)

  # Local Face to Local Vertex map
  LFToLV = flatten_tuples(((1,3), (2, 4), (1,2), (3,4)))
  for e = 1:nelems
    for lf = 1:4
      gf = EToF[lf, e]
      if FToE[1, gf] == 0
        @assert FToLF[1, gf] == 0
        FToE[1, gf] = e
        FToLF[1, gf] = lf
        EToO[lf, e] = true
        EToS[lf, e] = 1
      else
        @assert FToE[2, gf] == 0
        @assert FToLF[2, gf] == 0
        FToE[2, gf] = e
        FToLF[2, gf] = lf
        EToS[lf, e] = 2

        ne = FToE[1, gf]
        nf = FToLF[1, gf]

        nv = EToV[LFToLV[:,nf], ne]
        lv = EToV[LFToLV[:,lf], e]
        
        if nv == lv
          EToO[lf, e] = true
        elseif nv[end:-1:1] == lv
          EToO[lf, e] = false
        else
          error("problem with connectivity")
        end
      end
    end
  end
  (FToE, FToLF, EToO, EToS)
end
#}}}


# flatten tuples to arrays
flatten_tuples(x) = reshape(collect(Iterators.flatten(x)), length(x[1]),
                            length(x))

export create_metrics, get_operators, computetraction,   plot_blocks, plot_connectivity
export bdry_vec!
export rateandstate, newtbndv, bcstarts, connectivityarrays, flatten_tuples
export LocalGlobalOperators, SBPLocalOperator1, global_operator
