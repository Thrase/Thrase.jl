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



function create_metrics(Nr, Ns, exact_mu, exact_lambda,
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
  λ = exact_lambda(x, y)
  
  J = xr .* ys - xs .* yr
  # The Jacobian determinant of mapping from logical to physical.  (large if mapping from domain larger than 1. )
 
  @assert minimum(J) > 0

  JI = 1 ./ J

  rx =  ys ./ J
  sx = -yr ./ J
  ry = -xs ./ J
  sy =  xr ./ J

  # variable coefficient matrix components
  c11 = J .* (rx .* (λ + 2μ) .* rx .+ ry .* μ .* ry)
  c12 = J .* (rx .* (λ + 2μ) .* sx .+ ry .* μ .* sy)
  c13 = J .* (rx .* λ .*  ry + rx .* μ .* ry);
  c14 = J .* (rx .* λ .* sy .+ ry .* μ .* sx);

  c21 = J .* (rx .* (λ + 2μ) .* sx .+ ry .* μ .* sy);
  c22 = J .* (sx .* (λ + 2μ) .* sx .+ sy .* μ .* sy);
  c23 = J .* (ry .* λ .* sx .+ rx .* μ .* sy);
  c24 = J .* (sx .* λ .* sy .+ sx .* μ .* sy);

  c31 = J .* (rx .* μ .* ry .+ rx .* λ .* ry);
  c32 = J .* (rx .* μ .* sy .+ ry .* λ .* sx);
  c33 = J .* (rx .* μ .* rx .+ ry .* (λ + 2μ) .* ry);
  c34 = J .* (rx .* μ .* sx .+ ry .* (λ + 2μ) .* sy);

  c41 = J .* (ry .* μ .* sx .+ rx .* λ .* sy);
  c42 = J .* (sx .* μ .* sy .+ sx .* λ .* sy);
  c43 = J .* (rx .* μ .* sx .+ ry .* (λ + 2μ) .* sy);
  c44 = J .* (sx .* μ .* sx .+ sy .* (λ + 2μ) .* sy);


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
  c11 = c11, c12 = c12, c13 = c13, c14 = c14, c21 = c21, c22 = c22, c23 = c23, c24 = c24, 
  c31 = c31, c32 = c32, c33 = c33, c34 = c34, c41 = c41, c42 = c42, c43 = c43, c44 = c44,
  J=J,
  JI = JI,
  sJ = (sJ1, sJ2, sJ3, sJ4),
  nx = (nx1, nx2, nx3, nx4),
  ny = (ny1, ny2, ny3, ny4),
  rx = rx, ry = ry, sx = sx, sy = sy)
end


function get_pure_rr(p, Nr, Nrp, Nsp, crr)

  Np = Nrp * Nsp

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

  return Drr, Sr0, SrN
end

function get_pure_ss(p, Ns, Nsp, Nrp, css)

  Np = Nrp * Nsp

  # Set up the ss derivative matrix
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

  return Dss, Ss0, SsN
end


function locoperator(p, Nr, Ns, μ, λ, metrics=create_metrics(Nr,Ns,μ,λ),  
            LFToB = (BC_DIRICHLET, BC_DIRICHLET,
                     BC_DIRICHLET, BC_DIRICHLET);
                     c11 = metrics.c11, c12 = metrics.c12,c13 = metrics.c13, c14 = metrics.c14,
                     c21 = metrics.c21, c22 = metrics.c22,c23 = metrics.c23, c24 = metrics.c24,
                     c31 = metrics.c31, c32 = metrics.c32,c33 = metrics.c33, c34 = metrics.c34,
                     c41 = metrics.c41, c42 = metrics.c42,c43 = metrics.c43, c44 = metrics.c44)
  
  
  # {{{ read in parameters
  J = metrics.J


  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp
  # }}}

  # {{{ Create 1D operators in logical space #TODO: higher order accuracy for computing tractions. 
  (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
  (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))
  # Change boundary stencils:
  # (_, S0r, SNr, _, _, _) = diagonal_sbp_D2(p, Nr; xc = (-1,1))
  # (_, S0s, SNs, _, _, _) = diagonal_sbp_D2(p, Ns; xc = (-1,1))
  # Dr[1,:] = S0r[1,:]
  # Dr[end,:] = SNr[end,:]
  # Ds[1,:] = S0s[1,:]
  # Ds[end,:] = SNs[end,:]
  # }}}

  # {{{ 2D Identity matrices
  Ir = sparse(I, Nrp, Nrp)
  Is = sparse(I, Nsp, Nsp)
  # }}}

  # {{{ Set up the operators in logical space

  Drr_11, Sr0_11, SrN_11 = get_pure_rr(p, Nr, Nrp, Nsp, c11)
  Dss_22, Ss0_22, SsN_22 = get_pure_ss(p, Ns, Nrp, Nsp, c22)
  
  Drr_13, Sr0_13, SrN_13 = get_pure_rr(p, Nr, Nrp, Nsp, c13)
  Dss_24, Ss0_24, SsN_24 = get_pure_ss(p, Ns, Nrp, Nsp, c24)
  
  Drr_31, Sr0_31, SrN_31 = get_pure_rr(p, Nr, Nrp, Nsp, c31)
  Dss_42, Ss0_42, SsN_42 = get_pure_ss(p, Ns, Nrp, Nsp, c42)

  Drr_33, Sr0_33, SrN_33 = get_pure_rr(p, Nr, Nrp, Nsp, c33)
  Dss_44, Ss0_44, SsN_44 = get_pure_ss(p, Ns, Nrp, Nsp, c44)


  # {{{ Boundary point matrices and surface mass matrices 

  er0 = sparse([1  ], [1], [1], Nrp, 1)
  erN = sparse([Nrp], [1], [1], Nrp, 1)
  es0 = sparse([1  ], [1], [1], Nsp, 1)
  esN = sparse([Nsp], [1], [1], Nsp, 1)


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
  
 
  # Surface mass matrices
  H1 = Hs
  H2 = Hs
  H3 = Hr
  H4 = Hr

  # Surface quadrature matrices
  H1 = H2 = Hs 
  H3 = H4 = Hr

  H = (Hs, Hs, Hr, Hr)
  Hinv = HsI ⊗ HrI
  
  # h1 and h2 in A&D are in logical space
  h1 = Hr[1, 1] 
  h2 = Hr[Nrp, Nrp] 
  h3 = Hs[1, 1] 
  h4 = Hs[Nsp, Nsp] 


  H̃ = Hs ⊗ Hr            # Need this for convergence tests.

  # }}}

  # {{{ Jacobian
  J = spdiagm(0 => reshape(metrics.J, Nrp*Nsp))
  # }}}



  #{{{ Coefficent matrices on each of the four faces
  C11_1 = spdiagm(0 => c11[1, :])
  C12_1 = spdiagm(0 => c12[1, :])
  C13_1 = spdiagm(0 => c13[1, :])
  C14_1 = spdiagm(0 => c14[1, :])
  C21_1 = spdiagm(0 => c21[1, :])
  C22_1 = spdiagm(0 => c22[1, :])
  C23_1 = spdiagm(0 => c23[1, :])
  C24_1 = spdiagm(0 => c24[1, :])
  C31_1 = spdiagm(0 => c31[1, :])
  C32_1 = spdiagm(0 => c32[1, :])
  C33_1 = spdiagm(0 => c33[1, :])
  C34_1 = spdiagm(0 => c34[1, :])
  C41_1 = spdiagm(0 => c41[1, :])
  C42_1 = spdiagm(0 => c42[1, :])
  C43_1 = spdiagm(0 => c43[1, :])
  C44_1 = spdiagm(0 => c44[1, :])

  C11_2 = spdiagm(0 => c11[Nrp, :])
  C12_2 = spdiagm(0 => c12[Nrp, :])
  C13_2 = spdiagm(0 => c13[Nrp, :])
  C14_2 = spdiagm(0 => c14[Nrp, :])
  C21_2 = spdiagm(0 => c21[Nrp, :])
  C22_2 = spdiagm(0 => c22[Nrp, :])
  C23_2 = spdiagm(0 => c23[Nrp, :])
  C24_2 = spdiagm(0 => c24[Nrp, :])
  C31_2 = spdiagm(0 => c31[Nrp, :])
  C32_2 = spdiagm(0 => c32[Nrp, :])
  C33_2 = spdiagm(0 => c33[Nrp, :])
  C34_2 = spdiagm(0 => c34[Nrp, :])
  C41_2 = spdiagm(0 => c41[Nrp, :])
  C42_2 = spdiagm(0 => c42[Nrp, :])
  C43_2 = spdiagm(0 => c43[Nrp, :])
  C44_2 = spdiagm(0 => c44[Nrp, :])

  
  
  C11_3 = spdiagm(0 => c11[:, 1])
  C12_3 = spdiagm(0 => c12[:, 1])
  C13_3 = spdiagm(0 => c13[:, 1])
  C14_3 = spdiagm(0 => c14[:, 1])
  C21_3 = spdiagm(0 => c21[:, 1])
  C22_3 = spdiagm(0 => c22[:, 1])
  C23_3 = spdiagm(0 => c23[:, 1])
  C24_3 = spdiagm(0 => c24[:, 1])
  C31_3 = spdiagm(0 => c31[:, 1])
  C32_3 = spdiagm(0 => c32[:, 1])
  C33_3 = spdiagm(0 => c33[:, 1])
  C34_3 = spdiagm(0 => c34[:, 1])
  C41_3 = spdiagm(0 => c41[:, 1])
  C42_3 = spdiagm(0 => c42[:, 1])
  C43_3 = spdiagm(0 => c43[:, 1])
  C44_3 = spdiagm(0 => c44[:, 1])
  
  C11_4 = spdiagm(0 => c11[:, Nsp])
  C12_4 = spdiagm(0 => c12[:, Nsp])
  C13_4 = spdiagm(0 => c13[:, Nsp])
  C14_4 = spdiagm(0 => c14[:, Nsp])
  C21_4 = spdiagm(0 => c21[:, Nsp])
  C22_4 = spdiagm(0 => c22[:, Nsp])
  C23_4 = spdiagm(0 => c23[:, Nsp])
  C24_4 = spdiagm(0 => c24[:, Nsp])
  C31_4 = spdiagm(0 => c31[:, Nsp])
  C32_4 = spdiagm(0 => c32[:, Nsp])
  C33_4 = spdiagm(0 => c33[:, Nsp])
  C34_4 = spdiagm(0 => c34[:, Nsp])
  C41_4 = spdiagm(0 => c41[:, Nsp])
  C42_4 = spdiagm(0 => c42[:, Nsp])
  C43_4 = spdiagm(0 => c43[:, Nsp])
  C44_4 = spdiagm(0 => c44[:, Nsp])

  

  # }}}
 
  # {{{ discrete boundary traction operators in logical space on each of the 4 faces:
  T11_1 = -Sr0_11 - ((C12_1*Ds) ⊗ Ir)
  T12_1 = -Sr0_13 - ((C14_1*Ds) ⊗ Ir)
  T21_1 = -Sr0_31 - ((C32_1*Ds) ⊗ Ir)
  T22_1 = -Sr0_33 - ((C34_1*Ds) ⊗ Ir)

  T11_2 = SrN_11 + ((C12_2*Ds) ⊗ Ir)
  T12_2 = SrN_13 + ((C14_2*Ds) ⊗ Ir)
  T21_2 = SrN_31 + ((C32_2*Ds) ⊗ Ir)
  T22_2 = SrN_33 + ((C34_2*Ds) ⊗ Ir)

  T11_3 = -(Is ⊗ (C21_3*Dr)) - Ss0_22
  T12_3 = -(Is ⊗ (C32_3*Dr)) - Ss0_24
  T21_3 = -(Is ⊗ (C41_3*Dr)) - Ss0_42
  T22_3 = -(Is ⊗ (C43_3*Dr)) - Ss0_44

  T11_4 = (Is ⊗ (C21_4*Dr)) + SsN_22
  T12_4 = (Is ⊗ (C32_4*Dr)) + SsN_24
  T21_4 = (Is ⊗ (C41_4*Dr)) + SsN_42
  T22_4 = (Is ⊗ (C43_4*Dr)) + SsN_44
 

  T_1 = [[T11_1] [T12_1]; [T21_1] [T22_1]]
  T_2 = [[T11_2] [T12_2]; [T21_2] [T22_2]]
  T_3 = [[T11_3] [T12_3]; [T21_3] [T22_3]]
  T_4 = [[T11_4] [T12_4]; [T21_4] [T22_4]]

  
  T = (T_1, T_2, T_3, T_4)   


  TT_1 = [[T11_1'] [1 * T21_1']; [1 * T12_1'] [T22_1']]
  TT_2 = [[T11_2'] [1 * T21_2']; [1 * T12_2'] [T22_2']]
  TT_3 = [[T11_3'] [1 * T21_3']; [1 * T12_3'] [T22_3']]
  TT_4 = [[T11_4'] [1 * T21_4']; [1 * T12_4'] [T22_4']]

  TT = (1 .* TT_1, 1 .* TT_2, 1 .* TT_3, 1 .* TT_4)   
  # }}}
    


  # {{{ Dirichlet BC penalties

  β = 1
  d = 2 # dimension
  
  sJ1 = spdiagm(0 => metrics.sJ[1])
  sJ2 = spdiagm(0 => metrics.sJ[2])
  sJ3 = spdiagm(0 => metrics.sJ[3])
  sJ4 = spdiagm(0 => metrics.sJ[4])


 
  # Z' penalty parameters (for Dirichlet BC) in logical space: 
  sJZ11_1 = β * (d ./ h1) * ((sJ1 * C11_1) ⊗ Ir) # TODO check sJ in these
  sJZ12_1 = β * (d ./ h1) * ((sJ1 * C13_1) ⊗ Ir)
  sJZ21_1 = β * (d ./ h1) * ((sJ1 * C31_1) ⊗ Ir)
  sJZ22_1 = β * (d ./ h1) * ((sJ1 * C33_1) ⊗ Ir) 

  sJZ11_2 = β * (d ./ h2) * ((sJ2 * C11_2) ⊗ Ir)
  sJZ12_2 = β * (d ./ h2) * ((sJ2 *C13_2) ⊗ Ir)
  sJZ21_2 = β * (d ./ h2) * ((sJ2 *C31_2) ⊗ Ir)
  sJZ22_2 = β * (d ./ h2) * ((sJ2 *C33_2) ⊗ Ir) 

  sJZ11_3 = β * (d ./ h3) * (Is ⊗ (sJ3 * C22_3)) 
  sJZ12_3 = β * (d ./ h3) * (Is ⊗ (sJ3 * C24_3))  
  sJZ21_3 = β * (d ./ h3) * (Is ⊗ (sJ3 * C42_3)) 
  sJZ22_3 = β * (d ./ h3) * (Is ⊗ (sJ3 * C44_3)) 

  sJZ11_4 = β * (d ./ h4) * (Is ⊗ (sJ4 * C22_4)) 
  sJZ12_4 = β * (d ./ h4) * (Is ⊗ (sJ4 * C24_4))  
  sJZ21_4 = β * (d ./ h4) * (Is ⊗ (sJ4 * C42_4)) 
  sJZ22_4 = β * (d ./ h4) * (Is ⊗ (sJ4 * C44_4)) 
  


  sJZ_1 = [[sJZ11_1] [sJZ12_1]; [sJZ21_1] [sJZ22_1]]
  sJZ_2 = [[sJZ11_2] [sJZ12_2]; [sJZ21_2] [sJZ22_2]]
  sJZ_3 = [[sJZ11_3] [sJZ12_3]; [sJZ21_3] [sJZ22_3]]
  sJZ_4 = [[sJZ11_4] [sJZ12_4]; [sJZ21_4] [sJZ22_4]]


  sJZ = (sJZ_1, sJZ_2, sJZ_3, sJZ_4)   


  # Dirichlet faces: 

  dSAT_1_11 = Hinv * (T11_1 .- sJZ11_1)' * eRS[1] * H1 * eRS[1]'   #multiplies first variable and added to first equation
  dSAT_1_12 = Hinv * (T12_1 .- sJZ12_1)' * eRS[1] * H1 * eRS[1]'     # mulitplies first variable and added to second
  dSAT_1_21 = Hinv * (T21_1 .- sJZ21_1)' * eRS[1] * H1 * eRS[1]' # multiplies second variable and added to first
  dSAT_1_22 = Hinv * (T22_1 .- sJZ22_1)' * eRS[1] * H1 * eRS[1]'  # multiplies second variable and added to second equaiton

  dSAT_2_11 = Hinv * (T11_2 .- sJZ11_2)' * eRS[2] * H2 * eRS[2]'   
  dSAT_2_12 = Hinv * (T12_2 .- sJZ12_2)' * eRS[2] * H2 * eRS[2]'
  dSAT_2_21 = Hinv * (T21_2 .- sJZ21_2)' * eRS[2] * H2 * eRS[2]'
  dSAT_2_22 = Hinv * (T22_2 .- sJZ22_2)' * eRS[2] * H2 * eRS[2]'

  dSAT_3_11 = Hinv * (T11_3 .- sJZ11_3)' * eRS[3] * H3 * eRS[3]'   
  dSAT_3_12 = Hinv * (T12_3 .- sJZ12_3)' * eRS[3] * H3 * eRS[3]'
  dSAT_3_21 = Hinv * (T21_3 .- sJZ21_3)' * eRS[3] * H3 * eRS[3]'
  dSAT_3_22 = Hinv * (T22_3 .- sJZ22_3)' * eRS[3] * H3 * eRS[3]'

  dSAT_4_11 = Hinv * (T11_4 .- sJZ11_4)' * eRS[4] * H4 * eRS[4]'   
  dSAT_4_12 = Hinv * (T12_4 .- sJZ12_4)' * eRS[4] * H4 * eRS[4]'
  dSAT_4_21 = Hinv * (T21_4 .- sJZ21_4)' * eRS[4] * H4 * eRS[4]'
  dSAT_4_22 = Hinv * (T22_4 .- sJZ22_4)' * eRS[4] * H4 * eRS[4]'
  
  dSAT_1 = [[dSAT_1_11] [dSAT_1_12]; [dSAT_1_21] [dSAT_1_22]]
  dSAT_2 = [[dSAT_2_11] [dSAT_2_12]; [dSAT_2_21] [dSAT_2_22]]
  dSAT_3 = [[dSAT_3_11] [dSAT_3_12]; [dSAT_3_21] [dSAT_3_22]]
  dSAT_4 = [[dSAT_4_11] [dSAT_4_12]; [dSAT_4_21] [dSAT_4_22]]

  dSAT = (dSAT_1, dSAT_2, dSAT_3, dSAT_4)
  # }}}

  # {{{ Traction faces: 
  nSAT_1_11 = -Hinv * eRS[1] * H1 * eRS[1]' * T11_1   # see above. 
  nSAT_1_12 = -Hinv * eRS[1] * H1 * eRS[1]' * T12_1    
  nSAT_1_21 = -Hinv * eRS[1] * H1 * eRS[1]' * T21_1   
  nSAT_1_22 = -Hinv * eRS[1] * H1 * eRS[1]' * T22_1      
  
  nSAT_2_11 = -Hinv * eRS[2] * H2 * eRS[2]' * T11_2  
  nSAT_2_12 = -Hinv * eRS[2] * H2 * eRS[2]' * T12_2    
  nSAT_2_21 = -Hinv * eRS[2] * H2 * eRS[2]' * T21_2
  nSAT_2_22 = -Hinv * eRS[2] * H2 * eRS[2]' * T22_2   

  nSAT_3_11 = -Hinv * eRS[3] * H3 * eRS[3]' * T11_3 
  nSAT_3_12 = -Hinv * eRS[3] * H3 * eRS[3]' * T12_3    
  nSAT_3_21 = -Hinv * eRS[3] * H3 * eRS[3]' * T21_3
  nSAT_3_22 = -Hinv * eRS[3] * H3 * eRS[3]' * T22_3   

  nSAT_4_11 = -Hinv * eRS[4] * H4 * eRS[4]' * T11_4 
  nSAT_4_12 = -Hinv * eRS[4] * H4 * eRS[4]' * T12_4    
  nSAT_4_21 = -Hinv * eRS[4] * H4 * eRS[4]' * T21_4
  nSAT_4_22 = -Hinv * eRS[4] * H4 * eRS[4]' * T22_4 

  nSAT_1 = [[nSAT_1_11] [nSAT_1_12]; [nSAT_1_21] [nSAT_1_22]]
  nSAT_2 = [[nSAT_2_11] [nSAT_2_12]; [nSAT_2_21] [nSAT_2_22]]
  nSAT_3 = [[nSAT_3_11] [nSAT_3_12]; [nSAT_3_21] [nSAT_3_22]]
  nSAT_4 = [[nSAT_4_11] [nSAT_4_12]; [nSAT_4_21] [nSAT_4_22]]


  nSAT = (nSAT_1, nSAT_2, nSAT_3, nSAT_4)     
  # }}}


  # {{{ Interface Penalties #TODO Edit for plane strain
  # Z' penalty parameters (for interface Dirichlet conditions) in logical space:
  
  IsJZ_1_11 =  β * (d ./ (4 * h1)) * ((sJ1 * C11_1) ⊗ Ir) # TODO check sJ in these
  IsJZ_1_12 =  β * (d ./ (4 * h1)) * ((sJ1 * C13_1) ⊗ Ir)
  IsJZ_1_21 =  β * (d ./ (4 * h1)) * ((sJ1 * C31_1) ⊗ Ir)
  IsJZ_1_22 =  β * (d ./ (4 * h1)) * ((sJ1 * C33_1) ⊗ Ir) 

  IsJZ_2_11 =  β * (d ./ (4 * h2)) * ((sJ2 * C11_2) ⊗ Ir)
  IsJZ_2_12 =  β * (d ./ (4 * h2)) * ((sJ2 *C13_2) ⊗ Ir)
  IsJZ_2_21 =  β * (d ./ (4 * h2)) * ((sJ2 *C31_2) ⊗ Ir)
  IsJZ_2_22 =  β * (d ./ (4 * h2)) * ((sJ2 *C33_2) ⊗ Ir) 

  IsJZ_3_11 =  β * (d ./ (4 * h3)) * (Is ⊗ (sJ3 * C22_3)) 
  IsJZ_3_12 =  β * (d ./ (4 * h3)) * (Is ⊗ (sJ3 * C24_3))  
  IsJZ_3_21 =  β * (d ./ (4 * h3)) * (Is ⊗ (sJ3 * C42_3)) 
  IsJZ_3_22 =  β * (d ./ (4 * h3)) * (Is ⊗ (sJ3 * C44_3)) 

  IsJZ_4_11 =  β * (d ./ (4 * h4)) * (Is ⊗ (sJ4 * C22_4)) 
  IsJZ_4_12 =  β * (d ./ (4 * h4)) * (Is ⊗ (sJ4 * C24_4))  
  IsJZ_4_21 =  β * (d ./ (4 * h4)) * (Is ⊗ (sJ4 * C42_4)) 
  IsJZ_4_22 =  β * (d ./ (4 * h4)) * (Is ⊗ (sJ4 * C44_4)) 


  IsJZ_1 = [[IsJZ_1_11] [IsJZ_1_12]; [IsJZ_1_21] [IsJZ_1_22]]
  IsJZ_2 = [[IsJZ_2_11] [IsJZ_2_12]; [IsJZ_2_21] [IsJZ_2_22]]
  IsJZ_3 = [[IsJZ_3_11] [IsJZ_3_12]; [IsJZ_3_21] [IsJZ_3_22]]
  IsJZ_4 = [[IsJZ_4_11] [IsJZ_4_12]; [IsJZ_4_21] [IsJZ_4_22]]


  IsJZ = (IsJZ_1, IsJZ_2, IsJZ_3, IsJZ_4)
     


  
  # {{{ RHS VECTOR
  # boundary data operators for Dirichlet faces 
  dB1_11 = Hinv * (T11_1 .- sJZ11_1)' * eRS[1] * H1    # multiplies first variable data, added to first
  dB1_12 = Hinv * (T12_1 .- sJZ12_1)' * eRS[1] * H1    # mulitopleis first vrriable data, added to second
  dB1_21 = Hinv * (T21_1 .- sJZ21_1)' * eRS[1] * H1    # mulitiplies second variable, added to first
  dB1_22 = Hinv * (T22_1 .- sJZ22_1)' * eRS[1] * H1    # mulitplies second var, added to second
  
  dB2_11 = Hinv * (T11_2 .- sJZ11_2)' * eRS[2] * H2    # multiplies first variable data, added to first
  dB2_12 = Hinv * (T12_2 .- sJZ12_2)' * eRS[2] * H2    # mulitopleis first vrriable data, added to second
  dB2_21 = Hinv * (T21_2 .- sJZ21_2)' * eRS[2] * H2    # mulitiplies second variable, added to first
  dB2_22 = Hinv * (T22_2 .- sJZ22_2)' * eRS[2] * H2    # mulitplies second var, added to second
  
  dB3_11 = Hinv * (T11_3 .- sJZ11_3)' * eRS[3] * H3    # multiplies first variable data, added to first
  dB3_12 = Hinv * (T12_3 .- sJZ12_3)' * eRS[3] * H3    # mulitopleis first vrriable data, added to second
  dB3_21 = Hinv * (T21_3 .- sJZ21_3)' * eRS[3] * H3    # mulitiplies second variable, added to first
  dB3_22 = Hinv * (T22_3 .- sJZ22_3)' * eRS[3] * H3    # mulitplies second var, added to second
  
  dB4_11 = Hinv * (T11_4 .- sJZ11_4)' * eRS[4] * H4    # multiplies first variable data, added to first
  dB4_12 = Hinv * (T12_4 .- sJZ12_4)' * eRS[4] * H4    # mulitopleis first vrriable data, added to second
  dB4_21 = Hinv * (T21_4 .- sJZ21_4)' * eRS[4] * H4    # mulitiplies second variable, added to first
  dB4_22 = Hinv * (T22_4 .- sJZ22_4)' * eRS[4] * H4    # mulitplies second var, added to second
  
  
  
  
  dB1 = [[dB1_11] [dB1_12]; [dB1_21] [dB1_22]]
  dB2 = [[dB2_11] [dB2_12]; [dB2_21] [dB2_22]]
  dB3 = [[dB3_11] [dB3_12]; [dB3_21] [dB3_22]]
  dB4 = [[dB4_11] [dB4_12]; [dB4_21] [dB4_22]]

  dB = (dB1, dB2, dB3, dB4)


  # boundary data operators for Neumann faces 
  nB1 = -Hinv * eRS[1] * H1
  nB2 = -Hinv * eRS[2] * H2 
  nB3 = -Hinv * eRS[3] * H3 
  nB4 = -Hinv * eRS[4] * H4 
  nB = (nB1, nB2, nB3, nB4)


  Drs_C12 = (Is ⊗ Dr) * sparse(1:length(c12), 1:length(c12), view(c12, :)) * (Ds ⊗ Ir)
  Dsr_C21 = (Ds ⊗ Ir) * sparse(1:length(c21), 1:length(c21), view(c21, :)) * (Is ⊗ Dr)

  Drs_C14 = (Is ⊗ Dr) * sparse(1:length(c14), 1:length(c14), view(c14, :)) * (Ds ⊗ Ir)
  Dsr_C23 = (Ds ⊗ Ir) * sparse(1:length(c23), 1:length(c23), view(c23, :)) * (Is ⊗ Dr)

  Drs_C32 = (Is ⊗ Dr) * sparse(1:length(c32), 1:length(c32), view(c32, :)) * (Ds ⊗ Ir)
  Dsr_C41 = (Ds ⊗ Ir) * sparse(1:length(c41), 1:length(c41), view(c41, :)) * (Is ⊗ Dr)

  Drs_C34 = (Is ⊗ Dr) * sparse(1:length(c34), 1:length(c34), view(c34, :)) * (Ds ⊗ Ir)
  Dsr_C43 = (Ds ⊗ Ir) * sparse(1:length(c43), 1:length(c43), view(c43, :)) * (Is ⊗ Dr)


  # {{{ Build component operators for element local
  A11 = Drr_11 + Drs_C12 + Dsr_C21 + Dss_22;
  A12 = Drr_13 + Drs_C14 + Dsr_C23 + Dss_24;

  A21 = Drr_31 + Drs_C32 + Dsr_C41 + Dss_42;
  A22 = Drr_33 + Drs_C34 + Dsr_C43 + Dss_44;

  # }}} 

  # {{{ Build operators to connect interface (modify operator to handle the boundary conditions)
  bctype=(BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE)
  
  for lf = 1:4
    if LFToB[lf] == BC_NEUMANN
      A11 += nSAT[lf][1, 1]
      A12 += nSAT[lf][1, 2] 
      A21 += nSAT[lf][2, 1] 
      A22 += nSAT[lf][2, 2]

    elseif LFToB[lf] == BC_DIRICHLET
      A11 += dSAT[lf][1, 1]
      A12 += dSAT[lf][2, 1]# these seem transposed, but are not.
      A21 += dSAT[lf][1, 2]
      A22 += dSAT[lf][2, 2]
    elseif (LFToB[lf] == BC_LOCKED_INTERFACE) || (LFToB[lf] == BC_JUMP_INTERFACE) || (LFToB[lf] == RS_FAULT) || (LFToB[lf] == VP_FAULT)
      # this happens in global assembly
    else  
      print("invalid bc")
     end
  end
  bctype=(LFToB[1], LFToB[2], LFToB[3], LFToB[4])
  # }}}

  # Full operator
  A = (A11, A12, A21, A22)
 

  return (A, dB, nB, eRST, H̃, T, TT, eRS, J, coord = metrics.coord, IsJZ, Hinv, H,
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
  s = source(xloc[:], yloc[:], volargs...)
  ge[:, 1] += Jloc * s[:,1]
  ge[:, 2] += Jloc * s[:,2]
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
  ge[:, :] .= 0
  for lf = 1:4
    if LFToB[lf] == BC_DIRICHLET
      # fault (Dirichlet):
      gD = bc_Dirichlet(lf, xf[lf], yf[lf], bcargs...) # this should return a two column vector
      vf_1 = dB[lf][1, 1]*gD[:, 1] + dB[lf][2, 1]*gD[:, 2]
      vf_2 = dB[lf][1, 2]*gD[:, 1] + dB[lf][2, 2]*gD[:, 2]
      

    elseif LFToB[lf] == BC_NEUMANN
      gN = bc_Neumann(lf, xf[lf], yf[lf], nx[lf], ny[lf], bcargs...) # this should return a two column vector
      vf_1 = nB[lf] * (sJ[lf] .* gN[:, 1])
      vf_2 = nB[lf] * (sJ[lf] .* gN[:, 2])


    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here

    elseif LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
      
      f = EToF[lf, e]  # get global face number

      (em, ep) = FToE[:, f]  # find the two elements that share global face f; ep is on plus side
      (fm, fp) = FToLF[:, f] # Local face numbers. 

      if em == e  # if I'm currently working on the minus side:
        nf = fp   # then set other face to plus side
      else
        nf = fm   #otherwise set other face to minus side
      end
      


      sJZm11 = 1*(lop.IsJZ[lf][1, 1])' + (neighborZ[nf][1, 1])'
      sJZm12 = 1*(lop.IsJZ[lf][1, 2])' + (neighborZ[nf][1, 2])'
      sJZm21 = 1*(lop.IsJZ[lf][2, 1])' + (neighborZ[nf][2, 1])'
      sJZm22 = 1*(lop.IsJZ[lf][2, 2])' + (neighborZ[nf][2, 2])'
  

      # Jump in displacement

      gJ = in_jump(lf, xf[lf], yf[lf], bcargs...)  # this should return a two column vector containing slip components

      cc = 0
      A11 = -lop.Hinv * sJZm11 * lop.eRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][1, 1] * lop.eRS[lf] * lop.H[lf]
      A21 = -lop.Hinv * sJZm21 * lop.eRS[lf] * lop.H[lf] + cc * 0.5 * lop.Hinv * lop.TT[lf][1, 2] * lop.eRS[lf] * lop.H[lf]
      A12 = -lop.Hinv * sJZm12 * lop.eRS[lf] * lop.H[lf] + cc * 0.5 * lop.Hinv * lop.TT[lf][2, 1] * lop.eRS[lf] * lop.H[lf]
      A22 = -lop.Hinv * sJZm22 * lop.eRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][2, 2] * lop.eRS[lf] * lop.H[lf]
      
      # Jump in traction
      gT = in_tractionjump(lf, xf[lf], yf[lf], bcargs...) # this should return a two column vector
    
      B = -0.5 * lop.Hinv * lop.eRS[lf] * lop.H[lf]
    
      vf_1 = A11*gJ[:,1] + A21*gJ[:,2] + B*(sJ[lf] .* gT[:, 1]) 
      vf_2 = A12*gJ[:,1] + A22*gJ[:,2] + B*(sJ[lf] .* gT[:, 2]) 

    else
      error("invalid bc")
    end

    ge[:, 1] += vf_1 
    ge[:, 2] += vf_2 
  end

end




function loc_bdry_vec_v2!(ge, lop, neighborZ, LFToB, EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, bcargs = ())

  e = bcargs[1]
  dB = lop.dB 
  nB = lop.nB 



  (xf, yf) = lop.facecoord
  sJ = lop.sJ 
  nx = lop.nx
  ny = lop.ny
  ge[:, :] .= 0
  for lf = 1:4
    if LFToB[lf] == BC_DIRICHLET

      gD = bc_Dirichlet(lf, xf[lf], yf[lf], bcargs...) # this should return B[1]*g
   
      vf_1 = dB[lf][1, 1]*gD[1] + dB[lf][2, 1]*gD[2]
      vf_2 = dB[lf][1, 2]*gD[1] + dB[lf][2, 2]*gD[2]
    
      ge[:, 1] += vf_1 
      ge[:, 2] += vf_2 
 
    elseif LFToB[lf] == BC_NEUMANN
 
      gN = bc_Neumann(lf, xf[lf], yf[lf], nx[lf], ny[lf], bcargs...) #this should return B[3]*g3
   
      vf_1 = nB[lf] * (sJ[lf] .* gN[1])
      vf_2 = nB[lf] * (sJ[lf] .* gN[2])
      
      ge[:, 1] += vf_1 
      ge[:, 2] += vf_2 

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


      sJZm11 = 1*(lop.IsJZ[lf][1, 1])' + (neighborZ[nf][1, 1])'
      sJZm12 = 1*(lop.IsJZ[lf][1, 2])' + (neighborZ[nf][1, 2])'
      sJZm21 = 1*(lop.IsJZ[lf][2, 1])' + (neighborZ[nf][2, 1])'
      sJZm22 = 1*(lop.IsJZ[lf][2, 2])' + (neighborZ[nf][2, 2])'
  
      # Jump in displacement

      δ = in_jump(lf, xf[lf], yf[lf], bcargs...)  # this should return a two column vector containing slip components (parallel followed by normal)

      # Calculate jumps in horizontal (u) and vertical (w) displacement: 

      du =  δ[:, 1] .* abs.(ny[lf])  # horizontal jump
      dw =  δ[:, 1] .* abs.(nx[lf])  #vertical jump


      cc = 0 # TODO: address this with A&D
      A11 = -lop.Hinv * sJZm11 * lop.eRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][1, 1] * lop.eRS[lf] * lop.H[lf]
      A21 = -lop.Hinv * sJZm21 * lop.eRS[lf] * lop.H[lf] + cc * 0.5 * lop.Hinv * lop.TT[lf][1, 2] * lop.eRS[lf] * lop.H[lf]
      A12 = -lop.Hinv * sJZm12 * lop.eRS[lf] * lop.H[lf] + cc * 0.5 * lop.Hinv * lop.TT[lf][2, 1] * lop.eRS[lf] * lop.H[lf]
      A22 = -lop.Hinv * sJZm22 * lop.eRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][2, 2] * lop.eRS[lf] * lop.H[lf]
      
  3
      vf_1 = A11*du + A21*dw
      vf_2 = A12*du + A22*dw

      ge[:, 1] += vf_1 
      ge[:, 2] += vf_2 

    else
      error("invalid bc")
    end

   
   
   
  end

end



function computetraction(lop, lf, u1, u2)


  T11 = lop.T[lf][1, 1]
  T12 = lop.T[lf][1, 2]
  T21 = lop.T[lf][2, 1]
  T22 = lop.T[lf][2, 2]

  sJ = lop.sJ[lf]
  eRST = lop.eRST[lf]

  t1 = (eRST * (T11 * u1 + T12 * u2)) ./ sJ   
  t2 = (eRST * (T21 * u1 + T22 * u2)) ./ sJ  

  return (t1, t2)
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

  # Compute the total number of volume and jump (δ) points
  VNp = vstarts[nelems+1]-1
  #δNp = FToδstarts[nfaces+1]-1

  M11 = spzeros(VNp, VNp)
  M12 = spzeros(VNp, VNp)
  M21 = spzeros(VNp, VNp)
  M22 = spzeros(VNp, VNp)
  for e = 1:nelems
    M11[vstarts[e]:vstarts[e+1]-1, vstarts[e]:vstarts[e+1]-1] = lop[e].A[1]  # put local A along block diagonal
    M12[vstarts[e]:vstarts[e+1]-1, vstarts[e]:vstarts[e+1]-1] = lop[e].A[2]  # put local A along block diagonal
    M21[vstarts[e]:vstarts[e+1]-1, vstarts[e]:vstarts[e+1]-1] = lop[e].A[3]  # put local A along block diagonal
    M22[vstarts[e]:vstarts[e+1]-1, vstarts[e]:vstarts[e+1]-1] = lop[e].A[4]  # put local A along block diagonal

  end

  for f = 1:nfaces  # Loop over the global faces to put in off-diagonal interface ops.
    
    if (FToB[f] == BC_JUMP_INTERFACE) || (FToB[f] == BC_LOCKED_INTERFACE) || (FToB[f] == RS_FAULT) || (FToB[f] == VP_FAULT)  # check if face is an interface
   
      (em, ep) = FToE[:, f]  # find the two elements that share global face f.
      (fm, fp) = FToLF[:, f]  # find corresponding local faces corresponding to global face f.
     
      sJZ_11 = (lop[em].IsJZ[fm][1, 1])' + (lop[ep].IsJZ[fp][1, 1])'
      sJZ_12 = (lop[em].IsJZ[fm][1, 2])' + (lop[ep].IsJZ[fp][1, 2])'
      sJZ_21 = (lop[em].IsJZ[fm][2, 1])' + (lop[ep].IsJZ[fp][2, 1])'
      sJZ_22 = (lop[em].IsJZ[fm][2, 2])' + (lop[ep].IsJZ[fp][2, 2])'

      Tm_11 = lop[em].T[fm][1, 1]
      Tm_12 = lop[em].T[fm][1, 2]
      Tm_21 = lop[em].T[fm][2, 1]
      Tm_22 = lop[em].T[fm][2, 2]

      Tp_11 = lop[ep].T[fp][1, 1]
      Tp_12 = lop[ep].T[fp][1, 2]
      Tp_21 = lop[ep].T[fp][2, 1]
      Tp_22 = lop[ep].T[fp][2, 2]

      
      TTm_11 = lop[em].TT[fm][1, 1]
      TTm_12 = lop[em].TT[fm][1, 2]
      TTm_21 = lop[em].TT[fm][2, 1]
      TTm_22 = lop[em].TT[fm][2, 2]

      TTp_11 = lop[ep].TT[fp][1, 1]
      TTp_12 = lop[ep].TT[fp][1, 2]
      TTp_21 = lop[ep].TT[fp][2, 1]
      TTp_22 = lop[ep].TT[fp][2, 2]


        cc = 0 # TODO: Address with A&D 

      # Local effects, for both plus and minus sides of face f, they should match signs:
      M11[vstarts[em]:vstarts[em+1]-1, vstarts[em]:vstarts[em+1]-1] += (-lop[em].Hinv * (sJZ_11 - 1 * 0.5 * TTm_11) * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] + 
                                                                      -  0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] * Tm_11)
      M11[vstarts[ep]:vstarts[ep+1]-1, vstarts[ep]:vstarts[ep+1]-1] += (-lop[ep].Hinv * (sJZ_11 - 1 * 0.5 * TTp_11) * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] + 
                                                                      -  0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] * Tp_11)

      M12[vstarts[em]:vstarts[em+1]-1, vstarts[em]:vstarts[em+1]-1] += (-lop[em].Hinv * (sJZ_21 - cc * 0.5 * TTm_12) * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] + 
                                                                      -  0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] * Tm_12)
      M12[vstarts[ep]:vstarts[ep+1]-1, vstarts[ep]:vstarts[ep+1]-1] +=  (-lop[ep].Hinv * (sJZ_21 - cc * 0.5 * TTp_12) * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] + 
                                                                      -  0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] * Tp_12)

      M21[vstarts[em]:vstarts[em+1]-1, vstarts[em]:vstarts[em+1]-1] += (-lop[em].Hinv * (sJZ_12 - cc * 0.5 * TTm_21) * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] + 
                                                                      -  0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] * Tm_21)
      M21[vstarts[ep]:vstarts[ep+1]-1, vstarts[ep]:vstarts[ep+1]-1] += (-lop[ep].Hinv * (sJZ_12 - cc * 0.5 * TTp_21) * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] + 
                                                                      -  0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] * Tp_21)

      M22[vstarts[em]:vstarts[em+1]-1, vstarts[em]:vstarts[em+1]-1] += (-lop[em].Hinv * (sJZ_22 - 1 * 0.5 * TTm_22) * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] + 
                                                                      -  0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * lop[em].eRST[fm] * Tm_22)
      M22[vstarts[ep]:vstarts[ep+1]-1, vstarts[ep]:vstarts[ep+1]-1] += (-lop[ep].Hinv * (sJZ_22 - 1 * 0.5 * TTp_22) * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] + 
                                                                      -  0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * lop[ep].eRST[fp] * Tp_22)


      # If face orientation on plus side does not match minus side, then flip
      Np = (fm <= 2 ? Ns[em]+1 : Nr[em]+1) # Set Np to be either Nr or Ns, depending on face
      if EToO[fp, ep]                     # orientation matches
        R = sparse(I, Np, Np)         
      else  
        R = rotr90(sparse(I, Np, Np)) # orientation doesn't match
      end

  
      # account for opposite side; 
      M11[vstarts[em]:vstarts[em+1]-1, vstarts[ep]:vstarts[ep+1]-1] +=  (lop[em].Hinv * (sJZ_11 - 1 * 0.5 * TTm_11) * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] + 
                                                                      -  0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] * Tp_11)
      M11[vstarts[ep]:vstarts[ep+1]-1, vstarts[em]:vstarts[em+1]-1] += (lop[ep].Hinv * (sJZ_11 - 1 * 0.5 * TTp_11) * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] + 
                                                                      -  0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] * Tm_11)

      M12[vstarts[em]:vstarts[em+1]-1, vstarts[ep]:vstarts[ep+1]-1] += (lop[em].Hinv * (sJZ_21 - cc * 0.5 * TTm_12) * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] + 
                                                                      -  0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] * Tp_12)

      M12[vstarts[ep]:vstarts[ep+1]-1, vstarts[em]:vstarts[em+1]-1] += (lop[ep].Hinv * (sJZ_21 - cc * 0.5 * TTp_12) * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] + 
                                                                      -  0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] * Tm_12)

      
      M21[vstarts[em]:vstarts[em+1]-1, vstarts[ep]:vstarts[ep+1]-1] += (lop[em].Hinv * (sJZ_12 - cc * 0.5 * TTm_21) * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] + 
                                                                      -  0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] * Tp_21)
      M21[vstarts[ep]:vstarts[ep+1]-1, vstarts[em]:vstarts[em+1]-1] += (lop[ep].Hinv * (sJZ_12 - cc * 0.5 * TTp_21) * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] + 
                                                                      -  0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] * Tm_21)
         
      M22[vstarts[em]:vstarts[em+1]-1, vstarts[ep]:vstarts[ep+1]-1] +=  (lop[em].Hinv * (sJZ_22 - 1 * 0.5 * TTm_22) * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] + 
                                                                      -  0.5*lop[em].Hinv * lop[em].eRS[fm] * lop[em].H[fm] * R * lop[ep].eRST[fp] * Tp_22)
      
        
      M22[vstarts[ep]:vstarts[ep+1]-1, vstarts[em]:vstarts[em+1]-1] += (lop[ep].Hinv * (sJZ_22 - 1 * 0.5 * TTp_22) * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] + 
                                                                      -  0.5*lop[ep].Hinv * lop[ep].eRS[fp] * lop[ep].H[fp] * R * lop[em].eRST[fm] * Tm_22)
        
    end
  end

  M = [M11 M12; M21 M22]
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
export bdry_vec!, loc_bdry_vec_v2!
export rateandstate, newtbndv, bcstarts, connectivityarrays, flatten_tuples
export LocalGlobalOperators, SBPLocalOperator1, global_operator
