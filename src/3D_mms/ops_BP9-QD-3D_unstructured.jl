using SparseArrays
using LinearAlgebra
using SparseArrayKit
#using UnicodePlots

⊗(A,B) = kron(A, B)

const BC_DIRICHLET        = 1
const BC_NEUMANN          = 2
const EARTH_SURFACE       = 3
const BC_LOCKED_INTERFACE = 0
const BC_JUMP_INTERFACE   = 7


include("SBP_3D.jl")

function rotate_diagonal(M::AbstractMatrix)
    # diagind(M) gets the linear indices of the main diagonal
    # diag(M) extracts the diagonal elements
    N = copy(M)
    N[diagind(N)] .= reverse(diag(N))
    return N
end



function create_metrics(pm, Nq, Nr, Ns, λ, μ, 
                        xf=(q,r,s)->(q, ones(size(q)), zeros(size(r)), zeros(size(s))),
                        yf=(q,r,s)->(r, zeros(size(q)), ones(size(r)), zeros(size(s))),
                        zf=(q,r,s)->(s, zeros(size(q)), zeros(size(r)), ones(size(s))))
  
  Nqp = Nq + 1
  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nqp * Nrp * Nsp

  # Derivative operators for the metric terms
  @assert pm <= 8
  pp = pm == 6 ? 8 : pm

  q = range(-1, stop=1, length=Nqp)
  r = range(-1, stop=1, length=Nrp)
  s = range(-1, stop=1, length=Nsp)

  # Create the mesh
 

  q = (q * ones(1, Nrp)) .* reshape(ones(Nsp), 1, 1, :)
  r = (ones(Nqp) * r') .* reshape(ones(Nsp), 1, 1, :)
  s = (ones(Nqp) * ones(1, Nrp)) .* reshape(s, 1, 1, :) # possible source for error? 



  (x, xq, xr, xs) = xf(q, r, s)
  (y, yq, yr, ys) = yf(q, r, s)
  (z, zq, zr, zs) = zf(q, r, s)



  J = xq .* (yr .* zs - zr .* ys) - xr .* (yq .* zs - zq .* ys) + xs .* (yq .* zr - zq .* yr)
  @assert minimum(J) >= 1e-12

  qx =  (yr .* zs - ys .* zr) ./ J 
  qy = -(xr .* zs - xs .* zr) ./ J 
  qz =  (xr .* ys - xs .* yr) ./ J 
  rx =  -(yq .* zs - ys .* zq) ./ J
  ry =  (xq .* zs - xs .* zq) ./ J
  rz = -(xq .* ys - xs .* yq) ./ J
  sx =  (yq .* zr - yr .* zq) ./ J
  sy = -(xq .* zr - xr .* zq) ./ J
  sz =  (xq .* yr - xr .* yq) ./ J

  qrs = (qx, qy, qz, rx, ry, rz, sx, sy, sz)


  # variable coefficient matrix components
  #C = fill(Array{Float64, 3}(undef, 3, 4, 5), 3, 3)

  C = fill(zeros(Nqp, Nrp, Nsp), 3, 3, 3, 3)
  for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                C[i, j, k, l] = J .* (λ(x, y, z) .* F(j, i, qrs) .* F(l, k, qrs) + 
                                   μ(x, y, z) .* (F(1, i, qrs) .* delt(j, l) .* F(1, k, qrs) + F(2, i, qrs) .* delt(j, l) .* F(2, k, qrs) + F(3, i, qrs) .* delt(j, l) .* F(3, k, qrs)) +
                                   μ(x, y, z) .* (F(l, i, qrs) .* F(j, k, qrs)))
            end
        end
    end
  end

  #
  # Block surface matrices
  #

  (xf1, yf1, zf1) = (view(x, 1, :, :), view(y, 1, :, :), view(z, 1, :, :))
  nx1 =   ys[1, :, :] .* zr[1, :, :] - yr[1, :, :] .* zs[1, :, :]
  ny1 = -(xs[1, :, :] .* zr[1, :, :] - xr[1, :, :] .* zs[1, :, :])
  nz1 =   xs[1, :, :] .* yr[1, :, :] - xr[1, :, :] .* ys[1, :, :]
  sJ1 = hypot.(nx1, ny1, nz1)
  nx1 = nx1 ./ sJ1
  ny1 = ny1 ./ sJ1
  nz1 = nz1 ./ sJ1


  (xf2, yf2, zf2) = (view(x, Nqp, :, :), view(y, Nqp, :, :), view(z, Nqp, :, :))
  nx2 =    yr[end, :, :] .* zs[end, :, :] - ys[end, :, :] .* zr[end, :, :]
  ny2 =  -(xr[end, :, :] .* zs[end, :, :] - xs[end, :, :] .* zr[end, :, :])
  nz2 =    xr[end, :, :] .* ys[end, :, :] - xs[end, :, :] .* yr[end, :, :]
  sJ2 = hypot.(nx2, ny2, nz2)
  nx2 = nx2 ./ sJ2
  ny2 = ny2 ./ sJ2
  nz2 = nz2 ./ sJ2

  (xf3, yf3, zf3) = (view(x, :, 1, :), view(y, :, 1, :), view(z, :, 1, :))
  nx3 =   yq[:, 1, :] .* zs[:, 1, :] - ys[:, 1, :] .* zq[:, 1, :]
  ny3 = -(xq[:, 1, :] .* zs[:, 1, :] - xs[:, 1, :] .* zq[:, 1, :])
  nz3 =   xq[:, 1, :] .* ys[:, 1, :] - xs[:, 1, :] .* yq[:, 1, :]
  sJ3 = hypot.(nx3, ny3, nz3)
  nx3 = nx3 ./ sJ3
  ny3 = ny3 ./ sJ3
  nz3 = nz3 ./ sJ3


  (xf4, yf4, zf4) = (view(x, :, Nrp, :), view(y, :, Nrp, :), view(z, :, Nrp, :))
  nx4 =    ys[:, end, :] .* zq[:, end, :] - yq[:, end, :] .* zs[:, end, :]
  ny4 =  -(xs[:, end, :] .* zq[:, end, :] - xq[:, end, :] .* zs[:, end, :])
  nz4 =    xs[:, end, :] .* yq[:, end, :] - xq[:, end, :] .* ys[:, end, :]
  sJ4 = hypot.(nx4, ny4, nz4)
  nx4 = nx4 ./ sJ4
  ny4 = ny4 ./ sJ4
  nz4 = nz4 ./ sJ4

  (xf5, yf5, zf5) = (view(x, :, :, 1), view(y, :, :, 1), view(z, :, :, 1))
  nx5 =   yr[:, :, 1] .* zq[:, :, 1] - yq[:, :, 1] .* zr[:, :, 1]
  ny5 = -(xr[:, :, 1] .* zq[:, :, 1] - xq[:, :, 1] .* zr[:, :, 1])
  nz5 =  xr[:, :, 1] .* yq[:, :, 1] - xq[:, :, 1] .* yr[:, :, 1]
  sJ5 = hypot.(nx5, ny5, nz5)
  nx5 = nx5 ./ sJ5
  ny5 = ny5 ./ sJ5
  nz5 = nz5 ./ sJ5

  (xf6, yf6, zf6) = (view(x, :, :, Nsp), view(y, :, :, Nsp), view(z, :, :, Nsp))
  nx6 =   yq[:, :, end] .* zr[:, :, end] - yr[:, :, end] .* zq[:, :, end]
  ny6 = -(xq[:, :, end] .* zr[:, :, end] - xr[:, :, end] .* zq[:, :, end])
  nz6 =   xq[:, :, end] .* yr[:, :, end] - xr[:, :, end] .* yq[:, :, end]
 
  sJ6 = hypot.(nx6, ny6, nz6)
  nx6 = nx6 ./ sJ6
  ny6 = ny6 ./ sJ6
  nz6 = nz6 ./ sJ6


  (coord = (x,y,z),
   facecoord = ((xf1, xf2, xf3, xf4, xf5, xf6), (yf1, yf2, yf3, yf4, yf5, yf6), (zf1, zf2, zf3, zf4, zf5, zf6)),
   C = C,
   J=J,
   sJ = (sJ1[:], sJ2[:], sJ3[:], sJ4[:], sJ5[:], sJ6[:]),
   nx = (nx1[:], nx2[:], nx3[:], nx4[:], nx5[:], nx6[:]),
   ny = (ny1[:], ny2[:], ny3[:], ny4[:], ny5[:], ny6[:]),
   nz = (nz1[:], nz2[:], nz3[:], nz4[:], nz5[:], nz6[:]),
   qx = qx, qy = qy, qz = qz, 
   rx = rx, ry = ry, rz = rz,
   sx = sx, sy = sy, sz = sz)
end


function delt(i, j)
    if i == j 
        return 1
    else
        return 0
    end
end

function F(I, i, qrs)
    (qx, qy, qz, rx, ry, rz, sx, sy, sz) = qrs
    # F(I, i) = dxi/dXI, e.g. F(1, 1) = qx
    if i == 1
        if I == 1
            return qx
        elseif I == 2
            return qy
        else I == 3
            return qz
        end
    elseif i == 2
        if I == 1
            return rx
        elseif I == 2
            return ry
        else I == 3
            return rz
        end
    else i == 3
        if I == 1
            return sx
        elseif I == 2
            return sy
        else I == 3
            return sz
        end
    end
end



function locoperator(p, Nq, Nr, Ns, μ, λ, metrics=create_metrics(p, Nq, Nr, Ns, μ, λ),  
            LFToB = (BC_DIRICHLET, BC_DIRICHLET, BC_DIRICHLET,
                     BC_DIRICHLET, BC_DIRICHLET, BC_DIRICHLET);
                     C = metrics.C)
  

  # {{{ read in parameters

  
  C = metrics.C
  
    Nqp = Nq + 1
    Nrp = Nr + 1
    Nsp = Ns + 1
    Np = Nqp * Nrp * Nsp

       # {{{ Jacobian
        J = spdiagm(0 => reshape(metrics.J, Np))
      # }}}

    c = fill(Diagonal(zeros(Nqp*Nrp*Nsp)), 3, 3, 3, 3)

    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    c[i, j, k, l] = Diagonal(C[i, j, k, l][:])
                end
            end
        end
    end

  
  # }}}

  # {{{ Create 1D operators in logical space 
  
    (Dq, HqI, Hq, q) = diagonal_sbp_D1(p, Nq; xc = (-1,1))
    Qq = Hq * Dq
    QqT = sparse(transpose(Qq))

    (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
    Qr = Hr * Dr
    QrT = sparse(transpose(Qr))

    (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))
    Qs = Hs * Ds
    QsT = sparse(transpose(Qs))



  # {{{ 2D Identity matrices
  Iq = sparse(I, Nqp, Nqp)
  Ir = sparse(I, Nrp, Nrp)
  Is = sparse(I, Nsp, Nsp)
 
  # }}}

  # {{{ Set up the operators in logical space

    # create 1D operators:
    (D2q, S0q, SNq, _, _, _) = diagonal_sbp_D2(p, Nq; xc = (-1,1))
    (D2r, S0r, SNr, _, _, _) = diagonal_sbp_D2(p, Nr; xc = (-1,1))
    (D2s, S0s, SNs, _, _, _) = diagonal_sbp_D2(p, Ns; xc = (-1,1))

    D11 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D12 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D13 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D21 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D22 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D23 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D31 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D32 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D33 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)


     # Surface mass matrices
    #
    H1 = H2 = Hs ⊗ Hr
    H3 = H4 = Hs ⊗ Hq
    H5 = H6 = Hr ⊗ Hq
    H = (H1, H2, H3, H4, H5, H6)
    
    # Volume matrices
    Hinv = HsI ⊗ HrI ⊗ HqI

    H̃ = Hs ⊗ Hr ⊗ Hq           # Need this for convergence tests.

    # 3D first derivative operators:
    Dq3 = Is ⊗ Ir ⊗ Dq
    Dr3 = Is ⊗ Dr ⊗ Iq
    Ds3 = Ds ⊗ Ir ⊗ Iq

    # Face operators to reduce computations involing T11
    Sq = Is ⊗ Ir ⊗ (SNq+S0q)
    Sr = Is ⊗ (SNr+S0r) ⊗ Iq
    Ss = (SNs+S0s) ⊗ Ir ⊗ Iq
    

    for i = 1:3
        for j = 1:3
           (D11[i, j], _, _) = var_3D_D2q(p, Nqp, Nrp, Nsp, metrics.C[1, i, 1, j], HqI; xc = (-1, 1))
       
            D12[i, j] = Dq3 * c[1, i, 2, j] * Dr3  #c[1, i, 2, j] * (Is ⊗ Dr ⊗ Dq)
            D13[i, j] = Dq3 * c[1, i, 3, j] * Ds3  #c[1, i, 3, j] * (Ds ⊗ Ir ⊗ Dq)

            D21[i, j] = Dr3 * c[2, i, 1, j] * Dq3  #c[2, i, 1, j] * (Is ⊗ Dr ⊗ Dq)
           (D22[i, j], _, _) = var_3D_D2r(p, Nqp, Nrp, Nsp, metrics.C[2, i, 2, j], HrI; xc = (-1, 1))
            D23[i, j] = Dr3 * c[2, i, 3, j] * Ds3  #c[2, i, 3, j] * (Ds ⊗ Dr ⊗ Iq)

            D31[i, j] = Ds3 * c[3, i, 1, j] * Dq3  #c[3, i, 1, j] * (Ds ⊗ Ir ⊗ Dq)
            D32[i, j] = Ds3 * c[3, i, 2, j] * Dr3  #c[3, i, 2, j] * (Ds ⊗ Dr ⊗ Iq)
           (D33[i, j], _, _) = var_3D_D2s(p, Nqp, Nrp, Nsp, metrics.C[3, i, 3, j], HsI; xc = (-1, 1))
        end
    end

    


  # {{{ Boundary point matrices and surface mass matrices 

    eq0 = sparse([1  ], [1], [1], Nqp, 1)
    eqN = sparse([Nqp], [1], [1], Nqp, 1)
    er0 = sparse([1  ], [1], [1], Nrp, 1)
    erN = sparse([Nrp], [1], [1], Nrp, 1)
    es0 = sparse([1  ], [1], [1], Nsp, 1)
    esN = sparse([Nsp], [1], [1], Nsp, 1)
    e1 = Is ⊗ Ir ⊗ eq0
    e2 = Is ⊗ Ir ⊗ eqN
    e3 = Is ⊗ er0 ⊗ Iq
    e4 = Is ⊗ erN ⊗ Iq
    e5 = es0 ⊗ Ir ⊗ Iq
    e6 = esN ⊗ Ir ⊗ Iq

    eQRS = (e1, e2, e3, e4, e5, e6)
    eQRST = (e1', e2', e3', e4', e5', e6')


  # }}}

 
  # {{{ discrete boundary traction operators in logical space on each of the 6 faces:
 # Create traction operators on each face


    T11_1 = -1 * (c[1,1,1,1]*Sq + c[1,1,2,1]*Dr3 + c[1,1,3,1]*Ds3)
    T12_1 = -1 * (c[1,1,1,2]*Sq + c[1,1,2,2]*Dr3 + c[1,1,3,2]*Ds3) 
    T21_1 = -1 * (c[1,2,1,1]*Sq + c[1,2,2,1]*Dr3 + c[1,2,3,1]*Ds3)
    T13_1 = -1 * (c[1,1,1,3]*Sq + c[1,1,2,3]*Dr3 + c[1,1,3,3]*Ds3) 
    T31_1 = -1 * (c[1,3,1,1]*Sq + c[1,3,2,1]*Dr3 + c[1,3,3,1]*Ds3) 
    T22_1 = -1 * (c[1,2,1,2]*Sq + c[1,2,2,2]*Dr3 + c[1,2,3,2]*Ds3) 
    T23_1 = -1 * (c[1,2,1,3]*Sq + c[1,2,2,3]*Dr3 + c[1,2,3,3]*Ds3) 
    T32_1 = -1 * (c[1,3,1,2]*Sq + c[1,3,2,2]*Dr3 + c[1,3,3,2]*Ds3)
    T33_1 = -1 * (c[1,3,1,3]*Sq + c[1,3,2,3]*Dr3 + c[1,3,3,3]*Ds3) 
    
    # FACE 2
    # (Nq, Nr, Ns) = (1, 0, 0)
    T11_2 = -T11_1
    T12_2 = -T12_1
    T21_2 = -T21_1
    T13_2 = -T13_1
    T31_2 = -T31_1
    T22_2 = -T22_1
    T23_2 = -T23_1
    T32_2 = -T32_1
    T33_2 = -T33_1

    
    #   # FACE 3
    #     # (Nq, Nr, Ns) = (0, -1, 0)
    T11_3 = -1 * (c[2,1,1,1]*Dq3 + c[2,1,2,1]*Sr + c[2,1,3,1]*Ds3) 
    T12_3 = -1 * (c[2,1,1,2]*Dq3 + c[2,1,2,2]*Sr + c[2,1,3,2]*Ds3) 
    T21_3 = -1 * (c[2,2,1,1]*Dq3 + c[2,2,2,1]*Sr + c[2,2,3,1]*Ds3)
    T13_3 = -1 * (c[2,1,1,3]*Dq3 + c[2,1,2,3]*Sr + c[2,1,3,3]*Ds3) 
    T31_3 = -1 * (c[2,3,1,1]*Dq3 + c[2,3,2,1]*Sr + c[2,3,3,1]*Ds3) 
    T22_3 = -1 * (c[2,2,1,2]*Dq3 + c[2,2,2,2]*Sr + c[2,2,3,2]*Ds3) 
    T23_3 = -1 * (c[2,2,1,3]*Dq3 + c[2,2,2,3]*Sr + c[2,2,3,3]*Ds3) 
    T32_3 = -1 * (c[2,3,1,2]*Dq3 + c[2,3,2,2]*Sr + c[2,3,3,2]*Ds3)
    T33_3 = -1 * (c[2,3,1,3]*Dq3 + c[2,3,2,3]*Sr + c[2,3,3,3]*Ds3) 

    # FACE 4
    # (Nq, Nr, Ns) = (0, 1, 0)
    T11_4 = -T11_3
    T12_4 = -T12_3
    T21_4 = -T21_3
    T13_4 = -T13_3
    T31_4 = -T31_3
    T22_4 = -T22_3
    T23_4 = -T23_3
    T32_4 = -T32_3
    T33_4 = -T33_3


    # FACE 5
    # (Nq, Nr, Ns) = (0, 0, -1)
    T11_5 = -1 * (c[3,1,1,1]*Dq3 + c[3,1,2,1]*Dr3 + c[3,1,3,1]*Ss)
    T12_5 = -1 *  (c[3,1,1,2]*Dq3 + c[3,1,2,2]*Dr3 + c[3,1,3,2]*Ss)
    T21_5 = -1 * (c[3,2,1,1]*Dq3 + c[3,2,2,1]*Dr3 + c[3,2,3,1]*Ss)
    T13_5 = -1 * (c[3,1,1,3]*Dq3 + c[3,1,2,3]*Dr3 + c[3,1,3,3]*Ss)
    T31_5 = -1 * (c[3,3,1,1]*Dq3 + c[3,3,2,1]*Dr3 + c[3,3,3,1]*Ss)
    T22_5 = -1 * (c[3,2,1,2]*Dq3 + c[3,2,2,2]*Dr3 + c[3,2,3,2]*Ss)
    T23_5 = -1 * (c[3,2,1,3]*Dq3 + c[3,2,2,3]*Dr3 + c[3,2,3,3]*Ss)
    T32_5 = -1 * (c[3,3,1,2]*Dq3 + c[3,3,2,2]*Dr3 + c[3,3,3,2]*Ss)
    T33_5 = -1 * (c[3,3,1,3]*Dq3 + c[3,3,2,3]*Dr3 + c[3,3,3,3]*Ss)


    # FACE 6
    # (Nq, Nr, Ns) = (0, 0, 1)
    T11_6 = -T11_5 
    T12_6 = -T12_5 
    T21_6 = -T21_5 
    T13_6 = -T13_5 
    T31_6 = -T31_5 
    T22_6 = -T22_5 
    T23_6 = -T23_5 
    T32_6 = -T32_5 
    T33_6 = -T33_5 

  T_1 = [[T11_1] [T12_1] [T13_1]; [T21_1] [T22_1] [T23_1]; [T31_1] [T32_1] [T33_1]]
  T_2 = -T_1 #[[T11_2] [T12_2] [T13_2]; [T21_2] [T22_2] [T23_2]; [T31_2] [T32_2] [T33_2]]
  T_3 = [[T11_3] [T12_3] [T13_3]; [T21_3] [T22_3] [T23_3]; [T31_3] [T32_3] [T33_3]]
  T_4 = -T_3 #[[T11_4] [T12_4] [T13_4]; [T21_4] [T22_4] [T23_4]; [T31_4] [T32_4] [T33_4]]
  T_5 = [[T11_5] [T12_5] [T13_5]; [T21_5] [T22_5] [T23_5]; [T31_5] [T32_5] [T33_5]]
  T_6 = -T_5 #[[T11_6] [T12_6] [T13_6]; [T21_6] [T22_6] [T23_6]; [T31_6] [T32_6] [T33_6]]

  T = (T_1, T_2, T_3, T_4, T_5, T_6)   

  
   TT_1 = [[T11_1'] [T21_1'] [T31_1']; [T12_1'] [T22_1'] [T32_1']; [T13_1'] [T23_1'] [T33_1']]
   TT_2 = -TT_1 #[[T11_2'] [T21_2'] [T31_2']; [T12_2'] [T22_2'] [T32_2']; [T13_2'] [T23_2'] [T33_2']]
   TT_3 = [[T11_3'] [T21_3'] [T31_3']; [T12_3'] [T22_3'] [T32_3']; [T13_3'] [T23_3'] [T33_3']]
   TT_4 = -TT_3 #[[T11_4'] [T21_4'] [T31_4']; [T12_4'] [T22_4'] [T32_4']; [T13_4'] [T23_4'] [T33_4']]
   TT_5 = [[T11_5'] [T21_5'] [T31_5']; [T12_5'] [T22_5'] [T32_5']; [T13_5'] [T23_5'] [T33_5']]
   TT_6 = -TT_5 #[T11_6'] [T21_6'] [T31_6']; [T12_6'] [T22_6'] [T32_6']; [T13_6'] [T23_6'] [T33_6']]

  TT = (TT_1, TT_2, TT_3, TT_4, TT_5, TT_6)   
  # }}}
    


  # {{{ Dirichlet BC penalties
    beta = 1
    h1 = Hq[1,1] 
    h3 = Hq[1,1]
    h5 = Hs[1,1]
    d = 3 #dimension? 

     sJZ11_1 = sJZ11_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,1,1,1][:])
     sJZ12_1 = sJZ12_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,1,1,2][:])
     sJZ21_1 = sJZ21_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,2,1,1][:])
     sJZ13_1 = sJZ13_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,1,1,3][:])
     sJZ31_1 = sJZ31_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,3,1,1][:])
     sJZ22_1 = sJZ22_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,2,1,2][:])
     sJZ23_1 = sJZ23_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,2,1,3][:])
     sJZ32_1 = sJZ32_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,3,1,2][:])
     sJZ33_1 = sJZ33_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,3,1,3][:])

    
    #  sJZ11_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,1,1,1][:])
    #  sJZ12_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,1,1,2][:])
    #  sJZ21_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,2,1,1][:])
    #  sJZ13_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,1,1,3][:])
    #  sJZ31_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,3,1,1][:])
    #  sJZ22_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,2,1,2][:])
    #  sJZ23_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,2,1,3][:])
    #  sJZ32_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,3,1,2][:])
    #  sJZ33_2 = sparse(1:Np, 1:Np, (beta/h1) * d * C[1,3,1,3][:])

    

     # FACE 3 
    # (nq, nr, ns) = (0, -1, 0)
     sJZ11_3 = sJZ11_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,1,2,1][:])
     sJZ12_3 = sJZ12_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,1,2,2][:])
     sJZ21_3 = sJZ21_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,2,2,1][:])
     sJZ13_3 = sJZ13_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,1,2,3][:])
     sJZ31_3 = sJZ31_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,3,2,1][:])
     sJZ22_3 = sJZ22_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,2,2,2][:])
     sJZ23_3 = sJZ23_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,2,2,3][:])
     sJZ32_3 = sJZ32_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,3,2,2][:])
     sJZ33_3 = sJZ33_4 = sparse(1:Np, 1:Np, (beta/h3) * d * C[2,3,2,3][:])

   # FACE 4 
    # (nq, nr, ns) = (0, 1, 0)
    #  sJZ11_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,1,2,1][:])
    #  sJZ12_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,1,2,2][:])
    #  sJZ21_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,2,2,1][:])
    #  sJZ13_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,1,2,3][:])
    #  sJZ31_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,3,2,1][:])
    #  sJZ22_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,2,2,2][:])
    #  sJZ23_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,2,2,3][:])
    #  sJZ32_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,3,2,2][:])
    #  sJZ33_4 = sparse(1:Np, 1:Np, (beta/h1) * d * C[2,3,2,3][:])
    # FACE 5 
    # (nq, nr, ns) = (0, 0, -1)
     sJZ11_5 = sJZ11_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,1,3,1][:])
     sJZ12_5 = sJZ12_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,1,3,2][:])
     sJZ21_5 = sJZ21_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,2,3,1][:])
     sJZ13_5 = sJZ13_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,1,3,3][:])
     sJZ31_5 = sJZ31_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,3,3,1][:])
     sJZ22_5 = sJZ22_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,2,3,2][:])
     sJZ23_5 = sJZ23_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,2,3,3][:])
     sJZ32_5 = sJZ32_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,3,3,2][:])
     sJZ33_5 = sJZ33_6 = sparse(1:Np, 1:Np, (beta/h5) * d * C[3,3,3,3][:])

    # FACE 6 
    # (nq, nr, ns) = (0, 0, 1)
    #  sJZ11_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,1,3,1][:])
    #  sJZ12_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,1,3,2][:])
    #  sJZ21_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,2,3,1][:])
    #  sJZ13_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,1,3,3][:])
    #  sJZ31_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,3,3,1][:])
    #  sJZ22_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,2,3,2][:])
    #  sJZ23_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,2,3,3][:])
    #  sJZ32_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,3,3,2][:])
    #  sJZ33_6 = sparse(1:Np, 1:Np, (beta/h1) * d * C[3,3,3,3][:])

  sJZ_1 = [[sJZ11_1] [sJZ12_1] [sJZ13_1]; [sJZ21_1] [sJZ22_1] [sJZ23_1]; [sJZ31_1] [sJZ32_1] [sJZ33_1]]
  sJZ_2 = [[sJZ11_2] [sJZ12_2] [sJZ13_2]; [sJZ21_2] [sJZ22_2] [sJZ23_2]; [sJZ31_2] [sJZ32_2] [sJZ33_2]]
  sJZ_3 = [[sJZ11_3] [sJZ12_3] [sJZ13_3]; [sJZ21_3] [sJZ22_3] [sJZ23_3]; [sJZ31_3] [sJZ32_3] [sJZ33_3]]
  sJZ_4 = [[sJZ11_4] [sJZ12_4] [sJZ13_4]; [sJZ21_4] [sJZ22_4] [sJZ23_4]; [sJZ31_4] [sJZ32_4] [sJZ33_4]]
  sJZ_5 = [[sJZ11_5] [sJZ12_5] [sJZ13_5]; [sJZ21_5] [sJZ22_5] [sJZ23_5]; [sJZ31_5] [sJZ32_5] [sJZ33_5]]
  sJZ_6 = [[sJZ11_6] [sJZ12_6] [sJZ13_6]; [sJZ21_6] [sJZ22_6] [sJZ23_6]; [sJZ31_6] [sJZ32_6] [sJZ33_6]]

  sJZ = (sJZ_1, sJZ_2, sJZ_3, sJZ_4, sJZ_5, sJZ_6)   



 # {{{ RHS VECTOR
  # boundary data operators for Dirichlet faces 
  dB1_11 = Hinv * (T11_1 .- sJZ11_1)' * eQRS[1] * H1    # 
  dB1_12 = Hinv * (T12_1 .- sJZ12_1)' * eQRS[1] * H1    # 
  dB1_13 = Hinv * (T13_1 .- sJZ13_1)' * eQRS[1] * H1    # 
  dB1_21 = Hinv * (T21_1 .- sJZ21_1)' * eQRS[1] * H1    # 
  dB1_22 = Hinv * (T22_1 .- sJZ22_1)' * eQRS[1] * H1    # 
  dB1_23 = Hinv * (T23_1 .- sJZ23_1)' * eQRS[1] * H1    # 
  dB1_31 = Hinv * (T31_1 .- sJZ31_1)' * eQRS[1] * H1    # 
  dB1_32 = Hinv * (T32_1 .- sJZ32_1)' * eQRS[1] * H1    # 
  dB1_33 = Hinv * (T33_1 .- sJZ33_1)' * eQRS[1] * H1    # 

  dB2_11 = Hinv * (T11_2 .- sJZ11_2)' * eQRS[2] * H2    # 
  dB2_12 = Hinv * (T12_2 .- sJZ12_2)' * eQRS[2] * H2    # 
  dB2_13 = Hinv * (T13_2 .- sJZ13_2)' * eQRS[2] * H2    # 
  dB2_21 = Hinv * (T21_2 .- sJZ21_2)' * eQRS[2] * H2    # 
  dB2_22 = Hinv * (T22_2 .- sJZ22_2)' * eQRS[2] * H2    # 
  dB2_23 = Hinv * (T23_2 .- sJZ23_2)' * eQRS[2] * H2    # 
  dB2_31 = Hinv * (T31_2 .- sJZ31_2)' * eQRS[2] * H2    #
  dB2_32 = Hinv * (T32_2 .- sJZ32_2)' * eQRS[2] * H2    # 
  dB2_33 = Hinv * (T33_2 .- sJZ33_2)' * eQRS[2] * H2    # 


  dB3_11 = Hinv * (T11_3 .- sJZ11_3)' * eQRS[3] * H3    # 
  dB3_12 = Hinv * (T12_3 .- sJZ12_3)' * eQRS[3] * H3    # 
  dB3_13 = Hinv * (T13_3 .- sJZ13_3)' * eQRS[3] * H3    # 
  dB3_21 = Hinv * (T21_3 .- sJZ21_3)' * eQRS[3] * H3    # 
  dB3_22 = Hinv * (T22_3 .- sJZ22_3)' * eQRS[3] * H3    # 
  dB3_23 = Hinv * (T23_3 .- sJZ23_3)' * eQRS[3] * H3    # 
  dB3_31 = Hinv * (T31_3 .- sJZ31_3)' * eQRS[3] * H3    #
  dB3_32 = Hinv * (T32_3 .- sJZ32_3)' * eQRS[3] * H3    # 
  dB3_33 = Hinv * (T33_3 .- sJZ33_3)' * eQRS[3] * H3    # 

  dB4_11 = Hinv * (T11_4 .- sJZ11_4)' * eQRS[4] * H4    # 
  dB4_12 = Hinv * (T12_4 .- sJZ12_4)' * eQRS[4] * H4    # 
  dB4_13 = Hinv * (T13_4 .- sJZ13_4)' * eQRS[4] * H4    # 
  dB4_21 = Hinv * (T21_4 .- sJZ21_4)' * eQRS[4] * H4    # 
  dB4_22 = Hinv * (T22_4 .- sJZ22_4)' * eQRS[4] * H4    # 
  dB4_23 = Hinv * (T23_4 .- sJZ23_4)' * eQRS[4] * H4    # 
  dB4_31 = Hinv * (T31_4 .- sJZ31_4)' * eQRS[4] * H4    #
  dB4_32 = Hinv * (T32_4 .- sJZ32_4)' * eQRS[4] * H4    # 
  dB4_33 = Hinv * (T33_4 .- sJZ33_4)' * eQRS[4] * H4    # 

  dB5_11 = Hinv * (T11_5 .- sJZ11_5)' * eQRS[5] * H5    # 
  dB5_12 = Hinv * (T12_5 .- sJZ12_5)' * eQRS[5] * H5    # 
  dB5_13 = Hinv * (T13_5 .- sJZ13_5)' * eQRS[5] * H5    # 
  dB5_21 = Hinv * (T21_5 .- sJZ21_5)' * eQRS[5] * H5    # 
  dB5_22 = Hinv * (T22_5 .- sJZ22_5)' * eQRS[5] * H5    # 
  dB5_23 = Hinv * (T23_5 .- sJZ23_5)' * eQRS[5] * H5    # 
  dB5_31 = Hinv * (T31_5 .- sJZ31_5)' * eQRS[5] * H5    #
  dB5_32 = Hinv * (T32_5 .- sJZ32_5)' * eQRS[5] * H5    # 
  dB5_33 = Hinv * (T33_5 .- sJZ33_5)' * eQRS[5] * H5    # 

  dB6_11 = Hinv * (T11_6 .- sJZ11_6)' * eQRS[6] * H6    # 
  dB6_12 = Hinv * (T12_6 .- sJZ12_6)' * eQRS[6] * H6    # 
  dB6_13 = Hinv * (T13_6 .- sJZ13_6)' * eQRS[6] * H6    # 
  dB6_21 = Hinv * (T21_6 .- sJZ21_6)' * eQRS[6] * H6    # 
  dB6_22 = Hinv * (T22_6 .- sJZ22_6)' * eQRS[6] * H6    # 
  dB6_23 = Hinv * (T23_6 .- sJZ23_6)' * eQRS[6] * H6    # 
  dB6_31 = Hinv * (T31_6 .- sJZ31_6)' * eQRS[6] * H6    #
  dB6_32 = Hinv * (T32_6 .- sJZ32_6)' * eQRS[6] * H6    # 
  dB6_33 = Hinv * (T33_6 .- sJZ33_6)' * eQRS[6] * H6    # 


  dB1 = [[dB1_11] [dB1_12] [dB1_13]; [dB1_21] [dB1_22] [dB1_23]; [dB1_31] [dB1_32] [dB1_33]]
  dB2 = [[dB2_11] [dB2_12] [dB2_13]; [dB2_21] [dB2_22] [dB2_23]; [dB2_31] [dB2_32] [dB2_33]]
  dB3 = [[dB3_11] [dB3_12] [dB3_13]; [dB3_21] [dB3_22] [dB3_23]; [dB3_31] [dB3_32] [dB3_33]]
  dB4 = [[dB4_11] [dB4_12] [dB4_13]; [dB4_21] [dB4_22] [dB4_23]; [dB4_31] [dB4_32] [dB4_33]]
  dB5 = [[dB5_11] [dB5_12] [dB5_13]; [dB5_21] [dB5_22] [dB5_23]; [dB5_31] [dB5_32] [dB5_33]]
  dB6 = [[dB6_11] [dB6_12] [dB6_13]; [dB6_21] [dB6_22] [dB6_23]; [dB6_31] [dB6_32] [dB6_33]]


  dB = (dB1, dB2, dB3, dB4, dB5, dB6)


  # boundary data operators for Neumann faces 
  nB1 = -Hinv * eQRS[1] * H1
  nB2 = -Hinv * eQRS[2] * H2 
  nB3 = -Hinv * eQRS[3] * H3 
  nB4 = -Hinv * eQRS[4] * H4 
  nB5 = -Hinv * eQRS[5] * H5 
  nB6 = -Hinv * eQRS[6] * H6 
  nB = (nB1, nB2, nB3, nB4, nB5, nB6)


  

  # Dirichlet faces: dSAT_{f}_{jl} = dB{f}_{jl} * eQRS[f]'
  # (dB already computed above; reuse to skip recomputing Hinv*(T-sJZ)'*eQRS[f]*Hf)

  dSAT_1_11 = dB1_11 * eQRST[1]
  dSAT_1_12 = dB1_12 * eQRST[1]
  dSAT_1_13 = dB1_13 * eQRST[1]
  dSAT_1_21 = dB1_21 * eQRST[1]
  dSAT_1_22 = dB1_22 * eQRST[1]
  dSAT_1_23 = dB1_23 * eQRST[1]
  dSAT_1_31 = dB1_31 * eQRST[1]
  dSAT_1_32 = dB1_32 * eQRST[1]
  dSAT_1_33 = dB1_33 * eQRST[1]

  dSAT_2_11 = dB2_11 * eQRST[2]
  dSAT_2_12 = dB2_12 * eQRST[2]
  dSAT_2_13 = dB2_13 * eQRST[2]
  dSAT_2_21 = dB2_21 * eQRST[2]
  dSAT_2_22 = dB2_22 * eQRST[2]
  dSAT_2_23 = dB2_23 * eQRST[2]
  dSAT_2_31 = dB2_31 * eQRST[2]
  dSAT_2_32 = dB2_32 * eQRST[2]
  dSAT_2_33 = dB2_33 * eQRST[2]

  dSAT_3_11 = dB3_11 * eQRST[3]
  dSAT_3_12 = dB3_12 * eQRST[3]
  dSAT_3_13 = dB3_13 * eQRST[3]
  dSAT_3_21 = dB3_21 * eQRST[3]
  dSAT_3_22 = dB3_22 * eQRST[3]
  dSAT_3_23 = dB3_23 * eQRST[3]
  dSAT_3_31 = dB3_31 * eQRST[3]
  dSAT_3_32 = dB3_32 * eQRST[3]
  dSAT_3_33 = dB3_33 * eQRST[3]

  dSAT_4_11 = dB4_11 * eQRST[4]
  dSAT_4_12 = dB4_12 * eQRST[4]
  dSAT_4_13 = dB4_13 * eQRST[4]
  dSAT_4_21 = dB4_21 * eQRST[4]
  dSAT_4_22 = dB4_22 * eQRST[4]
  dSAT_4_23 = dB4_23 * eQRST[4]
  dSAT_4_31 = dB4_31 * eQRST[4]
  dSAT_4_32 = dB4_32 * eQRST[4]
  dSAT_4_33 = dB4_33 * eQRST[4]

  dSAT_5_11 = dB5_11 * eQRST[5]
  dSAT_5_12 = dB5_12 * eQRST[5]
  dSAT_5_13 = dB5_13 * eQRST[5]
  dSAT_5_21 = dB5_21 * eQRST[5]
  dSAT_5_22 = dB5_22 * eQRST[5]
  dSAT_5_23 = dB5_23 * eQRST[5]
  dSAT_5_31 = dB5_31 * eQRST[5]
  dSAT_5_32 = dB5_32 * eQRST[5]
  dSAT_5_33 = dB5_33 * eQRST[5]

  dSAT_6_11 = dB6_11 * eQRST[6]
  dSAT_6_12 = dB6_12 * eQRST[6]
  dSAT_6_13 = dB6_13 * eQRST[6]
  dSAT_6_21 = dB6_21 * eQRST[6]
  dSAT_6_22 = dB6_22 * eQRST[6]
  dSAT_6_23 = dB6_23 * eQRST[6]
  dSAT_6_31 = dB6_31 * eQRST[6]
  dSAT_6_32 = dB6_32 * eQRST[6]
  dSAT_6_33 = dB6_33 * eQRST[6]

# 
  dSAT_1 = [[dSAT_1_11] [dSAT_1_12] [dSAT_1_13]; [dSAT_1_21] [dSAT_1_22] [dSAT_1_23];  [dSAT_1_31] [dSAT_1_32] [dSAT_1_33]]
  dSAT_2 = [[dSAT_2_11] [dSAT_2_12] [dSAT_2_13]; [dSAT_2_21] [dSAT_2_22] [dSAT_2_23];  [dSAT_2_31] [dSAT_2_32] [dSAT_2_33]]
  dSAT_3 = [[dSAT_3_11] [dSAT_3_12] [dSAT_3_13]; [dSAT_3_21] [dSAT_3_22] [dSAT_3_23];  [dSAT_3_31] [dSAT_3_32] [dSAT_3_33]]
  dSAT_4 = [[dSAT_4_11] [dSAT_4_12] [dSAT_4_13]; [dSAT_4_21] [dSAT_4_22] [dSAT_4_23];  [dSAT_4_31] [dSAT_4_32] [dSAT_4_33]]
  dSAT_5 = [[dSAT_5_11] [dSAT_5_12] [dSAT_5_13]; [dSAT_5_21] [dSAT_5_22] [dSAT_5_23];  [dSAT_5_31] [dSAT_5_32] [dSAT_5_33]]
  dSAT_6 = [[dSAT_6_11] [dSAT_6_12] [dSAT_6_13]; [dSAT_6_21] [dSAT_6_22] [dSAT_6_23];  [dSAT_6_31] [dSAT_6_32] [dSAT_6_33]]


  dSAT = (dSAT_1, dSAT_2, dSAT_3, dSAT_4, dSAT_5, dSAT_6)

  # }}}

  # {{{ Traction faces:
  Hf = (H1, H2, H3, H4, H5, H6)
  HinvP = [Hinv * eQRS[f] * Hf[f] * eQRS[f]' for f in 1:6]

  nSAT_1_11 = -HinvP[1] * T11_1   # see above. 
  nSAT_1_12 = -HinvP[1] * T12_1   
  nSAT_1_13 = -HinvP[1] * T13_1  
  nSAT_1_21 = -HinvP[1] * T21_1   
  nSAT_1_22 = -HinvP[1] * T22_1     
  nSAT_1_23 = -HinvP[1] * T23_1   
  nSAT_1_31 = -HinvP[1] * T31_1   
  nSAT_1_32 = -HinvP[1] * T32_1     
  nSAT_1_33 = -HinvP[1] * T33_1  

  nSAT_2_11 = -HinvP[2] * T11_2   # see above. 
  nSAT_2_12 = -HinvP[2] * T12_2   
  nSAT_2_13 = -HinvP[2] * T13_2  
  nSAT_2_21 = -HinvP[2] * T21_2   
  nSAT_2_22 = -HinvP[2] * T22_2     
  nSAT_2_23 = -HinvP[2] * T23_2   
  nSAT_2_31 = -HinvP[2] * T31_2   
  nSAT_2_32 = -HinvP[2] * T32_2     
  nSAT_2_33 = -HinvP[2] * T33_2  
  
  nSAT_3_11 = -HinvP[3] * T11_3  # see above. 
  nSAT_3_12 = -HinvP[3] * T12_3  
  nSAT_3_13 = -HinvP[3] * T13_3 
  nSAT_3_21 = -HinvP[3] * T21_3  
  nSAT_3_22 = -HinvP[3] * T22_3    
  nSAT_3_23 = -HinvP[3] * T23_3  
  nSAT_3_31 = -HinvP[3] * T31_3  
  nSAT_3_32 = -HinvP[3] * T32_3    
  nSAT_3_33 = -HinvP[3] * T33_3 
  
  nSAT_4_11 = -HinvP[4] * T11_4
  nSAT_4_12 = -HinvP[4] * T12_4
  nSAT_4_13 = -HinvP[4] * T13_4
  nSAT_4_21 = -HinvP[4] * T21_4
  nSAT_4_22 = -HinvP[4] * T22_4
  nSAT_4_23 = -HinvP[4] * T23_4
  nSAT_4_31 = -HinvP[4] * T31_4
  nSAT_4_32 = -HinvP[4] * T32_4
  nSAT_4_33 = -HinvP[4] * T33_4

  nSAT_5_11 = -HinvP[5] * T11_5
  nSAT_5_12 = -HinvP[5] * T12_5
  nSAT_5_13 = -HinvP[5] * T13_5
  nSAT_5_21 = -HinvP[5] * T21_5
  nSAT_5_22 = -HinvP[5] * T22_5
  nSAT_5_23 = -HinvP[5] * T23_5
  nSAT_5_31 = -HinvP[5] * T31_5
  nSAT_5_32 = -HinvP[5] * T32_5
  nSAT_5_33 = -HinvP[5] * T33_5

  nSAT_6_11 = -HinvP[6] * T11_6
  nSAT_6_12 = -HinvP[6] * T12_6
  nSAT_6_13 = -HinvP[6] * T13_6
  nSAT_6_21 = -HinvP[6] * T21_6
  nSAT_6_22 = -HinvP[6] * T22_6
  nSAT_6_23 = -HinvP[6] * T23_6
  nSAT_6_31 = -HinvP[6] * T31_6
  nSAT_6_32 = -HinvP[6] * T32_6
  nSAT_6_33 = -HinvP[6] * T33_6
  

  nSAT_1 = [[nSAT_1_11] [nSAT_1_12] [nSAT_1_13]; [nSAT_1_21] [nSAT_1_22] [nSAT_1_23]; [nSAT_1_31] [nSAT_1_32] [nSAT_1_33]]
  nSAT_2 = [[nSAT_2_11] [nSAT_2_12] [nSAT_2_13]; [nSAT_2_21] [nSAT_2_22] [nSAT_2_23]; [nSAT_2_31] [nSAT_2_32] [nSAT_2_33]]
  nSAT_3 = [[nSAT_3_11] [nSAT_3_12] [nSAT_3_13]; [nSAT_3_21] [nSAT_3_22] [nSAT_3_23]; [nSAT_3_31] [nSAT_3_32] [nSAT_3_33]]
  nSAT_4 = [[nSAT_4_11] [nSAT_4_12] [nSAT_4_13]; [nSAT_4_21] [nSAT_4_22] [nSAT_4_23]; [nSAT_4_31] [nSAT_4_32] [nSAT_4_33]]
  nSAT_5 = [[nSAT_5_11] [nSAT_5_12] [nSAT_5_13]; [nSAT_5_21] [nSAT_5_22] [nSAT_5_23]; [nSAT_5_31] [nSAT_5_32] [nSAT_5_33]]
  nSAT_6 = [[nSAT_6_11] [nSAT_6_12] [nSAT_6_13]; [nSAT_6_21] [nSAT_6_22] [nSAT_6_23]; [nSAT_6_31] [nSAT_6_32] [nSAT_6_33]]


  nSAT = (nSAT_1, nSAT_2, nSAT_3, nSAT_4,  nSAT_5, nSAT_6)     
  # }}}

  # {{{ Interface Penalties 
  # Z' penalty parameteQRS (for interface Dirichlet conditions) in logical space:

  neighborZ = [similar(sJZ_1), similar(sJZ_2), similar(sJZ_3), similar(sJZ_4), similar(sJZ_4), similar(sJZ_5)]

     
 
  # {{{ Build component operators for element local
  A11 = D11[1, 1] .+ D12[1, 1] .+ D13[1, 1] .+ 
        D21[1, 1] .+ D22[1, 1] .+ D23[1, 1] .+ 
        D31[1, 1] .+ D32[1, 1] .+ D33[1, 1]

    A12 = D11[1, 2] .+ D12[1, 2] .+ D13[1, 2] .+ 
          D21[1, 2] .+ D22[1, 2] .+ D23[1, 2] .+ 
          D31[1, 2] .+ D32[1, 2] .+ D33[1, 2]

    A13 = D11[1, 3] .+ D12[1, 3] .+ D13[1, 3] .+ 
          D21[1, 3] .+ D22[1, 3] .+ D23[1, 3] .+ 
          D31[1, 3] .+ D32[1, 3] .+ D33[1, 3]

    A21 = D11[2, 1] .+ D12[2, 1] .+ D13[2, 1] .+ 
          D21[2, 1] .+ D22[2, 1] .+ D23[2, 1] .+ 
          D31[2, 1] .+ D32[2, 1] .+ D33[2, 1]

    A22 = D11[2, 2] .+ D12[2, 2] .+ D13[2, 2] .+ 
          D21[2, 2] .+ D22[2, 2] .+ D23[2, 2] .+ 
          D31[2, 2] .+ D32[2, 2] .+ D33[2, 2]      

    A23 = D11[2, 3] .+ D12[2, 3] .+ D13[2, 3] .+ 
          D21[2, 3] .+ D22[2, 3] .+ D23[2, 3] .+ 
          D31[2, 3] .+ D32[2, 3] .+ D33[2, 3]

    A31 = D11[3, 1] .+ D12[3, 1] .+ D13[3, 1] .+ 
          D21[3, 1] .+ D22[3, 1] .+ D23[3, 1] .+ 
          D31[3, 1] .+ D32[3, 1] .+ D33[3, 1]

    A32 = D11[3, 2] .+ D12[3, 2] .+ D13[3, 2] .+ 
          D21[3, 2] .+ D22[3, 2] .+ D23[3, 2] .+ 
          D31[3, 2] .+ D32[3, 2] .+ D33[3, 2]

    A33 = D11[3, 3] .+ D12[3, 3] .+ D13[3, 3] .+ 
          D21[3, 3] .+ D22[3, 3] .+ D23[3, 3] .+ 
          D31[3, 3] .+ D32[3, 3] .+ D33[3, 3]




  # }}} 

  # {{{ Build operators to connect interface (modify operator to handle the boundary conditions)
  bctype=(BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE)
  
  for lf = 1:6
    if (LFToB[lf] == BC_NEUMANN || LFToB[lf] == EARTH_SURFACE)
      A11 += nSAT[lf][1, 1]
      A12 += nSAT[lf][1, 2] 
      A13 += nSAT[lf][1, 3]
      A21 += nSAT[lf][2, 1] 
      A22 += nSAT[lf][2, 2]
      A23 += nSAT[lf][2, 3]
      A31 += nSAT[lf][3, 1] 
      A32 += nSAT[lf][3, 2]
      A33 += nSAT[lf][3, 3]

    elseif LFToB[lf] == BC_DIRICHLET
      A11 += dSAT[lf][1, 1]
      A12 += dSAT[lf][2, 1]
      A13 += dSAT[lf][3, 1]
      A21 += dSAT[lf][1, 2]
      A22 += dSAT[lf][2, 2]
      A23 += dSAT[lf][3, 2]
      A31 += dSAT[lf][1, 3]
      A32 += dSAT[lf][2, 3]
      A33 += dSAT[lf][3, 3]
    elseif (LFToB[lf] == BC_LOCKED_INTERFACE) || (LFToB[lf] == BC_JUMP_INTERFACE) || (LFToB[lf] == RS_FAULT) || (LFToB[lf] == VP_FAULT)
      # this happens in global assembly
    else  
      print("invalid bc")
     end
  end
  bctype=(LFToB[1], LFToB[2], LFToB[3], LFToB[4], LFToB[5], LFToB[6])
  # }}}

  # Full operator
  #A = [A11 A12 A13; A21 A22 A23; A31 A32 A33]
  A = (A11, A12, A13, A21, A22, A23, A31, A32, A33)
 

  return (A=A, dB=dB, nB=nB, eQRST=eQRST, H̃ = H̃,  H=H, Hinv = Hinv, T=T, TT=TT, eQRS=eQRS, 
          coord = metrics.coord, 
          J = J,  
          facecoord = metrics.facecoord,
          sJ = metrics.sJ,
          nx = metrics.nx,
          ny = metrics.ny,
          nz = metrics.nz, sJZ, neighborZ,  Nqp = Nqp, Nrp = Nrp, Nsp = Nsp)
end

#{{{ volsourcearray()
function locsourcearray!(ge, source, coord, J, volargs = ())

  (xloc, yloc, zloc) = coord
  Jloc = J 
  @show extrema(J)
  s = source(xloc[:], yloc[:], zloc[:], volargs...)

  ge[:, 1] += Jloc * s[:,1]
  ge[:, 2] += Jloc * s[:,2]
  ge[:, 3] += Jloc * s[:,3]
end

#}}}

function loc_bdry_vec!(ge, lop, LFToB, EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, in_tractionjump, bcargs = ())

  e = bcargs[1]
  
  dB = lop.dB 
  nB = lop.nB 
  

  (xf, yf, zf) = lop.facecoord
  sJ = lop.sJ 
  nx = lop.nx
  ny = lop.ny
  nz = lop.nz
  ge[:, :] .= 0
  for lf = 1:6
    if LFToB[lf] == BC_DIRICHLET
      # fault (Dirichlet):
      gD = bc_Dirichlet(lf, xf[lf], yf[lf], zf[lf], bcargs...) # this should return a three column vector
      vf_1 = dB[lf][1, 1]*gD[:, 1] + dB[lf][2, 1]*gD[:, 2] + dB[lf][3, 1]*gD[:, 3]
      vf_2 = dB[lf][1, 2]*gD[:, 1] + dB[lf][2, 2]*gD[:, 2] + dB[lf][3, 2]*gD[:, 3]
      vf_3 = dB[lf][1, 3]*gD[:, 1] + dB[lf][2, 3]*gD[:, 2] + dB[lf][3, 3]*gD[:, 3]
      

    elseif (LFToB[lf] == BC_NEUMANN || LFToB[lf] == EARTH_SURFACE)
      gN = bc_Neumann(lf, xf[lf], yf[lf], zf[lf], nx[lf], ny[lf], nz[lf], bcargs...) # this should return a three column vector
      vf_1 = nB[lf] * (sJ[lf] .* gN[:, 1])
      vf_2 = nB[lf] * (sJ[lf] .* gN[:, 2])
      vf_3 = nB[lf] * (sJ[lf] .* gN[:, 3])


    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here

    elseif LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
      
      f = EToF[lf, e]  # get global face number

      (em, ep) = FToE[:, f]  # find the two elements that share global face f; ep is on plus side
      (fm, fp) = FToLF[:, f] # Local face numbers. 


      sJZ11 = 0.5*(lop.sJZ[lf][1, 1])' + 0.5*(lop.neighborZ[lf][1, 1])'
      sJZ12 = 0.5*(lop.sJZ[lf][1, 2])' + 0.5*(lop.neighborZ[lf][1, 2])'
      sJZ13 = 0.5*(lop.sJZ[lf][1, 3])' + 0.5*(lop.neighborZ[lf][1, 3])'
      sJZ21 = 0.5*(lop.sJZ[lf][2, 1])' + 0.5*(lop.neighborZ[lf][2, 1])'
      sJZ22 = 0.5*(lop.sJZ[lf][2, 2])' + 0.5*(lop.neighborZ[lf][2, 2])'
      sJZ23 = 0.5*(lop.sJZ[lf][2, 3])' + 0.5*(lop.neighborZ[lf][2, 3])'
      sJZ31 = 0.5*(lop.sJZ[lf][3, 1])' + 0.5*(lop.neighborZ[lf][3, 1])'
      sJZ32 = 0.5*(lop.sJZ[lf][3, 2])' + 0.5*(lop.neighborZ[lf][3, 2])'
      sJZ33 = 0.5*(lop.sJZ[lf][3, 3])' + 0.5*(lop.neighborZ[lf][3, 3])'


   
      # Jump in displacement

      gJ = in_jump(lf, xf[lf], yf[lf], zf[lf], bcargs...)  # this should return a three column vector containing slip components

      A11 = -lop.Hinv * sJZ11 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][1, 1] * lop.eQRS[lf] * lop.H[lf]
      A21 = (-lop.Hinv * sJZ21 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][1, 2] * lop.eQRS[lf] * lop.H[lf])
      A31 = (-lop.Hinv * sJZ31 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][1, 3] * lop.eQRS[lf] * lop.H[lf])

      A12 = (-lop.Hinv * sJZ12 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][2, 1] * lop.eQRS[lf] * lop.H[lf])
      A22 = -lop.Hinv * sJZ22 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][2, 2] * lop.eQRS[lf] * lop.H[lf]
      A32 = (-lop.Hinv * sJZ32 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][2, 3] * lop.eQRS[lf] * lop.H[lf])

      A13 = (-lop.Hinv * sJZ13 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][3, 1] * lop.eQRS[lf] * lop.H[lf])
      A23 = -lop.Hinv * sJZ23 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][3, 2] * lop.eQRS[lf] * lop.H[lf]
      A33 = (-lop.Hinv * sJZ33 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][3, 3] * lop.eQRS[lf] * lop.H[lf])

      # Jump in traction
      gT = in_tractionjump(lf, xf[lf], yf[lf], zf[lf], bcargs...) # this should return a three column vector
    
      B = -0.5 * lop.Hinv * lop.eQRS[lf] * lop.H[lf]
    

      vf_1 = A11*gJ[:,1] + A21*gJ[:,2] + A31*gJ[:,3]  + B*(sJ[lf] .* gT[:, 1]) 
      vf_2 = A12*gJ[:,1] + A22*gJ[:,2] + A32*gJ[:,3] +  B*(sJ[lf] .* gT[:, 2])  
      vf_3 = A13*gJ[:,1] + A23*gJ[:,2] + A33*gJ[:,3] +  B*(sJ[lf] .* gT[:, 3])  

    

    else
      error("invalid bc")
    end

    ge[:, 1] += vf_1 
    ge[:, 2] += vf_2 
    ge[:, 3] += vf_3
  end

end


function simple_bdry_vec!(ge, lop, LFToB, bc_Dirichlet, bc_Neumann)
  
  dB = lop.dB 
  nB = lop.nB 
  

  (xf, yf, zf) = lop.facecoord
  sJ = lop.sJ 
  nx = lop.nx
  ny = lop.ny
  nz = lop.nz
  ge[:, :] .= 0
  for lf = 1:6
    if LFToB[lf] == BC_DIRICHLET
      # fault (Dirichlet):
      gD = bc_Dirichlet(lf, xf[lf], yf[lf], zf[lf]) # this should return a three column vector
      vf_1 = dB[lf][1, 1]*gD[:, 1] + dB[lf][2, 1]*gD[:, 2] + dB[lf][3, 1]*gD[:, 3]
      vf_2 = dB[lf][1, 2]*gD[:, 1] + dB[lf][2, 2]*gD[:, 2] + dB[lf][3, 2]*gD[:, 3]
      vf_3 = dB[lf][1, 3]*gD[:, 1] + dB[lf][2, 3]*gD[:, 2] + dB[lf][3, 3]*gD[:, 3]
      

    elseif (LFToB[lf] == BC_NEUMANN || LFToB[lf] == EARTH_SURFACE)
      gN = bc_Neumann(lf, xf[lf], yf[lf], zf[lf], nx[lf], ny[lf], nz[lf]) # this should return a three column vector

      vf_1 = nB[lf] * (sJ[lf] .* gN[:, 1])
      vf_2 = nB[lf] * (sJ[lf] .* gN[:, 2])
      vf_3 = nB[lf] * (sJ[lf] .* gN[:, 3])


    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here

    elseif LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
      
      f = EToF[lf, e]  # get global face number

      (em, ep) = FToE[:, f]  # find the two elements that share global face f; ep is on plus side
      (fm, fp) = FToLF[:, f] # Local face numbers. 


      sJZ11 = 0.5*(lop.sJZ[lf][1, 1])' + 0.5*(lop.neighborZ[lf][1, 1])'
      sJZ12 = 0.5*(lop.sJZ[lf][1, 2])' + 0.5*(lop.neighborZ[lf][1, 2])'
      sJZ13 = 0.5*(lop.sJZ[lf][1, 3])' + 0.5*(lop.neighborZ[lf][1, 3])'
      sJZ21 = 0.5*(lop.sJZ[lf][2, 1])' + 0.5*(lop.neighborZ[lf][2, 1])'
      sJZ22 = 0.5*(lop.sJZ[lf][2, 2])' + 0.5*(lop.neighborZ[lf][2, 2])'
      sJZ23 = 0.5*(lop.sJZ[lf][2, 3])' + 0.5*(lop.neighborZ[lf][2, 3])'
      sJZ31 = 0.5*(lop.sJZ[lf][3, 1])' + 0.5*(lop.neighborZ[lf][3, 1])'
      sJZ32 = 0.5*(lop.sJZ[lf][3, 2])' + 0.5*(lop.neighborZ[lf][3, 2])'
      sJZ33 = 0.5*(lop.sJZ[lf][3, 3])' + 0.5*(lop.neighborZ[lf][3, 3])'


   
      # Jump in displacement

      gJ = in_jump(lf, xf[lf], yf[lf], zf[lf], bcargs...)  # this should return a three column vector containing slip components

      A11 = -lop.Hinv * sJZ11 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][1, 1] * lop.eQRS[lf] * lop.H[lf]
      A21 = (-lop.Hinv * sJZ21 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][1, 2] * lop.eQRS[lf] * lop.H[lf])
      A31 = (-lop.Hinv * sJZ31 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][1, 3] * lop.eQRS[lf] * lop.H[lf])

      A12 = (-lop.Hinv * sJZ12 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][2, 1] * lop.eQRS[lf] * lop.H[lf])
      A22 = -lop.Hinv * sJZ22 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][2, 2] * lop.eQRS[lf] * lop.H[lf]
      A32 = (-lop.Hinv * sJZ32 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][2, 3] * lop.eQRS[lf] * lop.H[lf])

      A13 = (-lop.Hinv * sJZ13 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][3, 1] * lop.eQRS[lf] * lop.H[lf])
      A23 = -lop.Hinv * sJZ23 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][3, 2] * lop.eQRS[lf] * lop.H[lf]
      A33 = (-lop.Hinv * sJZ33 * lop.eQRS[lf] * lop.H[lf] + 0.5 * lop.Hinv * lop.TT[lf][3, 3] * lop.eQRS[lf] * lop.H[lf])

      # Jump in traction
      gT = in_tractionjump(lf, xf[lf], yf[lf], zf[lf], bcargs...) # this should return a three column vector
    
      B = -0.5 * lop.Hinv * lop.eQRS[lf] * lop.H[lf]
    

      vf_1 = A11*gJ[:,1] + A21*gJ[:,2] + A31*gJ[:,3]  + B*(sJ[lf] .* gT[:, 1]) 
      vf_2 = A12*gJ[:,1] + A22*gJ[:,2] + A32*gJ[:,3] +  B*(sJ[lf] .* gT[:, 2])  
      vf_3 = A13*gJ[:,1] + A23*gJ[:,2] + A33*gJ[:,3] +  B*(sJ[lf] .* gT[:, 3])  

    

    else
      error("invalid bc")
    end

    ge[:, 1] += vf_1 
    ge[:, 2] += vf_2 
    ge[:, 3] += vf_3
  end

end



function loc_bdry_vec_v2!(ge, lop, LFToB, EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, bcargs = ())

  e = bcargs[1]
  dB = lop.dB 
  nB = lop.nB 


  (xf, yf, zf) = lop.facecoord
  sJ = lop.sJ 



  nx = lop.nx
  ny = lop.ny
  nz = lop.nz
  ge[:, :] .= 0
  for lf = 1:6
    if LFToB[lf] == BC_DIRICHLET

      gD = bc_Dirichlet(lf, xf[lf], yf[lf], zf[lf], bcargs...) # this should return B[1]*g
   
      vf_1 = dB[lf][1, 1]*gD[1] + dB[lf][2, 1]*gD[2] + dB[lf][3, 1]*gD[3]
      vf_2 = dB[lf][1, 2]*gD[1] + dB[lf][2, 2]*gD[2] + dB[lf][3, 2]*gD[3]
      vf_3 = dB[lf][1, 3]*gD[1] + dB[lf][2, 3]*gD[2] + dB[lf][3, 3]*gD[3]
    
      ge[:, 1] += vf_1 
      ge[:, 2] += vf_2 
      ge[:, 3] += vf_3 
 
    elseif (LFToB[lf] == BC_NEUMANN || LFToB[lf] == EARTH_SURFACE)
 
      gN = bc_Neumann(lf, xf[lf], yf[lf], zf[lf], nx[lf], ny[lf], nz[lf], bcargs...) #this should return B[3]*g3
   
      vf_1 = nB[lf] * (sJ[lf] .* gN[1])
      vf_2 = nB[lf] * (sJ[lf] .* gN[2])
      vf_3 = nB[lf] * (sJ[lf] .* gN[3])
      
      ge[:, 1] += vf_1 
      ge[:, 2] += vf_2 
      ge[:, 3] += vf_3 

    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here

    elseif LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
     

    sJZ11 = 0.5*(lop.sJZ[lf][1, 1])' + 0.5*(lop.neighborZ[lf][1, 1])'
    sJZ12 = 0.5*(lop.sJZ[lf][1, 2])' + 0.5*(lop.neighborZ[lf][1, 2])'
    sJZ13 = 0.5*(lop.sJZ[lf][1, 3])' + 0.5*(lop.neighborZ[lf][1, 3])'
    sJZ21 = 0.5*(lop.sJZ[lf][2, 1])' + 0.5*(lop.neighborZ[lf][2, 1])'
    sJZ22 = 0.5*(lop.sJZ[lf][2, 2])' + 0.5*(lop.neighborZ[lf][2, 2])'
    sJZ23 = 0.5*(lop.sJZ[lf][2, 3])' + 0.5*(lop.neighborZ[lf][2, 3])'
    sJZ31 = 0.5*(lop.sJZ[lf][3, 1])' + 0.5*(lop.neighborZ[lf][3, 1])'
    sJZ32 = 0.5*(lop.sJZ[lf][3, 2])' + 0.5*(lop.neighborZ[lf][3, 2])'
    sJZ33 = 0.5*(lop.sJZ[lf][3, 3])' + 0.5*(lop.neighborZ[lf][3, 3])'

      # Jump in displacement

      δ = in_jump(lf, xf[lf], yf[lf], zf[lf], bcargs...)  # return a two column vector containing slip components (parallel followed by normal)
  
      du =  -δ[:, 1] .* abs.(ny[lf])   
      dv =   δ[:, 1] .* abs.(ny[lf])  
      dw =   δ[:, 1] .* abs.(nx[lf])  


      A11 = -lop.Hinv * sJZ11 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][1, 1] * lop.eQRS[lf] * lop.H[lf]
      A21 = -lop.Hinv * sJZ21 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][1, 2] * lop.eQRS[lf] * lop.H[lf]
      A31 = -lop.Hinv * sJZ31 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][1, 3] * lop.eQRS[lf] * lop.H[lf]

      A12 = -lop.Hinv * sJZ12 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][2, 1] * lop.eQRS[lf] * lop.H[lf]
      A22 = -lop.Hinv * sJZ22 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][2, 2] * lop.eQRS[lf] * lop.H[lf]
      A32 = -lop.Hinv * sJZ32 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][2, 3] * lop.eQRS[lf] * lop.H[lf]

      A13 = -lop.Hinv * sJZ13 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][3, 1] * lop.eQRS[lf] * lop.H[lf]
      A23 = -lop.Hinv * sJZ23 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][3, 2] * lop.eQRS[lf] * lop.H[lf]
      A33 = -lop.Hinv * sJZ33 * lop.eQRS[lf] * lop.H[lf] + 1 * 0.5 * lop.Hinv * lop.TT[lf][3, 3] * lop.eQRS[lf] * lop.H[lf]

      vf_1 = A11*du + A21*dv + A31*dw 
      vf_2 = A12*du + A22*dv + A32*dw
      vf_3 = A13*du + A23*dv + A33*dw

      ge[:, 1] += vf_1 
      ge[:, 2] += vf_2 
      ge[:, 3] += vf_3 

    else
      error("invalid bc")
    end

   
   
   
  end

end



function computetraction(lop, e1, e2, lf1, lf2, u, w, urng1, urng2, EToF, EToO, FToδstarts, δ; correction = true)


  T11_1 = lop[e1].T_high[lf1][1, 1]
  T12_1 = lop[e1].T_high[lf1][1, 2]
  T13_1 = lop[e1].T_high[lf1][1, 3]
  T21_1 = lop[e1].T_high[lf1][2, 1]
  T22_1 = lop[e1].T_high[lf1][2, 2]
  T23_1 = lop[e1].T_high[lf1][2, 3]
  T31_1 = lop[e1].T_high[lf1][3, 1]
  T32_1 = lop[e1].T_high[lf1][3, 2]
  T33_1 = lop[e1].T_high[lf1][3, 3]

  sJ1 = lop[e1].sJ[lf1]
  eQRST_1 = lop[e1].eQRST[lf1]
     
  t1 = (eQRST_1 * (T11_1 * u[urng1] + T12_1 * v[urng1] + T13_1 * w[urng1])) ./ sJ1
  t2 = (eQRST_1 * (T21_1 * u[urng1] + T22_1 * v[urng1] + T23_1 * w[urng1])) ./ sJ1
  t3 = (eQRST_1 * (T31_1 * u[urng1] + T32_1 * v[urng1] + T33_1 * w[urng1])) ./ sJ1

 
  T11_2 = lop[e2].T_high[lf2][1, 1]
  T12_2 = lop[e2].T_high[lf2][1, 2]
  T13_2 = lop[e2].T_high[lf2][1, 3]
  T21_2 = lop[e2].T_high[lf2][2, 1]
  T22_2 = lop[e2].T_high[lf2][2, 2]
  T23_2 = lop[e2].T_high[lf2][2, 3]
  T31_2 = lop[e2].T_high[lf2][3, 1]
  T32_2 = lop[e2].T_high[lf2][3, 2]
  T33_2 = lop[e2].T_high[lf2][3, 3]

  sJ2 = lop[e2].sJ[lf2]
  eQRST_2 = lop[e2].eQRST[lf2]
  
  t4 = (eQRST_2 * (T11_2 * u[urng2] + T12_2 * v[urng2] + T13_2 * w[urng2])) ./ sJ2
  t5 = (eQRST_2 * (T21_2 * u[urng2] + T22_2 * v[urng2] + T23_2 * w[urng2])) ./ sJ2
  t6 = (eQRST_2 * (T31_2 * u[urng2] + T32_2 * v[urng2] + T33_2 * w[urng2])) ./ sJ2

  # compute actual slip, by default this is with respect to minus side face orientation 
  f = EToF[lf1, e1]  # Get global face number
  δ1 = δ[FToδstarts[f]:(FToδstarts[f+1]-1), 1]  # TODO: CHECK FOR 3D
  δ2 = δ[FToδstarts[f]:(FToδstarts[f+1]-1), 2]  # WRT minus side face orientation
  δ3 = δ[FToδstarts[f]:(FToδstarts[f+1]-1), 3]  # WRT minus side face orientation

  # compute displacement jumps from linear solve with respect to volume orientation

  u1 = lop[e1].eQRST[lf1]*u[urng1]
  v1 = lop[e1].eQRST[lf1]*v[urng1]
  w1 = lop[e1].eQRST[lf1]*w[urng1] 
  u2 = lop[e2].eQRST[lf2]*u[urng2]
  v2 = lop[e2].eQRST[lf2]*v[urng2]
  w2 = lop[e2].eQRST[lf2]*w[urng2]

  du_minus = similar(u1)
  dv_minus = similar(v1)
  dw_minus = similar(w1)
  du_plus = similar(u2)
  dv_plus = similar(v2)
  dw_plus = similar(w2)

  δ1_minus = similar(δ1)
  δ1_plus = similar(δ1)
  δ2_minus = similar(δ2)
  δ2_plus = similar(δ2)
  δ3_minus = similar(δ3)
  δ3_plus = similar(δ3)

  if EToO[lf2, e2] # TODO:  CHECK for 3D
    δ1_minus .= δ1
    δ1_plus .= -δ1
    δ2_minus .= δ2
    δ2_plus .= -δ2
    δ3_minus .= δ3
    δ4_plus .= -δ3 # TODO: check for sign errors
    du_minus .= u2 - u1
    dv_minus .= v2 - v1
    dw_minus .= w2 - w1
    du_plus .=  u1 - u2
    dv_plus .=  v1 - v2
    dw_plus .=  w1 - w2
  else
    δ1_minus .= δ1# TODO: check for sign errors
    δ1_plus .= -δ1[end:-1:1]
    δ2_minus .= δ2# TODO: check for sign errors
    δ2_plus .= -δ2[end:-1:1]
    δ3_minus .= δ3# TODO: check for sign errors
    δ3_plus .= -δ3[end:-1:1]
    du_minus .= u2[end:-1:1] - u1
    dv_minus .= v2[end:-1:1] - v1
    dw_minus .= w2[end:-1:1] - w1
    du_plus .=  u1[end:-1:1] - u2
    dv_plus .=  v1[end:-1:1] - v2
    dw_plus .=  w1[end:-1:1] - w2
  end

      
  #interface penalty terms wrt minus side face orientation
  sJZ_11_minus = (lop[e1].sJZ[lf1][1, 1])' + (lop[e1].neighborZ[lf1][1, 1])' 
  sJZ_12_minus = (lop[e1].sJZ[lf1][1, 2])' + (lop[e1].neighborZ[lf1][1, 2])' 
  sJZ_13_minus = (lop[e1].sJZ[lf1][1, 3])' + (lop[e1].neighborZ[lf1][1, 3])' 
  sJZ_21_minus = (lop[e1].sJZ[lf1][2, 1])' + (lop[e1].neighborZ[lf1][2, 1])' 
  sJZ_22_minus = (lop[e1].sJZ[lf1][2, 2])' + (lop[e1].neighborZ[lf1][2, 2])' 
  sJZ_23_minus = (lop[e1].sJZ[lf1][2, 3])' + (lop[e1].neighborZ[lf1][2, 3])' 
  sJZ_31_minus = (lop[e1].sJZ[lf1][3, 1])' + (lop[e1].neighborZ[lf1][3, 1])' 
  sJZ_32_minus = (lop[e1].sJZ[lf1][3, 2])' + (lop[e1].neighborZ[lf1][3, 2])' 
  sJZ_33_minus = (lop[e1].sJZ[lf1][3, 3])' + (lop[e1].neighborZ[lf1][3, 3])' 

  #interface penalty terms wrt plus side face orientation
  sJZ_11_plus = (lop[e2].sJZ[lf2][1, 1])' + (lop[e2].neighborZ[lf2][1, 1])' 
  sJZ_12_plus = (lop[e2].sJZ[lf2][1, 2])' + (lop[e2].neighborZ[lf2][1, 2])' 
  sJZ_13_plus = (lop[e2].sJZ[lf2][1, 3])' + (lop[e2].neighborZ[lf2][1, 3])' 
  sJZ_21_plus = (lop[e2].sJZ[lf2][2, 1])' + (lop[e2].neighborZ[lf2][2, 1])' 
  sJZ_22_plus = (lop[e2].sJZ[lf2][2, 2])' + (lop[e2].neighborZ[lf2][2, 2])' 
  sJZ_23_plus = (lop[e2].sJZ[lf2][2, 3])' + (lop[e2].neighborZ[lf2][2, 3])' 
  sJZ_31_plus = (lop[e2].sJZ[lf2][3, 1])' + (lop[e2].neighborZ[lf2][3, 1])' 
  sJZ_32_plus = (lop[e2].sJZ[lf2][3, 2])' + (lop[e2].neighborZ[lf2][3, 2])' 
  sJZ_33_plus = (lop[e2].sJZ[lf2][3, 3])' + (lop[e2].neighborZ[lf2][3, 3])' 
  
  du_minus_actual = -δ1_minus .* abs.(lop[e1].ny[lf1])   # horizontal jump
  dv_minus_actual = -δ2_minus .* abs.(lop[e1].ny[lf1])   # horizontal jump # TODO these need to be checked for 3D
  dw_minus_actual =  δ3_minus .* abs.(lop[e1].nx[lf1])  #vertical jump
  du_plus_actual =  -δ1_plus  .* abs.(lop[e2].ny[lf2])   # horizontal jump
  dv_plus_actual =  -δ1_plus  .* abs.(lop[e2].ny[lf2])   # horizontal jump
  dw_plus_actual =   δ1_plus  .* abs.(lop[e2].nx[lf2])  #vertical jump

  # add corrections to traction wrt volume #TODO: check for 3D
  c1 = lop[e1].eQRST[lf1] * (sJZ_11_minus* lop[e1].eQRS[lf1]*(du_minus .- du_minus_actual) .+ sJZ_12_minus* lop[e1].eQRS[lf1]*(dv_minus - dv_minus_actual) .+ sJZ_13_minus* lop[e1].eQRS[lf1]*(dw_minus - dw_minus_actual))
  c2 = lop[e1].eQRST[lf1] * (sJZ_21_minus* lop[e1].eQRS[lf1]*(du_minus .- du_minus_actual) .+ sJZ_22_minus* lop[e1].eQRS[lf1]*(dv_minus - dv_minus_actual) .+ sJZ_23_minus* lop[e1].eQRS[lf1]*(dw_minus - dw_minus_actual))
  c3 = lop[e1].eQRST[lf1] * (sJZ_31_minus* lop[e1].eQRS[lf1]*(du_minus .- du_minus_actual) .+ sJZ_32_minus* lop[e1].eQRS[lf1]*(dv_minus - dv_minus_actual) .+ sJZ_33_minus* lop[e1].eQRS[lf1]*(dw_minus - dw_minus_actual))


  c4 = lop[e2].eQRST[lf2] * (sJZ_11_plus * lop[e2].eQRS[lf2]*(du_plus  .- du_plus_actual) .+  sJZ_12_plus * lop[e2].eQRS[lf2]*(dv_plus - dv_plus_actual) .+  sJZ_13_plus * lop[e2].eQRS[lf2]*(dw_plus - dw_plus_actual))
  c5 = lop[e2].eQRST[lf2] * (sJZ_21_plus * lop[e2].eQRS[lf2]*(du_plus  .- du_plus_actual) .+  sJZ_22_plus * lop[e2].eQRS[lf2]*(dv_plus - dv_plus_actual) .+  sJZ_23_plus * lop[e2].eQRS[lf2]*(dw_plus - dw_plus_actual))
  c6 = lop[e2].eQRST[lf2] * (sJZ_31_plus * lop[e2].eQRS[lf2]*(du_plus  .- du_plus_actual) .+  sJZ_32_plus * lop[e2].eQRS[lf2]*(dv_plus - dv_plus_actual) .+  sJZ_33_plus * lop[e2].eQRS[lf2]*(dw_plus - dw_plus_actual))

 
  # add in correction
  if correction
    t1 .+= c1 ./ sJ1
    t2 .+= c2 ./ sJ1
    t3 .+= c3 ./ sJ1
    t4 .+= c4 ./ sJ2
    t5 .+= c5 ./ sJ2
    t6 .+= c6 ./ sJ2
  end


  return (t1, t2, t3, t4, t5, t6)
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

#{{{ calculate_neighbors: determine neighbor penalty parameters
function calculate_neighbors(lop, e, FToB, EToF, FToE, FToLF, EToO)

  neighborZ = [similar(lop[e].sJZ[1]), similar(lop[e].sJZ[2]), similar(lop[e].sJZ[3]),
               similar(lop[e].sJZ[4]), similar(lop[e].sJZ[5]), similar(lop[e].sJZ[6])]

  LFToB = FToB[EToF[:, e]]
  for lf = 1:6
    if LFToB[lf] == BC_LOCKED_INTERFACE || LFToB[lf] == BC_JUMP_INTERFACE ||
       LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
      f  = EToF[lf, e]
      (em, ep) = FToE[:, f]
      (fm, fp) = FToLF[:, f]
      if em == e
        eo = ep; nf = fp
      else
        eo = em; nf = fm
      end

      # Face grid dimensions for the current element's face lf
      # (faces 1,2: r×s nodes; faces 3,4: q×s nodes; faces 5,6: q×r nodes)
      if lf <= 2
        n1, n2 = lop[e].Nrp, lop[e].Nsp
      elseif lf <= 4
        n1, n2 = lop[e].Nqp, lop[e].Nsp
      else
        n1, n2 = lop[e].Nqp, lop[e].Nrp
      end

      # EToO[fp, ep] is Tuple{Bool,Bool,Bool}: orientation of plus-side face
      # relative to global canonical (which matches the minus-side ordering).
      # (true,true,true) = perfectly aligned; any false = permutation needed.
      ori = EToO[fp, ep]

      for i in 1:3, j in 1:3
        # Neighbor's penalty projected onto its face nodes (Nface × Nface)
        Zf = lop[eo].eQRST[nf] * lop[eo].sJZ[nf][i, j] * lop[eo].eQRS[nf]

        # Reorder neighbor face nodes to match current element's face ordering.
        # ori[1]: axes same (true) or transposed (false)
        # ori[2]: first axis same direction (true) or reversed (false)
        # ori[3]: second axis same direction (true) or reversed (false)
        if ori != (true, true, true)
          Nf = size(Zf, 1)
          perm = zeros(Int, Nf)
          for jj in 1:n2, ii in 1:n1
            k = ii + n1*(jj - 1)
            if ori[1]   # axes not transposed
              ni = ori[2] ? ii : n1 - ii + 1
              nj = ori[3] ? jj : n2 - jj + 1
              perm[k] = ni + n1*(nj - 1)
            else         # axes transposed: (ii,jj) ↦ (jj,ii) in neighbor coords
              ni = ori[2] ? jj : n2 - jj + 1
              nj = ori[3] ? ii : n1 - ii + 1
              perm[k] = ni + n2*(nj - 1)
            end
          end
          Zf = em == e ? Zf[perm, perm] : Zf[invperm(perm), invperm(perm)]
        end

        # Embed into volume: Np×Np, correct for all 6 face types.
        # eQRS[lf] is Np×Nface; eQRST[lf] is Nface×Np.
        neighborZ[lf][i, j] = lop[e].eQRS[lf] * Zf * lop[e].eQRST[lf]
      end
    end
  end

  return neighborZ
end

# }}}

#{{{ global_operator: Build the global operator
function global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nq, Nr, Ns)
  
  nelems = length(lop)
  nfaces = length(FToB)

  VNp = vstarts[nelems+1]-1

  # COO storage for each of the 9 block matrices
  I11 = Int[]; J11 = Int[]; V11 = Float64[]
  I12 = Int[]; J12 = Int[]; V12 = Float64[]
  I13 = Int[]; J13 = Int[]; V13 = Float64[]
  I21 = Int[]; J21 = Int[]; V21 = Float64[]
  I22 = Int[]; J22 = Int[]; V22 = Float64[]
  I23 = Int[]; J23 = Int[]; V23 = Float64[]
  I31 = Int[]; J31 = Int[]; V31 = Float64[]
  I32 = Int[]; J32 = Int[]; V32 = Float64[]
  I33 = Int[]; J33 = Int[]; V33 = Float64[]

  # Append nonzeros from a sparse block into COO vectors, offset by (roff, coff)
  function coo_add!(Iv, Jv, Vv, block, roff, coff)
    II, JJ, VV = findnz(block)
    append!(Iv, II .+ roff)
    append!(Jv, JJ .+ coff)
    append!(Vv, VV)
  end

  # Diagonal blocks: local volume operators
  for e = 1:nelems
    off = vstarts[e] - 1
    coo_add!(I11, J11, V11, lop[e].A[1], off, off)
    coo_add!(I12, J12, V12, lop[e].A[2], off, off)
    coo_add!(I13, J13, V13, lop[e].A[3], off, off)
    coo_add!(I21, J21, V21, lop[e].A[4], off, off)
    coo_add!(I22, J22, V22, lop[e].A[5], off, off)
    coo_add!(I23, J23, V23, lop[e].A[6], off, off)
    coo_add!(I31, J31, V31, lop[e].A[7], off, off)
    coo_add!(I32, J32, V32, lop[e].A[8], off, off)
    coo_add!(I33, J33, V33, lop[e].A[9], off, off)
  end


  for f = 1:nfaces  # Loop over the global faces to put in off-diagonal interface ops.
    
    if (FToB[f] == BC_JUMP_INTERFACE) || (FToB[f] == BC_LOCKED_INTERFACE) || (FToB[f] == RS_FAULT) || (FToB[f] == VP_FAULT)  # check if face is an interface
   
      (em, ep) = FToE[:, f]  # find the two elements that share global face f.
      (fm, fp) = FToLF[:, f]  # find corresponding local faces corresponding to global face f.

      # If face orientation on plus side does not match minus side, then flip
      ori = EToO[fp, ep]
      if fm <= 2
        n1, n2 = Nr[em]+1, Ns[em]+1
      elseif fm <= 4
        n1, n2 = Nq[em]+1, Ns[em]+1
      else
        n1, n2 = Nq[em]+1, Nr[em]+1
      end
      Nf = n1 * n2
      perm = zeros(Int, Nf)
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
      R = sparse(1:Nf, perm, ones(Nf), Nf, Nf)
   

      # ------------------------------------------------------------------
      # Precompute shared suffix products (once per face, not per component)
      #   eHe_m    = eQRS_m * H_m * eQRST_m        (Nf×Nf, for em-em blocks)
      #   eHRe_mp  = eQRS_m * H_m * R * eQRST_p    (Nf×Nf, for em-ep blocks)
      #   HinveHe_m   = Hinv_m * eHe_m              (Nf×Nf, shared T-term in em-em)
      #   HinveHRe_mp = Hinv_m * eHRe_mp            (Nf×Nf, shared T-term in em-ep)
      # and p-side equivalents. Note: ep-em uses R' (R not self-inverse in 3D).
      # ------------------------------------------------------------------
      eHe_m    = lop[em].eQRS[fm] * lop[em].H[fm] * lop[em].eQRST[fm]
      eHRe_mp  = lop[em].eQRS[fm] * lop[em].H[fm] * R * lop[ep].eQRST[fp]
      HinveHe_m   = lop[em].Hinv * eHe_m
      HinveHRe_mp = lop[em].Hinv * eHRe_mp

      eHe_p    = lop[ep].eQRS[fp] * lop[ep].H[fp] * lop[ep].eQRST[fp]
      eHRe_pm  = lop[ep].eQRS[fp] * lop[ep].H[fp] * R' * lop[em].eQRST[fm]
      HinveHe_p   = lop[ep].Hinv * eHe_p
      HinveHRe_pm = lop[ep].Hinv * eHRe_pm

      # sJZ components (with neighborZ, transposed)
      sJZ_11m = 0.5*(lop[em].sJZ[fm][1, 1])' + 0.5*(lop[em].neighborZ[fm][1, 1])' 
      sJZ_12m = 0.5*(lop[em].sJZ[fm][1, 2])' + 0.5*(lop[em].neighborZ[fm][1, 2])' 
      sJZ_13m = 0.5*(lop[em].sJZ[fm][1, 3])' + 0.5*(lop[em].neighborZ[fm][1, 3])' 
      sJZ_21m = 0.5*(lop[em].sJZ[fm][2, 1])' + 0.5*(lop[em].neighborZ[fm][2, 1])' 
      sJZ_22m = 0.5*(lop[em].sJZ[fm][2, 2])' + 0.5*(lop[em].neighborZ[fm][2, 2])' 
      sJZ_23m = 0.5*(lop[em].sJZ[fm][2, 3])' + 0.5*(lop[em].neighborZ[fm][2, 3])' 
      sJZ_31m = 0.5*(lop[em].sJZ[fm][3, 1])' + 0.5*(lop[em].neighborZ[fm][3, 1])' 
      sJZ_32m = 0.5*(lop[em].sJZ[fm][3, 2])' + 0.5*(lop[em].neighborZ[fm][3, 2])' 
      sJZ_33m = 0.5*(lop[em].sJZ[fm][3, 3])' + 0.5*(lop[em].neighborZ[fm][3, 3])' 

      sJZ_11p = 0.5*(lop[ep].sJZ[fp][1, 1])' + 0.5*(lop[ep].neighborZ[fp][1, 1])' 
      sJZ_12p = 0.5*(lop[ep].sJZ[fp][1, 2])' + 0.5*(lop[ep].neighborZ[fp][1, 2])' 
      sJZ_13p = 0.5*(lop[ep].sJZ[fp][1, 3])' + 0.5*(lop[ep].neighborZ[fp][1, 3])' 
      sJZ_21p = 0.5*(lop[ep].sJZ[fp][2, 1])' + 0.5*(lop[ep].neighborZ[fp][2, 1])' 
      sJZ_22p = 0.5*(lop[ep].sJZ[fp][2, 2])' + 0.5*(lop[ep].neighborZ[fp][2, 2])' 
      sJZ_23p = 0.5*(lop[ep].sJZ[fp][2, 3])' + 0.5*(lop[ep].neighborZ[fp][2, 3])' 
      sJZ_31p = 0.5*(lop[ep].sJZ[fp][3, 1])' + 0.5*(lop[ep].neighborZ[fp][3, 1])' 
      sJZ_32p = 0.5*(lop[ep].sJZ[fp][3, 2])' + 0.5*(lop[ep].neighborZ[fp][3, 2])' 
      sJZ_33p = 0.5*(lop[ep].sJZ[fp][3, 3])' + 0.5*(lop[ep].neighborZ[fp][3, 3])' 

      Tm_11 = lop[em].T[fm][1, 1]; Tm_12 = lop[em].T[fm][1, 2]; Tm_13 = lop[em].T[fm][1, 3]
      Tm_21 = lop[em].T[fm][2, 1]; Tm_22 = lop[em].T[fm][2, 2]; Tm_23 = lop[em].T[fm][2, 3]
      Tm_31 = lop[em].T[fm][3, 1]; Tm_32 = lop[em].T[fm][3, 2]; Tm_33 = lop[em].T[fm][3, 3]

      Tp_11 = lop[ep].T[fp][1, 1]; Tp_12 = lop[ep].T[fp][1, 2]; Tp_13 = lop[ep].T[fp][1, 3]
      Tp_21 = lop[ep].T[fp][2, 1]; Tp_22 = lop[ep].T[fp][2, 2]; Tp_23 = lop[ep].T[fp][2, 3]
      Tp_31 = lop[ep].T[fp][3, 1]; Tp_32 = lop[ep].T[fp][3, 2]; Tp_33 = lop[ep].T[fp][3, 3]

      TTm_11 = lop[em].TT[fm][1, 1]; TTm_12 = lop[em].TT[fm][1, 2]; TTm_13 = lop[em].TT[fm][1, 3]
      TTm_21 = lop[em].TT[fm][2, 1]; TTm_22 = lop[em].TT[fm][2, 2]; TTm_23 = lop[em].TT[fm][2, 3]
      TTm_31 = lop[em].TT[fm][3, 1]; TTm_32 = lop[em].TT[fm][3, 2]; TTm_33 = lop[em].TT[fm][3, 3]

      TTp_11 = lop[ep].TT[fp][1, 1]; TTp_12 = lop[ep].TT[fp][1, 2]; TTp_13 = lop[ep].TT[fp][1, 3]
      TTp_21 = lop[ep].TT[fp][2, 1]; TTp_22 = lop[ep].TT[fp][2, 2]; TTp_23 = lop[ep].TT[fp][2, 3]
      TTp_31 = lop[ep].TT[fp][3, 1]; TTp_32 = lop[ep].TT[fp][3, 2]; TTp_33 = lop[ep].TT[fp][3, 3]

      rem = vstarts[em] - 1
      rep = vstarts[ep] - 1

      # ------------------------------------------------------------------
      # Each block now uses the precomputed suffix factors.
      # em-em:  -(Hinv_m * (sJZ - 0.5*TT)) * eHe_m    - 0.5 * HinveHe_m * Tm
      # ep-ep:  -(Hinv_p * (sJZ - 0.5*TT)) * eHe_p    - 0.5 * HinveHe_p * Tp
      # em-ep:  +(Hinv_m * (sJZ - 0.5*TT)) * eHRe_mp  - 0.5 * HinveHRe_mp * Tp
      # ep-em:  +(Hinv_p * (sJZ - 0.5*TT)) * eHRe_pm  - 0.5 * HinveHRe_pm * Tm
      # ------------------------------------------------------------------

      # M11
      HinvSZm_11 = lop[em].Hinv * (sJZ_11m - 0.5 * TTm_11)
      HinvSZp_11 = lop[ep].Hinv * (sJZ_11p - 0.5 * TTp_11)
      coo_add!(I11, J11, V11, -HinvSZm_11 * eHe_m   - 0.5 * HinveHe_m   * Tm_11, rem, rem)
      coo_add!(I11, J11, V11, -HinvSZp_11 * eHe_p   - 0.5 * HinveHe_p   * Tp_11, rep, rep)
      coo_add!(I11, J11, V11,  HinvSZm_11 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_11, rem, rep)
      coo_add!(I11, J11, V11,  HinvSZp_11 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_11, rep, rem)

      # M12
      HinvSZm_12 = lop[em].Hinv * (sJZ_21m - 0.5 * TTm_12)
      HinvSZp_12 = lop[ep].Hinv * (sJZ_21p - 0.5 * TTp_12)
      coo_add!(I12, J12, V12, -HinvSZm_12 * eHe_m   - 0.5 * HinveHe_m   * Tm_12, rem, rem)
      coo_add!(I12, J12, V12, -HinvSZp_12 * eHe_p   - 0.5 * HinveHe_p   * Tp_12, rep, rep)
      coo_add!(I12, J12, V12,  HinvSZm_12 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_12, rem, rep)
      coo_add!(I12, J12, V12,  HinvSZp_12 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_12, rep, rem)

      # M13
      HinvSZm_13 = lop[em].Hinv * (sJZ_31m - 0.5 * TTm_13)
      HinvSZp_13 = lop[ep].Hinv * (sJZ_31p - 0.5 * TTp_13)
      coo_add!(I13, J13, V13, -HinvSZm_13 * eHe_m   - 0.5 * HinveHe_m   * Tm_13, rem, rem)
      coo_add!(I13, J13, V13, -HinvSZp_13 * eHe_p   - 0.5 * HinveHe_p   * Tp_13, rep, rep)
      coo_add!(I13, J13, V13,  HinvSZm_13 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_13, rem, rep)
      coo_add!(I13, J13, V13,  HinvSZp_13 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_13, rep, rem)

      # M21
      HinvSZm_21 = lop[em].Hinv * (sJZ_12m - 0.5 * TTm_21)
      HinvSZp_21 = lop[ep].Hinv * (sJZ_12p - 0.5 * TTp_21)
      coo_add!(I21, J21, V21, -HinvSZm_21 * eHe_m   - 0.5 * HinveHe_m   * Tm_21, rem, rem)
      coo_add!(I21, J21, V21, -HinvSZp_21 * eHe_p   - 0.5 * HinveHe_p   * Tp_21, rep, rep)
      coo_add!(I21, J21, V21,  HinvSZm_21 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_21, rem, rep)
      coo_add!(I21, J21, V21,  HinvSZp_21 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_21, rep, rem)

      # M22
      HinvSZm_22 = lop[em].Hinv * (sJZ_22m - 0.5 * TTm_22)
      HinvSZp_22 = lop[ep].Hinv * (sJZ_22p - 0.5 * TTp_22)
      coo_add!(I22, J22, V22, -HinvSZm_22 * eHe_m   - 0.5 * HinveHe_m   * Tm_22, rem, rem)
      coo_add!(I22, J22, V22, -HinvSZp_22 * eHe_p   - 0.5 * HinveHe_p   * Tp_22, rep, rep)
      coo_add!(I22, J22, V22,  HinvSZm_22 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_22, rem, rep)
      coo_add!(I22, J22, V22,  HinvSZp_22 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_22, rep, rem)

      # M23  (note: em-em uses sJZ_32m/TTm_23; em-ep uses sJZ_23m/TTm_32 — preserved from original)
      HinvSZm_23  = lop[em].Hinv * (sJZ_32m - 0.5 * TTm_23)
      HinvSZp_23  = lop[ep].Hinv * (sJZ_32p - 0.5 * TTp_23)
      coo_add!(I23, J23, V23, -HinvSZm_23 * eHe_m   - 0.5 * HinveHe_m   * Tm_23, rem, rem)
      coo_add!(I23, J23, V23, -HinvSZp_23 * eHe_p   - 0.5 * HinveHe_p   * Tp_23, rep, rep)
      coo_add!(I23, J23, V23,  HinvSZm_23 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_23, rem, rep)
      coo_add!(I23, J23, V23,  HinvSZp_23 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_23, rep, rem)

      # M31
      HinvSZm_31 = lop[em].Hinv * (sJZ_13m - 0.5 * TTm_31)
      HinvSZp_31 = lop[ep].Hinv * (sJZ_13p - 0.5 * TTp_31)
      coo_add!(I31, J31, V31, -HinvSZm_31 * eHe_m   - 0.5 * HinveHe_m   * Tm_31, rem, rem)
      coo_add!(I31, J31, V31, -HinvSZp_31 * eHe_p   - 0.5 * HinveHe_p   * Tp_31, rep, rep)
      coo_add!(I31, J31, V31,  HinvSZm_31 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_31, rem, rep)
      coo_add!(I31, J31, V31,  HinvSZp_31 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_31, rep, rem)

      # M32  (em-em and em-ep both use sJZ_23m/TTm_32)
      HinvSZm_32 = lop[em].Hinv * (sJZ_23m - 0.5 * TTm_32)
      HinvSZp_32 = lop[ep].Hinv * (sJZ_23p - 0.5 * TTp_32)
      coo_add!(I32, J32, V32, -HinvSZm_32 * eHe_m   - 0.5 * HinveHe_m   * Tm_32, rem, rem)
      coo_add!(I32, J32, V32, -HinvSZp_32 * eHe_p   - 0.5 * HinveHe_p   * Tp_32, rep, rep)
      coo_add!(I32, J32, V32,  HinvSZm_32 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_32, rem, rep)
      coo_add!(I32, J32, V32,  HinvSZp_32 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_32, rep, rem)

      # M33
      HinvSZm_33 = lop[em].Hinv * (sJZ_33m - 0.5 * TTm_33)
      HinvSZp_33 = lop[ep].Hinv * (sJZ_33p - 0.5 * TTp_33)
      coo_add!(I33, J33, V33, -HinvSZm_33 * eHe_m   - 0.5 * HinveHe_m   * Tm_33, rem, rem)
      coo_add!(I33, J33, V33, -HinvSZp_33 * eHe_p   - 0.5 * HinveHe_p   * Tp_33, rep, rep)
      coo_add!(I33, J33, V33,  HinvSZm_33 * eHRe_mp - 0.5 * HinveHRe_mp * Tp_33, rem, rep)
      coo_add!(I33, J33, V33,  HinvSZp_33 * eHRe_pm - 0.5 * HinveHRe_pm * Tm_33, rep, rem)

    end
  end

 
  M11 = sparse(I11, J11, V11, VNp, VNp)
  M12 = sparse(I12, J12, V12, VNp, VNp)
  M13 = sparse(I13, J13, V13, VNp, VNp)
  M21 = sparse(I21, J21, V21, VNp, VNp)
  M22 = sparse(I22, J22, V22, VNp, VNp)
  M23 = sparse(I23, J23, V23, VNp, VNp)
  M31 = sparse(I31, J31, V31, VNp, VNp)
  M32 = sparse(I32, J32, V32, VNp, VNp)
  M33 = sparse(I33, J33, V33, VNp, VNp)

  M = [M11 M12 M13; M21 M22 M23; M31 M32 M33]
  return M
end
#}}}


function bcstarts(FToB, FToE, FToLF, bctype, Nq, Nr, Ns)
    # compute an array of indices to loop over vertices on element faces
    nfaces = length(FToB)   # get total number of faces
    bcstarts = Array{Int64, 1}(undef, nfaces + 1) # initialize array
    bcstarts[1] = 1                               # array starts at 1
    for f = 1:nfaces                              # loop over the global face numbeQRS
      if FToB[f] ∈ bctype                         # determine if face is of type bctype, e.g. BC_JUMP_INTERFACE
        e  = FToE[1,f]                            # if so, set e = first element shared by face (plus side)
        lf = FToLF[1,f]                           # get the local face number for first element
        
        if lf ∈ (1,2) 
          bcs = (Nr[e]+1)*(Ns[e]+1)
        elseif lf ∈ (3,4) 
          bcs = (Nq[e]+1)*(Ns[e]+1)
        else # i.e. lf ∈ (5,6)
          bcs = (Nq[e]+1)*(Nr[e]+1)
        end
        bcstarts[f+1] = bcstarts[f] + bcs

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
  EToO = Array{Tuple{Bool, Bool, Bool},2}(undef, 6,  nelems)
  EToS  = zeros(Int64, 6, nelems) 

  # Local Face to Local Vertex map

  LFToLV = flatten_tuples(((1,3,5,7), (2,4,6,8), (1,2,5,6), (3,4,7,8), (1,2,3,4), (5,6,7,8))) 
  for e = 1:nelems
    for lf = 1:6
      gf = EToF[lf, e]
      if FToE[1, gf] == 0
        @assert FToLF[1, gf] == 0
        FToE[1, gf] = e
        FToLF[1, gf] = lf
        EToO[lf, e] = (1,1,1) # default corresponding to (normal vectors oriented the same, vector for fast direction is pointing in same direction, vector in secondary direction pointing in same direction)
        EToS[lf, e] = 1
      else
        @assert FToE[2, gf] == 0
        @assert FToLF[2, gf] == 0
        FToE[2, gf] = e
        FToLF[2, gf] = lf
        EToS[lf, e] = 2

        ne = FToE[1, gf] # first (other) element to share this global face
        nf = FToLF[1, gf] # second element's local face corresponding to gf. 

        nv = EToV[LFToLV[:,nf], ne] # global vertices of first element sharing global face gf. 
        lv = EToV[LFToLV[:,lf], e] # global vertices of second element sharing global face gf.
        
    
        if nv == lv
          EToO[lf, e] = (1, 1, 1)
        elseif nv[[1;3;2;4]] == lv
          EToO[lf, e] = (0, 1, 1)

        elseif nv[[2;1;4;3]] == lv
          EToO[lf, e] = (1, 0, 1)
        elseif nv[[2;4;1;3]] == lv
          EToO[lf, e] = (0, 0, 1)

        elseif nv[[3;1;4;2]] == lv
          EToO[lf, e] = (0, 1, 0)
        elseif nv[[3;4;1;2]] == lv
          EToO[lf, e] = (1, 1, 0)

         elseif nv[[4;2;3;1]] == lv
          EToO[lf, e] = (0, 0, 0)
        elseif nv[[4;3;2;1]] == lv
          EToO[lf, e] = (1, 0, 0)
        
        else
          error("problem with connectivity")
        end
      end
    end
  end
  (FToE, FToLF, EToO, EToS)
end
#}}}



function compute_surface_displacements!(u1, u2, u3, u, v, w, nfaces, FToB, FToE, FToLF, vstarts, lop, FToFREEstarts)
    for f = 1:nfaces  # loop over faces
        if FToB[f] ∈ EARTH_SURFACE # determine if a free surface
            (e, _) = FToE[:, f] # get element number
            (lf, _) = FToLF[:, f] # get local face number
    
                u1e = u[vstarts[e]:vstarts[e+1]-1, :] # restrict displacement to element
                u2e = v[vstarts[e]:vstarts[e+1]-1, :] # restrict displacement to element
                u3e = w[vstarts[e]:vstarts[e+1]-1, :] # restrict displacement to element

                u1_surf = lop[e].eQRST[lf]*u1e   # restrict element displacement to surface
                u2_surf = lop[e].eQRST[lf]*u2e   # restrict element displacement to surface
                u3_surf = lop[e].eQRST[lf]*u3e   # restrict element displacement to surface

                Frng = FToFREEstarts[f]:(FToFREEstarts[f+1]-1)
              
                for n = 1:length(Frng)
                    u1[Frng[n]] = u1_surf[n]
                    u2[Frng[n]] = u2_surf[n]
                    u3[Frng[n]] = u3_surf[n]
                end
          
        end
    end
  end


# flatten tuples to arrays
flatten_tuples(x) = reshape(collect(Iterators.flatten(x)), length(x[1]),
                            length(x))

export create_metrics, get_operators, computetraction,   plot_blocks, plot_connectivity
export bdry_vec!, loc_bdry_vec_v2!, compute_surface_displacements!
export rateandstate, newtbndv, bcstarts, connectivityarrays, flatten_tuples
export LocalGlobalOperators, SBPLocalOperator1, global_operator
