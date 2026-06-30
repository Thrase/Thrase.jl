using SparseArrays
using LinearAlgebra


⊗(A,B) = kron(A, B)

const BC_DIRICHLET        = 1
const BC_NEUMANN          = 2
const BC_LOCKED_INTERFACE = 0
const BC_JUMP_INTERFACE   = 7

include("transfinite.jl")

function create_metrics_3D(pm, Nq, Nr, Ns, λ, μ, K, B_p,
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
  @show extrema(J)

  @assert minimum(J) >= 1e-12

  qx =  (yr .* zs - ys .* zr) ./ J 
  qy = -(xr .* zs - xs .* zr) ./ J 
  qz =  (xr .* ys - xs .* yr) ./ J 
  rx =  (yq .* zs - ys .* zq) ./ J
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
                C[i, j, k, l] = J .* (λ(x, y, z, B_p) .* F(j, i, qrs) .* F(l, k, qrs) + 
                                   μ(x, y, z, B_p) .* (F(1, i, qrs) .* delt(j, l) .* F(1, k, qrs) + F(2, i, qrs) .* delt(j, l) .* F(2, k, qrs) + F(3, i, qrs) .* delt(j, l) .* F(3, k, qrs)) +
                                   μ(x, y, z, B_p) .* (F(l, i, qrs) .* F(j, k, qrs)))
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


  (coord = (x,y, z),
   facecoord = ((xf1, xf2, xf3, xf4, xf5, xf6), (yf1, yf2, yf3, yf4, yf5, yf6), (zf1, zf2, zf3, zf4, zf5, zf6)),
   C = C,
   J=J,
   sJ = (sJ1, sJ2, sJ3, sJ4, sJ5, sJ6),
   nx = (nx1, nx2, nx3, nx4, nx5, nx6),
   ny = (ny1, ny2, ny3, ny4, ny5, ny6),
   nz = (nz1, nz2, nz3, nz4, nz5, nz6),
   qx = qx, qy = qy, qz = qz, 
   rx = rx, ry = ry, rz = rz,
   sx = sx, sy = sy, sz = sz)
end

function transforms_e(Lw, r̂, l)
    
    
    A = (Lw - Lw*r̂ - Lw)/(2*tanh((r̂-1)/l) + tanh(-2/l)*(r̂ - 1))
    b = (A*tanh(-2/l) + Lw)/2
    c = Lw - b
    xt = (q, r, s) -> (A .* tanh.((q .- 1) ./ l) .+ b .* q .+ c,
                   ((A .* sech.((q .- 1) ./ l).^2) ./ l) .+ b,
                   zeros(size(r)), zeros(size(r)))
    yt = (q, r, s) -> (A .* tanh.((r .- 1) ./ l) .+ b.*r .+ c,
                   zeros(size(q)),
                   ((A .* sech.((r .- 1) ./ l).^2) ./ l) .+ b, zeros(size(s)))
    zt = (q, r, s) -> (A .* tanh.((s .- 1) ./ l) .+ b.*s .+ c,
                   zeros(size(q)), zeros(size(r)),
                   ((A .* sech.((s .- 1) ./ l).^2) ./ l) .+ b)
    
    
    return xt, yt, zt
    
end

function constantgrid(x1, x2, y1, y2, z1, z2)
     xt = (q, r, s) -> ((x2-x1)*(q+1)/2 .+ x1 ,
                   (x2-x1)/2,
                   zeros(size(r)), zeros(size(s)))
    
    yt = (q, r, s) -> ((y2-y1)*(r+1)/2 .+ y1,
                   zeros(size(q)),
                   (x2-x1)/2, 
                   zeros(size(s)))

    zt = (q, r, s) -> ((z2-z1)*(s+1)/2 .+ z1,
                   zeros(size(q)),
                   zeros(size(r)),
                   (x2-x1)/2)


    return xt, yt, zt
end

function transforms_ne(Lw, el_x, el_y, el_z)
    
    xt = (q, r, s) -> (el_x .* tan.(atan((Lw)/el_x).* (0.5*q .+ 0.5)),
                   el_x .* sec.(atan((Lw)/el_x).* (0.5*q .+ 0.5)).^2 * atan((Lw)/el_x) * 0.5 ,
                   zeros(size(r)), zeros(size(s)))
    
    yt = (q, r, s) -> (el_y .* tan.(atan((Lw)/el_y).* (0.5*r .+ 0.5)),
                   zeros(size(q)),
                   el_y .* sec.(atan((Lw)/el_y).*(0.5*r .+ 0.5)) .^2 * atan((Lw)/el_y) * 0.5, 
                   zeros(size(s)))

    zt = (q, r, s) -> (el_z .* tan.(atan((Lw)/el_z).* (0.5*s .+ 0.5)),
                   zeros(size(q)),
                   zeros(size(r)),
                   el_z .* sec.(atan((Lw)/el_z).*(0.5*s .+ 0.5)) .^2 * atan((Lw)/el_z) * 0.5)


    return xt, yt, zt
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



function locoperator_3D(p, Nq, Nr, Ns, metrics=create_metrics_3D(pm, Nq, Nr, Ns, λ, μ, K, B_p), LFToB = (BC_DIRICHLET, BC_DIRICHLET,
                     BC_DIRICHLET, BC_DIRICHLET, BC_DIRICHLET, BC_DIRICHLET);
                     C = metrics.C)
    
    # {{{ read in parameters  
    C = metrics.C   
    J = metrics.J
                
    Nqp = Nq + 1
    Nrp = Nr + 1
    Nsp = Ns + 1
    Np = Nqp * Nrp * Nsp
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
    # }}}

    # {{{ 1D Identity matrices
    Iq = sparse(I, Nqp, Nqp)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)
    # }}}

    # {{{ Set up the operators in logical space
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


    for i = 1:3
        for j = 1:3
            #D11[i, j] = c[1, i, 1, j] * JI * (Is ⊗ Ir ⊗ D2q)
           (D11[i, j], _, _) = var_3D_D2q(p, Nqp, Nrp, Nsp, metrics.C[1, i, 1, j], HqI; xc = (-1, 1))
            D11[i, j] = JI * D11[i, j]
            D12[i, j] = c[1, i, 2, j] * JI * (Is ⊗ Dr ⊗ Dq)
            D13[i, j] = c[1, i, 3, j] * JI * (Ds ⊗ Ir ⊗ Dq)

            D21[i, j] = c[2, i, 1, j] * JI * (Is ⊗ Dr ⊗ Dq)
            #D22[i, j] = c[2, i, 2, j] * JI * (Is ⊗ D2r ⊗ Iq)
           (D22[i, j], _, _) = var_3D_D2r(p, Nqp, Nrp, Nsp, metrics.C[2, i, 2, j], HrI; xc = (-1, 1))
            D22[i, j] = JI * D22[i, j]
            D23[i, j] = c[2, i, 3, j] * JI * (Ds ⊗ Dr ⊗ Iq)

            D31[i, j] = c[3, i, 1, j] * JI * (Ds ⊗ Ir ⊗ Dq)
            D32[i, j] = c[3, i, 2, j] * JI * (Ds ⊗ Dr ⊗ Iq)
            #D33[i, j] = c[3, i, 3, j] * JI * (D2s ⊗ Ir ⊗ Iq)
           (D33[i, j], _, _) = var_3D_D2s(p, Nqp, Nrp, Nsp, metrics.C[3, i, 3, j], HsI; xc = (-1, 1))
            D33[i, j] = JI * D33[i, j]

        end
    end


    M11 = D11[1, 1] .+ D12[1, 1] .+ D13[1, 1] .+ 
          D21[1, 1] .+ D22[1, 1] .+ D23[1, 1] .+ 
          D31[1, 1] .+ D32[1, 1] .+ D33[1, 1]

    M12 = D11[1, 2] .+ D12[1, 2] .+ D13[1, 2] .+ 
          D21[1, 2] .+ D22[1, 2] .+ D23[1, 2] .+ 
          D31[1, 2] .+ D32[1, 2] .+ D33[1, 2]

    M13 = D11[1, 3] .+ D12[1, 3] .+ D13[1, 3] .+ 
          D21[1, 3] .+ D22[1, 3] .+ D23[1, 3] .+ 
          D31[1, 3] .+ D32[1, 3] .+ D33[1, 3]

    M21 = D11[2, 1] .+ D12[2, 1] .+ D13[2, 1] .+ 
          D21[2, 1] .+ D22[2, 1] .+ D23[2, 1] .+ 
          D31[2, 1] .+ D32[2, 1] .+ D33[2, 1]

    M22 = D11[2, 2] .+ D12[2, 2] .+ D13[2, 2] .+ 
          D21[2, 2] .+ D22[2, 2] .+ D23[2, 2] .+ 
          D31[2, 2] .+ D32[2, 2] .+ D33[2, 2]      

    M23 = D11[2, 3] .+ D12[2, 3] .+ D13[2, 3] .+ 
          D21[2, 3] .+ D22[2, 3] .+ D23[2, 3] .+ 
          D31[2, 3] .+ D32[2, 3] .+ D33[2, 3]

    M31 = D11[3, 1] .+ D12[3, 1] .+ D13[3, 1] .+ 
          D21[3, 1] .+ D22[3, 1] .+ D23[3, 1] .+ 
          D31[3, 1] .+ D32[3, 1] .+ D33[3, 1]

    M32 = D11[3, 2] .+ D12[3, 2] .+ D13[3, 2] .+ 
          D21[3, 2] .+ D22[3, 2] .+ D23[3, 2] .+ 
          D31[3, 2] .+ D32[3, 2] .+ D33[3, 2]

    M33 = D11[3, 3] .+ D12[3, 3] .+ D13[3, 3] .+ 
          D21[3, 3] .+ D22[3, 3] .+ D23[3, 3] .+ 
          D31[3, 3] .+ D32[3, 3] .+ D33[3, 3]
    # }}}
    

    # {{{ Boundary point matrices and surface mass matrices 
    
    H1 = H2 = Hs ⊗ Hr
    H1I = H2I = HsI ⊗ HrI

    H3 = H4 = Hs ⊗ Hq
    H3I = H4I = HsI ⊗ HqI

    H5 = H6 = Hr ⊗ Hq
    H5I = H6I = HrI ⊗ HqI

    # Volume matrices
    H = Hs ⊗ Hr ⊗ Hq
    HI = HsI ⊗ HrI ⊗ HqI

    # Create Face Restriction Operators (e1T, e2T, etc) 
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
    e1T = Is ⊗ Ir ⊗ eq0'
    e2T = Is ⊗ Ir ⊗ eqN'
    e3T = Is ⊗ er0' ⊗ Iq
    e4T = Is ⊗ erN' ⊗ Iq
    e5T = es0' ⊗ Ir ⊗ Iq
    e6T = esN' ⊗ Ir ⊗ Iq


    # }}}

    
    # {{{ Jacobian and Surface Jacobian operators
    JHI = HI * JI
    JIm = spdiagm(0 => JI[:])

     (sJ1, sJ2, sJ3, sJ4, sJ5, sJ6) = metrics.sJ
   

    nx = metrics.nx
    ny = metrics.ny
    nz = metrics.nz
    
    # Set up matrices holding normal vectors at faces:
    # FACE 1
    Nx1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx1 .= 0
    Nx1[1, :, :] = metrics.nx[1]
    Nx1 = spdiagm(0 => Nx1[:])
    
    Ny1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny1 .= 0
    Ny1[1, :, :] = metrics.ny[1]
    Ny1 = spdiagm(0 => Ny1[:])
    Nz1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz1 .= 0
    Nz1[1, :, :] = metrics.nz[1]
    Nz1 = spdiagm(0 => Nz1[:])

  
    # FACE 2
    Nx2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx2 .= 0
    Nx2[Nqp, :, :] = metrics.nx[2]
    Nx2 = spdiagm(0 => Nx2[:])
    Ny2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny2 .= 0
    Ny2[Nqp, :, :] = metrics.ny[2]
    Ny2 = spdiagm(0 => Ny2[:])
    Nz2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz2 .= 0
    Nz2[Nqp, :, :] = metrics.nz[2]
    Nz2 = spdiagm(0 => Nz2[:])

    # FACE 3
    Nx3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx3 .= 0
    Nx3[:, 1, :] = metrics.nx[3]
    Nx3 = spdiagm(0 => Nx3[:])
    Ny3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny3 .= 0
    Ny3[:, 1, :] = metrics.ny[3]
    Ny3 = spdiagm(0 => Ny3[:])
    Nz3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz3 .= 0
    Nz3[:, 1, :] = metrics.nz[3]
    Nz3 = spdiagm(0 => Nz3[:])

    # FACE 3
    Nx3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx3 .= 0
    Nx3[:, 1, :] = metrics.nx[3]
    Nx3 = spdiagm(0 => Nx3[:])
    Ny3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny3 .= 0
    Ny3[:, 1, :] = metrics.ny[3]
    Ny3 = spdiagm(0 => Ny3[:])
    Nz3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz3 .= 0
    Nz3[:, 1, :] = metrics.nz[3]
    Nz3 = spdiagm(0 => Nz3[:])

    # FACE 4
    Nx4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx4 .= 0
    Nx4[:, Nrp, :] = metrics.nx[4]
    Nx4 = spdiagm(0 => Nx4[:])
    Ny4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny4 .= 0
    Ny4[:, Nrp, :] = metrics.ny[4]
    Ny4 = spdiagm(0 => Ny4[:])
    Nz4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz4 .= 0
    Nz4[:, Nrp, :] = metrics.nz[4]
    Nz4 = spdiagm(0 => Nz4[:])
   
    # FACE 5
    Nx5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx5 .= 0
    Nx5[:, :, 1] = metrics.nx[5]
    Nx5 = spdiagm(0 => Nx5[:])
    Ny5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny5 .= 0
    Ny5[:, :, 1] = metrics.ny[5]
    Ny5 = spdiagm(0 => Ny5[:])
    Nz5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz5 .= 0
    Nz5[:, :, 1] = metrics.nz[5]
    Nz5 = spdiagm(0 => Nz5[:])
  
    # FACE 6
    Nx6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx6 .= 0
    Nx6[:, :, Nsp] = metrics.nx[6]
    Nx6 = spdiagm(0 => Nx6[:])
    Ny6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny6 .= 0
    Ny6[:, :, Nsp] = metrics.ny[6]
    Ny6 = spdiagm(0 => Ny6[:])
    Nz6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz6 .= 0
    Nz6[:, :, Nsp] = metrics.nz[6]
    Nz6 = spdiagm(0 => Nz6[:])
    
    
    # define Jacobian matrix operators sJ/J evaluated on faces:
    EsJ1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ1 .= 0
    EsJ1[1, :, :] = sJ1 ./ J[1, :, :]
    EsJ2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ2 .= 0
    EsJ2[end, :, :] = sJ2 ./ J[end, :, :]

    EsJ3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ3 .= 0
    EsJ3[:, 1, :] = sJ3 ./ J[:, 1, :]
    EsJ4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ4 .= 0
    EsJ4[:, end, :] = sJ4 ./ J[:, end, :]

    EsJ5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ5 .= 0
    EsJ5[:, :, 1] = sJ5 ./ J[:, :, 1]
    EsJ6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ6 .= 0
    EsJ6[:, :, end] = sJ6 ./ J[:, :, end]


    # Compute inverse surface jacobian operators:
    sJI1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI1 .= 0
    sJI1[1, :, :] = 1 ./ sJ1
    sJI1 =   spdiagm(0 =>  sJI1[:])

    sJI2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI2 .= 0
    sJI2[end, :, :] = 1 ./ sJ2
    sJI2 =   spdiagm(0 =>  sJI2[:])

    sJI3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI3 .= 0
    sJI3[:, 1, :] = 1 ./ sJ3
    sJI3 =   spdiagm(0 =>  sJI3[:])

    sJI4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI4 .= 0
    sJI4[:, end, :] = 1 ./ sJ4
    sJI4 =   spdiagm(0 =>  sJI4[:])

    sJI5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI5 .= 0
    sJI5[:, :, 1] = 1 ./ sJ5
    sJI5 =   spdiagm(0 =>  sJI5[:])

    sJI6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI6 .= 0
    sJI6[:, :, end] = 1 ./ sJ6
    sJI6 =   spdiagm(0 =>  sJI6[:])


    # Turn J, sJ's EsJ's into diagonal matrices
    JI =  spdiagm(0 => 1 ./ J[:])
    J =   spdiagm(0 => J[:])
    sJ1 = spdiagm(0 => sJ1[:])
    sJ2 = spdiagm(0 => sJ2[:])
    sJ3 = spdiagm(0 => sJ3[:])
    sJ4 = spdiagm(0 => sJ4[:])
    sJ5 = spdiagm(0 => sJ5[:])
    sJ6 = spdiagm(0 => sJ6[:])

    EsJ1 = spdiagm(0 => EsJ1[:])
    EsJ2 = spdiagm(0 => EsJ2[:])
    EsJ3 = spdiagm(0 => EsJ3[:])
    EsJ4 = spdiagm(0 => EsJ4[:])
    EsJ5 = spdiagm(0 => EsJ5[:])
    EsJ6 = spdiagm(0 => EsJ6[:])
    # }}}



    # {{{ discrete boundary traction operators in logical space:
    # Create 3D ops from 1D
    Dq3 = Is ⊗ Ir ⊗ Dq
    Dr3 = Is ⊗ Dr ⊗ Iq
    Ds3 = Ds ⊗ Ir ⊗ Iq

    # Face operators to reduce computations involing T11
    # S0q = Is ⊗ Ir ⊗ S0q
    # SNq = Is ⊗ Ir ⊗ SNq
    # S0r = Is ⊗ S0r ⊗ Iq
    # SNr = Is ⊗ SNr ⊗ Iq
    # S0s = S0s ⊗ Ir ⊗ Iq
    # SNs = SNs ⊗ Ir ⊗ Iq

    Sq = Is ⊗ Ir ⊗ (SNq+S0q)
    Sr = Is ⊗ (SNr+S0r) ⊗ Iq
    Ss = (SNs+S0s) ⊗ Ir ⊗ Iq
    # SNq = Is ⊗ Ir ⊗ SNq
    # S0r = Is ⊗ S0r ⊗ Iq
    # SNr = Is ⊗ SNr ⊗ Iq
    # S0s = S0s ⊗ Ir ⊗ Iq
    # SNs = SNs ⊗ Ir ⊗ Iq


   
  
    #{{{ Turn C into a diagonal matrix
    c = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3, 3, 3)

    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    c[i, j, k, l] = spdiagm(0 => C[i, j, k, l][:])
                end
            end
        end
    end
    
  
    T11_1 = (-sJI1) * (c[1,1,1,1]*Sq + c[1,1,2,1]*Dr3 + c[1,1,3,1]*Ds3)
    T12_1 = (-sJI1) * (c[1,1,1,2]*Sq + c[1,1,2,2]*Dr3 + c[1,1,3,2]*Ds3) 
    T21_1 = (-sJI1) * (c[1,2,1,1]*Sq + c[1,2,2,1]*Dr3 + c[1,2,3,1]*Ds3)
    T13_1 = (-sJI1) * (c[1,1,1,3]*Sq + c[1,1,2,3]*Dr3 + c[1,1,3,3]*Ds3) 
    T31_1 = (-sJI1) * (c[1,3,1,1]*Sq + c[1,3,2,1]*Dr3 + c[1,3,3,1]*Ds3) 
    T22_1 = (-sJI1) * (c[1,2,1,2]*Sq + c[1,2,2,2]*Dr3 + c[1,2,3,2]*Ds3) 
    T23_1 = (-sJI1) * (c[1,2,1,3]*Sq + c[1,2,2,3]*Dr3 + c[1,2,3,3]*Ds3) 
    T32_1 = (-sJI1) * (c[1,3,1,2]*Sq + c[1,3,2,2]*Dr3 + c[1,3,3,2]*Ds3)
    T33_1 = (-sJI1) * (c[1,3,1,3]*Sq + c[1,3,2,3]*Dr3 + c[1,3,3,3]*Ds3) 
    
    # FACE 2
    # (Nq, Nr, Ns) = (1, 0, 0)
    T11_2 = (sJI2) * (c[1,1,1,1]*Sq + c[1,1,2,1]*Dr3 + c[1,1,3,1]*Ds3) 
    T12_2 = (sJI2) * (c[1,1,1,2]*Sq + c[1,1,2,2]*Dr3 + c[1,1,3,2]*Ds3) 
    T21_2 = (sJI2) * (c[1,2,1,1]*Sq + c[1,2,2,1]*Dr3 + c[1,2,3,1]*Ds3) 
    T13_2 = (sJI2) * (c[1,1,1,3]*Sq + c[1,1,2,3]*Dr3 + c[1,1,3,3]*Ds3) 
    T31_2 = (sJI2) * (c[1,3,1,1]*Sq + c[1,3,2,1]*Dr3 + c[1,3,3,1]*Ds3) 
    T22_2 = (sJI2) * (c[1,2,1,2]*Sq + c[1,2,2,2]*Dr3 + c[1,2,3,2]*Ds3) 
    T23_2 = (sJI2) * (c[1,2,1,3]*Sq + c[1,2,2,3]*Dr3 + c[1,2,3,3]*Ds3) 
    T32_2 = (sJI2) * (c[1,3,1,2]*Sq + c[1,3,2,2]*Dr3 + c[1,3,3,2]*Ds3) 
    T33_2 = (sJI2) * (c[1,3,1,3]*Sq + c[1,3,2,3]*Dr3 + c[1,3,3,3]*Ds3) 

    
    #   # FACE 3
    #     # (Nq, Nr, Ns) = (0, -1, 0)
    T11_3 = (-sJI3) * (c[2,1,1,1]*Dq3 + c[2,1,2,1]*Sr + c[2,1,3,1]*Ds3) 
    T12_3 = (-sJI3) * (c[2,1,1,2]*Dq3 + c[2,1,2,2]*Sr + c[2,1,3,2]*Ds3) 
    T21_3 = (-sJI3) * (c[2,2,1,1]*Dq3 + c[2,2,2,1]*Sr + c[2,2,3,1]*Ds3)
    T13_3 = (-sJI3) * (c[2,1,1,3]*Dq3 + c[2,1,2,3]*Sr + c[2,1,3,3]*Ds3) 
    T31_3 = (-sJI3) * (c[2,3,1,1]*Dq3 + c[2,3,2,1]*Sr + c[2,3,3,1]*Ds3) 
    T22_3 = (-sJI3) * (c[2,2,1,2]*Dq3 + c[2,2,2,2]*Sr + c[2,2,3,2]*Ds3) 
    T23_3 = (-sJI3) * (c[2,2,1,3]*Dq3 + c[2,2,2,3]*Sr + c[2,2,3,3]*Ds3) 
    T32_3 = (-sJI3) * (c[2,3,1,2]*Dq3 + c[2,3,2,2]*Sr + c[2,3,3,2]*Ds3)
    T33_3 = (-sJI3) * (c[2,3,1,3]*Dq3 + c[2,3,2,3]*Sr + c[2,3,3,3]*Ds3) 

    # FACE 4
    # (Nq, Nr, Ns) = (0, 1, 0)
    T11_4 = (sJI4) * (c[2,1,1,1]*Dq3 + c[2,1,2,1]*Sr + c[2,1,3,1]*Ds3) 
    T12_4 = (sJI4) * (c[2,1,1,2]*Dq3 + c[2,1,2,2]*Sr + c[2,1,3,2]*Ds3)
    T21_4 = (sJI4) * (c[2,2,1,1]*Dq3 + c[2,2,2,1]*Sr + c[2,2,3,1]*Ds3)
    T13_4 = (sJI4) * (c[2,1,1,3]*Dq3 + c[2,1,2,3]*Sr + c[2,1,3,3]*Ds3)
    T31_4 = (sJI4) * (c[2,3,1,1]*Dq3 + c[2,3,2,1]*Sr + c[2,3,3,1]*Ds3) 
    T22_4 = (sJI4) * (c[2,2,1,2]*Dq3 + c[2,2,2,2]*Sr + c[2,2,3,2]*Ds3)
    T23_4 = (sJI4) * (c[2,2,1,3]*Dq3 + c[2,2,2,3]*Sr + c[2,2,3,3]*Ds3)
    T32_4 = (sJI4) * (c[2,3,1,2]*Dq3 + c[2,3,2,2]*Sr + c[2,3,3,2]*Ds3)
    T33_4 = (sJI4) * (c[2,3,1,3]*Dq3 + c[2,3,2,3]*Sr + c[2,3,3,3]*Ds3)


    # FACE 5
    # (Nq, Nr, Ns) = (0, 0, -1)
    T11_5 = (-sJI5) * (c[3,1,1,1]*Dq3 + c[3,1,2,1]*Dr3 + c[3,1,3,1]*Ss)
    T12_5 = (-sJI5) * (c[3,1,1,2]*Dq3 + c[3,1,2,2]*Dr3 + c[3,1,3,2]*Ss)
    T21_5 = (-sJI5) * (c[3,2,1,1]*Dq3 + c[3,2,2,1]*Dr3 + c[3,2,3,1]*Ss)
    T13_5 = (-sJI5) * (c[3,1,1,3]*Dq3 + c[3,1,2,3]*Dr3 + c[3,1,3,3]*Ss)
    T31_5 = (-sJI5) * (c[3,3,1,1]*Dq3 + c[3,3,2,1]*Dr3 + c[3,3,3,1]*Ss)
    T22_5 = (-sJI5) * (c[3,2,1,2]*Dq3 + c[3,2,2,2]*Dr3 + c[3,2,3,2]*Ss)
    T23_5 = (-sJI5) * (c[3,2,1,3]*Dq3 + c[3,2,2,3]*Dr3 + c[3,2,3,3]*Ss)
    T32_5 = (-sJI5) * (c[3,3,1,2]*Dq3 + c[3,3,2,2]*Dr3 + c[3,3,3,2]*Ss)
    T33_5 = (-sJI5) * (c[3,3,1,3]*Dq3 + c[3,3,2,3]*Dr3 + c[3,3,3,3]*Ss)


    # FACE 6
    # (Nq, Nr, Ns) = (0, 0, 1)
    T11_6 = sJI6 * (c[3,1,1,1]*Dq3 + c[3,1,2,1]*Dr3 + c[3,1,3,1]*Ss)
    T12_6 = sJI6 *  (c[3,1,1,2]*Dq3 + c[3,1,2,2]*Dr3 + c[3,1,3,2]*Ss)
    T21_6 = sJI6 *  (c[3,2,1,1]*Dq3 + c[3,2,2,1]*Dr3 + c[3,2,3,1]*Ss)
    T13_6 = sJI6 *  (c[3,1,1,3]*Dq3 + c[3,1,2,3]*Dr3 + c[3,1,3,3]*Ss)
    T31_6 = sJI6 *  (c[3,3,1,1]*Dq3 + c[3,3,2,1]*Dr3 + c[3,3,3,1]*Ss)
    T22_6 = sJI6 *  (c[3,2,1,2]*Dq3 + c[3,2,2,2]*Dr3 + c[3,2,3,2]*Ss)
    T23_6 = sJI6 *  (c[3,2,1,3]*Dq3 + c[3,2,2,3]*Dr3 + c[3,2,3,3]*Ss)
    T32_6 = sJI6 *  (c[3,3,1,2]*Dq3 + c[3,3,2,2]*Dr3 + c[3,3,3,2]*Ss)
    T33_6 = sJI6 *  (c[3,3,1,3]*Dq3 + c[3,3,2,3]*Dr3 + c[3,3,3,3]*Ss)
    # }}}
        
    # {{{ Dirichlet BC penalties
    β = 1
    h1 = Hq[1,1] #TODO: fix this
    d = 3
    g = 1

    Z11_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,1,1,1]*Nx1 + c[1,1,2,1]*Ny1 + c[1,1,3,1]*Nz1) + Ny1 * (c[2,1,1,1]*Nx1 + c[2,1,2,1]*Ny1 + c[2,1,3,1]*Nz1) + Nz1 * (c[3,1,1,1]*Nx1 + c[3,1,2,1]*Ny1 + c[3,1,3,1]*Nz1)) 
    Z12_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,1,1,2]*Nx1 + c[1,1,2,2]*Ny1 + c[1,1,3,2]*Nz1) + Ny1 * (c[2,1,1,2]*Nx1 + c[2,1,2,2]*Ny1 + c[2,1,3,2]*Nz1) + Nz1 * (c[3,1,1,2]*Nx1 + c[3,1,2,2]*Ny1 + c[3,1,3,2]*Nz1)) 
    Z21_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,2,1,1]*Nx1 + c[1,2,2,1]*Ny1 + c[1,2,3,1]*Nz1) + Ny1 * (c[2,2,1,1]*Nx1 + c[2,2,2,1]*Ny1 + c[2,2,3,1]*Nz1) + Nz1 * (c[3,2,1,1]*Nx1 + c[3,2,2,1]*Ny1 + c[3,2,3,1]*Nz1)) 
    Z13_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,1,1,3]*Nx1 + c[1,1,2,3]*Ny1 + c[1,1,3,3]*Nz1) + Ny1 * (c[2,1,1,3]*Nx1 + c[2,1,2,3]*Ny1 + c[2,1,3,3]*Nz1) + Nz1 * (c[3,1,1,3]*Nx1 + c[3,1,2,3]*Ny1 + c[3,1,3,3]*Nz1)) 
    Z31_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,3,1,1]*Nx1 + c[1,3,2,1]*Ny1 + c[1,3,3,1]*Nz1) + Ny1 * (c[2,3,1,1]*Nx1 + c[2,3,2,1]*Ny1 + c[2,3,3,1]*Nz1) + Nz1 * (c[3,3,1,1]*Nx1 + c[3,3,2,1]*Ny1 + c[3,3,3,1]*Nz1)) 
    Z22_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,2,1,2]*Nx1 + c[1,2,2,2]*Ny1 + c[1,2,3,2]*Nz1) + Ny1 * (c[2,2,1,2]*Nx1 + c[2,2,2,2]*Ny1 + c[2,2,3,2]*Nz1) + Nz1 * (c[3,2,1,2]*Nx1 + c[3,2,2,2]*Ny1 + c[3,2,3,2]*Nz1)) 
    Z23_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,2,1,3]*Nx1 + c[1,2,2,3]*Ny1 + c[1,2,3,3]*Nz1) + Ny1 * (c[2,2,1,3]*Nx1 + c[2,2,2,3]*Ny1 + c[2,2,3,3]*Nz1) + Nz1 * (c[3,2,1,3]*Nx1 + c[3,2,2,3]*Ny1 + c[3,2,3,3]*Nz1)) 
    Z32_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,3,1,2]*Nx1 + c[1,3,2,2]*Ny1 + c[1,3,3,2]*Nz1) + Ny1 * (c[2,3,1,2]*Nx1 + c[2,3,2,2]*Ny1 + c[2,3,3,2]*Nz1) + Nz1 * (c[3,3,1,2]*Nx1 + c[3,3,2,2]*Ny1 + c[3,3,3,2]*Nz1)) 
    Z33_1 = (β/h1) * d * (EsJ1) * (Nx1 * (c[1,3,1,3]*Nx1 + c[1,3,2,3]*Ny1 + c[1,3,3,3]*Nz1) + Ny1 * (c[2,3,1,3]*Nx1 + c[2,3,2,3]*Ny1 + c[2,3,3,3]*Nz1) + Nz1 * (c[3,3,1,3]*Nx1 + c[3,3,2,3]*Ny1 + c[3,3,3,3]*Nz1)) 

    
    Z11_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,1,1,1]*Nx2 + c[1,1,2,1]*Ny2 + c[1,1,3,1]*Nz2) + Ny2 * (c[2,1,1,1]*Nx2 + c[2,1,2,1]*Ny2 + c[2,1,3,1]*Nz2) + Nz2 * (c[3,1,1,1]*Nx2 + c[3,1,2,1]*Ny2 + c[3,1,3,1]*Nz2)) 
    Z12_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,1,1,2]*Nx2 + c[1,1,2,2]*Ny2 + c[1,1,3,2]*Nz2) + Ny2 * (c[2,1,1,2]*Nx2 + c[2,1,2,2]*Ny2 + c[2,1,3,2]*Nz2) + Nz2 * (c[3,1,1,2]*Nx2 + c[3,1,2,2]*Ny2 + c[3,1,3,2]*Nz2)) 
    Z21_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,2,1,1]*Nx2 + c[1,2,2,1]*Ny2 + c[1,2,3,1]*Nz2) + Ny2 * (c[2,2,1,1]*Nx2 + c[2,2,2,1]*Ny2 + c[2,2,3,1]*Nz2) + Nz2 * (c[3,2,1,1]*Nx2 + c[3,2,2,1]*Ny2 + c[3,2,3,1]*Nz2)) 
    Z13_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,1,1,3]*Nx2 + c[1,1,2,3]*Ny2 + c[1,1,3,3]*Nz2) + Ny2 * (c[2,1,1,3]*Nx2 + c[2,1,2,3]*Ny2 + c[2,1,3,3]*Nz2) + Nz2 * (c[3,1,1,3]*Nx2 + c[3,1,2,3]*Ny2 + c[3,1,3,3]*Nz2)) 
    Z31_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,3,1,1]*Nx2 + c[1,3,2,1]*Ny2 + c[1,3,3,1]*Nz2) + Ny2 * (c[2,3,1,1]*Nx2 + c[2,3,2,1]*Ny2 + c[2,3,3,1]*Nz2) + Nz2 * (c[3,3,1,1]*Nx2 + c[3,3,2,1]*Ny2 + c[3,3,3,1]*Nz2)) 
    Z22_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,2,1,2]*Nx2 + c[1,2,2,2]*Ny2 + c[1,2,3,2]*Nz2) + Ny2 * (c[2,2,1,2]*Nx2 + c[2,2,2,2]*Ny2 + c[2,2,3,2]*Nz2) + Nz2 * (c[3,2,1,2]*Nx2 + c[3,2,2,2]*Ny2 + c[3,2,3,2]*Nz2)) 
    Z23_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,2,1,3]*Nx2 + c[1,2,2,3]*Ny2 + c[1,2,3,3]*Nz2) + Ny2 * (c[2,2,1,3]*Nx2 + c[2,2,2,3]*Ny2 + c[2,2,3,3]*Nz2) + Nz2 * (c[3,2,1,3]*Nx2 + c[3,2,2,3]*Ny2 + c[3,2,3,3]*Nz2)) 
    Z32_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,3,1,2]*Nx2 + c[1,3,2,2]*Ny2 + c[1,3,3,2]*Nz2) + Ny2 * (c[2,3,1,2]*Nx2 + c[2,3,2,2]*Ny2 + c[2,3,3,2]*Nz2) + Nz2 * (c[3,3,1,2]*Nx2 + c[3,3,2,2]*Ny2 + c[3,3,3,2]*Nz2)) 
    Z33_2 = (β/h1) * d * (EsJ2) * (Nx2 * (c[1,3,1,3]*Nx2 + c[1,3,2,3]*Ny2 + c[1,3,3,3]*Nz2) + Ny2 * (c[2,3,1,3]*Nx2 + c[2,3,2,3]*Ny2 + c[2,3,3,3]*Nz2) + Nz2 * (c[3,3,1,3]*Nx2 + c[3,3,2,3]*Ny2 + c[3,3,3,3]*Nz2)) 
    

     # FACE 3 
    # (nq, nr, ns) = (0, -1, 0)
    Z11_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,1,1,1]*Nx3 + c[1,1,2,1]*Ny3 + c[1,1,3,1]*Nz3) + Ny3 * (c[2,1,1,1]*Nx3 + c[2,1,2,1]*Ny3 + c[2,1,3,1]*Nz3) + Nz3 * (c[3,1,1,1]*Nx3 + c[3,1,2,1]*Ny3 + c[3,1,3,1]*Nz3)) 
    Z12_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,1,1,2]*Nx3 + c[1,1,2,2]*Ny3 + c[1,1,3,2]*Nz3) + Ny3 * (c[2,1,1,2]*Nx3 + c[2,1,2,2]*Ny3 + c[2,1,3,2]*Nz3) + Nz3 * (c[3,1,1,2]*Nx3 + c[3,1,2,2]*Ny3 + c[3,1,3,2]*Nz3)) 
    Z21_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,2,1,1]*Nx3 + c[1,2,2,1]*Ny3 + c[1,2,3,1]*Nz3) + Ny3 * (c[2,2,1,1]*Nx3 + c[2,2,2,1]*Ny3 + c[2,2,3,1]*Nz3) + Nz3 * (c[3,2,1,1]*Nx3 + c[3,2,2,1]*Ny3 + c[3,2,3,1]*Nz3)) 
    Z13_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,1,1,3]*Nx3 + c[1,1,2,3]*Ny3 + c[1,1,3,3]*Nz3) + Ny3 * (c[2,1,1,3]*Nx3 + c[2,1,2,3]*Ny3 + c[2,1,3,3]*Nz3) + Nz3 * (c[3,1,1,3]*Nx3 + c[3,1,2,3]*Ny3 + c[3,1,3,3]*Nz3)) 
    Z31_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,3,1,1]*Nx3 + c[1,3,2,1]*Ny3 + c[1,3,3,1]*Nz3) + Ny3 * (c[2,3,1,1]*Nx3 + c[2,3,2,1]*Ny3 + c[2,3,3,1]*Nz3) + Nz3 * (c[3,3,1,1]*Nx3 + c[3,3,2,1]*Ny3 + c[3,3,3,1]*Nz3)) 
    Z22_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,2,1,2]*Nx3 + c[1,2,2,2]*Ny3 + c[1,2,3,2]*Nz3) + Ny3 * (c[2,2,1,2]*Nx3 + c[2,2,2,2]*Ny3 + c[2,2,3,2]*Nz3) + Nz3 * (c[3,2,1,2]*Nx3 + c[3,2,2,2]*Ny3 + c[3,2,3,2]*Nz3)) 
    Z23_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,2,1,3]*Nx3 + c[1,2,2,3]*Ny3 + c[1,2,3,3]*Nz3) + Ny3 * (c[2,2,1,3]*Nx3 + c[2,2,2,3]*Ny3 + c[2,2,3,3]*Nz3) + Nz3 * (c[3,2,1,3]*Nx3 + c[3,2,2,3]*Ny3 + c[3,2,3,3]*Nz3)) 
    Z32_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,3,1,2]*Nx3 + c[1,3,2,2]*Ny3 + c[1,3,3,2]*Nz3) + Ny3 * (c[2,3,1,2]*Nx3 + c[2,3,2,2]*Ny3 + c[2,3,3,2]*Nz3) + Nz3 * (c[3,3,1,2]*Nx3 + c[3,3,2,2]*Ny3 + c[3,3,3,2]*Nz3)) 
    Z33_3 = (β/h1) * d * (EsJ3) * (Nx3 * (c[1,3,1,3]*Nx3 + c[1,3,2,3]*Ny3 + c[1,3,3,3]*Nz3) + Ny3 * (c[2,3,1,3]*Nx3 + c[2,3,2,3]*Ny3 + c[2,3,3,3]*Nz3) + Nz3 * (c[3,3,1,3]*Nx3 + c[3,3,2,3]*Ny3 + c[3,3,3,3]*Nz3)) 

   # FACE 4 
    # (nq, nr, ns) = (0, 1, 0)
    Z11_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,1,1,1]*Nx4 + c[1,1,2,1]*Ny4 + c[1,1,3,1]*Nz4) + Ny4 * (c[2,1,1,1]*Nx4 + c[2,1,2,1]*Ny4 + c[2,1,3,1]*Nz4) + Nz4 * (c[3,1,1,1]*Nx4 + c[3,1,2,1]*Ny4 + c[3,1,3,1]*Nz4)) 
    Z12_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,1,1,2]*Nx4 + c[1,1,2,2]*Ny4 + c[1,1,3,2]*Nz4) + Ny4 * (c[2,1,1,2]*Nx4 + c[2,1,2,2]*Ny4 + c[2,1,3,2]*Nz4) + Nz4 * (c[3,1,1,2]*Nx4 + c[3,1,2,2]*Ny4 + c[3,1,3,2]*Nz4)) 
    Z21_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,2,1,1]*Nx4 + c[1,2,2,1]*Ny4 + c[1,2,3,1]*Nz4) + Ny4 * (c[2,2,1,1]*Nx4 + c[2,2,2,1]*Ny4 + c[2,2,3,1]*Nz4) + Nz4 * (c[3,2,1,1]*Nx4 + c[3,2,2,1]*Ny4 + c[3,2,3,1]*Nz4)) 
    Z13_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,1,1,3]*Nx4 + c[1,1,2,3]*Ny4 + c[1,1,3,3]*Nz4) + Ny4 * (c[2,1,1,3]*Nx4 + c[2,1,2,3]*Ny4 + c[2,1,3,3]*Nz4) + Nz4 * (c[3,1,1,3]*Nx4 + c[3,1,2,3]*Ny4 + c[3,1,3,3]*Nz4)) 
    Z31_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,3,1,1]*Nx4 + c[1,3,2,1]*Ny4 + c[1,3,3,1]*Nz4) + Ny4 * (c[2,3,1,1]*Nx4 + c[2,3,2,1]*Ny4 + c[2,3,3,1]*Nz4) + Nz4 * (c[3,3,1,1]*Nx4 + c[3,3,2,1]*Ny4 + c[3,3,3,1]*Nz4)) 
    Z22_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,2,1,2]*Nx4 + c[1,2,2,2]*Ny4 + c[1,2,3,2]*Nz4) + Ny4 * (c[2,2,1,2]*Nx4 + c[2,2,2,2]*Ny4 + c[2,2,3,2]*Nz4) + Nz4 * (c[3,2,1,2]*Nx4 + c[3,2,2,2]*Ny4 + c[3,2,3,2]*Nz4)) 
    Z23_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,2,1,3]*Nx4 + c[1,2,2,3]*Ny4 + c[1,2,3,3]*Nz4) + Ny4 * (c[2,2,1,3]*Nx4 + c[2,2,2,3]*Ny4 + c[2,2,3,3]*Nz4) + Nz4 * (c[3,2,1,3]*Nx4 + c[3,2,2,3]*Ny4 + c[3,2,3,3]*Nz4)) 
    Z32_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,3,1,2]*Nx4 + c[1,3,2,2]*Ny4 + c[1,3,3,2]*Nz4) + Ny4 * (c[2,3,1,2]*Nx4 + c[2,3,2,2]*Ny4 + c[2,3,3,2]*Nz4) + Nz4 * (c[3,3,1,2]*Nx4 + c[3,3,2,2]*Ny4 + c[3,3,3,2]*Nz4)) 
    Z33_4 = (β/h1) * d * (EsJ4) * (Nx4 * (c[1,3,1,3]*Nx4 + c[1,3,2,3]*Ny4 + c[1,3,3,3]*Nz4) + Ny4 * (c[2,3,1,3]*Nx4 + c[2,3,2,3]*Ny4 + c[2,3,3,3]*Nz4) + Nz4 * (c[3,3,1,3]*Nx4 + c[3,3,2,3]*Ny4 + c[3,3,3,3]*Nz4)) 

    # FACE 5 
    # (nq, nr, ns) = (0, 0, -1)
    Z11_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,1,1,1]*Nx5 + c[1,1,2,1]*Ny5 + c[1,1,3,1]*Nz5) + Ny5 * (c[2,1,1,1]*Nx5 + c[2,1,2,1]*Ny5 + c[2,1,3,1]*Nz5) + Nz5 * (c[3,1,1,1]*Nx5 + c[3,1,2,1]*Ny5 + c[3,1,3,1]*Nz5)) 
    Z12_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,1,1,2]*Nx5 + c[1,1,2,2]*Ny5 + c[1,1,3,2]*Nz5) + Ny5 * (c[2,1,1,2]*Nx5 + c[2,1,2,2]*Ny5 + c[2,1,3,2]*Nz5) + Nz5 * (c[3,1,1,2]*Nx5 + c[3,1,2,2]*Ny5 + c[3,1,3,2]*Nz5)) 
    Z21_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,2,1,1]*Nx5 + c[1,2,2,1]*Ny5 + c[1,2,3,1]*Nz5) + Ny5 * (c[2,2,1,1]*Nx5 + c[2,2,2,1]*Ny5 + c[2,2,3,1]*Nz5) + Nz5 * (c[3,2,1,1]*Nx5 + c[3,2,2,1]*Ny5 + c[3,2,3,1]*Nz5)) 
    Z13_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,1,1,3]*Nx5 + c[1,1,2,3]*Ny5 + c[1,1,3,3]*Nz5) + Ny5 * (c[2,1,1,3]*Nx5 + c[2,1,2,3]*Ny5 + c[2,1,3,3]*Nz5) + Nz5 * (c[3,1,1,3]*Nx5 + c[3,1,2,3]*Ny5 + c[3,1,3,3]*Nz5)) 
    Z31_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,3,1,1]*Nx5 + c[1,3,2,1]*Ny5 + c[1,3,3,1]*Nz5) + Ny5 * (c[2,3,1,1]*Nx5 + c[2,3,2,1]*Ny5 + c[2,3,3,1]*Nz5) + Nz5 * (c[3,3,1,1]*Nx5 + c[3,3,2,1]*Ny5 + c[3,3,3,1]*Nz5)) 
    Z22_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,2,1,2]*Nx5 + c[1,2,2,2]*Ny5 + c[1,2,3,2]*Nz5) + Ny5 * (c[2,2,1,2]*Nx5 + c[2,2,2,2]*Ny5 + c[2,2,3,2]*Nz5) + Nz5 * (c[3,2,1,2]*Nx5 + c[3,2,2,2]*Ny5 + c[3,2,3,2]*Nz5)) 
    Z23_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,2,1,3]*Nx5 + c[1,2,2,3]*Ny5 + c[1,2,3,3]*Nz5) + Ny5 * (c[2,2,1,3]*Nx5 + c[2,2,2,3]*Ny5 + c[2,2,3,3]*Nz5) + Nz5 * (c[3,2,1,3]*Nx5 + c[3,2,2,3]*Ny5 + c[3,2,3,3]*Nz5)) 
    Z32_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,3,1,2]*Nx5 + c[1,3,2,2]*Ny5 + c[1,3,3,2]*Nz5) + Ny5 * (c[2,3,1,2]*Nx5 + c[2,3,2,2]*Ny5 + c[2,3,3,2]*Nz5) + Nz5 * (c[3,3,1,2]*Nx5 + c[3,3,2,2]*Ny5 + c[3,3,3,2]*Nz5)) 
    Z33_5 = (β/h1) * d * (EsJ5) * (Nx5 * (c[1,3,1,3]*Nx5 + c[1,3,2,3]*Ny5 + c[1,3,3,3]*Nz5) + Ny5 * (c[2,3,1,3]*Nx5 + c[2,3,2,3]*Ny5 + c[2,3,3,3]*Nz5) + Nz5 * (c[3,3,1,3]*Nx5 + c[3,3,2,3]*Ny5 + c[3,3,3,3]*Nz5)) 

    # FACE 6 
    # (nq, nr, ns) = (0, 0, 1)
    Z11_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,1,1,1]*Nx6 + c[1,1,2,1]*Ny6 + c[1,1,3,1]*Nz6) + Ny6 * (c[2,1,1,1]*Nx6 + c[2,1,2,1]*Ny6 + c[2,1,3,1]*Nz6) + Nz6 * (c[3,1,1,1]*Nx6 + c[3,1,2,1]*Ny6 + c[3,1,3,1]*Nz6)) 
    Z12_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,1,1,2]*Nx6 + c[1,1,2,2]*Ny6 + c[1,1,3,2]*Nz6) + Ny6 * (c[2,1,1,2]*Nx6 + c[2,1,2,2]*Ny6 + c[2,1,3,2]*Nz6) + Nz6 * (c[3,1,1,2]*Nx6 + c[3,1,2,2]*Ny6 + c[3,1,3,2]*Nz6)) 
    Z21_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,2,1,1]*Nx6 + c[1,2,2,1]*Ny6 + c[1,2,3,1]*Nz6) + Ny6 * (c[2,2,1,1]*Nx6 + c[2,2,2,1]*Ny6 + c[2,2,3,1]*Nz6) + Nz6 * (c[3,2,1,1]*Nx6 + c[3,2,2,1]*Ny6 + c[3,2,3,1]*Nz6)) 
    Z13_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,1,1,3]*Nx6 + c[1,1,2,3]*Ny6 + c[1,1,3,3]*Nz6) + Ny6 * (c[2,1,1,3]*Nx6 + c[2,1,2,3]*Ny6 + c[2,1,3,3]*Nz6) + Nz6 * (c[3,1,1,3]*Nx6 + c[3,1,2,3]*Ny6 + c[3,1,3,3]*Nz6)) 
    Z31_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,3,1,1]*Nx6 + c[1,3,2,1]*Ny6 + c[1,3,3,1]*Nz6) + Ny6 * (c[2,3,1,1]*Nx6 + c[2,3,2,1]*Ny6 + c[2,3,3,1]*Nz6) + Nz6 * (c[3,3,1,1]*Nx6 + c[3,3,2,1]*Ny6 + c[3,3,3,1]*Nz6)) 
    Z22_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,2,1,2]*Nx6 + c[1,2,2,2]*Ny6 + c[1,2,3,2]*Nz6) + Ny6 * (c[2,2,1,2]*Nx6 + c[2,2,2,2]*Ny6 + c[2,2,3,2]*Nz6) + Nz6 * (c[3,2,1,2]*Nx6 + c[3,2,2,2]*Ny6 + c[3,2,3,2]*Nz6)) 
    Z23_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,2,1,3]*Nx6 + c[1,2,2,3]*Ny6 + c[1,2,3,3]*Nz6) + Ny6 * (c[2,2,1,3]*Nx6 + c[2,2,2,3]*Ny6 + c[2,2,3,3]*Nz6) + Nz6 * (c[3,2,1,3]*Nx6 + c[3,2,2,3]*Ny6 + c[3,2,3,3]*Nz6)) 
    Z32_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,3,1,2]*Nx6 + c[1,3,2,2]*Ny6 + c[1,3,3,2]*Nz6) + Ny6 * (c[2,3,1,2]*Nx6 + c[2,3,2,2]*Ny6 + c[2,3,3,2]*Nz6) + Nz6 * (c[3,3,1,2]*Nx6 + c[3,3,2,2]*Ny6 + c[3,3,3,2]*Nz6)) 
    Z33_6 = (β/h1) * d * (EsJ6) * (Nx6 * (c[1,3,1,3]*Nx6 + c[1,3,2,3]*Ny6 + c[1,3,3,3]*Nz6) + Ny6 * (c[2,3,1,3]*Nx6 + c[2,3,2,3]*Ny6 + c[2,3,3,3]*Nz6) + Nz6 * (c[3,3,1,3]*Nx6 + c[3,3,2,3]*Ny6 + c[3,3,3,3]*Nz6)) 


    # Dirichlet faces:
    dSAT_11 =  JHI * ((n*T11_1 .- Z11_1)'*e1*sJ1*H1*e1T + (n*T11_2 .- Z11_2)'*e2*sJ2*H2*e2T + (n*T11_3 .- Z11_3)'*e3*sJ3*H3*e3T + (n*T11_4 .- Z11_4)'*e4*sJ4*H4*e4T + (n*T11_5 .- Z11_5)'*e5*sJ5*H5*e5T + (n*T11_6 .- Z11_6)'*e6*sJ6*H6*e6T)
    dSAT_12 =  JHI * ((n*T21_1 .- Z21_1)'*e1*sJ1*H1*e1T + (n*T21_2 .- Z21_2)'*e2*sJ2*H2*e2T + (n*T21_3 .- Z21_3)'*e3*sJ3*H3*e3T + (n*T21_4 .- Z21_4)'*e4*sJ4*H4*e4T + (n*T21_5 .- Z21_5)'*e5*sJ5*H5*e5T + (n*T21_6 .- Z21_6)'*e6*sJ6*H6*e6T)    
    dSAT_13 =  JHI * ((n*T31_1 .- Z31_1)'*e1*sJ1*H1*e1T + (n*T31_2 .- Z31_2)'*e2*sJ2*H2*e2T + (n*T31_3 .- Z31_3)'*e3*sJ3*H3*e3T + (n*T31_4 .- Z31_4)'*e4*sJ4*H4*e4T + (n*T31_5 .- Z31_5)'*e5*sJ5*H5*e5T + (n*T31_6 .- Z31_6)'*e6*sJ6*H6*e6T) 
    dSAT_21 =  JHI * ((n*T12_1 .- Z12_1)'*e1*sJ1*H1*e1T + (n*T12_2 .- Z12_2)'*e2*sJ2*H2*e2T + (n*T12_3 .- Z12_3)'*e3*sJ3*H3*e3T + (n*T12_4 .- Z12_4)'*e4*sJ4*H4*e4T + (n*T12_5 .- Z12_5)'*e5*sJ5*H5*e5T + (n*T12_6 .- Z12_6)'*e6*sJ6*H6*e6T) 
    dSAT_22 =  JHI * ((n*T22_1 .- Z22_1)'*e1*sJ1*H1*e1T + (n*T22_2 .- Z22_2)'*e2*sJ2*H2*e2T + (n*T22_3 .- Z22_3)'*e3*sJ3*H3*e3T + (n*T22_4 .- Z22_4)'*e4*sJ4*H4*e4T + (n*T22_5 .- Z22_5)'*e5*sJ5*H5*e5T + (n*T22_6 .- Z22_6)'*e6*sJ6*H6*e6T) 
    dSAT_23 =  JHI * ((n*T32_1 .- Z32_1)'*e1*sJ1*H1*e1T + (n*T32_2 .- Z32_2)'*e2*sJ2*H2*e2T + (n*T32_3 .- Z32_3)'*e3*sJ3*H3*e3T + (n*T32_4 .- Z32_4)'*e4*sJ4*H4*e4T + (n*T32_5 .- Z32_5)'*e5*sJ5*H5*e5T + (n*T32_6 .- Z32_6)'*e6*sJ6*H6*e6T) 
    dSAT_31 = JHI * ((n*T13_1 .- Z13_1)'*e1*sJ1*H1*e1T + (n*T13_2 .- Z13_2)'*e2*sJ2*H2*e2T + (n*T13_3 .- Z13_3)'*e3*sJ3*H3*e3T + (n*T13_4 .- Z13_4)'*e4*sJ4*H4*e4T + (n*T13_5 .- Z13_5)'*e5*sJ5*H5*e5T + (n*T13_6 .- Z13_6)'*e6*sJ6*H6*e6T)
    dSAT_32 = JHI * ((n*T23_1 .- Z23_1)'*e1*sJ1*H1*e1T + (n*T23_2 .- Z23_2)'*e2*sJ2*H2*e2T + (n*T23_3 .- Z23_3)'*e3*sJ3*H3*e3T + (n*T23_4 .- Z23_4)'*e4*sJ4*H4*e4T + (n*T23_5 .- Z23_5)'*e5*sJ5*H5*e5T + (n*T23_6 .- Z23_6)'*e6*sJ6*H6*e6T)
    dSAT_33 = JHI * ((n*T33_1 .- Z33_1)'*e1*sJ1*H1*e1T + (n*T33_2 .- Z33_2)'*e2*sJ2*H2*e2T + (n*T33_3 .- Z33_3)'*e3*sJ3*H3*e3T + (n*T33_4 .- Z33_4)'*e4*sJ4*H4*e4T + (n*T33_5 .- Z33_5)'*e5*sJ5*H5*e5T + (n*T33_6 .- Z33_6)'*e6*sJ6*H6*e6T)  
    
    dSAT = (dSAT_11, dSAT_12, dSAT_13, dSAT_21, dSAT_22, dSAT_23, dSAT_31, dSAT_32, dSAT_33)
    # }}}
    
    # {{{ Traction faces:

    nSAT_11 = -JHI * (e1*sJ1*H1*e1'*T11_1 + e2*sJ2*H2*e2'*T11_2  + e3*sJ3*H3*e3'*T11_3 + e4*sJ4*H4*e4'*T11_4 + e5*sJ5*H5*e5'*T11_5 + e6*sJ6*H6*e6'*T11_6)
    nSAT_12 = -JHI * (e1*sJ1*H1*e1'*T11_1 + e2*sJ2*H2*e2'*T11_2  + e3*sJ3*H3*e3'*T12_3 + e4*sJ4*H4*e4'*T12_4 + e5*sJ5*H5*e5'*T12_5 + e6*sJ6*H6*e6'*T12_6)

    # }}}
    
    n = 1 #default to n = 1
    # Create SAT vectors - all DIRICHLET conditions
     
    # b11 = [JHI*(n*T11_1 .- Z11_1)'*e1*sJ1*H1, JHI*(n*T11_2 .- Z11_2)'*e2*sJ2*H2, JHI*(n*T11_3 .- Z11_3)'*e3*sJ3*H3, JHI*(n*T11_4 .- Z11_4)'*e4*sJ4*H4, JHI*(n*T11_5 .- Z11_5)'*e5*sJ5*H5, JHI*(n*T11_6 .- Z11_6)'*e6*sJ6*H6]
    # b12 = [JHI*(n*T21_1 .- Z21_1)'*e1*sJ1*H1, JHI*(n*T21_2 .- Z21_2)'*e2*sJ2*H2, JHI*(n*T21_3 .- Z21_3)'*e3*sJ3*H3, JHI*(n*T21_4 .- Z21_4)'*e4*sJ4*H4, JHI*(n*T21_5 .- Z21_5)'*e5*sJ5*H5, JHI*(n*T21_6 .- Z21_6)'*e6*sJ6*H6]
    # b13 = [JHI*(n*T31_1 .- Z31_1)'*e1*sJ1*H1, JHI*(n*T31_2 .- Z31_2)'*e2*sJ2*H2, JHI*(n*T31_3 .- Z31_3)'*e3*sJ3*H3, JHI*(n*T31_4 .- Z31_4)'*e4*sJ4*H4, JHI*(n*T31_5 .- Z31_5)'*e5*sJ5*H5, JHI*(n*T31_6 .- Z31_6)'*e6*sJ6*H6]

     
    # b21 = [JHI*(n*T12_1 .- Z12_1)'*e1*sJ1*H1, JHI*(n*T12_2 .- Z12_2)'*e2*sJ2*H2, JHI*(n*T12_3 .- Z12_3)'*e3*sJ3*H3, JHI*(n*T12_4 .- Z12_4)'*e4*sJ4*H4, JHI*(n*T12_5 .- Z12_5)'*e5*sJ5*H5, JHI*(n*T12_6 .- Z12_6)'*e6*sJ6*H6]
    # b22 = [JHI*(n*T22_1 .- Z22_1)'*e1*sJ1*H1, JHI*(n*T22_2 .- Z22_2)'*e2*sJ2*H2, JHI*(n*T22_3 .- Z22_3)'*e3*sJ3*H3, JHI*(n*T22_4 .- Z22_4)'*e4*sJ4*H4, JHI*(n*T22_5 .- Z22_5)'*e5*sJ5*H5, JHI*(n*T22_6 .- Z22_6)'*e6*sJ6*H6]
    # b23 = [JHI*(n*T32_1 .- Z32_1)'*e1*sJ1*H1, JHI*(n*T32_2 .- Z32_2)'*e2*sJ2*H2, JHI*(n*T32_3 .- Z32_3)'*e3*sJ3*H3, JHI*(n*T32_4 .- Z32_4)'*e4*sJ4*H4, JHI*(n*T32_5 .- Z32_5)'*e5*sJ5*H5, JHI*(n*T32_6 .- Z32_6)'*e6*sJ6*H6]

          
    # b31 = [JHI*(n*T13_1 .- Z13_1)'*e1*sJ1*H1, JHI*(n*T13_2 .- Z13_2)'*e2*sJ2*H2, JHI*(n*T13_3 .- Z13_3)'*e3*sJ3*H3, JHI*(n*T13_4 .- Z13_4)'*e4*sJ4*H4, JHI*(n*T13_5 .- Z13_5)'*e5*sJ5*H5, JHI*(n*T13_6 .- Z13_6)'*e6*sJ6*H6]
    # b32 = [JHI*(n*T23_1 .- Z23_1)'*e1*sJ1*H1, JHI*(n*T23_2 .- Z23_2)'*e2*sJ2*H2, JHI*(n*T23_3 .- Z23_3)'*e3*sJ3*H3, JHI*(n*T23_4 .- Z23_4)'*e4*sJ4*H4, JHI*(n*T23_5 .- Z23_5)'*e5*sJ5*H5, JHI*(n*T23_6 .- Z23_6)'*e6*sJ6*H6]
    # b33 = [JHI*(n*T33_1 .- Z33_1)'*e1*sJ1*H1, JHI*(n*T33_2 .- Z33_2)'*e2*sJ2*H2, JHI*(n*T33_3 .- Z33_3)'*e3*sJ3*H3, JHI*(n*T33_4 .- Z33_4)'*e4*sJ4*H4, JHI*(n*T33_5 .- Z33_5)'*e5*sJ5*H5, JHI*(n*T33_6 .- Z33_6)'*e6*sJ6*H6] 
         

    # Create SAT vectors - DIRICHLET conditions on faces 1 and 2, else traction

    - JHI * (e1*sJ1*H1*e1'*T11_1 + e2*sJ2*H2*e2'*T11_2  + 3*sJ3*H3*e3'*T11_3 + e4*sJ4*H4*e4'*T11_4 + e5*sJ5*H5*e5'*T11_5 + e6*sJ6*H6*e6'*T11_6)
     - JHI * (e3*sJ3*H3*e3'*T12_3 + e4*sJ4*H4*e4'*T12_4 + e5*sJ5*H5*e5'*T12_5 + e6*sJ6*H6*e6'*T12_6)
     - JHI * (e3*sJ3*H3*e3'*T13_3 + e4*sJ4*H4*e4'*T13_4 + e5*sJ5*H5*e5'*T13_5 + e6*sJ6*H6*e6'*T13_6)
    #S11 =  JHI * ((T11_1 .- Z11_1)'*e1*sJ1*H1*e1T + (T11_2 .- Z11_2)'*e2*sJ2*H2*e2T + (T11_3 .- Z11_3)'*e3*sJ3*H3*e3T + (T11_4 .- Z11_4)'*e4*sJ4*H4*e4T + (T11_5 .- Z11_5)'*e5*sJ5*H5*e5T + (T11_6 .- Z11_6)'*e6*sJ6*H6*e6T)
    #S12 =  JHI * ((T21_1 .- Z21_1)'*e1*sJ1*H1*e1T + (T21_2 .- Z21_2)'*e2*sJ2*H2*e2T + (T21_3 .- Z21_3)'*e3*sJ3*H3*e3T + (T21_4 .- Z21_4)'*e4*sJ4*H4*e4T + (T21_5 .- Z21_5)'*e5*sJ5*H5*e5T + (T21_6 .- Z21_6)'*e6*sJ6*H6*e6T)    
    #S13 =  JHI * ((T31_1 .- Z31_1)'*e1*sJ1*H1*e1T + (T31_2 .- Z31_2)'*e2*sJ2*H2*e2T + (T31_3 .- Z31_3)'*e3*sJ3*H3*e3T + (T31_4 .- Z31_4)'*e4*sJ4*H4*e4T + (T31_5 .- Z31_5)'*e5*sJ5*H5*e5T + (T31_6 .- Z31_6)'*e6*sJ6*H6*e6T)
    b11 = [JHI*(T11_1 .- Z11_1)'*e1*sJ1*H1, JHI*(T11_2 .- Z11_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b12 = [JHI*(T21_1 .- Z21_1)'*e1*sJ1*H1, JHI*(T21_2 .- Z21_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b13 = [JHI*(T31_1 .- Z31_1)'*e1*sJ1*H1, JHI*(T31_2 .- Z31_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]

    S21 =  JHI * ((T12_1 .- Z12_1)'*e1*sJ1*H1*e1T + (T12_2 .- Z12_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T21_3 + e4*sJ4*H4*e4'*T21_4 + e5*sJ5*H5*e5'*T21_5 + e6*sJ6*H6*e6'*T21_6)
    S22 =  JHI * ((T22_1 .- Z22_1)'*e1*sJ1*H1*e1T + (T22_2 .- Z22_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T22_3 + e4*sJ4*H4*e4'*T22_4 + e5*sJ5*H5*e5'*T22_5 + e6*sJ6*H6*e6'*T22_6)
    S23 =  JHI * ((T32_1 .- Z32_1)'*e1*sJ1*H1*e1T + (T32_2 .- Z32_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T23_3 + e4*sJ4*H4*e4'*T23_4 + e5*sJ5*H5*e5'*T23_5 + e6*sJ6*H6*e6'*T23_6)
    b21 = [JHI*(T12_1 .- Z12_1)'*e1*sJ1*H1, JHI*(T12_2 .- Z12_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b22 = [JHI*(T22_1 .- Z22_1)'*e1*sJ1*H1, JHI*(T22_2 .- Z22_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b23 = [JHI*(T32_1 .- Z32_1)'*e1*sJ1*H1, JHI*(T32_2 .- Z32_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]

    S31 = JHI * ((T13_1 .- Z13_1)'*e1*sJ1*H1*e1T + (T13_2 .- Z13_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T31_3 + e4*sJ4*H4*e4'*T31_4 + e5*sJ5*H5*e5'*T31_5 + e6*sJ6*H6*e6'*T31_6)
    S32 = JHI * ((T23_1 .- Z23_1)'*e1*sJ1*H1*e1T + (T23_2 .- Z23_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T32_3 + e4*sJ4*H4*e4'*T32_4 + e5*sJ5*H5*e5'*T32_5 + e6*sJ6*H6*e6'*T32_6)
    S33 = JHI * ((T33_1 .- Z33_1)'*e1*sJ1*H1*e1T + (T33_2 .- Z33_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T33_3 + e4*sJ4*H4*e4'*T33_4 + e5*sJ5*H5*e5'*T33_5 + e6*sJ6*H6*e6'*T33_6)        
    b31 = [JHI*(T13_1 .- Z13_1)'*e1*sJ1*H1, JHI*(T13_2 .- Z13_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b32 = [JHI*(T23_1 .- Z23_1)'*e1*sJ1*H1, JHI*(T23_2 .- Z23_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b33 = [JHI*(T33_1 .- Z33_1)'*e1*sJ1*H1, JHI*(T33_2 .- Z33_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)] 
         


    # Equation 1: J = 1
    # 0 = M11*u1 + M12*u2 + M13*u3 + f1 + SAT1
    # OR 
    # 0 = M11*u1 + M12*u2 + M13*u3 + f1 + S11*u1 + S12*u2 + S13*u3 - b1
     
    # Equation 2: J = 2
    # 0 = M21*u1 + M22*u2 + M23*u3 + f2 + SAT2
    # OR
    # 0 = M21*u1 + M22*u2 + M23*u3 + f2 + S21*u1 + S22*u2 + S23*u3 - b3

    # Equation 3: J = 3
    # 0 = M31*u1 + M32*u2 + M33*u3 + f3 + SAT3
    # OR
    # 0 = M31*u1 + M32*u2 + M33*u3 + f3 + S31*u1 + S32*u2 + S33*u3

    # AND ALL TOGETHER: MU = [B11*g1 + B12*g2 + B13*g3; B21*g1 + B22*g2 + B23*g3; B31*g1 + B32*g2 + B33*g3] + J*H*f  where
    M = [M11 M12 M13; M21 M22 M23; M31 M32 M33]
    S = [S11 S12 S13; S21 S22 S23; S31 S32 S33]
    
    A = M + S

    
    B = (b11, 1*b12, 1*b13, 1*b21, b22, 1*b23, 1*b31, 1*b32, b33)
   
    # and f = [f1; f2; f3]
    # where U = [u1; u2; u3]
    JH = J*H
    
    #return (M, B, JH, A, S, HqI, HrI, HsI)


    return (A, dB, nB, iB, eRST, tSAT, H̃, T, TT, eRS, J, coord = metrics.coord, IsJZ, Hinv, H,
    facecoord = metrics.facecoord, 
    sJ = metrics.sJ,
    nx = metrics.nx,
    ny = metrics.ny,
    sJZ)



end


function var_3D_D2q(p, Nqp, Nrp, Nsp, C, HIq; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    Iq = sparse(I, Nqp, Nqp)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)

    N = Nqp*Nrp*Nsp
    D2q = spzeros(N, N) # initialize
    S0q = spzeros(N, N) # initialize
    SNq = spzeros(N, N) # initialize
    for i = 1:Nrp
        for j = 1:Nsp
            B = C[:, i, j]# get coefficient on 1D line in q-direction
            (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nqp-1, B; xc = (-1,1))
            ej = spzeros(Nsp, 1)
            ej[j] = 1
            ei = spzeros(Nrp, 1)
            ei[i] = 1
            D2q += (ej ⊗ ei ⊗ Iq) * D2 * (ej' ⊗ ei' ⊗ Iq)
            S0q += (ej ⊗ ei ⊗ Iq) * S0 * (ej' ⊗ ei' ⊗ Iq)
            SNq += (ej ⊗ ei ⊗ Iq) * SN * (ej' ⊗ ei' ⊗ Iq)
        end
    end
    return D2q, S0q, SNq
end

function var_3D_D2r(p, Nqp, Nrp, Nsp, C, HIr; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    Iq = sparse(I, Nqp, Nqp)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)

    N = Nqp*Nrp*Nsp
    D2r = spzeros(N, N) # initialize
    S0r = spzeros(N, N) # initialize
    SNr = spzeros(N, N) # initialize
    for i = 1:Nqp
        for j = 1:Nsp
            B = C[i, :, j]# get coefficient on 1D line in r-direction
            (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nrp-1, B; xc = (-1,1))
            ej = spzeros(Nsp, 1)
            ej[j] = 1
            ei = spzeros(Nqp, 1)
            ei[i] = 1
            D2r += (ej ⊗ Ir ⊗ ei) * D2 * (ej' ⊗ Ir ⊗ ei')
            S0r += (ej ⊗ Ir ⊗ ei) * S0 * (ej' ⊗ Ir ⊗ ei')
            SNr += (ej ⊗ Ir ⊗ ei) * SN * (ej' ⊗ Ir ⊗ ei')
        end
    end
    return D2r, S0r, SNr
end

function var_3D_D2s(p, Nqp, Nrp, Nsp, C, HIq; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    Iq = sparse(I, Nqp, Nqp)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)

    N = Nqp*Nrp*Nsp
    D2s = spzeros(N, N) # initialize
    S0s = spzeros(N, N) # initialize
    SNs = spzeros(N, N) # initialize
    for i = 1:Nrp
        for j = 1:Nqp
            B = C[j, i, :]# get coefficient on 1D line in s-direction
            (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nsp-1, B; xc = (-1,1))
            ej = spzeros(Nqp, 1)
            ej[j] = 1
            ei = spzeros(Nrp, 1)
            ei[i] = 1
            D2s += (Is ⊗ ei ⊗ ej) * D2 * (Is ⊗ ei' ⊗ ej')
            S0s += (Is ⊗ ei ⊗ ej) * S0 * (Is ⊗ ei' ⊗ ej')
            SNs += (Is ⊗ ei ⊗ ej) * SN * (Is ⊗ ei' ⊗ ej')
        end
    end
    return D2s, S0s, SNs
end
