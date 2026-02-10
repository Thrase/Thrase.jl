using LinearAlgebra
using SparseArrays

function get_operators(p, Nx, Nz, μ; xc = (-1, 1), zc = (-1, 1))

  Nxp = Nx + 1
  Nzp = Nz + 1


# Derivative operators for the rest of the computation
  (Dx, HxI, Hx, x) = diagonal_sbp_D1(p, Nx; xc = xc)
  (Dz, HzI, Hz, z) = diagonal_sbp_D1(p, Nz; xc = zc)


# Compatible SBP
  (D2x, S0x, SNx, _, _, _) = diagonal_sbp_D2(p, Nx; xc = xc)
  (D2z, S0z, SNz, _, _, _) = diagonal_sbp_D2(p, Nz; xc = zc)

  ⊗(A,B) = kron(A, B)
  Ix = sparse(I, Nxp, Nxp)
  Iz = sparse(I, Nzp, Nzp)

  ex0 = sparse([1  ], [1], [1], Nxp, 1)
  exN = sparse([Nxp], [1], [1], Nxp, 1)
  ez0 = sparse([1  ], [1], [1], Nzp, 1)
  ezN = sparse([Nzp], [1], [1], Nzp, 1)

  ex0T = sparse([1], [1  ], [1], 1, Nxp)
  exNT = sparse([1], [Nxp], [1], 1, Nxp)
  ez0T = sparse([1], [1  ], [1], 1, Nzp)
  ezNT = sparse([1], [Nzp], [1], 1, Nzp)

  # Face lift operators (their transpose gives the face extraction operators)
  e1 = Iz ⊗ ex0
  e2 = Iz ⊗ exN
  e3 = ez0 ⊗ Ix
  e4 = ezN ⊗ Ix
  e = (e1, e2, e3, e4)

  #
  # Surface quadtrature matrices
  #
  H1 = Hz
  H1I = HzI

  H2 = Hz
  H2I = HzI

  H3 = Hx
  H3I = HxI

  H4 = Hx
  H4I = HxI

  Hinv = HzI ⊗ HxI
  # Penalty terms according to A&D
  β = 1
  d = 2                 # 2-dimensional
  h1 = Hx[1,1]          # TODO: or Hx[2,2]? or Hz?
  

  
  # penalty matrices:
  Z_1 = β * (d/h1) * μ *  (Iz ⊗ Ix) # TODO: are these correct?
  Z_2 = β * (d/h1) * μ *  (Iz ⊗ Ix)


  # discrete boundary traction operators:
  T_1 = -μ * (Iz ⊗ S0x)
  T_2 = μ * (Iz ⊗ SNx)
  T_3 = -μ * (S0z ⊗ Ix)
  T_4 = μ * (SNz ⊗ Ix)  # TODO: are these correct?  Is it okay to just include boundary operation?

  T = (T_1, T_2, T_3, T_4)


  # Dirichlet on faces 1 and 2:
  SAT_1 = Hinv * (T_1 .- Z_1)' * e1 * H1 * e1'   #TODO: fix these with A&D paper
  SAT_2 = Hinv * (T_2 .- Z_2)' * e2 * H2 * e2'   #TODO: fix these with A&D paper

  # Neumann on faces 3 and 4:
  SAT_3 = -Hinv * e3 * H3 * e3' * T_3        #TODO: fix these with A&D paper
  SAT_4 = -Hinv * e4 * H4 * e4' * T_4        #TODO: fix these with A&D paper

  # boundary data operators for Dirichlet faces 1 and 2: 
  B1 = Hinv * (T_1 .- Z_1)' * e1 * H1   #TODO: fix this with A&D paper
  B2 = Hinv * (T_2 .- Z_2)' * e2 * H2   #TODO: fix this with A&D paper

  # boundary data operators for Neumann faces 3 and 4: 
  B3 = -Hinv * e3 * H3  #TODO: fix this with A&D paper
  B4 = -Hinv * e4 * H4  #TODO: fix this with A&D paper

  D2 = μ * ((Iz ⊗ D2x) + (D2z ⊗ Ix))
  A = D2 + SAT_1 + SAT_2 + SAT_3 + SAT_4  # Creation of LHS matrix

  B = (B1, B2, B3, B4)  # These will multiply the boundary data when forming the RHS.

  H̃ = Hz ⊗ Hx            # Need this for convergence tests.




  return (A, B, H̃, T, e)
end

function bdry_vec_strip!(g, B, slip_data, remote_data, free_surface_data)

    
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
  g[:] += B[3] * vf

  # FACE 4 (Neumann):
  gN = free_surface_data
  vf = gN
  g[:] += B[4] * vf 
  
end


function computetraction_stripped(T, u, e)
    e1 = e[1]
    return (-e1' * T[1]* u) 
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

  export get_operators, bdry_vec_strip!, computetraction_stripped
  export rateandstate, newtbndv