
include("diagonal_sbp.jl")

using LinearAlgebra

function get_operators(p, Nx, Nz, μ; xc = (-1, 1), zc = (-1, 1))

  
  Nxp = Nx + 1
  Nzp = Nz + 1

  # "coefficient" matrices 
  cxx = μ * ones(Nxp, Nzp)
  czz = μ * ones(Nxp, Nzp)

# Derivative operators for the rest of the computation
  (Dx, HxI, Hx, x) = diagonal_sbp_D1(p, Nx; xc = xc)
  Qx = Hx * Dx
  QxT = sparse(transpose(Qx))

  (Dz, HzI, Hz, z) = diagonal_sbp_D1(p, Nz; xc = xc)
  Qz = Hz * Dz
  QzT = sparse(transpose(Qz))

  ⊗(A,B) = kron(A, B)
  Ix = sparse(I, Nxp, Nxp)
  Iz = sparse(I, Nzp, Nzp)
  (_, S0x, SNx, _, _, Ax, _) = variable_diagonal_sbp_D2(p, Nx, μ * ones(Nxp); xc = xc)
  (_, S0z, SNz, _, _, Az, _) = variable_diagonal_sbp_D2(p, Nz, μ * ones(Nzp); xc = zc)


  Ex0 = sparse([1], [1], [1], Nxp, Nxp)
  ExN = sparse([Nxp], [Nxp], [1], Nxp, Nxp)
  Ez0 = sparse([1], [1], [1], Nzp, Nzp)
  EzN = sparse([Nzp], [Nzp], [1], Nzp, Nzp)

  ex0 = sparse([1  ], [1], [1], Nxp, 1)
  exN = sparse([Nxp], [1], [1], Nxp, 1)
  ez0 = sparse([1  ], [1], [1], Nzp, 1)
  ezN = sparse([Nzp], [1], [1], Nzp, 1)

  ex0T = sparse([1], [1  ], [1], 1, Nxp)
  exNT = sparse([1], [Nxp], [1], 1, Nxp)
  ez0T = sparse([1], [1  ], [1], 1, Nzp)
  ezNT = sparse([1], [Nzp], [1], 1, Nzp)


  #
  # Surface mass matrices
  #
  H1 = Hz
  H1I = HzI

  H2 = Hz
  H2I = HzI

  H3 = Hx
  H3I = HxI

  H4 = Hx
  H4I = HxI

  #
  # Penalty terms
  #
  if p == 2
      l = 2
      β = 0.363636363
      α = 1 / 2
  elseif p == 4
      l = 4
      β = 0.2505765857
      α = 17 / 48
  elseif p == 6
      l = 7
      β = 0.1878687080
      α = 13649 / 43200
  else
      error("unknown order")
  end

  ψmin = reshape((cxx + czz - sqrt.((cxx - czz).^2)) / 2, Nxp, Nzp)
  @assert minimum(ψmin) > 0

  hx = 2 / Nx
  hz = 2 / Nz

  ψ1 = ψmin[  1, :]
  ψ2 = ψmin[Nxp, :]
  ψ3 = ψmin[:,   1]
  ψ4 = ψmin[:, Nzp]
  for k = 2:l
      ψ1 = min.(ψ1, ψmin[k, :])
      ψ2 = min.(ψ2, ψmin[Nxp+1-k, :])
      ψ3 = min.(ψ3, ψmin[:, k])
      ψ4 = min.(ψ4, ψmin[:, Nzp+1-k])
  end
  τscale = 2

  τ1 = (2τscale / hx) * (cxx[  1, :].^2 / β  / α) ./ ψ1
  τ2 = (2τscale / hx) * (cxx[Nxp, :].^2 / β  / α) ./ ψ2
  τ3 = (2τscale / hz) * (czz[:,   1].^2 / β / α) ./ ψ3
  τ4 = (2τscale / hz) * (czz[:, Nzp].^2 / β  / α) ./ ψ4

  τ1 = sparse(1:Nzp, 1:Nzp, τ1)
  τ2 = sparse(1:Nzp, 1:Nzp, τ2)
  τ3 = sparse(1:Nxp, 1:Nxp, τ3)
  τ4 = sparse(1:Nxp, 1:Nxp, τ4)

  Ãxx = Hz ⊗ Ax
  Ãzz = Az ⊗ Hx
  Ã = Ãxx + Ãzz 

  Sx0 = Hz ⊗ S0x
  SxN = Hz ⊗ SNx
  Sz0 = S0z ⊗ Hx
  SzN = SNz ⊗ Hx

  Sx0T = Hz ⊗ sparse(transpose(S0x))
  SxNT = Hz ⊗ sparse(transpose(SNx))
  Sz0T = sparse(transpose(S0z)) ⊗ Hx
  SzNT = sparse(transpose(SNz)) ⊗ Hx

  C̃1 =  (Sx0 + Sx0T) + ((τ1 * H1) ⊗ Ex0)
  C̃2 = -(SxN + SxNT) + ((τ2 * H2) ⊗ ExN)
  C̃3 =  (Sz0 + Sz0T) + (Ez0 ⊗ (τ3 * H3))
  C̃4 = -(SzN + SzNT) + (EzN ⊗ (τ4 * H4))

  # TODO: Fix minus sign (reverse of the paper)
  G1 = -(Iz ⊗ ex0T) * Sx0 
  G2 = +(Iz ⊗ exNT) * SxN 
  G3 = -(ez0T ⊗ Ix) * Sz0
  G4 = +(ezNT ⊗ Ix) * SzN

  F1 = G1' - ((τ1 * H1) ⊗ ex0)
  F2 = G2' - ((τ2 * H2) ⊗ exN)
  F3 = G3' - (ez0 ⊗ (τ3 * H3))
  F4 = G4' - (ezN ⊗ (τ4 * H4))


  HfI_F1T = H1I * G1 - (τ1 ⊗ ex0')
  HfI_F2T = H2I * G2 - (τ2 ⊗ exN')
  HfI_F3T = H3I * G3 - (ez0' ⊗ τ3)
  HfI_F4T = H4I * G4 - (ezN' ⊗ τ4)


  M̃ = Ã + C̃1 + C̃2 + C̃3 + C̃4

  H̃ = Hz ⊗ Hx
  # Modify the operator to handle the boundary conditions
  F = (F1, F2, F3, F4)
  τ = (τ1, τ2, τ3, τ4)
  HfI = (H1I, H2I, H3I, H4I)

  # Modify operators for the BC (Neumann only only faces 3 and 4)
  M̃ -= F[3] * (Diagonal(1 ./ (diag(τ[3]))) * HfI[3]) * F[3]'
  M̃ -= F[4] * (Diagonal(1 ./ (diag(τ[4]))) * HfI[4]) * F[4]'

  HfI_FT = (HfI_F1T, HfI_F2T, HfI_F3T, HfI_F4T)
  return (M̃ , F, τ, H̃, HfI_FT)
end

function bdry_vec_strip!(g, F, τ, x, z, slip_data, remote_data, free_surface_data, Lx, Lz)

    
  g[:] .= 0

  # fault (Dirichlet):
  vf = slip_data
  g[:] -= F[1] * vf

  # FACE 2 (Dirichlet):
  vf = remote_data
  g[:] -= F[2] * vf

  # FACE 3 (Neumann):
  gN = free_surface_data
  vf = gN  ./ diag(τ[3])
  g[:] -= F[3] * vf

  # FACE 4 (Neumann):
  gN = free_surface_data
  vf = gN  ./ diag(τ[4])
  g[:] -= F[4] * vf 
  

end

function bdry_vec_mod!(g, F, τ, x, z, bc_Dirichlet, bc_Neumann, Lx, Lz)

    
    g[:] .= 0

    # FACE 1 (Dirichlet):
    vf = bc_Dirichlet(1, 0, z)
    g[:] -= F[1] * vf

    # FACE 2 (Dirichlet):
    vf = bc_Dirichlet(2, Lx, z)
    g[:] -= F[2] * vf

    # FACE 3 (Neumann):
    gN = bc_Neumann(3, x, 0, 0, -1)
    vf = gN  ./ diag(τ[3])
    g[:] -= F[3] * vf

    # FACE 4 (Neumann):
    gN = bc_Neumann(4, x, Lz, 0, 1)
    vf = gN  ./ diag(τ[4])
    g[:] -= F[4] * vf 
    

end

function computetraction_stripped(HfI_FT, τ, u, δ)
    HfI_FT = HfI_FT[1]
    τf = τ[1]

    return (HfI_FT * u + τf * (δ .- δ / 2)) 
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

  export get_operators, bdry_vec_mod!, bdry_vec_strip!, computetraction_stripped
  export rateandstate, newtbndv