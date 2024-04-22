
# using OrdinaryDiffEq
# using DiffEqCallbacks
# using Printf
# using Plots
# using DelimitedFiles
 using SpecialFunctions
# EXAMINE LINES ~[166], [266], AND [316] TO SWITCH BETWEEN AGING LAW AND SLIP LAW

function odefun_BP6(dψV, ψδ, p, t)
  
  # unpack structure of parameters:
  Vp = p.Vp
  M = p.M
  u = p.u
  Δτ = p.Δτ
  g = p.g
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn_0 = p.σn_0
  η = p.η
  RSV0 = p.RSV0
  τz0 = p.τz0
  RSD_RS = p.RSD_RS
  RSf0 = p.RSf0
  δNp = p.δNp
  δlf = p.δlf
  Nz = p.Ns
  F = p.F
  τf = p.τf
  τ = p.τ
  x = p.x 
  z = p.z
  zf = p.zf
  HfI_FT = p.HfI_FT
  P = p.P
  metrics = p.metrics
  q_0 = p.q_0
  q = p.q
  t_off = p.t_off
  α = p.α
  β = p.β
  φ = p.φ
  k = p.k
  η_visc = p.η_visc
  sJ = p.sJ


  # show current time:
  @show t ./ 31556926 

  # unpack vector into ψ and δ:
  ψ  = @view ψδ[        (1:δNp) ]
  δ  = ψδ[ δNp .+ (1:Nz+1) ]

  # update functions for calculating boundary data:
  function bc_Dirichlet(lf, x, z)
    if lf == 1      # left
      ans = δ ./ 2
    elseif lf == 2  # right
      ans = zeros(size(z))
    else            # top or bottom
      ans = zeros(size(x))
    end
    return ans
  end
  bc_Neumann   = (lf, x, z, nx, nz) -> zeros(size(x))
  bdry_vec_mod_BP6!(g, F, τ, x, z, bc_Dirichlet, bc_Neumann, metrics)
  # solve for displacements everywhere in domain
  u[:] = M \ g

  # set up rates of change for  state and slip and intialize to 0:
  dψ = @view dψV[       (1:δNp) ]
  V  = @view dψV[ δNp .+ (1:Nz+1)]
  dψ .= 0
  V  .= 0


  # Update the fault data
  Δτ .= 0
  lf1 = 1 # recall that face 1 (lf1) is face 1

  # define pressure and necessary functions

  function G(z, t)
      if t > 0
        f = @. exp(-z^2/(4*α*t))
        g = @. erfc(abs(z)/sqrt(4*α*t))
        return sqrt(t) * (f / sqrt(π) - (abs(z)/sqrt(4*α*t)) * g)
      else
        return 0
      end
  end

  function heaviside(t)
    if t > 0
      return 1
    elseif t < 0
      return 0
    else
      return 0.5
    end
  end 

  function Pfun(z, t)
    p = @. (q_0/(β*φ*sqrt(α))) * (G(z, t) * heaviside(t) - G(z, t-t_off) * heaviside(t-t_off))
    return p
  end 

  function dG_dz(z, t)
      if t > 0
        f = @. exp(-z^2/(4*α*t))
        df_dz = @. f * (-2*z/(4*α*t))
        g = @. erfc(abs(z)/sqrt(4*α*t))
        dg_dz = @. (-2/sqrt(π)) * f * sign(z) / sqrt(4*α*t)
        return sqrt(t) * (df_dz / sqrt(π) - (sign(z)/sqrt(4*α*t))*g - (abs(z)/sqrt(4*α*t))*dg_dz)
      else
        return 0
      end
  end 

  function dP_dz(z, t)
    dPdz = @. (q_0/(β*φ*sqrt(α))) * (dG_dz(z, t) * heaviside(t) - dG_dz(z, t-t_off) * heaviside(t-t_off))
    return dPdz
  end


  P .= Pfun(zf, t)
  q .= -1e3 * (k/η_visc) * dP_dz(zf, t)   # (scale by 1e3 to get m/s units)) similar to how τ is stored   # where is this used?? - it would be used if we solved heat equation non-analytically

  # Update effective normal stress
  σn_everywhere = σn_0 .- P

  # compute change in shear stress:
  Δτ .= -computetraction_stripped_BP6(HfI_FT, sJ, τ, lf1, u, δ)
  τf .= Δτ  .+ τz0

  # Loop over fault indices and update rates of change:
    # δNP = number of points along rate and state portion
  for n = 1:Nz+1
    if δlf <= n <= δlf+δNp-1  # "if n is an index inside RS portion"
      ψn = ψ[n-δlf+1]  # loop over all points on fault and if n is in RS portion do this
    end
    an = RSa

    τf[n] = Δτ[n] + τz0 # store for writing out to file.
    
    σn = σn_everywhere[n]

    τn = τf[n]
    VR = abs(τn / η)
    VL = -VR
    if δlf <= n <= δlf+δNp-1
      Vn = V[n]
      obj_rs(V) = rateandstate(V, ψn, σn, τn, η, an, RSV0)
      (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-9,
                                 atolx = 1e-9, rtolx = 1e-9)



      V[n] = Vn  # 
       dψ[n-δlf+1] = (RSb * RSV0 / RSD_RS) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0)  # aging law
      #dψ[n-δlf+1] = (-RSb * Vn / RSD_RS) * log((Vn / RSV0) * exp((ψn - RSf0) / RSb) )  # slip law
    end
   
   
  end

  V[1:δlf-1] .= 0 
  V[δlf+δNp:end] .= 0


  nothing
end



function moment_density_rate(V, μshear, dz)
  # trapezoidal Riemann Sum
  return  sum(μshear * 0.5 * (V[1:end-1] + V[2:end]) * dz) * 1e12
end

  

export odefun_BP6