const year_seconds = 31556926
const sim_years = 1500
global const ctr = Ref{Int64}(1) 

using DifferentialEquations
using Printf

using DelimitedFiles


function odefun(dψV, ψδ, p, t)
  
  Vp = p.Vp
  M = p.M
  u = p.u
  Δτ = p.Δτ
  τf = p.τf
  g = p.g
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn = p.σn
  η = p.η
  RSV0 = p.RSV0
  τz0 = p.τz0
  RSDc = p.RSDc
  RSf0 = p.RSf0
  δNp = p.δNp
  N = p.N
  F = p.F
  τ = p.τ
  x = p.x 
  z = p.z
  HfI_FT = p.HfI_FT
  Lx = p.Lx
  Lz = p.Lz

  current_time = t ./ 31556926
  print("TIME [YRS] = $(current_time).\n")

  ψ  = @view ψδ[        (1:δNp) ]
  δ  = ψδ[ δNp .+ (1:N+1) ]
  
  bc_Dirichlet = (lf, x, z) -> (2-lf)*(δ ./ 2) + (lf-1)*fill(t .* Vp./2, size(z))
  bc_Neumann   = (lf, x, z, nx, nz) -> zeros(size(x))
  
  bdry_vec_mod!(g, F, τ, x, z, bc_Dirichlet, bc_Neumann, Lx, Lz)

  # solve for displacements everywhere in domain
  u[:] = M \ g

  # set up rates of change for  state and slip
  dψ = @view dψV[       (1:δNp) ]
  V  = @view dψV[ δNp .+ (1:N+1)]

  dψ .= 0 # initialize values to 0
  V  .= 0 # initialize values to 0

  # Update the fault data
  Δτ .= 0
  lf1 = 1  # fault is at face 1

  Δτ .= -computetraction_stripped(HfI_FT, τ, lf1, u, δ)
  τf .= τz0 .+ Δτ
  

  for n = 1:δNp
    ψn = ψ[n]
    an = RSa[n]

    τn = Δτ[n] + τz0
  
    VR = abs(τn / η)
    VL = -VR
    Vn = V[n]
    obj_rs(V) = rateandstate(V, ψn, σn, τn, η, an, RSV0)
    (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-9,
                                 atolx = 1e-9, rtolx = 1e-9)


    
    V[n] = Vn

    dψ[n] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0)
   
  end
 
  V[δNp+1:N+1] .= Vp


  nothing
end



export odefun