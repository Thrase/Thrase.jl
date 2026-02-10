
using DifferentialEquations
using Printf
using DelimitedFiles


function odefun(dψV, ψδ, p, t)
  
  Vp = p.Vp
  A = p.A
  sJ = p.sJ
  u = p.u
  Δτ = p.Δτ
  τf = p.τf
  b = p.b
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn = p.σn
  η = p.η
  RSV0 = p.RSV0
  τ0 = p.τ0
  RSDc = p.RSDc
  RSf0 = p.RSf0
  δNp = p.δNp
  Ns = p.Ns
  Nr = p.Nr
  N = p.N
  B = p.B
  T = p.T
  e = p.e


  current_time = t ./ 31556926
  print("TIME [YRS] = $(current_time).\n")

  ψ  = @view ψδ[        (1:δNp) ]
  δ  = ψδ[ δNp .+ (1:Ns+1) ]
  

  # fill in initial boundary data into b
  bdry_vec!(b, B, δ ./ 2, (t .* Vp./2)*ones(Ns+1), zeros(Nr+1), zeros(Nr+1), sJ)


  # solve for displacements everywhere in domain
  u[:] = A \ b

  # set up rates of change for  state and slip
  dψ = @view dψV[       (1:δNp) ]
  V  = @view dψV[ δNp .+ (1:Ns+1)]

  dψ .= 0 # initialize values to 0
  V  .= 0 # initialize values to 0

  # Update the fault data
  Δτ .= 0

  Δτ .= computetraction_stripped(T, u, e, sJ)
  τf .= τ0 .+ Δτ
  

  for n = 1:δNp
    ψn = ψ[n]
    an = RSa[n]

    τn = Δτ[n] + τ0
  
    VR = abs(τn / η)
    VL = -VR
    Vn = V[n]
    obj_rs(V) = rateandstate(V, ψn, σn, τn, η, an, RSV0)
    (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-9,
                                 atolx = 1e-9, rtolx = 1e-9)


    
    V[n] = Vn

    dψ[n] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0)
   
  end
 
  V[δNp+1:Ns+1] .= Vp


  nothing
end


export odefun