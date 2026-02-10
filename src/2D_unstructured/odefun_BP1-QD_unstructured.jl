
using DifferentialEquations
using Printf
using DelimitedFiles


function odefun(dψV, ψδ, p, t)
  
  RS_FAULT = 7
  VP_FAULT = 8

  reject_step = p.reject_step
  Vp = p.Vp
  A = p.A
  lop = p.lop
  neighborZ = p.neighborZ
  EToF = p.EToF
  EToS = p.EToS
  FToE = p.FToE
  FToLF = p.FToLF
  EToO = p.EToO
  FToB = p.FToB
  FToδstarts = p.FToδstarts
  b = p.b
  u = p.u
  τ = p.τ
  Δτ = p.Δτ
  vstarts = p.vstarts
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn = p.σn
  η = p.η
  RSV0 = p.RSV0
  τ0 = p.τ0
  RSDc = p.RSDc
  RSf0 = p.RSf0
  Lx = p.Lx

  

  nelems = length(lop)
  nfaces = size(FToE, 2)


  creep(x, y, t) = (x ./ Lx) .* (Vp/2) .* t
  bc_Dirichlet = (lf, x, y, e, δ, t) -> creep(x,y,t)
  bc_Neumann   = (lf, x, y, nx, ny, e, δ, t) -> zeros(size(x))
  
  in_jump      = (lf, x, y, e, δ, t) -> begin
    f = EToF[lf, e]
    if EToS[lf, e] == 1
      if EToO[lf, e]
        return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else
        return -δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
      end
    else
      if EToO[lf, e]
        return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else
        return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
      end
    end
  end

  




  if reject_step[1]
    return
  end

  current_time = t ./ 31556926
  #print("TIME [YRS] = $(current_time).\n")

  δNp = div(length(ψδ), 2)  # TODO: fix this to not do friction at depth

  ψ  = @view ψδ[        (1:δNp) ]
  δ  = ψδ[  δNp .+ (1:δNp)]

  
  b .= 0
  # fill in boundary data into b
  for e = 1:nelems
    loc_bdry_vec_v2!((@view b[vstarts[e]:vstarts[e+1]-1]), lop[e], neighborZ[e], FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, (e, δ, t))
  end


  # solve for displacements everywhere in domain
  u[:] = A \ b

  # set up rates of change for state and slip
  dψ = @view dψV[       (1:δNp) ]
  V  = @view dψV[ δNp .+ (1:δNp) ]

  dψ .= 0 # initialize values to 0
  V  .= 0 # initialize values to 0

  # Update Δτ (the shear stress due to quasi-static deformation):
  Δτ .= 0
  τ .= 0
  show_val = false
  for f = 1:nfaces
    
    # TODO: fix this to not do friction at depth
    if FToB[f] == RS_FAULT
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      urng1 = vstarts[e1]:(vstarts[e1+1]-1)
      urng2 = vstarts[e2]:(vstarts[e2+1]-1)

      y1 = computetraction(lop[e1], lf1, u[urng1])
   
      y2 = computetraction(lop[e2], lf2, u[urng2])

      #@views Δτ[δrng] .= -(y1 .+ y2) ./ 2
      @views Δτ[δrng] .= y1
    
      for n = 1:length(δrng)
        δn = δrng[n]
        τ[δn] = Δτ[δn] + τ0[δn]
   
        if isnan(τ[δn])
          reject_step[1] = true
          return
        end
        VR = abs(τ[δn] / η)
        VL = -VR
        Vn = V[n]
        (Vnew, ~, iter) = newtbndv((V) -> rateandstate(V, ψ[δn], σn,
                                                       τ[δn], η, RSa[δn],
                                                       RSV0),
                                   VL, VR, Vn; atolx=1e-12, rtolx=1e-12,
                                   ftol=1e-12)
        if show_val
          show_val = false
          @show (ψ[δn], σn[δn], τ[δn], η, RSa[δn], RSV0)
        end
        if isnan(Vnew) || iter < 0
          # @show (VL, VR, V[δn], Vnew, Tz[n], η, RSa[δn], RSV0)
          println("V reject")
          Vnew = 1e10
          reject_step[1] = true
          return
          #error()
        end
    
        V[δn] = Vnew #TODO: fix to not do friction at depth

        dψ[δn] = (RSb * RSV0 / RSDc) * (exp((RSf0-ψ[δn]) / RSb) - abs(V[δn])
                                        / RSV0)
        # if !isfinite(dψ[δn])
        #   println("ψ reject")
        #   dψ[δn] = 0
        #   reject_step[1] = true
        #   return
        # end

      end
      
    elseif FToB[f] == VP_FAULT
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      urng1 = vstarts[e1]:(vstarts[e1+1]-1)
      urng2 = vstarts[e2]:(vstarts[e2+1]-1)

      y1 = computetraction(lop[e1], lf1, u[urng1])
      y2 = computetraction(lop[e2], lf2, u[urng2])
     

      #@views Δτ[δrng] .= -(y1 .+ y2) ./ 2
      @views Δτ[δrng] .= y1
      # TODO: update τ here.
      (~, ~, ~, ~, ~, ~, ~, ~,~, ~, ~, ~, ~, ~, ~,~, ~, nx, ~, ~) = lop[e1]
      for δn = FToδstarts[f]:(FToδstarts[f+1]-1)
        V[δn] = sign(nx[lf1][1]) * Vp  
        #@show sign(nx[lf1][1])
        #V[δn] = Vp
      end
    end
    #@show norm(Δτ)
    
    
  end

  nothing
end



export odefun