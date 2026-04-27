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
  EToDomain = p.EToDomain
  b = p.b
  u = p.u
  w = p.w
  τ = p.τ
  σ = p.σ
  Δτ = p.Δτ
  Δσ = p.Δσ
  vstarts = p.vstarts
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn = p.σn
  η = p.η
  RSV0 = p.RSV0
  τ0 = p.τ0
  σ0 = p.σ0
  RSDc = p.RSDc
  RSf0 = p.RSf0
  #fault_nodes = p.fault_nodes
  psi = p.psi

  

  nelems = length(lop)
  nfaces = size(FToE, 2)



    # Set boundary data functions:
    function creep(x,y,t, e, EToDomain)
        if EToDomain[e] == 1 # left/top hand side of fault
            return [(Vp/2) .* t * cosd(psi) .+ 0 .* x .+ 0 .* y, (Vp/2) .* t * sind(psi).+ 0 .* x .+ 0 .* y]
        elseif EToDomain[e] == 2
            return [-(Vp/2) .* t * cosd(psi).+ 0 .* x .+ 0 .* y, -(Vp/2) .* t * sind(psi).+ 0 .* x .+ 0 .* y]
        else
            error("shouldn't get here")
        end
            
    end


  bc_Dirichlet = (lf, x, y, e, δ, t, EToDomain) -> creep(x,y,t,e,EToDomain)
  bc_Neumann   = (lf, x, y, nx, ny, e, δ, t, EToDomain) -> [zeros(size(x)), zeros(size(x))]

  
  # define the interface jump in displacement function
    in_jump      = (lf, x, y, e, δ, t, EToDomain) -> begin
      f = EToF[lf, e]  # Get global face number
      if EToS[lf, e] == 1  # check if face is on minus side
        if EToO[lf, e]     # check if on minus side, A&D define δ as minus side minus plus side
          return δ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else
          error("shouldn't get here")  # this is because "correct" orientation is always true of a face on the minus side. 
        end
      else                  # face on plus side, add minus sign
        if EToO[lf, e]      # check if orientation is correct
          return  -δ[FToδstarts[f]:(FToδstarts[f+1]-1), :]
        else                # if orientation is reversed, reverse the data
          return  -δ[(FToδstarts[f+1]-1):-1:FToδstarts[f], :]
        end
      end
    end
  
  if reject_step[1]
    @show "oh great"
    return
  end

  current_time = t ./ 31556926
  #print("TIME [YRS] = $(current_time).\n")

  δNp = div(length(ψδ), 3)  # TODO: fix this to not do friction at depth

  ψ  = @view ψδ[        (1:δNp) ]
  δ1  = ψδ[  δNp .+ (1:δNp)]    # fault parallel jump
  δ2  = ψδ[  2*δNp .+ (1:δNp)]  # fault perpendicular jump
  δ = [δ1 δ2]

  
  b .= 0
  # fill in boundary data into b
  for e = 1:nelems
        loc_bdry_vec_v2!((@view b[vstarts[e]:vstarts[e+1]-1, :]), lop[e], neighborZ[e], FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, (e, δ, t, EToDomain))
  end

  # solve for displacements everywhere in domain
  U = A \ b[:]
  VNp = Integer(length(U)/2)
  u[:] = U[1:VNp]
  w[:] = U[VNp+1:2*VNp]

  # set up rates of change for state and slip
  dψ =  @view dψV[       (1:δNp) ]
  V1  = @view dψV[ δNp .+ (1:δNp) ]   # fault parallel rate
  V2  = @view dψV[ 2*δNp .+ (1:δNp) ] # fault perpendicular rate

  dψ .= 0 # initialize values to 0
  V1  .= 0 # initialize values to 0
  V2  .= 0 # initialize values to 0 

  # Update Δτ and Δσ (the shear/normal stress change due to quasi-static deformation):
  Δτ .= 0
  τ .= 0
  Δσ .= 0
  σ .= 0

  show_val = false
  for f = 1:nfaces
   
    # TODO: fix this to not do friction at depth (not a bug just doing more work than necessary)
    if FToB[f] == RS_FAULT
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      urng1 = vstarts[e1]:(vstarts[e1+1]-1)
      urng2 = vstarts[e2]:(vstarts[e2+1]-1)

      nxf1 = lop[e1].nx[lf1]
      nyf1 = lop[e1].ny[lf1]
      sJ1 = lop[e1].sJ[lf1]

      nxf2 = lop[e2].nx[lf2]
      nyf2 = lop[e2].ny[lf2]
      sJ2 = lop[e2].sJ[lf2]

    

      t1, t2 = computetraction(lop[e1], lf1, u[urng1], w[urng1])
      t3, t4 = computetraction(lop[e2], lf2, u[urng2], w[urng2]) 
 
        Δτ1 = (-nyf1) .* t1 .+ nxf1 .* t2   
        Δτ2 = (-nyf2) .* t3 .+ nxf2 .* t4
        Δσ1 = nxf1 .* t1 .+ nyf1 .* t2   
        Δσ2 = nxf2 .* t3 .+ nyf2 .* t4
      

        @views Δτ[δrng] .= (Δτ1 .+ Δτ2) ./ 2 
        @views Δσ[δrng] .= (Δσ1 .+ Δσ2) ./ 2 
    
      for n = 1:length(δrng)
        δn = δrng[n]
        τ[δn] = Δτ[δn] + τ0[δn]
        σ[δn] = Δσ[δn] + σ0[δn]

        if isnan(τ[δn])
          @show "uh oh"
          reject_step[1] = true
          return
        end
        VR = abs(τ[δn] / η)
        VL = -VR
        Vn = V1[n]
        (Vnew, ~, iter) = newtbndv((V) -> rateandstate(V, ψ[δn], σ[δn],
                                                       τ[δn], η, RSa[δn],
                                                       RSV0),
                                   VL, VR, Vn; atolx=1e-12, rtolx=1e-12,
                                   ftol=1e-12)
       
        if show_val
          show_val = false
          #@show (ψ[δn], σ[δn], τ[δn], σ[δn], η, RSa[δn], RSV0)
        end
        if isnan(Vnew) || iter < 0
         #@show (VL, VR, V[δn], Vnew, τ[δn], σ[δn], η, RSa[δn], RSV0)
          println("V reject")
          Vnew = 1e10
          reject_step[1] = true
          return
          #error()
        end
    
  
        V1[δn] = Vnew #TODO: fix to not do friction at depth
        dψ[δn] = (RSb * RSV0 / RSDc) * (exp((RSf0-ψ[δn]) / RSb) - abs(V1[δn])
                                        / RSV0)
         if !isfinite(dψ[δn])
           println("ψ reject")
           dψ[δn] = 0
           reject_step[1] = true
           return
         end

      end
      
      
    elseif FToB[f] == VP_FAULT
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]

      nxf1 = lop[e1].nx[lf1]
      nyf1 = lop[e1].ny[lf1]
      sJ1 = lop[e1].sJ[lf1]

      nxf2 = lop[e2].nx[lf2]
      nyf2 = lop[e2].ny[lf2]
      sJ2 = lop[e2].sJ[lf2]


      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      urng1 = vstarts[e1]:(vstarts[e1+1]-1)
      urng2 = vstarts[e2]:(vstarts[e2+1]-1)


      t1, t2 = computetraction(lop[e1], lf1, u[urng1], w[urng1])
      t3, t4 = computetraction(lop[e2], lf2, u[urng2], w[urng2])
    
 
        Δτ1 = (-nyf1) .* t1 .+ nxf1 .* t2   
        Δτ2 = (-nyf2) .* t3 .+ nxf2 .* t4
        Δσ1 = nxf1 .* t1 .+ nyf1 .* t2   
        Δσ2 = nxf2 .* t3 .+ nyf2 .* t4
      
       @views Δτ[δrng] .= (Δτ1 .+ Δτ2) ./ 2 
       @views Δσ[δrng] .= (Δσ1 .+ Δσ2) ./ 2 
 
      for δn = FToδstarts[f]:(FToδstarts[f+1]-1)
  
        τ[δn] = Δτ[δn] + τ0[δn]
        σ[δn] = Δσ[δn] + σ0[δn]
        V1[δn] = Vp  

      end
   
    end
    
    
         

    
  end

  nothing
end



export odefun