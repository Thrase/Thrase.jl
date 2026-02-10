# On an unstructured mesh of quads:
# Solves the SEAS benchmark problem BP1-QD:
# https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf

using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using DelimitedFiles
using Dates

const year_seconds = 31556926
global const ctr = Ref{Int64}(1) 

const RS_FAULT = 7
const VP_FAULT = 8

include("ops_BP1-QD_unstructured.jl")
include("odefun_BP1-QD_unstructured.jl")
include("../utils_2D.jl")

function antiplane_solve(Nr, Ns, SBPp, verts, EToV, FToB, EToF, FToE, FToLF, EToO, EToS, nfaces, EToDomain)

  nelems = length(Nr)
  
  ################################## COORDINATE TRANSFORM ###################################
  #
  # Build the local volume operators
  #

  # Create an empty dictionary to store the operators;
  # index via element number, return structure containing local ops.
  OPTYPE = typeof(locoperator(2, 16, 16, exact_mu))
  lop = Dict{Int64, OPTYPE}()
  
  # Indices corresponding to vertices as you move through blocks/elements
  vstarts = Array{Int64, 1}(undef, nelems + 1)
  vstarts[1] = 1   # start at 1. 

  # Loop over blocks/elements and create local operators:
  for e = 1:nelems

      Np = (Nr[e]+1)*(Ns[e]+1)  # total number of volume points on each element
      vstarts[e+1] = vstarts[e] + Np # fill in vstarts

      # DOMAIN: ordered pairs that define the physical domain
      (x1,x2,x3,x4) = verts[1, EToV[:, e]] 
      (y1,y2,y3,y4) = verts[2, EToV[:, e]]

      # Details of transfinite interpolation: straight edges:
      # x-parameterization on each face
      ex = [(α) -> x1 * (1 .- α) / 2 + x3 * (1 .+ α) / 2,
      (α) -> x2 * (1 .- α) / 2 + x4 * (1 .+ α) / 2,
      (α) -> x1 * (1 .- α) / 2 + x2 * (1 .+ α) / 2,
      (α) -> x3 * (1 .- α) / 2 + x4 * (1 .+ α) / 2]

      # derivative of x-parameterization
      exα = [(α) -> -x1 / 2 + x3 / 2,
      (α) -> -x2 / 2 + x4 / 2,
      (α) -> -x1 / 2 + x2 / 2,
      (α) -> -x3 / 2 + x4 / 2]
      
      # y-parameterization on each face
      ey = [(α) -> y1 * (1 .- α) / 2 + y3 * (1 .+ α) / 2,
      (α) -> y2 * (1 .- α) / 2 + y4 * (1 .+ α) / 2,
      (α) -> y1 * (1 .- α) / 2 + y2 * (1 .+ α) / 2,
      (α) -> y3 * (1 .- α) / 2 + y4 * (1 .+ α) / 2]

      # derivative of y-parameterization
      eyα = [(α) -> -y1 / 2 + y3 / 2,
      (α) -> -y2 / 2 + y4 / 2,
      (α) -> -y1 / 2 + y2 / 2,
      (α) -> -y3 / 2 + y4 / 2]

      # Coordinate transformation with transfinite interpolation
      xt(r,s) = transfinite_blend(ex[1], ex[2], ex[3], ex[4], exα[1], exα[2], exα[3], exα[4], r, s)
      zt(r,s) = transfinite_blend(ey[1], ey[2], ey[3], ey[4], eyα[1], eyα[2], eyα[3], eyα[4], r, s)

      # Uncomment if coordinate transform user-specified, for example:
      #xt(r,s) = (2 * r, 2 * ones(size(r)), 0 * ones(size(s)))
      #zt(r,s) = (2 * s, 0 * ones(size(r)), 2 * ones(size(s)))

      metrics = create_metrics(Nr[e], Ns[e], exact_mu, xt, zt) # create coordinate transform


      ###################################################################### 
      # create local finite difference operators on computational domain:
      lop[e] = locoperator(SBPp, Nr[e], Ns[e], exact_mu, metrics, FToB[EToF[:, e]]) 
  end


  # Assemble the global volume operator and factor:
  M = global_operator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
  M = lu(M)

  # Get a unique array indices for the faces corresponding to the fault/jump interface
  FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT), Nr, Ns)
 
  # Compute the total number of volume and fault/jump (δ) points
  VNp = vstarts[nelems+1]-1
  δNp = FToδstarts[nfaces+1]-1
  
  # define the data for the jump in displacement on fault interface
  δ = zeros(δNp)
  for f = 1:nfaces
    if FToB[f] ∈ (RS_FAULT, VP_FAULT)
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]
      (xf, yf) = lop[e1].facecoord
      @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] =
      exact(xf[lf1], yf[lf1], e2, EToDomain) - exact(xf[lf1], yf[lf1], e1, EToDomain)
      
    end
  end

  # define the data for the jump in traction data on fault interface
  θ = zeros(δNp)
  for f = 1:nfaces
      if FToB[f] ∈ (RS_FAULT, VP_FAULT)
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]
      (xf, yf) = lop[e1].facecoord
      
      tf1 = lop[e1].nx[lf1] .* exact_mu(xf[lf1], yf[lf1]) .* exact_x(xf[lf1], yf[lf1], e1, EToDomain) + lop[e1].ny[lf1] .* exact_mu(xf[lf1], yf[lf1]) .* exact_z(xf[lf1], yf[lf1], e1, EToDomain)
      tf2 = lop[e2].nx[lf2] .* exact_mu(xf[lf1], yf[lf1]) .* exact_x(xf[lf1], yf[lf1], e2, EToDomain) + lop[e2].ny[lf2] .* exact_mu(xf[lf1], yf[lf1]) .* exact_z(xf[lf1], yf[lf1], e2, EToDomain)
      
      @views θ[FToδstarts[f]:(FToδstarts[f+1]-1)] = tf1 + tf2
      end
  end

  # define the interface jump in displacement function
  in_jump      = (lf, x, y, e, δ, θ) -> begin
    f = EToF[lf, e]  # Get global face number
    if EToS[lf, e] == 1  # check if face is on minus side
      if EToO[lf, e]     # check is correct orientation
        return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else
        error("shouldn't get here")  # this is because "correct" orientation is always true of a face on the minus side. 
      end
    else                  # face on plus side
      if EToO[lf, e]      # check if orientation is correct
        return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else                # if orientation is reversed, reverse the data
        return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
      end
    end
  end

  # define the interface jump in traction fucntion
  in_tractionjump      = (lf, x, y, e, δ, θ) -> begin
    f = EToF[lf, e]  # Get global face number
    if EToS[lf, e] == 1  # check if face is on minus side
      if EToO[lf, e]     # check is correct orientation
        return θ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else
        error("shouldn't get here")  # this is bcasue "correct" orientation is always true of a face on teh minus side. 
      end
    else                  # face on plus side
      if EToO[lf, e]      # check if orientation is correct
        return  θ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else                # if orientation is reversed, reverse the data
        return  θ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
      end
    end
  end

  # Set boundary data:
  bc_Dirichlet = (lf, x, y, e, δ, θ) -> exact(x, y, e, EToDomain)
  bc_Neumann   = (lf, x, y, nx, ny, e, δ, θ) -> (nx .* exact_mu(x, y) .* exact_x(x, y, e, EToDomain)
                                              + ny .* exact_mu(x,y) .* exact_z(x, y, e, EToDomain))

  # initialize RHS array that stores boundary data (linear system will be Au = b, where b = B*g)
  b = zeros(VNp)
  uexact = zeros(VNp)
  Jf = zeros(VNp)

  
  # # modify b to incorporate BC.
  for e = 1:nelems
      LFToB = FToB[EToF[:,e]]

      neighborZ = [similar(lop[e].IsJZ[1]), similar(lop[e].IsJZ[2]), similar(lop[e].IsJZ[3]), similar(lop[e].IsJZ[4])]
          
      # loop over the faces to get the neighboring penalty parameters:
      for lf = 1:4   
          if LFToB[lf] == BC_JUMP_INTERFACE || LFToB[lf] == RS_FAULT || LFToB[lf] == VP_FAULT
              f = EToF[lf, e]  # get global face number
              (em, ep) = FToE[:, f]  # find the two elements that share global face f.
              (fm, fp) = FToLF[:, f]

              if em == e
                  eo = ep 
                  nf = fp
              else
                  eo = em 
                  nf = fm
              end

              neighborZ[nf] .= lop[eo].IsJZ[nf] # need to store the Z from the other element's local face
              
          else  
          end
      
      end
  
 
  loc_bdry_vec!((@view b[vstarts[e]:vstarts[e+1]-1]), lop[e], neighborZ, FToB[EToF[:,e]], EToF, FToE, FToLF, bc_Dirichlet, bc_Neumann, in_jump, in_tractionjump, (e, δ, θ))
  uexact[vstarts[e]:vstarts[e+1]-1] = exact(lop[e].coord[1], lop[e].coord[2], e, EToDomain)
  locsourcearray!((@view Jf[vstarts[e]:vstarts[e+1]-1]), source, lop[e].coord, lop[e].J, (e, EToDomain))

  end
  
  # solve linear system with a backsolve (factorization) to obtain displacement.

  u = M \ (b-Jf) 
  
  

  return (u, uexact, lop, vstarts)

end

function main()


  ### input parameters
  (pth, stride_space, stride_time, xc, zc, Nx, Nz,
  sim_years, Vp, ρ, cs, σn, RSamin, RSamax, RSb, RSDc,
  RSf0, RSV0, RSVinit, RSH1,RSH2, RSWf, SBPp) = read_params(localARGS[1])

  try
      mkdir(pth)
  catch
      # folder already exists so make a new one.
      pth = pth*string(now())*"/"
      mkdir(pth)
  end

  year_seconds = 31556926
  μ = cs^2 * ρ 
  μshear = cs^2 * ρ
  η = μshear / (2 * cs)


   #FAULT = 7  # jump interface has code 7
  #VP_FAULT = 8  # shouldnt' need this. Stilll a jump interface, just with diff data.
  #(verts, EToV, EToF,
   #FToB, EToDomain) = read_inp_2d(joinpath(@__DIR__, "../../meshes/BP1_v1.inp"))
  verts = [-1 0 1 -1 0 1 -1 0 1;
            -1 -1 -1 0 0 0 1 1 1];
  EToV = [1 2 4 5;
          2 3 5 6;
          4 5 7 8;
          5 6 8 9];

   EToF = [1 2 4 5;
          2 3 5 6;
          7 10 8 11;
          8 11 9 12]; 
          
    FToB = [1; 7; 1; 1; 7; 1; 2; 0; 2; 2; 0; 2]

    EToDomain = [1; 1; 1; 1]
  # EToV defines the element by its vertices
  # EToF defines element by its four faces, in global face number
  # FToB defines whether face is Dirichlet (1), Neumann (2), interior jump (7)
  #      or just an interior interface (0)
  # EToDomain is 1 if element is inside circle; 2 otherwise

  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
  @show (nelems, nfaces)

  # Plot the connectivity
  #plot_connectivity(verts, EToV)#- this needs work.

  Nr = Ns = fill(ceil(Int, 3e3 / 50), nelems) # number of grid points in each direction on each element

  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)




  OPTYPE = typeof(locoperator(2, 8, 8))
  lop = Dict{Int64, OPTYPE}()

  # Loop over blocks and create local operators
  for e = 1:nelems
    # Get the element corners
    (x1, x2, x3, x4) = verts[1, EToV[:, e]]
    (y1, y2, y3, y4) = verts[2, EToV[:, e]]

    xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
    yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

    metrics = create_metrics(Nr[e], Ns[e], μ, xt, yt)

    # Build local operators
    lop[e] = locoperator(SBPp, Nr[e], Ns[e], metrics, FToB[EToF[:, e]])
  end

  plot_blocks(lop) 


  #
  # Assemble global operators
  #
  (M, FbarT, D, vstarts, FToλstarts) =
    LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                         (x) -> cholesky(Symmetric(x)))

  locfactors = M.F

  # Get a unique array indexes for the face to jumps map
  FToδstarts = bcstarts(FToB, FToE, FToLF, FAULT, Nr, Ns)

  # Compute the number of volume and jump (δ) points
  VNp = vstarts[nelems+1]-1
  δNp = FToδstarts[nfaces+1]-1

  
  # Assemble fault variables/data

  (u, g) = (zeros(VNp), zeros(VNp))
  RSa = zeros(δNp)
  for f = 1:nfaces
    if FToB[f] ∈ (RS_FAULT, VP_FAULT)
      (e1, _) = FToE[:, f]
      (lf1, _) = FToLF[:, f]
      xf = lop[e1].facecoord[1][lf1]
      yf = lop[e1].facecoord[2][lf1]
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      for n = 1:length(δrng)
        RSa[δrng[n]] = RSamin - (RSamin - RSamax) *
          min(1, max(0, (RSH1 + yf[n])/(RSH1 - RSH2)))
      end
    end
  end

  τz0 = fill(σn * RSamax * asinh(RSVinit / (2 * RSV0) *
                                 exp.((RSf0 + RSb * log.(RSV0 / RSVinit)) /
                                      RSamax)) + η * RSVinit,
             δNp)

  θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
      sinh.((τz0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)
  ψ0 = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)

  for f = 1:nfaces
    if FToB[f] == RS_FAULT
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]
      nx = lop[e1].nx
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      for n = 1:length(δrng)
        τz0[δrng[n]] = sign(nx[lf1][n])*abs(τz0[δrng[n]])
      end
    end
  end

  τ = zeros(δNp)
  ψδ = zeros(2δNp)
  ψδ[1:δNp] .= ψ0

  odeparam = (reject_step = [false],
              Vp=Vp,
              lop=lop,
              EToF=EToF,
              EToS=EToS,
              FToE=FToE,
              FToLF=FToLF,
              EToO=EToO,
              FToB=FToB,
              FToλstarts=FToλstarts,
              FToδstarts=FToδstarts,
              gδ=gδ,
              λ=λ,
              bλ=bλ,
              u=u,
              τ=τ,
              g=g,
              vstarts=vstarts,
              BF=BF,
              FbarT=FbarT,
              locfactors=locfactors,
              μshear=μshear,
              RSa=RSa,
              RSb=RSb,
              σn=σn,
              η=η,
              RSV0=RSV0,
              τz0=τz0,
              RSDc=RSDc,
              RSf0=RSf0,
             )
  dψV = zeros(2δNp)
  tspan = (0, sim_years * year_seconds)
  prob = ODEProblem(odefun, ψδ, tspan, odeparam)
  function stepcheck(_, p, _)
    if p.reject_step[1]
      p.reject_step[1] = false
      println("reject")
      return true
    end
    return false
  end
  stations_locations = [0 0
                        0 -2.5
                        0 -5
                        0 -7.5
                        0 -10
                        0 -12.5
                        0 -15
                        0 -17.5
                        0 -20
                        0 -25
                        0 -30
                       ]
  stations = setupfaultstations(stations_locations, lop, FToB, FToE, FToLF,
                                (RS_FAULT, VP_FAULT))
  cb = SavingCallback((ψδ, t, i)->savefaultstation(ψδ, t, i, stations,
                                                   FToδstarts, odeparam,
                                                   "BP1_N_$(Nr[1])_", 10year_seconds),
                      SavedValues(Float64, Float64))
  sol = solve(prob, Tsit5(); isoutofdomain=stepcheck, dt=year_seconds,
              atol = 1e-6, rtol = 1e-3, save_everystep=false,
              internalnorm=(x, _)->norm(x, Inf), callback=cb)

  return (sol, z, δNp, pth)

end


(S, z, δNp, pth) = main();




# example of how to plot slip contours (uncomment if desired):
# plot_slip(pth*"slip.dat")

# examples of how ot plot times series of shear stress:
# plot_fault_time_series("slip", pth*"fltst_strk000.txt")
# plot_fault_time_series("slip_rate", pth*"fltst_strk000.txt")
# plot_fault_time_series("shear_stress", pth*"fltst_strk+10.txt")
# plot_fault_time_series("state", pth*"fltst_strk+25.txt")