# On an unstructured mesh of quads:
# Solves the SEAS benchmark problem BP1-QD:
# https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf

using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using DelimitedFiles

#function main()


  ### input parameters
  (pth, stride_space, stride_time, xc, zc, Nx, Nz,
  sim_years, Vp, ρ, cs, σn, RSamin, RSamax, RSb, RSDc,
  RSf0, RSV0, RSVinit, RSH1,RSH2, RSWf, SBPp) = read_params(localARGS[1])

  try
      mkdir(pth)
  catch
      # folder already exists and data will be overwritten.
  end

  year_seconds = 31556926
  μ = cs^2 * ρ 
  μshear = cs^2 * ρ
  η = μshear / (2 * cs)


  RS_FAULT = 7
  VP_FAULT = 8
  (verts, EToV, EToF,
   FToB, EToDomain) = read_inp_2d(joinpath(@__DIR__, "meshes/BP1_v1.inp"))
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

  # Plot the original connectivity before mesh warping
  # plot_connectivity(verts, EToV)#- this needs work.

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

 # plot_blocks(lop) # needs work


  #
  # Do some assemble of the global volume operators
  #
  (M, FbarT, D, vstarts, FToλstarts) =
    LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                         (x) -> cholesky(Symmetric(x)))

  locfactors = M.F

  # Get a unique array indexes for the face to jumps map
  FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT, VP_FAULT), Nr, Ns)

  # Compute the number of volume, trace (λ), and jump (δ) points
  VNp = vstarts[nelems+1]-1
  λNp = FToλstarts[nfaces+1]-1
  δNp = FToδstarts[nfaces+1]-1

  # Build the (sparse) λ matrix using the schur complement and factor
  B = assembleλmatrix(FToλstarts, vstarts, EToF, FToB, locfactors, D, FbarT)
  BF = cholesky(Symmetric(B))

  # Assemble fault variables/data
  (bλ, λ, gδ) = (zeros(λNp), zeros(λNp), zeros(λNp))
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


end

nothing
=#