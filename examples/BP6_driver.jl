
# Solves the SEAS benchmark problem BP6
# https://strike.scec.org/cvws/seas/download/SEAS_BP6_Nov18.pdf


using Thrase
using LinearAlgebra
using DifferentialEquations
using DelimitedFiles

function main()
    
    ### input parameters
    (pth, stride_space, stride_time, xc, zc, Hx, Hz, Nr, Ns, dx, dz, el_r, el_s,
    sim_years, Vp, ρ, cs, σn_0, RSa, RSb, RSD_RS,
    RSf0, RSV0, RSVinit, RSLf, τ_init, q_0, 
    t_off, α, β, φ, k, η_visc, state_law, SBPp) = read_params_BP6(localARGS[1])

    Lz = zc[2]
    Lx = xc[2]
    lz = Lz - RSLf 
    μ = cs^2 * ρ 
    μshear = cs^2 * ρ 
    η = μshear / (2 * cs)  
    year_seconds = 31556926
    try
        mkdir(pth)
    catch
        # folder already exists and data will be overwritten.
    end

    ################################## COORDINATE TRANSFORM ###################################
    # Physical Domain: (x, z) in (0, Lx) x (-Lz, Lz)
    if dz < 2*Hz/Ns
      #for bp6-qd, need dz to be greater than 2Hz/Ns or will get errors
      print("need more grid points or increase dz\n")
      return
    end
    r_star = 1*Hx/(dx*Nr) # b/t 0 and 1; r_star = 1 is the limit of constant grid spacing
    s_star = 2*Hz/(dz*Ns) # b/t 0 and 1; r_star = 1 is the limit of constant grid spacing
    check_physical_domain(r_star, s_star, el_r, el_s, Nr, Ns)
    (Ax, Bx, Cx) = get_tanh_params(0, Lx, Hx, r_star, el_r)
    xt = (r,s) -> (Ax .* tanh.((((r .+ 1)/2) .- 1) ./ el_r) .+ Bx .* ((r .+ 1) ./ 2) .+ Cx, (Ax/(2*el_r)) .* (1 .- tanh.((((r .+ 1) ./ 2) .- 1)/el_r) .^2 ) .+ Bx/2 ,zeros(size(s)))
    (Az, Bz, Cz) = get_tanh_params(0, Lz, Hz, s_star, el_s)
    zt = (r,s) -> ((s .>= 0) .* (Az .* tanh.((s .- 1) ./ el_s) .+ Bz .* s .+ Cz) .+ (s .< 0) .* (-Az .* tanh.((-s .- 1) ./ el_s) .+ Bz .* s .- Cz) , zeros(size(r)), (s .>= 0) .* ((Az/el_s) .* (1 .- (tanh.((s .- 1) ./el_s)) .^ 2) .+ Bz) .+ (s .< 0) .* ((Az/el_s) .* (1 .- (tanh.((0 .- s .- 1) ./ el_s)) .^ 2) .+ Bz))
    metrics = create_metrics_BP6(SBPp, Nr, Ns, xt, zt) # create coordinate transform
    
    ###################################################################### 
    # create finite difference operators on computational domain:
    (M̃, F, HfI_FT, HfI_G, coord, facecoord, JH, sJ, nx, ny, Hf, HfI, τ) = get_operators_BP6(SBPp, Nr, Ns, μ, Lx, Lz; metrics)

    # and factor main matrix with Cholesky
    M = cholesky(Symmetric(M̃))    

    # initialize time and vector g that stores boundary data:
    t = 0
    g = zeros((Nr+1) * (Ns+1))

    # initialize slip:
    δ = zeros(Ns+1)

    # Initialize Pressure
    P = zeros(Ns+1)

    # Initialize Darcy Velocity
    q = zeros(Ns+1)
  
    # initialize boundary data on 4 corners of domain:
    # Dirichlet on []
    # bc_Dirichlet = (lf, x, z) -> (2-lf) .* (0 * x .+ 0 .* z) + (lf-1) .* (0 .* x .+ 0 .* z)
    function bc_Dirichlet(lf, x, z)
      if lf == 1      # left
        ans = zeros(size(z))
      elseif lf == 2  # right
        ans = zeros(size(z))
      elseif lf == 3           # top or bottom
        ans = zeros(size(x))
      elseif lf == 4
        ans = zeros(size(x))
      end
      return ans
    end


    (x, z) = metrics.coord
    
   
    # traction free (Neumann) []
    bc_Neumann   = (lf, x, z, nx, nz) -> zeros(size(x))
    # modify boundary data vector g with data from four sides:
  
    bdry_vec_mod_BP6!(g, F, τ, x, z, bc_Dirichlet, bc_Neumann, metrics)


    # Solve the linear system to get initial displacement vector u:
    u = M \ g # "\" refers to a direct solve of the system Mu = g  (i.e. kind of like an inverse)
  
    
    # create the computational vector zf that is fault coordinates:
    zf = z[1,:]   # i.e. just first row of z-variable here
  
  # Specifying RS domain

    (m, δlf) = findmin(abs.(zf .+ lz))          # index of z right at beginning of frictional domain 
  
    #z_below = zf[1:δlf]                  # below the frictional domain (z < 0)
    z_fric = zf[δlf  : Ns - δlf + 2]  # in the RS frictional domain (-lf ≤ z ≤ +lf)
    #z_above = zf[Ns - δlf + 2 : end]     # above the frictional domain (z > 0)
    δNp = length(z_fric)
  

    # initialize shear stress change:
    Δτ = zeros(Ns+1)

  
    # Calculate prestress:
    τz0 = τ_init + η*RSVinit

    # initialize the state variable θ:
    θ = ((RSD_RS ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
      sinh.((τ_init ./ (RSa .* σn_0) ))) .- RSf0 ./ RSb)) * ones(δNp)

   
    # calculate scaled state variable ψ (for numerical purposes):
    ψ0 = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSD_RS)

    # initialize vector sent to ODE [ψ, δ]:
    ψδ = zeros(δNp + Ns + 1)  # because length(ψ) = δNp = length(z_fric), length(δ) = Nz+1
    # fill vector:  
    ψδ[1:δNp] .= ψ0
    ψδ[δNp+1:δNp + Ns + 1] .= δ


    # Next part of code related to code output
    function find_station_index(stations, grid_points)
      numstations = length(stations)
      station_ind = zeros(numstations)
      for i in range(1, stop=numstations)
        station_ind[i] = argmin(abs.(grid_points .- stations[i]))
        station_ind[i]
      end
      return Integer.(station_ind)
    end

    # Stations to record on-fault time series:
    stations = [-1.5, 0, 0.5, 1, 1.5, 2.5, 3.5, 5, 7.5] # physical stations (km)

    station_indices = find_station_index(stations, zf)

    # Fault locations to record slip evolution output:
    flt_loc = -10:0.1:10 # physical stations every 100 m (units of km)
    flt_loc_indices = find_station_index(flt_loc, zf)

    u_old = copy(u)

  # set up parameters sent to right hand side of ODE:
  odeparam = (reject_step = [false],
              Vp=Vp,
              M = M,
              u=u,
              u_old=u_old,
              Δτ = Δτ,
              g = g,
              μshear=μshear,
              RSa=RSa,
              RSb=RSb,
              σn_0=σn_0,
              η=η,
              RSV0=RSV0,
              τz0=τz0,
              RSD_RS=RSD_RS,
              RSf0=RSf0,
              δNp = δNp,
              δlf = δlf,
              Ns = Ns,
              F = F,
              x = x ,
              z = z,
              zf = zf,
              τf = τz0*ones(Ns+1),
              τ = τ,
              HfI_FT = HfI_FT,
              save_stride_fields = stride_time, # save every save_stride_fields time steps
              P = P, 
              metrics = metrics,
              q_0 = q_0,
              q = q,
              t_off = t_off,
              α = α,
              β = β,
              φ = φ,
              k = k,
              η_visc = η_visc,
              sJ = sJ
             )


  tspan = (0, sim_years * year_seconds)

  # Set up the ODE problem to be solved, with rhs defined by odefun, 
  # on timespan tspan with initial condition ψδ and parameters odeparam:
  prob = ODEProblem(odefun_BP6, ψδ, tspan, odeparam)
  
  # this function gets called within rhs to enforce additional stability: 
  function stepcheck(_, p, _)
    if p.reject_step[1]
      p.reject_step[1] = false
      println("reject")
      return true
    end
    return false
  end

  #ODEresults = ODE_results([], [], [], Dict(i => [] for i = 1:length(stations)))


  # cb_fun gets called after every successful time step computed in ODE solver
  cb_fun = SavingCallback((ψδ, t, i) -> write_to_file_BP6(pth, ψδ, t, i, zf, flt_loc, flt_loc_indices,stations, station_indices, odeparam, μshear, dz, "BP1_", 0.1 * year_seconds), SavedValues(Float64, Float64))


  # Make text files to store on-fault time series and slip data,
  # Also initialize with initial data:
  create_text_files_BP6(pth, flt_loc, flt_loc_indices, stations, station_indices, 0, RSVinit, δ, τz0, θ, δlf, P, q, μshear)


  # Solve the ODE problem "prob" with Tsit5 (a Runge-Kutta method):
  sol = solve(prob, Tsit5(); dt=0.01,
              abstol = 1e-5, reltol = 1e-5, save_everystep=true, gamma = 0.5,
              internalnorm=(x, _)->norm(x, Inf), callback=cb_fun)

 
  return (sol, zf, δNp, pth)
end






# call main function: 
(S, zf, δNp, pth) = main()

# plot times series of shear stress:
plot_fault_time_series("slip", pth*"fltst_strk+10.txt")
# plot_fault_time_series("V", pth*"fltst_strk+10.txt")
# plot_fault_time_series("shear_stress", pth*"fltst_strk+10.txt")
# plot_fault_time_series("darcy_vel", pth*"fltst_strk+10.txt")
# plot_fault_time_series("state", pth*"fltst_strk+25.txt")