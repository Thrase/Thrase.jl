function create_text_files(pth, flt_loc, flt_loc_indices, stations, station_strings, station_indices, t, RSVinit, δ, τz0, θ)


    path_to_slip = pth * "slip.dat"
    # slip.dat is a file that stores time, max(V) and slip at all the stations:
    open(path_to_slip, "w") do io
      write(io,"0.0 0.0 ")
        for i in 1:length(flt_loc)
          write(io,"$(flt_loc[i]) ")
        end
          write(io,"\n")
      end
    
    #write out initial data into devol.txt:
    vv = Array{Float64}(undef, 1, 2+length(flt_loc))
      vv[1] = t
      vv[2] = log10(RSVinit)
      vv[3:end] = δ[flt_loc_indices]
      open(path_to_slip, "a") do io
          writedlm(io, vv)
      end
  
    # write out initial data into station files:
  
    # fltst_dpXXX.txt is a file that stores time and time-series of slip, log10(slip_rate), 
    # shear_stress and log10(state) at depth of z = XXX km, where XXX is each of the fault station depths.
    # First we write out initial data into each fltst_dpXXX.txt:
  
    for n = 1:length(station_strings)
      XXX = pth * "fltst_strk"*station_strings[n]*".txt"
      ww = Array{Float64}(undef, 1, 5)
      ww[1] = t
      ww[2] = δ[station_indices[n]]
      ww[3] = log10(RSVinit)
      ww[4] = τz0
      ww[5] = log10(θ[station_indices[n]])  # 
      open(XXX, "w") do io
        write(io, "# problem=SEAS Benchmark BP1-QD\n")  # 
        write(io, "# code=Thrase\n")
        write(io, "# modeler=B. A. Erickson\n")
        write(io, "# date=2023/01/09\n")
        write(io, "# element size=xx m\n")
        write(io, "# location=on fault, z = "*string(parse(Int64, station_strings[n])/10)*" km\n")
        write(io, "# Lz = 40 km\n")
        write(io, "t slip slip_rate shear_stress state\n")
  
        writedlm(io, ww)
      end
    end
  
  end
  
  function write_to_file(pth, ψδ, t, i, zf, flt_loc, flt_loc_indices, station_strings, station_indices, p, base_name="", tdump=100)
    
    path_to_slip = pth * "slip.dat"
    Vmax = 0.0
  
    if isdefined(i,:fsallast) 
      δNp = p.δNp
      Nz = p.N
      dψV = i.fsallast
      dψ = @view dψV[1:δNp]
      V = @view dψV[δNp .+ (1:Nz+1)]
      Vmax = maximum(abs.(extrema(V)))
      δ = @view ψδ[δNp .+ (1:Nz+1)]
      ψ = @view ψδ[1:δNp]
      τf = p.τf
    
   
      θ = (p.RSDc * exp.((ψ .- p.RSf0) ./ p.RSb)) / p.RSV0  # Invert ψ for θ.
    
      if mod(ctr[], p.save_stride_fields) == 0 || t == (p.sim_years ./ 31556926)
        vv = Array{Float64}(undef, 1, 2+length(flt_loc))
        vv[1] = t
        vv[2] = log10(Vmax)
        vv[3:end] = δ[flt_loc_indices]
        open(path_to_slip, "a") do io
          writedlm(io, vv)
        end
  
        for i = 1:length(station_indices)
          ww = Array{Float64}(undef, 1, 5)
          ww[1] = t
          ww[2] = δ[station_indices[i]]
          ww[3] = log10(V[station_indices[i]])
          ww[4] = τf[station_indices[i]]
          ww[5] = log10(θ[station_indices[i]])
  
          XXX = pth * "fltst_strk"*station_strings[i]*".txt"
          open(XXX, "a") do io
              writedlm(io, ww)
          end
        end
        
      end
    
      global ctr[] += 1
    end
  
    Vmax
  end
  

  function read_params_BP1_CT(f_name)
    f = open(f_name, "r")
    tmp_params = []
    while ! eof(f)
        s = readline(f)
        if s[1] != '#'
            push!(tmp_params, split(s, '=')[2])
            flush(stdout)
        end
    end
    close(f)
  
   
    #(pth, stride_space, stride_time, Lx, Lz, Hx, Hz, Nr, Ns, dx, 
   #dz, el_r, el_s, sim_years, Vp, ρ, cs, σn_0, RSamin, RSamax, 
    #RSb, RSDc, RSf0, RSV0, RSVinit, SBPp) 
  
      params = Vector{Any}(undef, 29)
      params[1] = strip(tmp_params[1])
      
      params[2] = parse(Int64, tmp_params[2])
      params[3] = parse(Int64, tmp_params[3])
      params[4] = parse(Float64, tmp_params[4])
      params[5] = parse(Float64, tmp_params[5])
      params[6] = parse(Float64, tmp_params[6])
      params[7] = parse(Float64, tmp_params[7])
      params[8] = parse(Int64, tmp_params[8])
      params[9] = parse(Int64, tmp_params[9])
      params[10] = parse(Float64, tmp_params[10])
  
      params[11] = parse(Float64, tmp_params[11])
      params[12] = parse(Float64, tmp_params[12])
      params[13] = parse(Float64, tmp_params[13])
      params[14] = parse(Float64, tmp_params[14])
  
      for i = 15:28
        params[i] = parse(Float64, tmp_params[i])
      end
      params[29] = parse(Int64, tmp_params[29])
    
    return params
  end
  
  
  
  function read_params_BP6(f_name)
    f = open(f_name, "r")
    tmp_params = []
    while ! eof(f)
        s = readline(f)
        if s[1] != '#'
            push!(tmp_params, split(s, '=')[2])
            flush(stdout)
        end
    end
    close(f)
  
   
    #(pth, stride_space, stride_time, xc, zc, Hx, Hz, Nr, Ns, dx, dz, el_r, el_s,
      #sim_years, Vp, ρ, cs, σn_0, RSa, RSb, RSD_RS,
      #RSf0, RSV0, RSVinit, RSLf, lz, μshear, τ_init, η, q_0, 
      #t_off, α, β, φ, k, η_visc, state_law, SBPp) = read_params(localARGS[1])
  
      params = Vector{Any}(undef, 35)
      params[1] = strip(tmp_params[1])
      
      params[2] = parse(Int64, tmp_params[2])
      params[3] = parse(Int64, tmp_params[3])
      params[4] = (parse(Float64, tmp_params[4]), parse(Float64, tmp_params[5]))
      params[5] = (parse(Float64, tmp_params[6]), parse(Float64, tmp_params[7]))
      params[6] = parse(Int64, tmp_params[8])
      params[7] = parse(Int64, tmp_params[9])
      params[8] = parse(Int64, tmp_params[10])
      params[9] = parse(Int64, tmp_params[11])
      for i = 12:length(tmp_params)-2
        params[i-2] = parse(Float64, tmp_params[i])
      end
      params[34] = tmp_params[36]
  
      params[35] = parse(Int64, tmp_params[37])
  
    return params
  end

  

function write_to_file_BP6(pth, ψδ, t, i, zf,flt_loc, flt_loc_indices, stations, station_indices, p, μshear, dz, base_name="", tdump=100)
  
    path_to_global = pth * "global.dat"
    path_to_slip = pth * "slip.dat"
  
    Vmax = 0.0
  
    if isdefined(i,:fsallast) 
      δNp = p.δNp
      Nz = p.Ns
      dψV = i.fsallast
      dψ = @view dψV[1:δNp]
      V = @view dψV[δNp .+ (1:Nz+1)]
      Vmax = maximum(abs.(extrema(V)))
      δ = @view ψδ[δNp .+ (1:Nz+1)]
      ψ = @view ψδ[1:δNp]
      τf = p.τf
      P = p.P
      q = p.q
      δlf = p.δlf
   
      θ = (p.RSD_RS * exp.((ψ .- p.RSf0) ./ p.RSb)) / p.RSV0  # Invert ψ for θ.
  
      # data for global.dat file
      uu = Array{Float64}(undef, 1, 3)
      uu[1] = t
      uu[2] = log10(Vmax)
      uu[3] = moment_density_rate(V, μshear, dz)
      open(path_to_global, "a") do io
        writedlm(io, uu)
      end
      
      if mod(ctr[], p.save_stride_fields) == 0 || t == (sim_years ./ 31556926)
        vv = Array{Float64}(undef, 1, 2+length(flt_loc))
        vv[1] = t
        vv[2] = log10(Vmax)
        vv[3:end] = δ[flt_loc_indices]
        open(path_to_slip, "a") do io
          writedlm(io, vv)
        end
  
        stations = ["-15", "+00", "+05", "+10", "+15", "+25", "+35", "+50", "+75"]
        
        for i = 1:length(station_indices)
          ww = Array{Float64}(undef, 1, 7)
          ww[1] = t
          ww[2] = δ[station_indices[i]]
          ww[3] = log10(V[station_indices[i]])
          ww[4] = τf[station_indices[i]]
          ww[5] = P[station_indices[i]]
          ww[6] = q[station_indices[i]]
          ww[7] = log10(θ[station_indices[i]-δlf+1])
  
          XXX = pth * "fltst_strk"*stations[i]*".txt"
          open(XXX, "a") do io
              writedlm(io, ww)
          end
        end
      end
    
        
    
      global ctr[] += 1
      @show ctr[]
    
  
    end
       Vmax
  
    
    
  end
  
  
      
  function create_text_files_BP6(pth, flt_loc, flt_loc_indices, stations, station_indices, t, RSVinit, δ, τz0, θ, δlf, P, q, μshear)
  
    path_to_global = pth * "global.dat"
    path_to_slip = pth * "slip.dat"
    # global.dat includes time series of maximum amplitude of slip rates, and moment density rates
    uu = Array{Float64}(undef, 1, 3)
    uu[1] = t
    uu[2] = log10(RSVinit)   # V = V_init everywhere
    uu[3] = μshear * RSVinit * 40 * 1e12 # constants come out of integral, int(dz) = length of RS domain = 40 km
    open(path_to_global, "w") do io
      # write(io, "# problem=SEAS Benchmark BP6-A\n")  # aging law
      write(io, "# problem=SEAS Benchmark BP6-S\n")  # slip law
      write(io, "# code=Thrase\n")
      write(io, "# modeler=J. Marcum\n")
      write(io, "# date=2022/10/20\n")
      write(io, "# element size=100 m\n")
      write(io, "# location=frictional domain\n")
      write(io, "# Column #1 = Time (s)\n")
      write(io, "# Column #2 = Max slip rate (log10 m/s)\n")
      write(io, "# Column #3 = Moment density rate (N/s)\n")
      write(io, "t max_slip_rate moment_rate\n")
      writedlm(io, uu)
    end
    
    # slip.dat is a file that stores time, max(V) and slip at all the stations:
    open(path_to_slip, "w") do io
      write(io,"0.0 0.0 ")
        for i in 1:length(flt_loc)
          write(io,"$(flt_loc[i]) ")
        end
          write(io,"\n")
      end
    
    #write out initial data into devol.txt:
    vv = Array{Float64}(undef, 1, 2+length(flt_loc))
      vv[1] = t
      vv[2] = log10(RSVinit)
      vv[3:end] = δ[flt_loc_indices]
      open(path_to_slip, "a") do io
          writedlm(io, vv)
      end
  
    # write out initial data into station files:
  
    # fltst_dpXXX.txt is a file that stores time and time-series of slip, log10(slip_rate), 
    # shear_stress and log10(state) at depth of z = XXX km, where XXX is each of the fault station depths.
    # First we write out initial data into each fltst_dpXXX.txt:
  
    stations = ["-15", "+00", "+05", "+10", "+15", "+25", "+35", "+50", "+75"]
    for n = 1:length(stations)
      XXX = pth * "fltst_strk"*stations[n]*".txt"
      ww = Array{Float64}(undef, 1, 7)
      ww[1] = t
      ww[2] = δ[station_indices[n]]
      ww[3] = log10(RSVinit)
      ww[4] = τz0
      ww[5] = P[station_indices[n]]
      ww[6] = q[station_indices[n]]
      ww[7] = log10(θ[station_indices[n]-δlf+1])  # subtract off number of points outside RS region?
      open(XXX, "w") do io
        # write(io, "# problem=SEAS Benchmark BP6-A\n")  # aging law
        write(io, "# problem=SEAS Benchmark BP6-S\n")  # slip law
        write(io, "# code=Thrase\n")
        write(io, "# modeler=J. Marcum\n")
        write(io, "# date=2022/10/16\n")
        write(io, "# element size=100 m\n")
        write(io, "# location=on fault, z = "*string(parse(Int64, stations[n])/10)*" km\n")
        write(io, "# Lz = 40 km\n")
        write(io, "t slip slip_rate shear_stress pore_pressure darcy_vel state\n")
  
        writedlm(io, ww)
      end
    end
  
  end

  # }}}

# function plot_connectivity(verts, EToV)
  
#   Lx = extrema(verts[1,:])
#   Lx = (floor(Int, Lx[1]), ceil(Int, Lx[2]))
#   Ly = extrema(verts[2,:])
#   Ly = (floor(Int, Ly[1]), ceil(Int, Ly[2]))
#   plt = Plot(BrailleCanvas(80, 40,
#                            origin_x = Lx[1], origin_y = Ly[1],
#                            width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))

#   #annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
#   #annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
#   #annotate!(plt, :bl, string(Lx[1]), color = :light_black)
#   #annotate!(plt, :br, string(Lx[2]), color = :light_black)
#   for e = 1:size(EToV, 2)
#     (v1, v2, v3, v4) = EToV[1:4, e]
#     x = verts[1, [v1 v2 v4 v3 v1]][:]
#     y = verts[2, [v1 v2 v4 v3 v1]][:]
#     lineplot!(plt, x, y)
#   end
#   title!(plt, "connectivity")
#   display(plt)
# end
