using Plots
using SparseArrays
using LinearAlgebra
using DelimitedFiles
using DifferentialEquations
      
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
  
    if mod(ctr[], p.save_stride_fields) == 0 || t == (sim_years ./ 31556926)
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






function interp1(xpt, ypt, x)

      knots = (xpt,) 
      itp = interpolate(knots, ypt, Gridded(Linear()))
      itp[x]  # endpoints of x must be between xpt[1] and xpt[end]
end


# find_ind() differentiates b/t phases by defining
# interseismic when max slip rate < 10^-3 m/s
# mv is maximum slip rate (log10 m/s) 
function find_ind(mv)
  ind = [1]
  int = 1
  cos = 0
  for i = 2:length(mv)
    if mv[i] > -3 && int == 1 && cos == 0
      append!(ind, i);
      int = 0;
      cos = 1;
    end
  
    if mv[i] < -3 && int == 0 && cos == 1
      append!(ind, i-1)
      int = 1
      cos = 0
    end
  end


  ind = append!(ind, length(mv));  #tack on for plotting any part of an incomplete coseismic/interseismic phase
  
  return ind
end

# plot_slip will plot slip contours from devol.txt - every 5 years in blue during interseismic, 
# every 1 second in red during coseismic
function plot_slip(filename)

  grid = readdlm(filename, Float64)
  sz = size(grid)
  flt_loc = grid[1,3:end]
  T = grid[2:sz[1],1]
  maxV = grid[2:end, 2]
  slip = grid[2:sz[1], 3:sz[2]]
  N = size(slip)[2]


  ind = find_ind(maxV);        #finds indices for inter/co-seismic phases
  interval = [5*31556926 1]   #plot every 5 years and every 1 second
  
  ct = 0   #this counts the number of events


  #Assumes an initial interseismic period
  #This for-loop only plots completed phases
  for i = 1:2:length(ind)-2
    
    T1 = T[ind[i]]:interval[1]:T[ind[i+1]];

    W1 = interp1(T,slip[:,1],T1)';
    
    for j = 2:N 
      w1 = interp1(T,slip[:,j],T1)';
      W1 = [W1; w1]
    end

    if i == 1
      plot(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
    else
      plot!(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
    end

   
    T1 = T[ind[i+1]]:interval[2]:T[ind[i+2]];


    W1 = interp1(T,slip[:,1],T1)';
    for j = 2:N 
      w1 = interp1(T,slip[:,j],T1)';
      W1 = [W1; w1]
    end

    plot!(W1, -flt_loc, linecolor = :red, legend = false) #interseismic phase

    ct = ct+1;
  end

  
  # plot remainder of an incomplete interseismic period:
  i = length(ind)-1;
  T1 = T[ind[i]]:interval[1]:T[ind[i+1]];
  W1 = interp1(T,slip[:,1],T1)';
      
      for j = 2:N 
        w1 = interp1(T,slip[:,j],T1)';
        W1 = [W1; w1]
      end
      if i == 1
        plot(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
      else
        plot!(W1, -flt_loc, linecolor = :blue, legend = false) #interseismic phase
      end

      xlabel!("Cumulative Slip (m)")
      ylabel!("Depth (km)")
end


function plot_global(field, filename)

  @show filename
  grid = readdlm(filename)#, Float64)  # some elements cannot be parsed as numbers, 
                                       # a heterogeneous array of numbers and strings is returned.
  sz = size(grid)
  
  # indexing `grid` starting at row 9 to skip the header info
  T = grid[11:end, 1]  # Get time.
  T = T ./ 31556926 # convert to years.
 @show field
  if field == "maxV"
    y = grid[11:sz[1],2]
    plot(T, y)
  elseif field == "moment_rate"
    y = grid[11:sz[1],3]
    plot(T, y)
  else
    print("field not recognized")
  end
  gui()
    #return nothing
end



# plot_fault_time_series will plot field "field" from "filename".
# "field" has to be one of "slip", "V", "shear_stress", "state"
function plot_fault_time_series(field, filename)

  @show filename
  grid = readdlm(filename)#, Float64)  # some elements cannot be parsed as numbers, 
                                       # a heterogeneous array of numbers and strings is returned.
  sz = size(grid)
  
  # indexing `grid` starting at row 9 to skip the header info
  T = grid[9:end, 1]  # Get time.
  T = T ./ 31556926 # convert to years.
 @show field
  if field == "slip"
    y = grid[9:sz[1],2]
    plot(T, y)
    ylabel!("slip [m]")
  elseif field == "slip_rate"
    y = grid[9:sz[1],3]
    plot(T, y)
    ylabel!("slip rate [m/s]")
  elseif field == "shear_stress"
    y = grid[9:sz[1],4]
    plot(T, y)
    ylabel!("shear stress [MPa]")
  elseif field == "state"
    y = grid[9:sz[1],7]
    plot(T, y)
    ylabel!("state")
  else
    print("field not recognized")
  end
  xlabel!("time [yr]")
  gui()
    #return nothing
end

# Function for reading in numerical parameters for basin simulations
function read_params(f_name)
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

  # (pth, stride_space, stride_time, x1, x2, z1, z2, 
  #  Nx, Nz, sim_years, Vp, 
  #  ρ, cs, σn, RSamin, RSamax, 
  #  RSb, RSDc, RSf0,
  #  RSV0, RSVinit, RSH1, RSH2, RSWf, SBPp) = tmp_params

    params = Vector{Any}(undef, 23)
    params[1] = strip(tmp_params[1])
    
    params[2] = parse(Int64, tmp_params[2])
    params[3] = parse(Int64, tmp_params[3])
    params[4] = (parse(Float64, tmp_params[4]), parse(Float64, tmp_params[5]))
    params[5] = (parse(Float64, tmp_params[6]), parse(Float64, tmp_params[7]))
    params[6] = parse(Int64, tmp_params[8])
    params[7] = parse(Int64, tmp_params[9])
    for i = 10:length(tmp_params)-1
      params[i-2] = parse(Float64, tmp_params[i])
    end

    params[23] = parse(Int64, tmp_params[25])

  return params
end

  export read_params, plot_slip, plot_fault_time_series
  export find_station_index
  export stepcheck, create_text_files, write_to_file