#using PyPlot 
using SparseArrays
using LinearAlgebra
using DelimitedFiles
using DifferentialEquations
using Interpolations
#using UnicodePlots

function interp1(xpt, ypt, x)

  knots = (xpt,) 
  itp = interpolate(knots, ypt, Gridded(Linear()))
  itp(x)  # endpoints of x must be between xpt[1] and xpt[end]
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

# Function for reading in numerical parameters 
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
    params = Vector{Any}(undef, 18)
    for i = 1:length(tmp_params)-1
      params[i] = parse(Float64, tmp_params[i])
    end
    params[18] = parse(Int64, tmp_params[18])
  return params
end


# Function for reading in numerical parameters on structured grid problems
function read_params_structured(f_name)
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
  params = Vector{Any}(undef, 25)
  params[1] = tmp_params[1]
  for i = 2:6
      params[i] = parse(Int64, tmp_params[i])
    end
   
    for i = 7:24
      params[i] = parse(Float64, tmp_params[i])
    end
    params[25] = parse(Int64, tmp_params[25])
  return params
end


# animate_slip will plot slip profiles against depth for every time step computed:
function animate_slip(S, δNp, zf, stride_time)

  m = length(zf)
  no_time_steps = size(S.t)
  slip_final = S.u[end][end]

  for i = 1:stride_time:no_time_steps[1]

    slip_t = S.u[i][δNp+1:end] # slip at time t
    #pyplot()
    display(plot(slip_t, -zf, xtickfont=font(18),
    ytickfont=font(18),
    guidefont=font(18),
    legendfont=font(18), ylabel = "Depth (km)", xlabel = "Slip (m)", xlims = (0, slip_final)))
    sleep(0.1)
  end

  #nothing
end



# {{{ Constructor for inp files
function read_inp_2d(T, S, filename::String; bc_map=1:10000)
  # {{{ Read in the file
  f = try
    open(filename)
  catch
    error("InpRead cannot open \"$filename\" ")
  end
  lines = readlines(f)
  close(f)
  # }}}

  # {{{ Read in nodes
  str = "NSET=ALLNODES"
  linenum = SeekToSubstring(lines, str);
  linenum > 0 || error("did not find: $str")
  num_nodes = 0
  for l = linenum+1:length(lines)
    occursin(r"^\s*[0-9]*\s*,.*", lines[l]) ? num_nodes+=1 : break
  end
  Vx = fill(S(NaN), num_nodes)
  Vy = fill(S(NaN), num_nodes)
  Vz = fill(S(NaN), num_nodes)
  for l = linenum .+ (1:num_nodes)
    node_data = split(lines[l], r"\s|,", keepempty=false)
    (node_num, node_x, node_y, node_z) = try
      (parse(T, node_data[1]),
       parse(S, node_data[2]),
       parse(S, node_data[3]),
       parse(S, node_data[4]))
    catch
      error("cannot parse line $l: \"$(lines[l])\" ")
    end

    Vx[node_num] = node_x
    Vy[node_num] = node_y
    Vz[node_num] = node_z
  end
  # }}}

  # {{{ Read in Elements
  str = "ELEMENT"
  linenum = SeekToSubstring(lines, str);
  num_elm = 0
  while linenum > 0
    for l = linenum .+ (1:length(lines))
      occursin(r"^\s*[0-9]*\s*,.*", lines[l]) ? num_elm+=1 : break
    end
    linenum = SeekToSubstring(lines, str; first=linenum+1)
  end
  num_elm > 0 || error("did not find any element")

  EToV = fill(T(0), 4, num_elm)
  EToBlock = fill(T(0), num_elm)
  linenum = SeekToSubstring(lines, str);
  while linenum > 0
    foo = split(lines[linenum], r"[^0-9]", keepempty=false)
    B = parse(T, foo[end])
    for l = linenum .+ (1:num_elm)
      elm_data = split(lines[l], r"\s|,", keepempty=false)
      # read into z-order
      (elm_num, elm_v1, elm_v2, elm_v4, elm_v3) = try
        (parse(T, elm_data[1]),
         parse(T, elm_data[2]),
        parse(T, elm_data[3]),
        parse(T, elm_data[4]),
        parse(T, elm_data[5]))
      catch
        break
      end
      EToV[:, elm_num] = [elm_v1, elm_v2, elm_v3, elm_v4]
      EToBlock[elm_num] = B
    end
    linenum = SeekToSubstring(lines, str; first=linenum+1)
  end
  # }}}

  # {{{ Determine connectivity
  EToF = fill(T(0), 4, num_elm)

  VsToF = Dict{Tuple{Int64, Int64}, Int64}()
  numfaces = 0
  for e = 1:num_elm
    for lf = 1:4
      if lf == 1
        Vs = (EToV[1, e], EToV[3, e])
      elseif lf == 2
        Vs = (EToV[2, e], EToV[4, e])
      elseif lf == 3
        Vs = (EToV[1, e], EToV[2, e])
      elseif lf == 4
        Vs = (EToV[3, e], EToV[4, e])
      end
      if Vs[1] > Vs[2]
        Vs = (Vs[2], Vs[1])
      end
      if haskey(VsToF, Vs)
        EToF[lf, e] = VsToF[Vs]
      else
        numfaces = numfaces + 1
        EToF[lf, e] = VsToF[Vs] = numfaces
      end
    end
  end
  #}}}

  # {{{ Read in side set info
  FToB = Array{T, 1}(undef, numfaces)
  fill!(FToB, BC_LOCKED_INTERFACE)
  linenum = SeekToSubstring(lines, "\\*ELSET")
  inp_to_zorder = [3,  2, 4, 1]
  while linenum > 0
    foo = split(lines[linenum], r"[^0-9]", keepempty=false)
    (bc, face) = try
      (parse(T, foo[1]),
       parse(T, foo[2]))
    catch
      error("cannot parse line $linenum: \"$(lines[linenum])\" ")
    end
    bc = bc_map[bc]
    face = inp_to_zorder[face]
    for l = linenum+1:length(lines)
      if !occursin(r"^\s*[0-9]+", lines[l])
        break
      end
      elms = split(lines[l], r"\s|,", keepempty=false)
      for elm in elms
        elm = try
          parse(T, elm)
        catch
          error("cannot parse line $linenum: \"$(lines[l])\" ")
        end
        if bc == 3
          bc = BC_LOCKED_INTERFACE
        end
        FToB[EToF[face, elm]] = bc
        @assert (bc == BC_DIRICHLET || bc == BC_NEUMANN ||
                 bc == BC_LOCKED_INTERFACE || bc >= BC_JUMP_INTERFACE)
      end
    end
    linenum = SeekToSubstring(lines, "\\*ELSET"; first=linenum+1)
  end
  # }}}

  ([Vx Vy]', EToV, EToF, FToB, EToBlock)
end
read_inp_2d(filename;kw...) = read_inp_2d(Int64, Float64, filename;kw...)

function SeekToSubstring(lines, substring; first=1)
  for l = first:length(lines)
    if occursin(Regex(".*$(substring).*"), lines[l])
      return l
    end
  end
  return -1
end


function plot_solution(u, nelems, vstarts, Nr, Ns, lop)

  ue = reshape(u[vstarts[1]:vstarts[2]-1], Nr[1]+1, Ns[1]+1)
  p = surface(lop[1].coord[1], lop[1].coord[2], ue, camera = (15, 35))
  for e = 2:nelems
    ue = reshape(u[vstarts[e]:vstarts[e+1]-1], Nr[e]+1, Ns[e]+1)
    surface!(lop[e].coord[1], lop[e].coord[2], ue, camera = (15, 35))
  end
  display(p)
end
    
function better_plot_solution(u, nelems, vstarts, Nr, Ns, lop)
  fig = figure()
  ax = fig.add_subplot(projection="3d")
  ue = reshape(u[vstarts[1]:vstarts[2]-1], Nr[1]+1, Ns[1]+1)
  ax.plot_surface(lop[1].coord[1], lop[1].coord[2], ue)
   for e = 2:nelems
     ue = reshape(u[vstarts[e]:vstarts[e+1]-1], Nr[e]+1, Ns[e]+1)
     ax.plot_surface(lop[e].coord[1], lop[e].coord[2], ue)
   end
  plt.show()
end


function setupfaultstations(locations, lop, FToB, FToE, FToLF, faults)
  T = eltype(locations)
  @assert size(locations, 2) == 2
  numstations = size(locations, 1)
  station_ind = zeros(Int64, numstations)
  station_face = zeros(Int64, numstations)
  nfaces = size(FToE, 2)
  for s = 1:numstations
    xs = locations[s, 1]
    ys = locations[s, 2]
    station_ind[s] = 0
    station_face[s] = 0
    d = typemax(T)
    for f = 1:nfaces
      if FToB[f] ∈ faults
        (e1, _) = FToE[:, f]
        (lf1, _) = FToLF[:, f]
        xf = lop[e1].facecoord[1][lf1]
        yf = lop[e1].facecoord[2][lf1]
        dA = hypot.(xf .- xs, yf .- ys)
        n = argmin(dA)
        if dA[n] < d
          station_ind[s] = n
          d = dA[n]
          station_face[s] = f
        end
      end
    end

    (e1, _) = FToE[:, station_face[s]]
    (lf1, _) = FToLF[:, station_face[s]]
    xf = lop[e1].facecoord[1][lf1][station_ind[s]]
    yf = lop[e1].facecoord[2][lf1][station_ind[s]]
  end
  return (ind=station_ind, face=station_face,
          xs = locations[:, 1],
          ys = locations[:, 2],
          tnext=zeros(T, 1),
          tdump=zeros(T, 1),
          t=Array{T, 1}(undef, 0),
          data=ntuple(n->(V=Array{T, 1}(undef, 0),
                          τ=Array{T, 1}(undef, 0),
                          θ=Array{T, 1}(undef, 0),
                          δ=Array{T, 1}(undef, 0)),
                      numstations))
end

function savefaultstation(ψδ, t, i, stations, FToδstarts, p, base_name="",
                         tdump=100)
  Vmax = 0.0
  
  if isdefined(i, :fsallast)
    δNp = div(length(ψδ), 2)
    dψV = i.fsallast
    V  = @view dψV[δNp .+ (1:δNp) ]
  
    #Vmax = maximum(abs.(extrema(V)))
    (Vmax, idx) = findmax(abs.(extrema(V)))
    dp = p.fault_nodes[idx]

    tlast = length(stations.t) > 0 ? stations.t[end] : -2year_seconds
    tnext = tlast + (Vmax > 1e-3 ? 0.1 : year_seconds)
    if (t >= tnext)
      ψ  = @view ψδ[        (1:δNp) ]
      δ  = @view ψδ[ δNp .+ (1:δNp) ]
      dψ = @view dψV[       (1:δNp) ]
      tlast = tnext
      #@show (t/year_seconds, Vmax)
      push!(stations.t, t)
      numstations = length(stations.data)
      for s = 1:numstations
        f = stations.face[s]
        n = stations.ind[s] + FToδstarts[f] - 1
        push!(stations.data[s].V, V[n])
        push!(stations.data[s].θ, (p.RSDc * exp((ψ[n] - p.RSf0) / p.RSb) /
                                   p.RSV0))
        push!(stations.data[s].δ, δ[n])
        push!(stations.data[s].τ, p.τ[n] - p.η * V[n])
      end
      println("took a step")
      #println(stations.t)
      if length(stations.t) == 1 ||
        ceil((stations.t[end] / tdump)) > ceil((stations.t[end-1] / tdump))
        println("saving data with basename = $base_name")
        for s = 1:numstations
          open("$(base_name)$(stations.xs[s])_$(stations.ys[s]).dat", "w") do f
            write(f, "t slip slip_rate shear_stress state\n")
            t = stations.t
            δ = stations.data[s].δ
            V = stations.data[s].V
            θ = stations.data[s].θ
            τ = stations.data[s].τ
            for n = 1:length(t)
              write(f, "$(t[n]) $(δ[n]) $(log10(abs(V[n]))) $(τ[n]) $(log10(θ[n]))\n")
            end
          end
        end
      end
    end
  end
  Vmax
end


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



  export setupfaultstations, savefaultstation
  export read_params, read_params_structured, plot_slip, plot_fault_time_series, find_ind 
  export find_station_index, read_inp_2d
  export stepcheck
  export animate_slip, interp1
  export better_plot_solution, plot_solution
  export create_text_files, write_to_file
  export read_params_BP1_CT, read_params_BP6, write_to_file_BP6, create_text_files_BP6