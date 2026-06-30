using PyPlot 
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

  grid = readdlm(filename, Float64, skipstart=1)
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
      plt.plot(W1, flt_loc, color = "blue", ls = "none", marker = ".", markersize = "0.5"); #interseismic phase
    else
      plt.plot(W1, flt_loc, color = "blue", ls = "none", marker = ".", markersize = "0.5"); #interseismic phase
    end

   
    T1 = T[ind[i+1]]:interval[2]:T[ind[i+2]];


    W1 = interp1(T,slip[:,1],T1)';
    for j = 2:N 
      w1 = interp1(T,slip[:,j],T1)';
      W1 = [W1; w1]
    end

    plt.plot(W1, flt_loc, color = "red", ls = "none", marker = ".", markersize = "0.5"); #coseismic phase

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
        plt.plot(W1, flt_loc, color = "blue", ls = "none", marker = ".", markersize = "0.5") #interseismic phase
      else
        plt.plot(W1, flt_loc, color = "blue", ls = "none", marker = ".", markersize = "0.5") #interseismic phase
      end

      plt.xlabel("Cumulative Slip (m)");
      plt.ylabel("Depth (km)");
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
    plt.plot(T, y)
  elseif field == "moment_rate"
    y = grid[11:sz[1],3]
    plt.plot(T, y)
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
    plt.plot(T, y)
    ylabel("slip [m]")
  elseif field == "slip_rate"
    y = grid[9:sz[1],3]
    plt.plot(T, y)
    ylabel("slip rate [m/s]")
  elseif field == "shear_stress"
    y = grid[9:sz[1],4]
    plt.plot(T, y)
    ylabel("shear stress [MPa]")
  elseif field == "state"
    y = grid[9:sz[1],7]
    plt.plot(T, y)
    ylabel("state")
  else
    print("field not recognized")
  end
  xlabel("time [yr]")
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
    params = Vector{Any}(undef, 20)
  params[1] = tmp_params[1]
  params[2] = tmp_params[2]
    for i = 3:length(tmp_params)-1
      params[i] = parse(Float64, tmp_params[i])
    end
    params[20] = parse(Int64, tmp_params[20])
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

function all_perms(t)
    # get the 8 possible permutations of t 
    (t1, t2, t3, t4) = t #unpack t

    a1 = (t1, t2, t3, t4)
    a2 = (t1, t3, t2, t4)
    a3 = (t2, t1, t4, t3)
    a4 = (t2, t4, t1, t3)
    a5 = (t3, t1, t4, t2)
    a6 = (t3, t4, t1, t2)
    a7 = (t4, t3, t2, t1)
    a8 = (t4, t2, t3, t1)

    my_list = [a1, a2, a3, a4, a5, a6, a7, a8]
    return my_list
end
  
function haskey_any(myDict::Dict, keys_to_check) 
    # change key if permutation already exists
    num_keys = length(keys_to_check)
  
    for i = 1:num_keys
        if haskey(myDict, keys_to_check[i])
             my_key = keys_to_check[i]
            return my_key, true
        end
    end
    my_key = ()
    return my_key, false
end





# {{{ Constructor for inp files
function read_inp_3d(T, S, filename::String; bc_map=1:10000)
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

  EToV = fill(T(0), 8, num_elm)
  EToBlock = fill(T(0), num_elm)
  linenum = SeekToSubstring(lines, str);
  while linenum > 0
    foo = split(lines[linenum], r"[^0-9]", keepempty=false)
    B = parse(T, foo[end])
    for l = linenum .+ (1:num_elm)
      elm_data = split(lines[l], r"\s|,", keepempty=false)
      # read into z-order
      (elm_num, elm_v1, elm_v2, elm_v4, elm_v3, elm_v5, elm_v6, elm_v8, elm_v7) = try #change bak
        (parse(T, elm_data[1]),
         parse(T, elm_data[2]),
         parse(T, elm_data[3]),
         parse(T, elm_data[4]),
         parse(T, elm_data[5]),
         parse(T, elm_data[6]),
         parse(T, elm_data[7]),
         parse(T, elm_data[8]),
         parse(T, elm_data[9]))
      catch
        break
      end
      EToV[:, elm_num] = [elm_v1, elm_v2, elm_v3, elm_v4, elm_v5, elm_v6, elm_v7, elm_v8]
      @show EToV
      EToBlock[elm_num] = B
    end
    linenum = SeekToSubstring(lines, str; first=linenum+1)
  end
  # }}}

  # {{{ Determine connectivity
  EToF = fill(T(0), 6, num_elm)

  VsToF = Dict{Tuple{Int64, Int64, Int64, Int64}, Int64}()
  numfaces = 0
  for e = 1:num_elm
    for lf = 1:6
      if lf == 1
        Vs = (EToV[1, e], EToV[5, e], EToV[3, e], EToV[7, e])
      elseif lf == 2
        Vs = (EToV[2, e], EToV[6, e], EToV[4, e], EToV[8, e])
      elseif lf == 3
        Vs = (EToV[1, e], EToV[2, e], EToV[3, e], EToV[4, e]) 
      elseif lf == 4
       Vs = (EToV[5, e], EToV[6, e], EToV[7, e], EToV[8, e])
      elseif lf == 5
        Vs = (EToV[1, e], EToV[2, e], EToV[5, e], EToV[6, e])
      elseif lf == 6
        Vs = (EToV[3, e], EToV[4, e], EToV[7, e], EToV[8, e])
      end

    Vsperms = all_perms(Vs)
    newkey, mybool = haskey_any(VsToF, Vsperms)
    if mybool
        Vs = newkey  # overwrite Vs is face has already been added; orientation won't be correct
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
  @show numfaces
  @show FToB
  fill!(FToB, BC_LOCKED_INTERFACE)
  linenum = SeekToSubstring(lines, "\\*ELSET")
  
  inp_to_zorder = [3, 4, 5, 2, 6, 1] # because .inp file face ordering different 
  while linenum > 0
    foo = split(lines[linenum], r"[^0-9]", keepempty=false)
    @show foo
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

  ([Vx Vy Vz]', EToV, EToF, FToB, EToBlock)
end
read_inp_3d(filename;kw...) = read_inp_3d(Int64, Float64, filename;kw...)

function SeekToSubstring(lines, substring; first=1)
  for l = first:length(lines)
    if occursin(Regex(".*$(substring).*"), lines[l])
      return l
    end
  end
  return -1
end

function normalize(t2, t1)

    n = [0, 0, 0, 0]
    for i = 1:4
        value_to_find = t2[i]
        idx = findfirst(==(value_to_find), t1)
        n[i] = idx
    end

    return n
end

function orient(t1::Tuple, t2::Tuple)
    # determine face permutation
    # return 1 if faces have same orientation

    (A1, B1, C1, D1) = (1, 2, 3, 4)  
    (A2, B2, C2, D2) = normalize(t2, t1)  # convert t2 wrt t1

    @show (A1, B1, C1, D1)
    @show (A2, B2, C2, D2)

    if (A1, B1, C1, D1) == (A2, B2, C2, D2)
        return 1  # No permuation
    elseif (A1, B1, C1, D1) == (B2, A2, D2, C2)
        return  2  # flip 2143 to 1234 (H)
    elseif (A1, B1, C1, D1) == (C2, D2, A2, B2)
        return 3   # flip 3412 to 1234 (V)
    elseif (A1, B1, C1, D1) == (D2, C2, B2, A2)
        return 4    # permute 4321 to 1234 (VH)
    elseif (A1, B1, C1, D1) == (A2, C2, B2, D2)
        return 5    # permute 1324 to 1234 (RH)
    elseif (A1, B1, C1, D1) == (B2, D2, A2, C2)
        return 6    # permute 2413 to 1234 (-R)
    elseif (A1, B1, C1, D1) == (C2, A2, D2, B2)
        return 7    # permute 3142 to 1234 (+R)
    elseif (A1, B1, C1, D1) == (D2, B2, C2, A2)
        return 8    # permute 4231 to 1234 (RV)
    else
        error("problem with connectivity")
    end 
end



#{{{ connectivityarrays
function connectivityarrays_3d(EToV, EToF)
  # number of elements
  nelems = size(EToV, 2)
  nfaces = maximum(maximum(EToF))

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number (elements that share the i'th column/face)
  # FToLF: Unique Global Face to Element local face number (ith column stores local face numbers of elements in FToE[i] sharing global face i)
  # EToO : Element to Unique Global Faces Orientation  (the i'th column, gf row of this stores gfo,  whether the i'th element's local face (gf) is oriented the 
            #same way as the how the adjacent element stores it). If 0, items need to be rotated
  # EToS : Element to Unique Global Face Side (the i'th column of this stores whether an element face is on the
    #        plus side or minus side of the global face)

  FToE  = zeros(Int64, 2, nfaces)
  FToLF = zeros(Int64, 2, nfaces)
  EToO  = Array{Int64,2}(undef, 6, nelems)
  EToS  = zeros(Int64, 6, nelems)


  # Local Face to Local Vertex map
  LFToLV = flatten_tuples(((1, 5, 3, 7), (2, 6, 4, 8), (1,2,3, 4), (5, 6, 7, 8), (1, 2, 5, 6), (3, 4, 7, 8)))
  for e = 1:nelems
    for lf = 1:6
      gf = EToF[lf, e]
      if FToE[1, gf] == 0
        @assert FToLF[1, gf] == 0
        FToE[1, gf] = e
        FToLF[1, gf] = lf
        EToO[lf, e] = 1
        EToS[lf, e] = 1
      else
        @assert FToE[2, gf] == 0
        @assert FToLF[2, gf] == 0
        FToE[2, gf] = e
        FToLF[2, gf] = lf
        EToS[lf, e] = 2

        ne = FToE[1, gf]  #get other element number that shares global face gf
        nf = FToLF[1, gf] # get the corresponding local face number of element ne for gf 

        nv = EToV[LFToLV[:,nf], ne] # get global vertex numbers nv from local 
        lv = EToV[LFToLV[:,lf], e] # check to see if new one (lv) matches older (nv)

        ov = orient(Tuple(nv), Tuple(lv))
        EToO[lf, e] = ov
 
      end
    end
  end
  (FToE, FToLF, EToO, EToS)
end
#}}}


#{{{ connectivityarrays
function connectivityarrays_2d(EToV, EToF)
  # number of elements
  nelems = size(EToV, 2)
  nfaces = maximum(maximum(EToF))

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side

  FToE  = zeros(Int64, 2, nfaces)
  FToLF = zeros(Int64, 2, nfaces)
  EToO  = Array{Bool,2}(undef, 4, nelems)
  EToS  = zeros(Int64, 4, nelems)

  # Local Face to Local Vertex map
  LFToLV = flatten_tuples(((1,3), (2, 4), (1,2), (3,4)))
  for e = 1:nelems
    for lf = 1:4
      gf = EToF[lf, e]
      if FToE[1, gf] == 0
        @assert FToLF[1, gf] == 0
        FToE[1, gf] = e
        FToLF[1, gf] = lf
        EToO[lf, e] = true
        EToS[lf, e] = 1
      else
        @assert FToE[2, gf] == 0
        @assert FToLF[2, gf] == 0
        FToE[2, gf] = e
        FToLF[2, gf] = lf
        EToS[lf, e] = 2

        ne = FToE[1, gf]  #get other element number that shares global face gf
        nf = FToLF[1, gf] # get the corresponding local face number of element ne for gf 

        nv = EToV[LFToLV[:,nf], ne] # 
        lv = EToV[LFToLV[:,lf], e]
        
        if nv == lv
          EToO[lf, e] = true
        elseif nv[end:-1:1] == lv
          EToO[lf, e] = false
        else
          error("problem with connectivity")
        end
      end
    end
  end
  (FToE, FToLF, EToO, EToS)
end
#}}}



# flatten tuples to arrays
flatten_tuples(x) = reshape(collect(Iterators.flatten(x)), length(x[1]),
                            length(x))


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


function setupfaultcoord(lop, FToB, FToE, FToLF, faults)
  T = Float64
  fault_coords = Matrix{Float64}(undef, 0, 2) # intialize with length 0; will be appended

  nfaces = size(FToE, 2) # total number of faces
  

  for f = 1:nfaces # loop over all faces
      if FToB[f] ∈ faults # determine if face is a fault
        (e1, _) = FToE[:, f] # get element on minus side of fault
        (lf1, _) = FToLF[:, f] # get element's local face
        xf = lop[e1].facecoord[1][lf1] # get physical coordinates of element's local face
        yf = lop[e1].facecoord[2][lf1]
        fault_coords = [fault_coords; [xf yf]]
      end
  end


  return (coords = fault_coords, Vmax = Array{T, 1}(undef, 0), slip = Vector{Vector{T}}(undef, 0), 
          t=Array{T, 1}(undef, 0))
end




function setupfaultstations(locations, lop, FToB, FToE, FToLF, faults)
  T = eltype(locations) # should be Float64 
  @assert size(locations, 2) == 2
  numstations = size(locations, 1) # total number of stations
  station_ind = zeros(Int64, numstations)
  station_face = zeros(Int64, numstations)
  nfaces = size(FToE, 2) # total number of faces
  for s = 1:numstations
    xs = locations[s, 1]
    ys = locations[s, 2] #location specified
    station_ind[s] = 0
    station_face[s] = 0
    d = typemax(T)
    for f = 1:nfaces # loop over all faces
      if FToB[f] ∈ faults # determine if face is a fault
        (e1, _) = FToE[:, f] # get element on minus side of fault
        (lf1, _) = FToLF[:, f] # get element's local face
        xf = lop[e1].facecoord[1][lf1] # get physical coordinates of element's local face
        yf = lop[e1].facecoord[2][lf1]
        dA = hypot.(xf .- xs, yf .- ys)
        n = argmin(dA) # find the index corresponding to the physical point on fault closest to station (this will only be a local index for the station). 
        if dA[n] < d
          station_ind[s] = n # store the local index of station
          d = dA[n]
          station_face[s] = f # store the face
        end
      end
    end

    (e1, _) = FToE[:, station_face[s]] # find element on minus side of station face
    (lf1, _) = FToLF[:, station_face[s]] # find local face
    xf = lop[e1].facecoord[1][lf1][station_ind[s]] # single coordinate corresponding to station
    yf = lop[e1].facecoord[2][lf1][station_ind[s]] # single coordinate corresponding to station
  end
  return (ind=station_ind, face=station_face,
          xs = locations[:, 1],
          ys = locations[:, 2],
          tnext=zeros(T, 1),
          tdump=zeros(T, 1),
          t=Array{T, 1}(undef, 0),
          data=ntuple(n->(V=Array{T, 1}(undef, 0), #initialize as length 0, since will be appended
                          τ=Array{T, 1}(undef, 0),
                          θ=Array{T, 1}(undef, 0),
                          δ=Array{T, 1}(undef, 0)),
                      numstations))
end

function savedatafields(ψδ, t, i, stations, fault, FToδstarts, p, base_name="", slipbase_name = "",
                         tdump=100)
  Vmax = 0.0
  T = Float64

  if isdefined(i, :fsallast)
    δNp = div(length(ψδ), 2)
    dψV = i.fsallast
    
    V  = @view dψV[δNp .+ (1:δNp) ]
    
  
    Vmax = maximum(abs.(extrema(V)))
    Vmax = Vmax[1]

    tlast = length(stations.t) > 0 ? stations.t[end] : -2year_seconds # if stations.t not empty, set tlast = stations.t[end]; otherwise set tlast = -2year_seconds
    tnext = tlast + (Vmax > 1e-3 ? 0.1 : year_seconds) # if Vmax > 1e-3 then add 0.1 to tlast; otherwise add year_seconds

    # the following is going to append lists, but not write them to file yet. 
    if (t >= tnext) #check if the time step taken is bigger that 0.1 s (coseismic) or 1 year (aseismic). If it is, then save data by appending lists.
      ψ  = @view ψδ[        (1:δNp) ]
      δ  = @view ψδ[ δNp .+ (1:δNp) ]
      dψ = @view dψV[       (1:δNp) ]

      tlast = tnext
      #@show (t/year_seconds, Vmax)
      
      push!(fault.t, t)
      push!(fault.slip, δ)
      push!(fault.Vmax, Vmax)

      push!(stations.t, t)
      numstations = length(stations.data)
      for s = 1:numstations
        f = stations.face[s] # get face station s is on. 
        n = stations.ind[s] + FToδstarts[f] - 1 # calculate global station index.
        push!(stations.data[s].V, V[n]) # append slip rate at station s
        push!(stations.data[s].θ, (p.RSDc * exp((ψ[n] - p.RSf0) / p.RSb) /
                                   p.RSV0))
        push!(stations.data[s].δ, δ[n])
        push!(stations.data[s].τ, p.τ[n] - p.η * V[n])
      end
      println("took a step")
     
      if t == 0
        println("saving data with basename = $base_name")

        open("$(slipbase_name).dat", "w") do f # This will overwrite the file everytime (it's not appending it!)
          write(f, "z t max_slip_rate slip\n")

          z = fault.coords[:, 2]
         # sort!(z)
     
          
          t = fault.t 
          Vmaxvec = fault.Vmax
          slip = fault.slip[:]
         
          for i = 1:2
            write(f, "$(0) ")
          end
          for k = 1:length(z)
            write(f, "$(z[k]) ")
          end
          write(f, "\n")
          for n = 1:length(t)
            write(f, "$(t[n]) $(log10(abs(Vmaxvec[n]))) ")
            for k = 1:length(slip[1])
              write(f, "$(slip[n][k]) ")
            end
            write(f, "\n")
          end
        end

        for s = 1:numstations
          open("$(base_name)$(stations.xs[s])_$(stations.ys[s]).dat", "w") do f # This will overwrite the file everytime (it's not appending it!)
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
        #NOW reset stations and fault to length 0 arrays:
        empty!(stations.t)
        for s = 1:numstations
          empty!(stations.data[s].V)
          empty!(stations.data[s].δ)
          empty!(stations.data[s].τ)
          empty!(stations.data[s].θ)
        end
        empty!(fault.Vmax)
        empty!(fault.slip)
        empty!(fault.t)
        


      
      elseif length(stations.t) > 1 &&
        ceil((stations.t[end] / tdump)) > ceil((stations.t[end-1] / tdump)) # if length(t) == 1 or if tdump = 10years have passed, dump/append data disk and reset lists to length 0
        println("saving data with basename = $base_name")

        open("$(slipbase_name).dat", "a") do f # This will overwrite the file everytime (it's not appending it!)
          
          t = fault.t 
          Vmaxvec = fault.Vmax
          slip = fault.slip[:]
         
          for n = 1:length(t)
            write(f, "$(t[n]) $(log10(abs(Vmaxvec[n]))) ")
            for k = 1:length(slip[1])
              write(f, "$(slip[n][k]) ")
            end
            write(f, "\n")
          end
        end

        for s = 1:numstations
          open("$(base_name)$(stations.xs[s])_$(stations.ys[s]).dat", "a") do f # This will overwrite the file everytime (it's not appending it!)
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

        # reset stations and fault to length 0 arrays:
        empty!(stations.t)
        for s = 1:numstations
          empty!(stations.data[s].V)
          empty!(stations.data[s].δ)
          empty!(stations.data[s].τ)
          empty!(stations.data[s].θ)
        end
        empty!(fault.Vmax)
        empty!(fault.slip)
        empty!(fault.t)
        
        
    
      else
      end
    end
  end

  Vmax
end





function savefaultstation(ψδ, t, i, stations, FToδstarts, p, base_name="",
                         tdump=100)
  Vmax = 0.0

  if isdefined(i, :fsallast)
    δNp = div(length(ψδ), 2)
    dψV = i.fsallast
    V  = @view dψV[δNp .+ (1:δNp) ]
    
  
    Vmax = maximum(abs.(extrema(V)))
    #(Vmax, idx) = findmax(abs.(extrema(V)))
    #dp = p.fault_nodes[idx]

    tlast = length(stations.t) > 0 ? stations.t[end] : -2year_seconds # if stations.t not empty, set tlast = stations.t[end]; otherwise set tlast = -2year_seconds
    tnext = tlast + (Vmax > 1e-3 ? 0.1 : year_seconds) # if Vmax > 1e-3 then add 0.1 to tlast; otherwise add year_seconds
    # the following is going to append lists, but not write them to file yet. 
    if (t >= tnext)
      ψ  = @view ψδ[        (1:δNp) ]
      δ  = @view ψδ[ δNp .+ (1:δNp) ]
      dψ = @view dψV[       (1:δNp) ]
      tlast = tnext
      #@show (t/year_seconds, Vmax)
      push!(stations.t, t)
      
      numstations = length(stations.data)
      for s = 1:numstations
        f = stations.face[s] # get face station s is on. 
        n = stations.ind[s] + FToδstarts[f] - 1 # calculate global station index.
        push!(stations.data[s].V, V[n]) # append slip rate at station s
        push!(stations.data[s].θ, (p.RSDc * exp((ψ[n] - p.RSf0) / p.RSb) /
                                   p.RSV0))
        push!(stations.data[s].δ, δ[n])
        push!(stations.data[s].τ, p.τ[n] - p.η * V[n])
      end
      println("took a step")
      #println(stations.t)
      # if length(t) == 1 or if tdump = 10years have passed, dump data to file
      if length(stations.t) == 1 ||
        ceil((stations.t[end] / tdump)) > ceil((stations.t[end-1] / tdump))
        println("saving data with basename = $base_name")
        for s = 1:numstations
          open("$(base_name)$(stations.xs[s])_$(stations.ys[s]).dat", "w") do f # This will overwrite the file everytime (it's not appending it!)
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



function create_metrics_3D(pm, Nq, Nr, Ns, λ, μ, K, B_p,
                        xf=(q,r,s)->(q, ones(size(q)), zeros(size(r)), zeros(size(s))),
                        yf=(q,r,s)->(r, zeros(size(q)), ones(size(r)), zeros(size(s))),
                        zf=(q,r,s)->(s, zeros(size(q)), zeros(size(r)), ones(size(s))))
  
  Nqp = Nq + 1
  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nqp * Nrp * Nsp

  # Derivative operators for the metric terms
  @assert pm <= 8
  pp = pm == 6 ? 8 : pm

  q = range(-1, stop=1, length=Nqp)
  r = range(-1, stop=1, length=Nrp)
  s = range(-1, stop=1, length=Nsp)

  # Create the mesh

  q = (q * ones(1, Nrp)) .* reshape(ones(Nsp), 1, 1, :)
  r = (ones(Nqp) * r') .* reshape(ones(Nsp), 1, 1, :)
  s = (ones(Nqp) * ones(1, Nrp)) .* reshape(s, 1, 1, :) # possible source for error? 
  
  (x, xq, xr, xs) = xf(q, r, s)

  (y, yq, yr, ys) = yf(q, r, s)
  (z, zq, zr, zs) = zf(q, r, s)

  J = xq .* (yr .* zs - zr .* ys) - xr .* (yq .* zs - zq .* ys) + xs .* (yq .* zr - zq .* yr)
  @show extrema(J)

  @assert minimum(J) >= 1e-12

  qx =  (yr .* zs - ys .* zr) ./ J 
  qy = -(xr .* zs - xs .* zr) ./ J 
  qz =  (xr .* ys - xs .* yr) ./ J 
  rx =  (yq .* zs - ys .* zq) ./ J
  ry =  (xq .* zs - xs .* zq) ./ J
  rz = -(xq .* ys - xs .* yq) ./ J
  sx =  (yq .* zr - yr .* zq) ./ J
  sy = -(xq .* zr - xr .* zq) ./ J
  sz =  (xq .* yr - xr .* yq) ./ J
  qrs = (qx, qy, qz, rx, ry, rz, sx, sy, sz)

 
  # variable coefficient matrix components
  #C = fill(Array{Float64, 3}(undef, 3, 4, 5), 3, 3)
  C = fill(zeros(Nqp, Nrp, Nsp), 3, 3, 3, 3)
  for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                C[i, j, k, l] = J .* (λ(x, y, z, B_p) .* F(j, i, qrs) .* F(l, k, qrs) + 
                                   μ(x, y, z, B_p) .* (F(1, i, qrs) .* delt(j, l) .* F(1, k, qrs) + F(2, i, qrs) .* delt(j, l) .* F(2, k, qrs) + F(3, i, qrs) .* delt(j, l) .* F(3, k, qrs)) +
                                   μ(x, y, z, B_p) .* (F(l, i, qrs) .* F(j, k, qrs)))
            end
        end
    end
  end


  #
  # Block surface matrices
  #
  (xf1, yf1, zf1) = (view(x, 1, :, :), view(y, 1, :, :), view(z, 1, :, :))
  nx1 =   ys[1, :, :] .* zr[1, :, :] - yr[1, :, :] .* zs[1, :, :]
  ny1 = -(xs[1, :, :] .* zr[1, :, :] - xr[1, :, :] .* zs[1, :, :])
  nz1 =   xs[1, :, :] .* yr[1, :, :] - xr[1, :, :] .* ys[1, :, :]
  sJ1 = hypot.(nx1, ny1, nz1)
  nx1 = nx1 ./ sJ1
  ny1 = ny1 ./ sJ1
  nz1 = nz1 ./ sJ1


  (xf2, yf2, zf2) = (view(x, Nqp, :, :), view(y, Nqp, :, :), view(z, Nqp, :, :))
  nx2 =    yr[end, :, :] .* zs[end, :, :] - ys[end, :, :] .* zr[end, :, :]
  ny2 =  -(xr[end, :, :] .* zs[end, :, :] - xs[end, :, :] .* zr[end, :, :])
  nz2 =    xr[end, :, :] .* ys[end, :, :] - xs[end, :, :] .* yr[end, :, :]
  sJ2 = hypot.(nx2, ny2, nz2)
  nx2 = nx2 ./ sJ2
  ny2 = ny2 ./ sJ2
  nz2 = nz2 ./ sJ2

  (xf3, yf3, zf3) = (view(x, :, 1, :), view(y, :, 1, :), view(z, :, 1, :))
  nx3 =   yq[:, 1, :] .* zs[:, 1, :] - ys[:, 1, :] .* zq[:, 1, :]
  ny3 = -(xq[:, 1, :] .* zs[:, 1, :] - xs[:, 1, :] .* zq[:, 1, :])
  nz3 =   xq[:, 1, :] .* ys[:, 1, :] - xs[:, 1, :] .* yq[:, 1, :]
  sJ3 = hypot.(nx3, ny3, nz3)
  nx3 = nx3 ./ sJ3
  ny3 = ny3 ./ sJ3
  nz3 = nz3 ./ sJ3


  (xf4, yf4, zf4) = (view(x, :, Nrp, :), view(y, :, Nrp, :), view(z, :, Nrp, :))
  nx4 =    ys[:, end, :] .* zq[:, end, :] - yq[:, end, :] .* zs[:, end, :]
  ny4 =  -(xs[:, end, :] .* zq[:, end, :] - xq[:, end, :] .* zs[:, end, :])
  nz4 =    xs[:, end, :] .* yq[:, end, :] - xq[:, end, :] .* ys[:, end, :]
  sJ4 = hypot.(nx4, ny4, nz4)
  nx4 = nx4 ./ sJ4
  ny4 = ny4 ./ sJ4
  nz4 = nz4 ./ sJ4

  (xf5, yf5, zf5) = (view(x, :, :, 1), view(y, :, :, 1), view(z, :, :, 1))
  nx5 =   yr[:, :, 1] .* zq[:, :, 1] - yq[:, :, 1] .* zr[:, :, 1]
  ny5 = -(xr[:, :, 1] .* zq[:, :, 1] - xq[:, :, 1] .* zr[:, :, 1])
  nz5 =  xr[:, :, 1] .* yq[:, :, 1] - xq[:, :, 1] .* yr[:, :, 1]
  sJ5 = hypot.(nx5, ny5, nz5)
  nx5 = nx5 ./ sJ5
  ny5 = ny5 ./ sJ5
  nz5 = nz5 ./ sJ5

  (xf6, yf6, zf6) = (view(x, :, :, Nsp), view(y, :, :, Nsp), view(z, :, :, Nsp))
  nx6 =   yq[:, :, end] .* zr[:, :, end] - yr[:, :, end] .* zq[:, :, end]
  ny6 = -(xq[:, :, end] .* zr[:, :, end] - xr[:, :, end] .* zq[:, :, end])
  nz6 =   xq[:, :, end] .* yr[:, :, end] - xr[:, :, end] .* yq[:, :, end]
 
  sJ6 = hypot.(nx6, ny6, nz6)
  nx6 = nx6 ./ sJ6
  ny6 = ny6 ./ sJ6
  nz6 = nz6 ./ sJ6


  (coord = (x,y, z),
   facecoord = ((xf1, xf2, xf3, xf4, xf5, xf6), (yf1, yf2, yf3, yf4, yf5, yf6), (zf1, zf2, zf3, zf4, zf5, zf6)),
   C = C,
   J=J,
   sJ = (sJ1, sJ2, sJ3, sJ4, sJ5, sJ6),
   nx = (nx1, nx2, nx3, nx4, nx5, nx6),
   ny = (ny1, ny2, ny3, ny4, ny5, ny6),
   nz = (nz1, nz2, nz3, nz4, nz5, nz6),
   qx = qx, qy = qy, qz = qz, 
   rx = rx, ry = ry, rz = rz,
   sx = sx, sy = sy, sz = sz)
end


function transforms_ne(Lw, el_x, el_y, el_z)
    
    xt = (q, r, s) -> (el_x .* tan.(atan((Lw)/el_x).* (0.5*q .+ 0.5)),
                   el_x .* sec.(atan((Lw)/el_x).* (0.5*q .+ 0.5)).^2 * atan((Lw)/el_x) * 0.5 ,
                   zeros(size(r)), zeros(size(s)))
    
    yt = (q, r, s) -> (el_y .* tan.(atan((Lw)/el_y).* (0.5*r .+ 0.5)),
                   zeros(size(q)),
                   el_y .* sec.(atan((Lw)/el_y).*(0.5*r .+ 0.5)) .^2 * atan((Lw)/el_y) * 0.5, 
                   zeros(size(s)))

    zt = (q, r, s) -> (el_z .* tan.(atan((Lw)/el_z).* (0.5*s .+ 0.5)),
                   zeros(size(q)),
                   zeros(size(r)),
                   el_z .* sec.(atan((Lw)/el_z).*(0.5*s .+ 0.5)) .^2 * atan((Lw)/el_z) * 0.5)


    return xt, yt, zt
end

function delt(i, j)
    if i == j 
        return 1
    else
        return 0
    end
end

function F(I, i, qrs)
    (qx, qy, qz, rx, ry, rz, sx, sy, sz) = qrs
    # F(I, i) = dxi/dXI, e.g. F(1, 1) = qx
    if i == 1
        if I == 1
            return qx
        elseif I == 2
            return qy
        else I == 3
            return qz
        end
    elseif i == 2
        if I == 1
            return rx
        elseif I == 2
            return ry
        else I == 3
            return rz
        end
    else i == 3
        if I == 1
            return sx
        elseif I == 2
            return sy
        else I == 3
            return sz
        end
    end
end



function locoperator_3D(p, Nq, Nr, Ns, metrics, C)
    
    Nqp = Nq + 1
    Nrp = Nr + 1
    Nsp = Ns + 1
    Np = Nqp * Nrp * Nsp

    (sJ1, sJ2, sJ3, sJ4, sJ5, sJ6) = metrics.sJ
    J = metrics.J
    C = metrics.C

    nx = metrics.nx
    ny = metrics.ny
    nz = metrics.nz
 
    # FACE 1
    Nx1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx1 .= 0
    Nx1[1, :, :] = metrics.nx[1]
    Nx1 = spdiagm(0 => Nx1[:])
    
    Ny1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny1 .= 0
    Ny1[1, :, :] = metrics.ny[1]
    Ny1 = spdiagm(0 => Ny1[:])
    Nz1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz1 .= 0
    Nz1[1, :, :] = metrics.nz[1]
    Nz1 = spdiagm(0 => Nz1[:])

  
    # FACE 2
    Nx2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx2 .= 0
    Nx2[Nqp, :, :] = metrics.nx[2]
    Nx2 = spdiagm(0 => Nx2[:])
    Ny2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny2 .= 0
    Ny2[Nqp, :, :] = metrics.ny[2]
    Ny2 = spdiagm(0 => Ny2[:])
    Nz2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz2 .= 0
    Nz2[Nqp, :, :] = metrics.nz[2]
    Nz2 = spdiagm(0 => Nz2[:])

    # FACE 3
    Nx3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx3 .= 0
    Nx3[:, 1, :] = metrics.nx[3]
    Nx3 = spdiagm(0 => Nx3[:])
    Ny3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny3 .= 0
    Ny3[:, 1, :] = metrics.ny[3]
    Ny3 = spdiagm(0 => Ny3[:])
    Nz3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz3 .= 0
    Nz3[:, 1, :] = metrics.nz[3]
    Nz3 = spdiagm(0 => Nz3[:])

        # FACE 3
    Nx3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx3 .= 0
    Nx3[:, 1, :] = metrics.nx[3]
    Nx3 = spdiagm(0 => Nx3[:])
    Ny3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny3 .= 0
    Ny3[:, 1, :] = metrics.ny[3]
    Ny3 = spdiagm(0 => Ny3[:])
    Nz3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz3 .= 0
    Nz3[:, 1, :] = metrics.nz[3]
    Nz3 = spdiagm(0 => Nz3[:])

    # FACE 4
    Nx4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx4 .= 0
    Nx4[:, Nrp, :] = metrics.nx[4]
    Nx4 = spdiagm(0 => Nx4[:])
    Ny4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny4 .= 0
    Ny4[:, Nrp, :] = metrics.ny[4]
    Ny4 = spdiagm(0 => Ny4[:])
    Nz4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz4 .= 0
    Nz4[:, Nrp, :] = metrics.nz[4]
    Nz4 = spdiagm(0 => Nz4[:])
   
    # FACE 5
    Nx5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx5 .= 0
    Nx5[:, :, 1] = metrics.nx[5]
    Nx5 = spdiagm(0 => Nx5[:])
    Ny5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny5 .= 0
    Ny5[:, :, 1] = metrics.ny[5]
    Ny5 = spdiagm(0 => Ny5[:])
    Nz5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz5 .= 0
    Nz5[:, :, 1] = metrics.nz[5]
    Nz5 = spdiagm(0 => Nz5[:])
  
    # FACE 6
    Nx6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nx6 .= 0
    Nx6[:, :, Nsp] = metrics.nx[6]
    Nx6 = spdiagm(0 => Nx6[:])
    Ny6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Ny6 .= 0
    Ny6[:, :, Nsp] = metrics.ny[6]
    Ny6 = spdiagm(0 => Ny6[:])
    Nz6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    Nz6 .= 0
    Nz6[:, :, Nsp] = metrics.nz[6]
    Nz6 = spdiagm(0 => Nz6[:])
    # define Jacobian matrix evaluated on faces: (are these same as surface J??)

    EsJ1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ1 .= 0
    EsJ1[1, :, :] = sJ1 ./ J[1, :, :]
    EsJ2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ2 .= 0
    EsJ2[end, :, :] = sJ2 ./ J[end, :, :]

    EsJ3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ3 .= 0
    EsJ3[:, 1, :] = sJ3 ./ J[:, 1, :]
    EsJ4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ4 .= 0
    EsJ4[:, end, :] = sJ4 ./ J[:, end, :]

    EsJ5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ5 .= 0
    EsJ5[:, :, 1] = sJ5 ./ J[:, :, 1]
    EsJ6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    EsJ6 .= 0
    EsJ6[:, :, end] = sJ6 ./ J[:, :, end]


    sJI1 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI1 .= 0
    sJI1[1, :, :] = 1 ./ sJ1
    sJI1 =   spdiagm(0 =>  sJI1[:])

    sJI2 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI2 .= 0
    sJI2[end, :, :] = 1 ./ sJ2
    sJI2 =   spdiagm(0 =>  sJI2[:])

    sJI3 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI3 .= 0
    sJI3[:, 1, :] = 1 ./ sJ3
    sJI3 =   spdiagm(0 =>  sJI3[:])

    sJI4 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI4 .= 0
    sJI4[:, end, :] = 1 ./ sJ4
    sJI4 =   spdiagm(0 =>  sJI4[:])

    sJI5 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI5 .= 0
    sJI5[:, :, 1] = 1 ./ sJ5
    sJI5 =   spdiagm(0 =>  sJI5[:])

    sJI6 = SparseArray{Float64, 3}(undef, (Nqp, Nrp, Nsp)) 
    sJI6 .= 0
    sJI6[:, :, end] = 1 ./ sJ6
    sJI6 =   spdiagm(0 =>  sJI6[:])


    # Turn J and sJ's into diagonal matrices
    JI =  spdiagm(0 => 1 ./ J[:])
    J =   spdiagm(0 => J[:])
    sJ1 = spdiagm(0 => sJ1[:])
    sJ2 = spdiagm(0 => sJ2[:])
    sJ3 = spdiagm(0 => sJ3[:])
    sJ4 = spdiagm(0 => sJ4[:])
    sJ5 = spdiagm(0 => sJ5[:])
    sJ6 = spdiagm(0 => sJ6[:])

    EsJ1 = spdiagm(0 => EsJ1[:])
    EsJ2 = spdiagm(0 => EsJ2[:])
    EsJ3 = spdiagm(0 => EsJ3[:])
    EsJ4 = spdiagm(0 => EsJ4[:])
    EsJ5 = spdiagm(0 => EsJ5[:])
    EsJ6 = spdiagm(0 => EsJ6[:])

  

    c = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3, 3, 3)
    @show sizeof(c)
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    c[i, j, k, l] = spdiagm(0 => C[i, j, k, l][:])
                end
            end
        end
    end
    @show sizeof(c)

    # First derivative operators:
    (Dq, HqI, Hq, q) = diagonal_sbp_D1(p, Nq; xc = (-1,1))
    Qq = Hq * Dq
    QqT = sparse(transpose(Qq))

    (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
    Qr = Hr * Dr
    QrT = sparse(transpose(Qr))

    (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))
    Qs = Hs * Ds
    QsT = sparse(transpose(Qs))

    # Identity matrices for the comuptation
    Iq = sparse(I, Nqp, Nqp)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)

    #Variable Coefficient Pure Second Derivative Operators
    #(D2q, BSq, _, _, q) = variable_diagonal_sbp_D2(p, Nq, B; xc = (-1,1))
    #(D2r, BSr, _, _, r) = variable_diagonal_sbp_D2(p, Nr, B; xc = (-1,1))
    #(D2s, BSs, _, _, s) = variable_diagonal_sbp_D2(p, Ns, B; xc = (-1,1))
    (D2q, S0q, SNq, _, _, _) = diagonal_sbp_D2(p, Nq; xc = (-1,1))
    (D2r, S0r, SNr, _, _, _) = diagonal_sbp_D2(p, Nr; xc = (-1,1))
    (D2s, S0s, SNs, _, _, _) = diagonal_sbp_D2(p, Ns; xc = (-1,1))

    D11 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D12 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D13 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D21 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D22 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D23 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D31 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D32 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)
    D33 = fill(spzeros(Nqp*Nrp*Nsp, Nqp*Nrp*Nsp), 3, 3)


    for i = 1:3
        for j = 1:3
            #D11[i, j] = c[1, i, 1, j] * JI * (Is ⊗ Ir ⊗ D2q)
           (D11[i, j], _, _) = var_3D_D2q(p, Nqp, Nrp, Nsp, metrics.C[1, i, 1, j], HqI; xc = (-1, 1))
            D11[i, j] = JI * D11[i, j]
            D12[i, j] = c[1, i, 2, j] * JI * (Is ⊗ Dr ⊗ Dq)
            D13[i, j] = c[1, i, 3, j] * JI * (Ds ⊗ Ir ⊗ Dq)

            D21[i, j] = c[2, i, 1, j] * JI * (Is ⊗ Dr ⊗ Dq)
            #D22[i, j] = c[2, i, 2, j] * JI * (Is ⊗ D2r ⊗ Iq)
           (D22[i, j], _, _) = var_3D_D2r(p, Nqp, Nrp, Nsp, metrics.C[2, i, 2, j], HrI; xc = (-1, 1))
            D22[i, j] = JI * D22[i, j]
            D23[i, j] = c[2, i, 3, j] * JI * (Ds ⊗ Dr ⊗ Iq)

            D31[i, j] = c[3, i, 1, j] * JI * (Ds ⊗ Ir ⊗ Dq)
            D32[i, j] = c[3, i, 2, j] * JI * (Ds ⊗ Dr ⊗ Iq)
            #D33[i, j] = c[3, i, 3, j] * JI * (D2s ⊗ Ir ⊗ Iq)
           (D33[i, j], _, _) = var_3D_D2s(p, Nqp, Nrp, Nsp, metrics.C[3, i, 3, j], HsI; xc = (-1, 1))
            D33[i, j] = JI * D33[i, j]

        end
    end

    @show sizeof(D33)

    A11 = D11[1, 1] .+ D12[1, 1] .+ D13[1, 1] .+ 
          D21[1, 1] .+ D22[1, 1] .+ D23[1, 1] .+ 
          D31[1, 1] .+ D32[1, 1] .+ D33[1, 1]

    A12 = D11[1, 2] .+ D12[1, 2] .+ D13[1, 2] .+ 
          D21[1, 2] .+ D22[1, 2] .+ D23[1, 2] .+ 
          D31[1, 2] .+ D32[1, 2] .+ D33[1, 2]

    A13 = D11[1, 3] .+ D12[1, 3] .+ D13[1, 3] .+ 
          D21[1, 3] .+ D22[1, 3] .+ D23[1, 3] .+ 
          D31[1, 3] .+ D32[1, 3] .+ D33[1, 3]

    A21 = D11[2, 1] .+ D12[2, 1] .+ D13[2, 1] .+ 
          D21[2, 1] .+ D22[2, 1] .+ D23[2, 1] .+ 
          D31[2, 1] .+ D32[2, 1] .+ D33[2, 1]

    A22 = D11[2, 2] .+ D12[2, 2] .+ D13[2, 2] .+ 
          D21[2, 2] .+ D22[2, 2] .+ D23[2, 2] .+ 
          D31[2, 2] .+ D32[2, 2] .+ D33[2, 2]      

    A23 = D11[2, 3] .+ D12[2, 3] .+ D13[2, 3] .+ 
          D21[2, 3] .+ D22[2, 3] .+ D23[2, 3] .+ 
          D31[2, 3] .+ D32[2, 3] .+ D33[2, 3]

    A31 = D11[3, 1] .+ D12[3, 1] .+ D13[3, 1] .+ 
          D21[3, 1] .+ D22[3, 1] .+ D23[3, 1] .+ 
          D31[3, 1] .+ D32[3, 1] .+ D33[3, 1]

    A32 = D11[3, 2] .+ D12[3, 2] .+ D13[3, 2] .+ 
          D21[3, 2] .+ D22[3, 2] .+ D23[3, 2] .+ 
          D31[3, 2] .+ D32[3, 2] .+ D33[3, 2]

    A33 = D11[3, 3] .+ D12[3, 3] .+ D13[3, 3] .+ 
          D21[3, 3] .+ D22[3, 3] .+ D23[3, 3] .+ 
          D31[3, 3] .+ D32[3, 3] .+ D33[3, 3]

    

     # Surface mass matrices
    #
    H1 = H2 = Hs ⊗ Hr
    H1I = H2I = HsI ⊗ HrI

    H3 = H4 = Hs ⊗ Hq
    H3I = H4I = HsI ⊗ HqI

    H5 = H6 = Hr ⊗ Hq
    H5I = H6I = HrI ⊗ HqI

    # Volume matrices
    H = Hs ⊗ Hr ⊗ Hq
    HI = HsI ⊗ HrI ⊗ HqI
    JHI = HI * JI

    JIm = spdiagm(0 => JI[:])
    # Create 3D ops from 1D
    Dq3 = Is ⊗ Ir ⊗ Dq
    Dr3 = Is ⊗ Dr ⊗ Iq
    Ds3 = Ds ⊗ Ir ⊗ Iq

    # Face operators to reduce computations involing T11
    # S0q = Is ⊗ Ir ⊗ S0q
    # SNq = Is ⊗ Ir ⊗ SNq
    # S0r = Is ⊗ S0r ⊗ Iq
    # SNr = Is ⊗ SNr ⊗ Iq
    # S0s = S0s ⊗ Ir ⊗ Iq
    # SNs = SNs ⊗ Ir ⊗ Iq

    Sq = Is ⊗ Ir ⊗ (SNq+S0q)
    Sr = Is ⊗ (SNr+S0r) ⊗ Iq
    Ss = (SNs+S0s) ⊗ Ir ⊗ Iq
    # SNq = Is ⊗ Ir ⊗ SNq
    # S0r = Is ⊗ S0r ⊗ Iq
    # SNr = Is ⊗ SNr ⊗ Iq
    # S0s = S0s ⊗ Ir ⊗ Iq
    # SNs = SNs ⊗ Ir ⊗ Iq
    
# Create traction operators on each face
# factor for turning T's on/off, default to c = 1

j = 1

# FACE 1
# (nq, nr, ns) = (-1, 0, 0)
@show size(J)
@show size(sJ1)
@show size(1 ./ sJ1)
@show size(c[1,1,1,1]*Sq)
@show size(sJI1)
T11_1 = (-sJI1) * (c[1,1,1,1]*Sq + c[1,1,2,1]*Dr3 + c[1,1,3,1]*Ds3)
T12_1 = (-sJI1) * (c[1,1,1,2]*Sq + c[1,1,2,2]*Dr3 + c[1,1,3,2]*Ds3) 
T21_1 = (-sJI1) * (c[1,2,1,1]*Sq + c[1,2,2,1]*Dr3 + c[1,2,3,1]*Ds3)
T13_1 = (-sJI1) * (c[1,1,1,3]*Sq + c[1,1,2,3]*Dr3 + c[1,1,3,3]*Ds3) 
T31_1 = (-sJI1) * (c[1,3,1,1]*Sq + c[1,3,2,1]*Dr3 + c[1,3,3,1]*Ds3) 
T22_1 = (-sJI1) * (c[1,2,1,2]*Sq + c[1,2,2,2]*Dr3 + c[1,2,3,2]*Ds3) 
T23_1 = (-sJI1) * (c[1,2,1,3]*Sq + c[1,2,2,3]*Dr3 + c[1,2,3,3]*Ds3) 
T32_1 = (-sJI1) * (c[1,3,1,2]*Sq + c[1,3,2,2]*Dr3 + c[1,3,3,2]*Ds3)
T33_1 = (-sJI1) * (c[1,3,1,3]*Sq + c[1,3,2,3]*Dr3 + c[1,3,3,3]*Ds3) 
   
# FACE 2
    # (Nq, Nr, Ns) = (1, 0, 0)
    T11_2 = (sJI2) * (c[1,1,1,1]*Sq + c[1,1,2,1]*Dr3 + c[1,1,3,1]*Ds3) 
    T12_2 = (sJI2) * (c[1,1,1,2]*Sq + c[1,1,2,2]*Dr3 + c[1,1,3,2]*Ds3) 
    T21_2 = (sJI2) * (c[1,2,1,1]*Sq + c[1,2,2,1]*Dr3 + c[1,2,3,1]*Ds3) 
    T13_2 = (sJI2) * (c[1,1,1,3]*Sq + c[1,1,2,3]*Dr3 + c[1,1,3,3]*Ds3) 
    T31_2 = (sJI2) * (c[1,3,1,1]*Sq + c[1,3,2,1]*Dr3 + c[1,3,3,1]*Ds3) 
    T22_2 = (sJI2) * (c[1,2,1,2]*Sq + c[1,2,2,2]*Dr3 + c[1,2,3,2]*Ds3) 
    T23_2 = (sJI2) * (c[1,2,1,3]*Sq + c[1,2,2,3]*Dr3 + c[1,2,3,3]*Ds3) 
    T32_2 = (sJI2) * (c[1,3,1,2]*Sq + c[1,3,2,2]*Dr3 + c[1,3,3,2]*Ds3) 
    T33_2 = (sJI2) * (c[1,3,1,3]*Sq + c[1,3,2,3]*Dr3 + c[1,3,3,3]*Ds3) 



    
#   # FACE 3
#     # (Nq, Nr, Ns) = (0, -1, 0)
T11_3 = (-sJI3) * (c[2,1,1,1]*Dq3 + c[2,1,2,1]*Sr + c[2,1,3,1]*Ds3) 
T12_3 = (-sJI3) * (c[2,1,1,2]*Dq3 + c[2,1,2,2]*Sr + c[2,1,3,2]*Ds3) 
T21_3 = (-sJI3) * (c[2,2,1,1]*Dq3 + c[2,2,2,1]*Sr + c[2,2,3,1]*Ds3)
T13_3 = (-sJI3) * (c[2,1,1,3]*Dq3 + c[2,1,2,3]*Sr + c[2,1,3,3]*Ds3) 
T31_3 = (-sJI3) * (c[2,3,1,1]*Dq3 + c[2,3,2,1]*Sr + c[2,3,3,1]*Ds3) 
T22_3 = (-sJI3) * (c[2,2,1,2]*Dq3 + c[2,2,2,2]*Sr + c[2,2,3,2]*Ds3) 
T23_3 = (-sJI3) * (c[2,2,1,3]*Dq3 + c[2,2,2,3]*Sr + c[2,2,3,3]*Ds3) 
T32_3 = (-sJI3) * (c[2,3,1,2]*Dq3 + c[2,3,2,2]*Sr + c[2,3,3,2]*Ds3)
T33_3 = (-sJI3) * (c[2,3,1,3]*Dq3 + c[2,3,2,3]*Sr + c[2,3,3,3]*Ds3) 

    # FACE 4
    # (Nq, Nr, Ns) = (0, 1, 0)
    T11_4 = (sJI4) * (c[2,1,1,1]*Dq3 + c[2,1,2,1]*Sr + c[2,1,3,1]*Ds3) 
    T12_4 = (sJI4) * (c[2,1,1,2]*Dq3 + c[2,1,2,2]*Sr + c[2,1,3,2]*Ds3)
    T21_4 = (sJI4) * (c[2,2,1,1]*Dq3 + c[2,2,2,1]*Sr + c[2,2,3,1]*Ds3)
    T13_4 = (sJI4) * (c[2,1,1,3]*Dq3 + c[2,1,2,3]*Sr + c[2,1,3,3]*Ds3)
    T31_4 = (sJI4) * (c[2,3,1,1]*Dq3 + c[2,3,2,1]*Sr + c[2,3,3,1]*Ds3) 
    T22_4 = (sJI4) * (c[2,2,1,2]*Dq3 + c[2,2,2,2]*Sr + c[2,2,3,2]*Ds3)
    T23_4 = (sJI4) * (c[2,2,1,3]*Dq3 + c[2,2,2,3]*Sr + c[2,2,3,3]*Ds3)
    T32_4 = (sJI4) * (c[2,3,1,2]*Dq3 + c[2,3,2,2]*Sr + c[2,3,3,2]*Ds3)
    T33_4 = (sJI4) * (c[2,3,1,3]*Dq3 + c[2,3,2,3]*Sr + c[2,3,3,3]*Ds3)



    # FACE 5
    # (Nq, Nr, Ns) = (0, 0, -1)
    T11_5 = (-sJI5) * (c[3,1,1,1]*Dq3 + c[3,1,2,1]*Dr3 + c[3,1,3,1]*Ss)
    T12_5 = (-sJI5) * (c[3,1,1,2]*Dq3 + c[3,1,2,2]*Dr3 + c[3,1,3,2]*Ss)
    T21_5 = (-sJI5) * (c[3,2,1,1]*Dq3 + c[3,2,2,1]*Dr3 + c[3,2,3,1]*Ss)
    T13_5 = (-sJI5) * (c[3,1,1,3]*Dq3 + c[3,1,2,3]*Dr3 + c[3,1,3,3]*Ss)
    T31_5 = (-sJI5) * (c[3,3,1,1]*Dq3 + c[3,3,2,1]*Dr3 + c[3,3,3,1]*Ss)
    T22_5 = (-sJI5) * (c[3,2,1,2]*Dq3 + c[3,2,2,2]*Dr3 + c[3,2,3,2]*Ss)
    T23_5 = (-sJI5) * (c[3,2,1,3]*Dq3 + c[3,2,2,3]*Dr3 + c[3,2,3,3]*Ss)
    T32_5 = (-sJI5) * (c[3,3,1,2]*Dq3 + c[3,3,2,2]*Dr3 + c[3,3,3,2]*Ss)
    T33_5 = (-sJI5) * (c[3,3,1,3]*Dq3 + c[3,3,2,3]*Dr3 + c[3,3,3,3]*Ss)


    # FACE 6
    # (Nq, Nr, Ns) = (0, 0, 1)
    T11_6 = sJI6 * (c[3,1,1,1]*Dq3 + c[3,1,2,1]*Dr3 + c[3,1,3,1]*Ss)
    T12_6 = sJI6 *  (c[3,1,1,2]*Dq3 + c[3,1,2,2]*Dr3 + c[3,1,3,2]*Ss)
    T21_6 = sJI6 *  (c[3,2,1,1]*Dq3 + c[3,2,2,1]*Dr3 + c[3,2,3,1]*Ss)
    T13_6 = sJI6 *  (c[3,1,1,3]*Dq3 + c[3,1,2,3]*Dr3 + c[3,1,3,3]*Ss)
    T31_6 = sJI6 *  (c[3,3,1,1]*Dq3 + c[3,3,2,1]*Dr3 + c[3,3,3,1]*Ss)
    T22_6 = sJI6 *  (c[3,2,1,2]*Dq3 + c[3,2,2,2]*Dr3 + c[3,2,3,2]*Ss)
    T23_6 = sJI6 *  (c[3,2,1,3]*Dq3 + c[3,2,2,3]*Dr3 + c[3,2,3,3]*Ss)
    T32_6 = sJI6 *  (c[3,3,1,2]*Dq3 + c[3,3,2,2]*Dr3 + c[3,3,3,2]*Ss)
    T33_6 = sJI6 *  (c[3,3,1,3]*Dq3 + c[3,3,2,3]*Dr3 + c[3,3,3,3]*Ss)


        
    beta = 1
    h1 = Hq[1,1] #TODO: fix this
    d = 2 #dimension? 
    g = 1


    Z11_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,1,1,1]*Nx1 + c[1,1,2,1]*Ny1 + c[1,1,3,1]*Nz1) + Ny1 * (c[2,1,1,1]*Nx1 + c[2,1,2,1]*Ny1 + c[2,1,3,1]*Nz1) + Nz1 * (c[3,1,1,1]*Nx1 + c[3,1,2,1]*Ny1 + c[3,1,3,1]*Nz1)) 
    Z12_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,1,1,2]*Nx1 + c[1,1,2,2]*Ny1 + c[1,1,3,2]*Nz1) + Ny1 * (c[2,1,1,2]*Nx1 + c[2,1,2,2]*Ny1 + c[2,1,3,2]*Nz1) + Nz1 * (c[3,1,1,2]*Nx1 + c[3,1,2,2]*Ny1 + c[3,1,3,2]*Nz1)) 
    Z21_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,2,1,1]*Nx1 + c[1,2,2,1]*Ny1 + c[1,2,3,1]*Nz1) + Ny1 * (c[2,2,1,1]*Nx1 + c[2,2,2,1]*Ny1 + c[2,2,3,1]*Nz1) + Nz1 * (c[3,2,1,1]*Nx1 + c[3,2,2,1]*Ny1 + c[3,2,3,1]*Nz1)) 
    Z13_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,1,1,3]*Nx1 + c[1,1,2,3]*Ny1 + c[1,1,3,3]*Nz1) + Ny1 * (c[2,1,1,3]*Nx1 + c[2,1,2,3]*Ny1 + c[2,1,3,3]*Nz1) + Nz1 * (c[3,1,1,3]*Nx1 + c[3,1,2,3]*Ny1 + c[3,1,3,3]*Nz1)) 
    Z31_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,3,1,1]*Nx1 + c[1,3,2,1]*Ny1 + c[1,3,3,1]*Nz1) + Ny1 * (c[2,3,1,1]*Nx1 + c[2,3,2,1]*Ny1 + c[2,3,3,1]*Nz1) + Nz1 * (c[3,3,1,1]*Nx1 + c[3,3,2,1]*Ny1 + c[3,3,3,1]*Nz1)) 
    Z22_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,2,1,2]*Nx1 + c[1,2,2,2]*Ny1 + c[1,2,3,2]*Nz1) + Ny1 * (c[2,2,1,2]*Nx1 + c[2,2,2,2]*Ny1 + c[2,2,3,2]*Nz1) + Nz1 * (c[3,2,1,2]*Nx1 + c[3,2,2,2]*Ny1 + c[3,2,3,2]*Nz1)) 
    Z23_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,2,1,3]*Nx1 + c[1,2,2,3]*Ny1 + c[1,2,3,3]*Nz1) + Ny1 * (c[2,2,1,3]*Nx1 + c[2,2,2,3]*Ny1 + c[2,2,3,3]*Nz1) + Nz1 * (c[3,2,1,3]*Nx1 + c[3,2,2,3]*Ny1 + c[3,2,3,3]*Nz1)) 
    Z32_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,3,1,2]*Nx1 + c[1,3,2,2]*Ny1 + c[1,3,3,2]*Nz1) + Ny1 * (c[2,3,1,2]*Nx1 + c[2,3,2,2]*Ny1 + c[2,3,3,2]*Nz1) + Nz1 * (c[3,3,1,2]*Nx1 + c[3,3,2,2]*Ny1 + c[3,3,3,2]*Nz1)) 
    Z33_1 = (beta/h1) * d * (EsJ1) * (Nx1 * (c[1,3,1,3]*Nx1 + c[1,3,2,3]*Ny1 + c[1,3,3,3]*Nz1) + Ny1 * (c[2,3,1,3]*Nx1 + c[2,3,2,3]*Ny1 + c[2,3,3,3]*Nz1) + Nz1 * (c[3,3,1,3]*Nx1 + c[3,3,2,3]*Ny1 + c[3,3,3,3]*Nz1)) 

    
    Z11_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,1,1,1]*Nx2 + c[1,1,2,1]*Ny2 + c[1,1,3,1]*Nz2) + Ny2 * (c[2,1,1,1]*Nx2 + c[2,1,2,1]*Ny2 + c[2,1,3,1]*Nz2) + Nz2 * (c[3,1,1,1]*Nx2 + c[3,1,2,1]*Ny2 + c[3,1,3,1]*Nz2)) 
    Z12_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,1,1,2]*Nx2 + c[1,1,2,2]*Ny2 + c[1,1,3,2]*Nz2) + Ny2 * (c[2,1,1,2]*Nx2 + c[2,1,2,2]*Ny2 + c[2,1,3,2]*Nz2) + Nz2 * (c[3,1,1,2]*Nx2 + c[3,1,2,2]*Ny2 + c[3,1,3,2]*Nz2)) 
    Z21_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,2,1,1]*Nx2 + c[1,2,2,1]*Ny2 + c[1,2,3,1]*Nz2) + Ny2 * (c[2,2,1,1]*Nx2 + c[2,2,2,1]*Ny2 + c[2,2,3,1]*Nz2) + Nz2 * (c[3,2,1,1]*Nx2 + c[3,2,2,1]*Ny2 + c[3,2,3,1]*Nz2)) 
    Z13_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,1,1,3]*Nx2 + c[1,1,2,3]*Ny2 + c[1,1,3,3]*Nz2) + Ny2 * (c[2,1,1,3]*Nx2 + c[2,1,2,3]*Ny2 + c[2,1,3,3]*Nz2) + Nz2 * (c[3,1,1,3]*Nx2 + c[3,1,2,3]*Ny2 + c[3,1,3,3]*Nz2)) 
    Z31_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,3,1,1]*Nx2 + c[1,3,2,1]*Ny2 + c[1,3,3,1]*Nz2) + Ny2 * (c[2,3,1,1]*Nx2 + c[2,3,2,1]*Ny2 + c[2,3,3,1]*Nz2) + Nz2 * (c[3,3,1,1]*Nx2 + c[3,3,2,1]*Ny2 + c[3,3,3,1]*Nz2)) 
    Z22_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,2,1,2]*Nx2 + c[1,2,2,2]*Ny2 + c[1,2,3,2]*Nz2) + Ny2 * (c[2,2,1,2]*Nx2 + c[2,2,2,2]*Ny2 + c[2,2,3,2]*Nz2) + Nz2 * (c[3,2,1,2]*Nx2 + c[3,2,2,2]*Ny2 + c[3,2,3,2]*Nz2)) 
    Z23_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,2,1,3]*Nx2 + c[1,2,2,3]*Ny2 + c[1,2,3,3]*Nz2) + Ny2 * (c[2,2,1,3]*Nx2 + c[2,2,2,3]*Ny2 + c[2,2,3,3]*Nz2) + Nz2 * (c[3,2,1,3]*Nx2 + c[3,2,2,3]*Ny2 + c[3,2,3,3]*Nz2)) 
    Z32_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,3,1,2]*Nx2 + c[1,3,2,2]*Ny2 + c[1,3,3,2]*Nz2) + Ny2 * (c[2,3,1,2]*Nx2 + c[2,3,2,2]*Ny2 + c[2,3,3,2]*Nz2) + Nz2 * (c[3,3,1,2]*Nx2 + c[3,3,2,2]*Ny2 + c[3,3,3,2]*Nz2)) 
    Z33_2 = (beta/h1) * d * (EsJ2) * (Nx2 * (c[1,3,1,3]*Nx2 + c[1,3,2,3]*Ny2 + c[1,3,3,3]*Nz2) + Ny2 * (c[2,3,1,3]*Nx2 + c[2,3,2,3]*Ny2 + c[2,3,3,3]*Nz2) + Nz2 * (c[3,3,1,3]*Nx2 + c[3,3,2,3]*Ny2 + c[3,3,3,3]*Nz2)) 
    

     # FACE 3 
    # (nq, nr, ns) = (0, -1, 0)
    Z11_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,1,1,1]*Nx3 + c[1,1,2,1]*Ny3 + c[1,1,3,1]*Nz3) + Ny3 * (c[2,1,1,1]*Nx3 + c[2,1,2,1]*Ny3 + c[2,1,3,1]*Nz3) + Nz3 * (c[3,1,1,1]*Nx3 + c[3,1,2,1]*Ny3 + c[3,1,3,1]*Nz3)) 
    Z12_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,1,1,2]*Nx3 + c[1,1,2,2]*Ny3 + c[1,1,3,2]*Nz3) + Ny3 * (c[2,1,1,2]*Nx3 + c[2,1,2,2]*Ny3 + c[2,1,3,2]*Nz3) + Nz3 * (c[3,1,1,2]*Nx3 + c[3,1,2,2]*Ny3 + c[3,1,3,2]*Nz3)) 
    Z21_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,2,1,1]*Nx3 + c[1,2,2,1]*Ny3 + c[1,2,3,1]*Nz3) + Ny3 * (c[2,2,1,1]*Nx3 + c[2,2,2,1]*Ny3 + c[2,2,3,1]*Nz3) + Nz3 * (c[3,2,1,1]*Nx3 + c[3,2,2,1]*Ny3 + c[3,2,3,1]*Nz3)) 
    Z13_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,1,1,3]*Nx3 + c[1,1,2,3]*Ny3 + c[1,1,3,3]*Nz3) + Ny3 * (c[2,1,1,3]*Nx3 + c[2,1,2,3]*Ny3 + c[2,1,3,3]*Nz3) + Nz3 * (c[3,1,1,3]*Nx3 + c[3,1,2,3]*Ny3 + c[3,1,3,3]*Nz3)) 
    Z31_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,3,1,1]*Nx3 + c[1,3,2,1]*Ny3 + c[1,3,3,1]*Nz3) + Ny3 * (c[2,3,1,1]*Nx3 + c[2,3,2,1]*Ny3 + c[2,3,3,1]*Nz3) + Nz3 * (c[3,3,1,1]*Nx3 + c[3,3,2,1]*Ny3 + c[3,3,3,1]*Nz3)) 
    Z22_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,2,1,2]*Nx3 + c[1,2,2,2]*Ny3 + c[1,2,3,2]*Nz3) + Ny3 * (c[2,2,1,2]*Nx3 + c[2,2,2,2]*Ny3 + c[2,2,3,2]*Nz3) + Nz3 * (c[3,2,1,2]*Nx3 + c[3,2,2,2]*Ny3 + c[3,2,3,2]*Nz3)) 
    Z23_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,2,1,3]*Nx3 + c[1,2,2,3]*Ny3 + c[1,2,3,3]*Nz3) + Ny3 * (c[2,2,1,3]*Nx3 + c[2,2,2,3]*Ny3 + c[2,2,3,3]*Nz3) + Nz3 * (c[3,2,1,3]*Nx3 + c[3,2,2,3]*Ny3 + c[3,2,3,3]*Nz3)) 
    Z32_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,3,1,2]*Nx3 + c[1,3,2,2]*Ny3 + c[1,3,3,2]*Nz3) + Ny3 * (c[2,3,1,2]*Nx3 + c[2,3,2,2]*Ny3 + c[2,3,3,2]*Nz3) + Nz3 * (c[3,3,1,2]*Nx3 + c[3,3,2,2]*Ny3 + c[3,3,3,2]*Nz3)) 
    Z33_3 = (beta/h1) * d * (EsJ3) * (Nx3 * (c[1,3,1,3]*Nx3 + c[1,3,2,3]*Ny3 + c[1,3,3,3]*Nz3) + Ny3 * (c[2,3,1,3]*Nx3 + c[2,3,2,3]*Ny3 + c[2,3,3,3]*Nz3) + Nz3 * (c[3,3,1,3]*Nx3 + c[3,3,2,3]*Ny3 + c[3,3,3,3]*Nz3)) 

   # FACE 4 
    # (nq, nr, ns) = (0, 1, 0)
    Z11_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,1,1,1]*Nx4 + c[1,1,2,1]*Ny4 + c[1,1,3,1]*Nz4) + Ny4 * (c[2,1,1,1]*Nx4 + c[2,1,2,1]*Ny4 + c[2,1,3,1]*Nz4) + Nz4 * (c[3,1,1,1]*Nx4 + c[3,1,2,1]*Ny4 + c[3,1,3,1]*Nz4)) 
    Z12_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,1,1,2]*Nx4 + c[1,1,2,2]*Ny4 + c[1,1,3,2]*Nz4) + Ny4 * (c[2,1,1,2]*Nx4 + c[2,1,2,2]*Ny4 + c[2,1,3,2]*Nz4) + Nz4 * (c[3,1,1,2]*Nx4 + c[3,1,2,2]*Ny4 + c[3,1,3,2]*Nz4)) 
    Z21_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,2,1,1]*Nx4 + c[1,2,2,1]*Ny4 + c[1,2,3,1]*Nz4) + Ny4 * (c[2,2,1,1]*Nx4 + c[2,2,2,1]*Ny4 + c[2,2,3,1]*Nz4) + Nz4 * (c[3,2,1,1]*Nx4 + c[3,2,2,1]*Ny4 + c[3,2,3,1]*Nz4)) 
    Z13_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,1,1,3]*Nx4 + c[1,1,2,3]*Ny4 + c[1,1,3,3]*Nz4) + Ny4 * (c[2,1,1,3]*Nx4 + c[2,1,2,3]*Ny4 + c[2,1,3,3]*Nz4) + Nz4 * (c[3,1,1,3]*Nx4 + c[3,1,2,3]*Ny4 + c[3,1,3,3]*Nz4)) 
    Z31_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,3,1,1]*Nx4 + c[1,3,2,1]*Ny4 + c[1,3,3,1]*Nz4) + Ny4 * (c[2,3,1,1]*Nx4 + c[2,3,2,1]*Ny4 + c[2,3,3,1]*Nz4) + Nz4 * (c[3,3,1,1]*Nx4 + c[3,3,2,1]*Ny4 + c[3,3,3,1]*Nz4)) 
    Z22_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,2,1,2]*Nx4 + c[1,2,2,2]*Ny4 + c[1,2,3,2]*Nz4) + Ny4 * (c[2,2,1,2]*Nx4 + c[2,2,2,2]*Ny4 + c[2,2,3,2]*Nz4) + Nz4 * (c[3,2,1,2]*Nx4 + c[3,2,2,2]*Ny4 + c[3,2,3,2]*Nz4)) 
    Z23_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,2,1,3]*Nx4 + c[1,2,2,3]*Ny4 + c[1,2,3,3]*Nz4) + Ny4 * (c[2,2,1,3]*Nx4 + c[2,2,2,3]*Ny4 + c[2,2,3,3]*Nz4) + Nz4 * (c[3,2,1,3]*Nx4 + c[3,2,2,3]*Ny4 + c[3,2,3,3]*Nz4)) 
    Z32_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,3,1,2]*Nx4 + c[1,3,2,2]*Ny4 + c[1,3,3,2]*Nz4) + Ny4 * (c[2,3,1,2]*Nx4 + c[2,3,2,2]*Ny4 + c[2,3,3,2]*Nz4) + Nz4 * (c[3,3,1,2]*Nx4 + c[3,3,2,2]*Ny4 + c[3,3,3,2]*Nz4)) 
    Z33_4 = (beta/h1) * d * (EsJ4) * (Nx4 * (c[1,3,1,3]*Nx4 + c[1,3,2,3]*Ny4 + c[1,3,3,3]*Nz4) + Ny4 * (c[2,3,1,3]*Nx4 + c[2,3,2,3]*Ny4 + c[2,3,3,3]*Nz4) + Nz4 * (c[3,3,1,3]*Nx4 + c[3,3,2,3]*Ny4 + c[3,3,3,3]*Nz4)) 

    # FACE 5 
    # (nq, nr, ns) = (0, 0, -1)
    Z11_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,1,1,1]*Nx5 + c[1,1,2,1]*Ny5 + c[1,1,3,1]*Nz5) + Ny5 * (c[2,1,1,1]*Nx5 + c[2,1,2,1]*Ny5 + c[2,1,3,1]*Nz5) + Nz5 * (c[3,1,1,1]*Nx5 + c[3,1,2,1]*Ny5 + c[3,1,3,1]*Nz5)) 
    Z12_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,1,1,2]*Nx5 + c[1,1,2,2]*Ny5 + c[1,1,3,2]*Nz5) + Ny5 * (c[2,1,1,2]*Nx5 + c[2,1,2,2]*Ny5 + c[2,1,3,2]*Nz5) + Nz5 * (c[3,1,1,2]*Nx5 + c[3,1,2,2]*Ny5 + c[3,1,3,2]*Nz5)) 
    Z21_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,2,1,1]*Nx5 + c[1,2,2,1]*Ny5 + c[1,2,3,1]*Nz5) + Ny5 * (c[2,2,1,1]*Nx5 + c[2,2,2,1]*Ny5 + c[2,2,3,1]*Nz5) + Nz5 * (c[3,2,1,1]*Nx5 + c[3,2,2,1]*Ny5 + c[3,2,3,1]*Nz5)) 
    Z13_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,1,1,3]*Nx5 + c[1,1,2,3]*Ny5 + c[1,1,3,3]*Nz5) + Ny5 * (c[2,1,1,3]*Nx5 + c[2,1,2,3]*Ny5 + c[2,1,3,3]*Nz5) + Nz5 * (c[3,1,1,3]*Nx5 + c[3,1,2,3]*Ny5 + c[3,1,3,3]*Nz5)) 
    Z31_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,3,1,1]*Nx5 + c[1,3,2,1]*Ny5 + c[1,3,3,1]*Nz5) + Ny5 * (c[2,3,1,1]*Nx5 + c[2,3,2,1]*Ny5 + c[2,3,3,1]*Nz5) + Nz5 * (c[3,3,1,1]*Nx5 + c[3,3,2,1]*Ny5 + c[3,3,3,1]*Nz5)) 
    Z22_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,2,1,2]*Nx5 + c[1,2,2,2]*Ny5 + c[1,2,3,2]*Nz5) + Ny5 * (c[2,2,1,2]*Nx5 + c[2,2,2,2]*Ny5 + c[2,2,3,2]*Nz5) + Nz5 * (c[3,2,1,2]*Nx5 + c[3,2,2,2]*Ny5 + c[3,2,3,2]*Nz5)) 
    Z23_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,2,1,3]*Nx5 + c[1,2,2,3]*Ny5 + c[1,2,3,3]*Nz5) + Ny5 * (c[2,2,1,3]*Nx5 + c[2,2,2,3]*Ny5 + c[2,2,3,3]*Nz5) + Nz5 * (c[3,2,1,3]*Nx5 + c[3,2,2,3]*Ny5 + c[3,2,3,3]*Nz5)) 
    Z32_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,3,1,2]*Nx5 + c[1,3,2,2]*Ny5 + c[1,3,3,2]*Nz5) + Ny5 * (c[2,3,1,2]*Nx5 + c[2,3,2,2]*Ny5 + c[2,3,3,2]*Nz5) + Nz5 * (c[3,3,1,2]*Nx5 + c[3,3,2,2]*Ny5 + c[3,3,3,2]*Nz5)) 
    Z33_5 = (beta/h1) * d * (EsJ5) * (Nx5 * (c[1,3,1,3]*Nx5 + c[1,3,2,3]*Ny5 + c[1,3,3,3]*Nz5) + Ny5 * (c[2,3,1,3]*Nx5 + c[2,3,2,3]*Ny5 + c[2,3,3,3]*Nz5) + Nz5 * (c[3,3,1,3]*Nx5 + c[3,3,2,3]*Ny5 + c[3,3,3,3]*Nz5)) 

    # FACE 6 
    # (nq, nr, ns) = (0, 0, 1)
    Z11_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,1,1,1]*Nx6 + c[1,1,2,1]*Ny6 + c[1,1,3,1]*Nz6) + Ny6 * (c[2,1,1,1]*Nx6 + c[2,1,2,1]*Ny6 + c[2,1,3,1]*Nz6) + Nz6 * (c[3,1,1,1]*Nx6 + c[3,1,2,1]*Ny6 + c[3,1,3,1]*Nz6)) 
    Z12_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,1,1,2]*Nx6 + c[1,1,2,2]*Ny6 + c[1,1,3,2]*Nz6) + Ny6 * (c[2,1,1,2]*Nx6 + c[2,1,2,2]*Ny6 + c[2,1,3,2]*Nz6) + Nz6 * (c[3,1,1,2]*Nx6 + c[3,1,2,2]*Ny6 + c[3,1,3,2]*Nz6)) 
    Z21_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,2,1,1]*Nx6 + c[1,2,2,1]*Ny6 + c[1,2,3,1]*Nz6) + Ny6 * (c[2,2,1,1]*Nx6 + c[2,2,2,1]*Ny6 + c[2,2,3,1]*Nz6) + Nz6 * (c[3,2,1,1]*Nx6 + c[3,2,2,1]*Ny6 + c[3,2,3,1]*Nz6)) 
    Z13_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,1,1,3]*Nx6 + c[1,1,2,3]*Ny6 + c[1,1,3,3]*Nz6) + Ny6 * (c[2,1,1,3]*Nx6 + c[2,1,2,3]*Ny6 + c[2,1,3,3]*Nz6) + Nz6 * (c[3,1,1,3]*Nx6 + c[3,1,2,3]*Ny6 + c[3,1,3,3]*Nz6)) 
    Z31_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,3,1,1]*Nx6 + c[1,3,2,1]*Ny6 + c[1,3,3,1]*Nz6) + Ny6 * (c[2,3,1,1]*Nx6 + c[2,3,2,1]*Ny6 + c[2,3,3,1]*Nz6) + Nz6 * (c[3,3,1,1]*Nx6 + c[3,3,2,1]*Ny6 + c[3,3,3,1]*Nz6)) 
    Z22_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,2,1,2]*Nx6 + c[1,2,2,2]*Ny6 + c[1,2,3,2]*Nz6) + Ny6 * (c[2,2,1,2]*Nx6 + c[2,2,2,2]*Ny6 + c[2,2,3,2]*Nz6) + Nz6 * (c[3,2,1,2]*Nx6 + c[3,2,2,2]*Ny6 + c[3,2,3,2]*Nz6)) 
    Z23_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,2,1,3]*Nx6 + c[1,2,2,3]*Ny6 + c[1,2,3,3]*Nz6) + Ny6 * (c[2,2,1,3]*Nx6 + c[2,2,2,3]*Ny6 + c[2,2,3,3]*Nz6) + Nz6 * (c[3,2,1,3]*Nx6 + c[3,2,2,3]*Ny6 + c[3,2,3,3]*Nz6)) 
    Z32_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,3,1,2]*Nx6 + c[1,3,2,2]*Ny6 + c[1,3,3,2]*Nz6) + Ny6 * (c[2,3,1,2]*Nx6 + c[2,3,2,2]*Ny6 + c[2,3,3,2]*Nz6) + Nz6 * (c[3,3,1,2]*Nx6 + c[3,3,2,2]*Ny6 + c[3,3,3,2]*Nz6)) 
    Z33_6 = (beta/h1) * d * (EsJ6) * (Nx6 * (c[1,3,1,3]*Nx6 + c[1,3,2,3]*Ny6 + c[1,3,3,3]*Nz6) + Ny6 * (c[2,3,1,3]*Nx6 + c[2,3,2,3]*Ny6 + c[2,3,3,3]*Nz6) + Nz6 * (c[3,3,1,3]*Nx6 + c[3,3,2,3]*Ny6 + c[3,3,3,3]*Nz6)) 




    # Create Face Restriction Operators (e1T, e2T, etc) - these need to be checked
    eq0 = sparse([1  ], [1], [1], Nqp, 1)
    eqN = sparse([Nqp], [1], [1], Nqp, 1)
    er0 = sparse([1  ], [1], [1], Nrp, 1)
    erN = sparse([Nrp], [1], [1], Nrp, 1)
    es0 = sparse([1  ], [1], [1], Nsp, 1)
    esN = sparse([Nsp], [1], [1], Nsp, 1)
    e1 = Is ⊗ Ir ⊗ eq0
    e2 = Is ⊗ Ir ⊗ eqN
    e3 = Is ⊗ er0 ⊗ Iq
    e4 = Is ⊗ erN ⊗ Iq
    e5 = es0 ⊗ Ir ⊗ Iq
    e6 = esN ⊗ Ir ⊗ Iq
    e1T = Is ⊗ Ir ⊗ eq0'
    e2T = Is ⊗ Ir ⊗ eqN'
    e3T = Is ⊗ er0' ⊗ Iq
    e4T = Is ⊗ erN' ⊗ Iq
    e5T = es0' ⊗ Ir ⊗ Iq
    e6T = esN' ⊗ Ir ⊗ Iq

    n = 1 #default to n = 1
    # Create SAT vectors - all DIRICHLET conditions
    # S11 =  JHI * ((n*T11_1 .- Z11_1)'*e1*sJ1*H1*e1T + (n*T11_2 .- Z11_2)'*e2*sJ2*H2*e2T + (n*T11_3 .- Z11_3)'*e3*sJ3*H3*e3T + (n*T11_4 .- Z11_4)'*e4*sJ4*H4*e4T + (n*T11_5 .- Z11_5)'*e5*sJ5*H5*e5T + (n*T11_6 .- Z11_6)'*e6*sJ6*H6*e6T)
    # S12 =  JHI * ((n*T21_1 .- Z21_1)'*e1*sJ1*H1*e1T + (n*T21_2 .- Z21_2)'*e2*sJ2*H2*e2T + (n*T21_3 .- Z21_3)'*e3*sJ3*H3*e3T + (n*T21_4 .- Z21_4)'*e4*sJ4*H4*e4T + (n*T21_5 .- Z21_5)'*e5*sJ5*H5*e5T + (n*T21_6 .- Z21_6)'*e6*sJ6*H6*e6T)    
    # S13 =  JHI * ((n*T31_1 .- Z31_1)'*e1*sJ1*H1*e1T + (n*T31_2 .- Z31_2)'*e2*sJ2*H2*e2T + (n*T31_3 .- Z31_3)'*e3*sJ3*H3*e3T + (n*T31_4 .- Z31_4)'*e4*sJ4*H4*e4T + (n*T31_5 .- Z31_5)'*e5*sJ5*H5*e5T + (n*T31_6 .- Z31_6)'*e6*sJ6*H6*e6T)
    # b11 = [JHI*(n*T11_1 .- Z11_1)'*e1*sJ1*H1, JHI*(n*T11_2 .- Z11_2)'*e2*sJ2*H2, JHI*(n*T11_3 .- Z11_3)'*e3*sJ3*H3, JHI*(n*T11_4 .- Z11_4)'*e4*sJ4*H4, JHI*(n*T11_5 .- Z11_5)'*e5*sJ5*H5, JHI*(n*T11_6 .- Z11_6)'*e6*sJ6*H6]
    # b12 = [JHI*(n*T21_1 .- Z21_1)'*e1*sJ1*H1, JHI*(n*T21_2 .- Z21_2)'*e2*sJ2*H2, JHI*(n*T21_3 .- Z21_3)'*e3*sJ3*H3, JHI*(n*T21_4 .- Z21_4)'*e4*sJ4*H4, JHI*(n*T21_5 .- Z21_5)'*e5*sJ5*H5, JHI*(n*T21_6 .- Z21_6)'*e6*sJ6*H6]
    # b13 = [JHI*(n*T31_1 .- Z31_1)'*e1*sJ1*H1, JHI*(n*T31_2 .- Z31_2)'*e2*sJ2*H2, JHI*(n*T31_3 .- Z31_3)'*e3*sJ3*H3, JHI*(n*T31_4 .- Z31_4)'*e4*sJ4*H4, JHI*(n*T31_5 .- Z31_5)'*e5*sJ5*H5, JHI*(n*T31_6 .- Z31_6)'*e6*sJ6*H6]

    # S21 =  JHI * ((n*T12_1 .- Z12_1)'*e1*sJ1*H1*e1T + (n*T12_2 .- Z12_2)'*e2*sJ2*H2*e2T + (n*T12_3 .- Z12_3)'*e3*sJ3*H3*e3T + (n*T12_4 .- Z12_4)'*e4*sJ4*H4*e4T + (n*T12_5 .- Z12_5)'*e5*sJ5*H5*e5T + (n*T12_6 .- Z12_6)'*e6*sJ6*H6*e6T) 
    # S22 =  JHI * ((n*T22_1 .- Z22_1)'*e1*sJ1*H1*e1T + (n*T22_2 .- Z22_2)'*e2*sJ2*H2*e2T + (n*T22_3 .- Z22_3)'*e3*sJ3*H3*e3T + (n*T22_4 .- Z22_4)'*e4*sJ4*H4*e4T + (n*T22_5 .- Z22_5)'*e5*sJ5*H5*e5T + (n*T22_6 .- Z22_6)'*e6*sJ6*H6*e6T) 
    # S23 =  JHI * ((n*T32_1 .- Z32_1)'*e1*sJ1*H1*e1T + (n*T32_2 .- Z32_2)'*e2*sJ2*H2*e2T + (n*T32_3 .- Z32_3)'*e3*sJ3*H3*e3T + (n*T32_4 .- Z32_4)'*e4*sJ4*H4*e4T + (n*T32_5 .- Z32_5)'*e5*sJ5*H5*e5T + (n*T32_6 .- Z32_6)'*e6*sJ6*H6*e6T) 
    # b21 = [JHI*(n*T12_1 .- Z12_1)'*e1*sJ1*H1, JHI*(n*T12_2 .- Z12_2)'*e2*sJ2*H2, JHI*(n*T12_3 .- Z12_3)'*e3*sJ3*H3, JHI*(n*T12_4 .- Z12_4)'*e4*sJ4*H4, JHI*(n*T12_5 .- Z12_5)'*e5*sJ5*H5, JHI*(n*T12_6 .- Z12_6)'*e6*sJ6*H6]
    # b22 = [JHI*(n*T22_1 .- Z22_1)'*e1*sJ1*H1, JHI*(n*T22_2 .- Z22_2)'*e2*sJ2*H2, JHI*(n*T22_3 .- Z22_3)'*e3*sJ3*H3, JHI*(n*T22_4 .- Z22_4)'*e4*sJ4*H4, JHI*(n*T22_5 .- Z22_5)'*e5*sJ5*H5, JHI*(n*T22_6 .- Z22_6)'*e6*sJ6*H6]
    # b23 = [JHI*(n*T32_1 .- Z32_1)'*e1*sJ1*H1, JHI*(n*T32_2 .- Z32_2)'*e2*sJ2*H2, JHI*(n*T32_3 .- Z32_3)'*e3*sJ3*H3, JHI*(n*T32_4 .- Z32_4)'*e4*sJ4*H4, JHI*(n*T32_5 .- Z32_5)'*e5*sJ5*H5, JHI*(n*T32_6 .- Z32_6)'*e6*sJ6*H6]

    # S31 = JHI * ((n*T13_1 .- Z13_1)'*e1*sJ1*H1*e1T + (n*T13_2 .- Z13_2)'*e2*sJ2*H2*e2T + (n*T13_3 .- Z13_3)'*e3*sJ3*H3*e3T + (n*T13_4 .- Z13_4)'*e4*sJ4*H4*e4T + (n*T13_5 .- Z13_5)'*e5*sJ5*H5*e5T + (n*T13_6 .- Z13_6)'*e6*sJ6*H6*e6T)
    # S32 = JHI * ((n*T23_1 .- Z23_1)'*e1*sJ1*H1*e1T + (n*T23_2 .- Z23_2)'*e2*sJ2*H2*e2T + (n*T23_3 .- Z23_3)'*e3*sJ3*H3*e3T + (n*T23_4 .- Z23_4)'*e4*sJ4*H4*e4T + (n*T23_5 .- Z23_5)'*e5*sJ5*H5*e5T + (n*T23_6 .- Z23_6)'*e6*sJ6*H6*e6T)
    # S33 = JHI * ((n*T33_1 .- Z33_1)'*e1*sJ1*H1*e1T + (n*T33_2 .- Z33_2)'*e2*sJ2*H2*e2T + (n*T33_3 .- Z33_3)'*e3*sJ3*H3*e3T + (n*T33_4 .- Z33_4)'*e4*sJ4*H4*e4T + (n*T33_5 .- Z33_5)'*e5*sJ5*H5*e5T + (n*T33_6 .- Z33_6)'*e6*sJ6*H6*e6T)        
    # b31 = [JHI*(n*T13_1 .- Z13_1)'*e1*sJ1*H1, JHI*(n*T13_2 .- Z13_2)'*e2*sJ2*H2, JHI*(n*T13_3 .- Z13_3)'*e3*sJ3*H3, JHI*(n*T13_4 .- Z13_4)'*e4*sJ4*H4, JHI*(n*T13_5 .- Z13_5)'*e5*sJ5*H5, JHI*(n*T13_6 .- Z13_6)'*e6*sJ6*H6]
    # b32 = [JHI*(n*T23_1 .- Z23_1)'*e1*sJ1*H1, JHI*(n*T23_2 .- Z23_2)'*e2*sJ2*H2, JHI*(n*T23_3 .- Z23_3)'*e3*sJ3*H3, JHI*(n*T23_4 .- Z23_4)'*e4*sJ4*H4, JHI*(n*T23_5 .- Z23_5)'*e5*sJ5*H5, JHI*(n*T23_6 .- Z23_6)'*e6*sJ6*H6]
    # b33 = [JHI*(n*T33_1 .- Z33_1)'*e1*sJ1*H1, JHI*(n*T33_2 .- Z33_2)'*e2*sJ2*H2, JHI*(n*T33_3 .- Z33_3)'*e3*sJ3*H3, JHI*(n*T33_4 .- Z33_4)'*e4*sJ4*H4, JHI*(n*T33_5 .- Z33_5)'*e5*sJ5*H5, JHI*(n*T33_6 .- Z33_6)'*e6*sJ6*H6] 
         

    # Create SAT vectors - DIRICHLET conditions on faces 1 and 2, else traction

    S11 = JHI * ((T11_1 .- Z11_1)'*e1*sJ1*H1*e1T + (T11_2 .- Z11_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T11_3 + e4*sJ4*H4*e4'*T11_4 + e5*sJ5*H5*e5'*T11_5 + e6*sJ6*H6*e6'*T11_6)
    S12 = JHI * ((T21_1 .- Z21_1)'*e1*sJ1*H1*e1T + (T21_2 .- Z21_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T12_3 + e4*sJ4*H4*e4'*T12_4 + e5*sJ5*H5*e5'*T12_5 + e6*sJ6*H6*e6'*T12_6)
    S13 = JHI * ((T31_1 .- Z31_1)'*e1*sJ1*H1*e1T + (T31_2 .- Z31_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T13_3 + e4*sJ4*H4*e4'*T13_4 + e5*sJ5*H5*e5'*T13_5 + e6*sJ6*H6*e6'*T13_6)
    #S11 =  JHI * ((T11_1 .- Z11_1)'*e1*sJ1*H1*e1T + (T11_2 .- Z11_2)'*e2*sJ2*H2*e2T + (T11_3 .- Z11_3)'*e3*sJ3*H3*e3T + (T11_4 .- Z11_4)'*e4*sJ4*H4*e4T + (T11_5 .- Z11_5)'*e5*sJ5*H5*e5T + (T11_6 .- Z11_6)'*e6*sJ6*H6*e6T)
    #S12 =  JHI * ((T21_1 .- Z21_1)'*e1*sJ1*H1*e1T + (T21_2 .- Z21_2)'*e2*sJ2*H2*e2T + (T21_3 .- Z21_3)'*e3*sJ3*H3*e3T + (T21_4 .- Z21_4)'*e4*sJ4*H4*e4T + (T21_5 .- Z21_5)'*e5*sJ5*H5*e5T + (T21_6 .- Z21_6)'*e6*sJ6*H6*e6T)    
    #S13 =  JHI * ((T31_1 .- Z31_1)'*e1*sJ1*H1*e1T + (T31_2 .- Z31_2)'*e2*sJ2*H2*e2T + (T31_3 .- Z31_3)'*e3*sJ3*H3*e3T + (T31_4 .- Z31_4)'*e4*sJ4*H4*e4T + (T31_5 .- Z31_5)'*e5*sJ5*H5*e5T + (T31_6 .- Z31_6)'*e6*sJ6*H6*e6T)
    b11 = [JHI*(T11_1 .- Z11_1)'*e1*sJ1*H1, JHI*(T11_2 .- Z11_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b12 = [JHI*(T21_1 .- Z21_1)'*e1*sJ1*H1, JHI*(T21_2 .- Z21_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b13 = [JHI*(T31_1 .- Z31_1)'*e1*sJ1*H1, JHI*(T31_2 .- Z31_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]

    S21 =  JHI * ((T12_1 .- Z12_1)'*e1*sJ1*H1*e1T + (T12_2 .- Z12_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T21_3 + e4*sJ4*H4*e4'*T21_4 + e5*sJ5*H5*e5'*T21_5 + e6*sJ6*H6*e6'*T21_6)
    S22 =  JHI * ((T22_1 .- Z22_1)'*e1*sJ1*H1*e1T + (T22_2 .- Z22_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T22_3 + e4*sJ4*H4*e4'*T22_4 + e5*sJ5*H5*e5'*T22_5 + e6*sJ6*H6*e6'*T22_6)
    S23 =  JHI * ((T32_1 .- Z32_1)'*e1*sJ1*H1*e1T + (T32_2 .- Z32_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T23_3 + e4*sJ4*H4*e4'*T23_4 + e5*sJ5*H5*e5'*T23_5 + e6*sJ6*H6*e6'*T23_6)
    b21 = [JHI*(T12_1 .- Z12_1)'*e1*sJ1*H1, JHI*(T12_2 .- Z12_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b22 = [JHI*(T22_1 .- Z22_1)'*e1*sJ1*H1, JHI*(T22_2 .- Z22_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b23 = [JHI*(T32_1 .- Z32_1)'*e1*sJ1*H1, JHI*(T32_2 .- Z32_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]

    S31 = JHI * ((T13_1 .- Z13_1)'*e1*sJ1*H1*e1T + (T13_2 .- Z13_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T31_3 + e4*sJ4*H4*e4'*T31_4 + e5*sJ5*H5*e5'*T31_5 + e6*sJ6*H6*e6'*T31_6)
    S32 = JHI * ((T23_1 .- Z23_1)'*e1*sJ1*H1*e1T + (T23_2 .- Z23_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T32_3 + e4*sJ4*H4*e4'*T32_4 + e5*sJ5*H5*e5'*T32_5 + e6*sJ6*H6*e6'*T32_6)
    S33 = JHI * ((T33_1 .- Z33_1)'*e1*sJ1*H1*e1T + (T33_2 .- Z33_2)'*e2*sJ2*H2*e2T) - JHI * (e3*sJ3*H3*e3'*T33_3 + e4*sJ4*H4*e4'*T33_4 + e5*sJ5*H5*e5'*T33_5 + e6*sJ6*H6*e6'*T33_6)        
    b31 = [JHI*(T13_1 .- Z13_1)'*e1*sJ1*H1, JHI*(T13_2 .- Z13_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b32 = [JHI*(T23_1 .- Z23_1)'*e1*sJ1*H1, JHI*(T23_2 .- Z23_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)]
    b33 = [JHI*(T33_1 .- Z33_1)'*e1*sJ1*H1, JHI*(T33_2 .- Z33_2)'*e2*sJ2*H2, -JHI*(e3*sJ3*H3), -JHI*(e4*sJ4*H4), -JHI*(e5*sJ5*H5), -JHI*(e6*sJ6*H6)] 
         


    # Equation 1: J = 1
    # 0 = A11*u1 + A12*u2 + A13*u3 + f1 + SAT1
    # OR 
    # 0 = A11*u1 + A12*u2 + A13*u3 + f1 + S11*u1 + S12*u2 + S13*u3 - b1
     
    # Equation 2: J = 2
    # 0 = A21*u1 + A22*u2 + A23*u3 + f2 + SAT2
    # OR
    # 0 = A21*u1 + A22*u2 + A23*u3 + f2 + S21*u1 + S22*u2 + S23*u3 - b3

    # Equation 3: J = 3
    # 0 = A31*u1 + A32*u2 + A33*u3 + f3 + SAT3
    # OR
    # 0 = A31*u1 + A32*u2 + A33*u3 + f3 + S31*u1 + S32*u2 + S33*u3

    # AND ALL TOGETHER: MU = [B11*g1 + B12*g2 + B13*g3; B21*g1 + B22*g2 + B23*g3; B31*g1 + B32*g2 + B33*g3] + J*H*f  where
    A = [A11 A12 A13; A21 A22 A23; A31 A32 A33]
    S = [S11 S12 S13; S21 S22 S23; S31 S32 S33]
    
    M = A + S
    @show sizeof(A)
    @show sizeof(S)
    
    B = (b11, 1*b12, 1*b13, 1*b21, b22, 1*b23, 1*b31, 1*b32, b33)
   
    # and f = [f1; f2; f3]
    # where U = [u1; u2; u3]
    JH = J*H
    
    return (M, B, JH, A, S, HqI, HrI, HsI)


end


function var_3D_D2q(p, Nqp, Nrp, Nsp, C, HIq; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    Iq = sparse(I, Nqp, Nqp)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)

    N = Nqp*Nrp*Nsp
    D2q = spzeros(N, N) # initialize
    S0q = spzeros(N, N) # initialize
    SNq = spzeros(N, N) # initialize
    for i = 1:Nrp
        for j = 1:Nsp
            B = C[:, i, j]# get coefficient on 1D line in q-direction
            (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nqp-1, B; xc = (-1,1))
            ej = spzeros(Nsp, 1)
            ej[j] = 1
            ei = spzeros(Nrp, 1)
            ei[i] = 1
            D2q += (ej ⊗ ei ⊗ Iq) * D2 * (ej' ⊗ ei' ⊗ Iq)
            S0q += (ej ⊗ ei ⊗ Iq) * S0 * (ej' ⊗ ei' ⊗ Iq)
            SNq += (ej ⊗ ei ⊗ Iq) * SN * (ej' ⊗ ei' ⊗ Iq)
        end
    end
    return D2q, S0q, SNq
end

function var_3D_D2r(p, Nqp, Nrp, Nsp, C, HIr; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    Iq = sparse(I, Nqp, Nqp)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)

    N = Nqp*Nrp*Nsp
    D2r = spzeros(N, N) # initialize
    S0r = spzeros(N, N) # initialize
    SNr = spzeros(N, N) # initialize
    for i = 1:Nqp
        for j = 1:Nsp
            B = C[i, :, j]# get coefficient on 1D line in r-direction
            (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nrp-1, B; xc = (-1,1))
            ej = spzeros(Nsp, 1)
            ej[j] = 1
            ei = spzeros(Nqp, 1)
            ei[i] = 1
            D2r += (ej ⊗ Ir ⊗ ei) * D2 * (ej' ⊗ Ir ⊗ ei')
            S0r += (ej ⊗ Ir ⊗ ei) * S0 * (ej' ⊗ Ir ⊗ ei')
            SNr += (ej ⊗ Ir ⊗ ei) * SN * (ej' ⊗ Ir ⊗ ei')
        end
    end
    return D2r, S0r, SNr
end

function var_3D_D2s(p, Nqp, Nrp, Nsp, C, HIq; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    Iq = sparse(I, Nqp, Nqp)
    Ir = sparse(I, Nrp, Nrp)
    Is = sparse(I, Nsp, Nsp)

    N = Nqp*Nrp*Nsp
    D2s = spzeros(N, N) # initialize
    S0s = spzeros(N, N) # initialize
    SNs = spzeros(N, N) # initialize
    for i = 1:Nrp
        for j = 1:Nqp
            B = C[j, i, :]# get coefficient on 1D line in s-direction
            (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nsp-1, B; xc = (-1,1))
            ej = spzeros(Nqp, 1)
            ej[j] = 1
            ei = spzeros(Nrp, 1)
            ei[i] = 1
            D2s += (Is ⊗ ei ⊗ ej) * D2 * (Is ⊗ ei' ⊗ ej')
            S0s += (Is ⊗ ei ⊗ ej) * S0 * (Is ⊗ ei' ⊗ ej')
            SNs += (Is ⊗ ei ⊗ ej) * SN * (Is ⊗ ei' ⊗ ej')
        end
    end
    return D2s, S0s, SNs
end


  export setupfaultstations, savefaultstation
  export read_params, read_params_structured, plot_slip, plot_fault_time_series, find_ind 
  export find_station_index, read_inp_3d
  export stepcheck
  export animate_slip, interp1
  export better_plot_solution, plot_solution
  export create_text_files, write_to_file
  export read_params_BP1_CT, read_params_BP6, write_to_file_BP6, create_text_files_BP6
  export create_metrics_3D, transforms_ne, locoperator_3D