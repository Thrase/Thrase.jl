# In development for reading unstructured meshes. 


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
  

read_inp_2d("square_circle.inp")
