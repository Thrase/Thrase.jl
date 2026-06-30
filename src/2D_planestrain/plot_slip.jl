
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

# plot_slip will plot slip contours every 5 years in blue during interseismic, 
# every 1 second in red during coseismic
function plot_slip(filename; headerlines=1, vert_limits = (-40, 0))

  grid = readdlm(filename, Float64, skipstart=headerlines)
  sz = size(grid)
  flt_loc = grid[1,3:end]
  if flt_loc[end] > 0
    flt_loc = -flt_loc
  end

  T = grid[2:sz[1],1]
  maxV = grid[2:end, 2]
  slip = grid[2:sz[1], 3:sz[2]]
  N = size(slip)[2]

  maxV = maxV[625:end]
  T = T[625:end]
  slip = slip[625:end,:]
  slip0 = 0 .* slip[1,:]

  ind = find_ind(maxV);        #finds indices for inter/co-seismic phases
  interval = [5*31556926 1]   #plot every 5 years and every 1 second
  ct = 0   #this counts the number of events


  #Assumes an initial interseismic period
  #This for-loop only plots completed phases
  for i = 1:2:length(ind)-2
    
    T1 = T[ind[i]]:interval[1]:T[ind[i+1]];

    W1 = interp1(T,slip[:,1],T1)';
    
    for j = 2:N 
      w1 = interp1(T,slip[:,j] .- slip0[j],T1)';

      W1 = [W1; w1]
    end

    if i == 1
      plt.plot(W1, flt_loc, color = "blue")#, ls = "none", marker = ".", markersize = "0.5"); #interseismic phase
    else
      plt.plot(W1, flt_loc, color = "blue")#, ls = "none", marker = ".", markersize = "0.5"); #interseismic phase
    end

   
    T1 = T[ind[i+1]]:interval[2]:T[ind[i+2]];


    W1 = interp1(T,slip[:,1],T1)';
    for j = 2:N 
      w1 = interp1(T,slip[:,j] .- slip0[j],T1)';
      W1 = [W1; w1]
    end

    plt.plot(W1, flt_loc, color = "red")#, ls = "none", marker = ".", markersize = "0.5"); #coseismic phase

    ct = ct+1;
  end

  
  # plot remainder of an incomplete interseismic period:
  i = length(ind)-1;
  T1 = T[ind[i]]:interval[1]:T[ind[i+1]];
  W1 = interp1(T,slip[:,1],T1)';
      
      for j = 2:N 
      w1 = interp1(T,slip[:,j] .- slip0[j],T1)';
        W1 = [W1; w1]
      end
      if i == 1
        plt.plot(W1, flt_loc, color = "blue", ls = "none", marker = ".", markersize = "0.5") #interseismic phase
      else
        plt.plot(W1, flt_loc, color = "blue", ls = "none", marker = ".", markersize = "0.5") #interseismic phase
      end

      plt.xlabel("Cumulative Slip (m)");
      plt.ylabel("Distance Down Dip (km)");
      plt.ylim(vert_limits)
end
