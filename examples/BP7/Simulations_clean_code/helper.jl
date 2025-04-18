# Helper functions for BP7-QD and BP7-FD simulations
using DifferentialEquations
using Printf
using Plots
using CUDA
using DelimitedFiles
using JSON

function _unpack(x::NamedTuple)
    kk = keys(x)
    vv = values(x)
    i = 1
    for k in kk
        @eval $k = $vv[$i]
        i += 1
    end
end


macro unpack_namedtuple(arg)
    quote
        _unpack($arg)
    end |> esc
end


struct coefficients
    # Table 1: Parameter values used in this benchmark problem 
    #        # variable names                                   # recommended values                
    ρ       # density                                           2670 kg/m³
    cs      # shear wave speed                                  3.464 km/s
    μ       # shear modulus                                     32.04 GPa
    ν       # Poisson's ratio                                   0.25
    a0      # rate-and-state parameter                          0.004
    amax    # rate-and-state parameter                          0.016
    b0      # rate-and-state parameter                          0.01
    σn      # effective normal stress                           25 MPa
    DRS     # characteristic slip distance                      0.50 / 0.53 mm
    Vp      # plate rate                                        10⁻⁹ m/s
    VL      # imposed fault slip rate                           10⁻⁹ m/s
    Vinit   # initial slip rate                                 10⁻⁹ m/s
    V0      # reference slip rate                               10⁻⁶ m/s 
    f0      # reference friction coefficient                    0.6
    RVW     # radius of uniform VW region                       200 m
    Wf      # half-width of rate-and-state fault                400 m
    lf      # half-length of rate-and-state fault               400 m
    Δz      # suggested cell size                               1000m
    tf      # final simulation time                             10 years 
    Δτ0     # max amplitude of nucleation stress perturbation   2 MPa 
    Rnuc    # radius of the nucleation stress perturbation      150 m
    T       # Duration of the nucleation stress perturbation    1s
    y2      # hypocenter for nucleation stress perturbation    -50
    y3      # hypocenter for nucleation stress perturbation     -50
end

# default constructor
coefficients() = coefficients(
    2670,                   # ρ
    3.464,                  # cs
    32.04,                   # μ
    0.25,                   # ν
    0.004,                  # a0
    0.016,                   # amax
    0.01,                   # b0            value for b in this problem
    25,                     # σn
    0.50,                   # DRS
    1E-9,                   # Vp
    1E-9,                    # VL
    1E-9,                   # Vinit
    1E-6,                   # V0
    0.6,                    # f0
    200,                    # RVW
    400,                    # Wf
    400,                    # lf
    10,                     # Δz in meter, 
    10,                     # tf
    2,                      # Δτ0
    150,                    # Rnuc
    1,                      # T 
    -50,                    # y2
    -50                     # y3
)

# initial state variable over the entire fault
# x2 
function θ0_func(Nx, Ny, BP7_coeff::coefficients)
    return fill(BP7_coeff.DRS / BP7_coeff.Vinit, Nx * Ny)
end

# for BP7, b is set to be the constant value b0
# using b in the variables of the functions below 
# to write more modular code easier to maintain

function f_func(V, θ, a, b, BP7_coeff::coefficients)
    return a * asinh(
        V / (2 * BP7_coeff.V0) * exp((BP7_coeff.f0 + b * ln(BP7_coeff.V0 / BP7_coeff.DRS)) / a) # is b the value b0 in coefficients?
    )
end

# TODO BP7 is homogeneous and isotropic, a_func and a_func_region is replaced
# function a_func(x2, x3, BP7_coeff::coefficients)
#     if (BP7_coeff.hs + BP7_coeff.ht ≤ x3 ≤ BP7_coeff.hs + BP7_coeff.ht + BP7_coeff.H) && (abs(x2) ≤ BP7_coeff.DRS / 2)
#         return BP7_coeff.a0
#     elseif (0 ≤ x3 ≤ BP7_coeff.hs) || (BP7_coeff.hs + 2 * BP7_coeff.ht + BP7_coeff.H ≤ x3 ≤ BP7_coeff.Wf) || (BP7_coeff.DRS / 2 + BP7_coeff.ht ≤ abs(x2) ≤ BP7_coeff.DRSf / 2)
#         return BP7_coeff.amax
#     else
#         r = max(abs(x3 - BP7_coeff.hs - BP7_coeff.ht - BP7_coeff.H / 2) - BP7_coeff.H / 2, abs(x2) - BP7_coeff.DRS / 2) / BP7_coeff.ht
#         return BP7_coeff.a0 + r * (BP7_coeff.amax - BP7_coeff.a0)
#     end
# end


function G1_func(r, Rnuc)
    if r < Rnuc
        return exp(r^2 / (r^2 - Rnuc^2))
    else
        return 0
    end
end

function G2_func(t, T)
    if t < T
        return exp((t - T)^2 / (t*(t - 2 * T)))
    else
        return 1
    end
end

# auxiliary function to determine which region a belongs to
# function a_func_region(x2, x3, BP7_coeff::coefficients)
#     if (BP7_coeff.hs + BP7_coeff.ht ≤ x3 ≤ BP7_coeff.hs + BP7_coeff.ht + BP7_coeff.H) && (abs(x2) ≤ BP7_coeff.DRS / 2)
#         return 0
#     elseif (0 ≤ x3 ≤ BP7_coeff.hs) || (BP7_coeff.hs + 2 * BP7_coeff.ht + BP7_coeff.H ≤ x3 ≤ BP7_coeff.Wf) || (BP7_coeff.DRS / 2 + BP7_coeff.ht ≤ abs(x2) ≤ BP7_coeff.DRSf / 2)
#         return 2
#     else
#         r = max(abs(x3 - BP7_coeff.hs - BP7_coeff.ht - BP7_coeff.H / 2) - BP7_coeff.H / 2, abs(x2) - BP7_coeff.DRS / 2) / BP7_coeff.ht
#         return 1
#     end
# end

# For BP7-QD, the scalar pre-stress τ⁰ is chosen as the steady-state stress
function τ0_QD_func(a, b, η, BP7_coeff::coefficients)
    return BP7_coeff.σn * a * asinh(BP7_coeff.Vinit / (2 * BP7_coeff.V0) * exp((BP7_coeff.f0 + b) / a)) + η * BP7_coeff.Vinit
end

# For BP7-QD, the scalar pre-stress τ⁰ is chosen as the steady-state stress
function τ0_FD_func(a, b, BP7_coeff::coefficients)
    return BP7_coeff.σn * asinh(BP7_coeff.Vinit / (2 * BP7_coeff.V0) * exp((BP7_coeff.f0 + b) / a))
end

# a higher pre-stress along the x2-direction
function τi0_FD_func(a, b, δ, τ, BP7_coeff::coefficients)
    return BP7_coeff.σn * asinh(BP7_coeff.Vinit / (2 * BP7_coeff.V0) * exp((BP7_coeff.f0 + b) / a)) + δ * τ
end


# quasi-static process zone
function Λ0_func(C, b, μ, BP7_coeff::coefficients)
    return C * μ * BP7_coeff.DRS / (b * BP7_coeff.σn)
end

# nucleation zone
function h_func(a, b, μ, BP7_coeff::coefficients)
    return π / 2 * (μ * b * BP7_coeff.DRS) / ((b - a)^2 * BP7_coeff.σn)^2
end

# fault strength
function F_func(f, Vbold, BP7_coeff::coefficients)
    return BP7_coeff.σn * f * Vbold / norm(V)
end


# boundary functions
# Dirichlet boundary conditions
function bc_Dirichlet(face, Vp, δ, x, t)
    if face == 1 # check this
        return (δ ./ 2)
    elseif face == 2 # check the 
        return fill(t * Vp / 2, size(x))
    end
end

# Neumann boundary conditions
function bc_Neumann(face, Vp, δ, x, t)
    return zeros(size(x))
end


# update boundary conditions
function boundary_update!(RHS, bc_Dirichlet, bc_Neumann)
    # TODO
end



function rateandstate(V2, V3, psi, σn, τ2, τ3, η, a, V0)

    V = sqrt(V2^2 + V3^2)
    dV_dV2 = 0.5*(V2^2 + V3^2)^(-0.5) .* 2 .* V2
    dV_dV3 = 0.5*(V2^2 + V3^2)^(-0.5) .* 2 .* V3
    
    Y = (1 ./ (2 .* V0)) .* exp.(psi ./ a)
    f = a .* asinh.(V .* Y)  # compute friction coefficient 
    df_dV2  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* (dV_dV2 .* Y)  # derivative wrt V_2
    df_dV3  = a .* (1 ./ sqrt.(1 + (V .* Y).^2)) .* (dV_dV3 .* Y)  # derivative wrt V_2

    g1 = σn .* f .* V2 / V   + η .* V2 - τ2
    g2 = σn .* f .* V3 / V   + η .* V3 - τ3
   

    A2 = V2/V
    A3 = V3/V

    dA2_dV2 = (V - V2*dV_dV2)/V^2
    dA2_dV3 = (-V2*dV_dV3)/V^2
    dA3_dV2 = (-V3*dV_dV2)/V^2
    dA3_dV3 = (V - V3*dV_dV3)/V^2

    dg1_dV2 = σn .* (df_dV2 .* V2 / V + f .* dA2_dV2) + η
    dg1_dV3 = σn .* (df_dV3 .* V2 / V + f .* dA2_dV3)

    dg2_dV2 = σn .* (df_dV2 .* V3 / V + f .* dA3_dV2) 
    dg2_dV3 = σn .* (df_dV3 .* V3 / V + f .* dA3_dV3) + η 

    return (g1, g2, dg1_dV2, dg1_dV3, dg2_dV2, dg2_dV3)
end
  

function newtbndv(func, x, y; ftol = 1e-12, maxiter = 500, 
                    atolx = 1e-4, rtolx = 1e-4)

    (f, g, dfx, dfy, dgx, dgy) = func(x, y)
    for iter = 1:maxiter

        z = [x; y] 
        (f, g, dfx, dfy, dgx, dgy) = func(x, y)
      
        J = [dfx dfy; dgx dgy] 
        dx, dy = -J\[f; g]
      
        x = x + dx
        y = y + dy
  
       if abs(f) < ftol && abs(dx) < atolx + rtolx * (abs(dx) + abs(x)) && abs(g) < ftol && abs(dy) < atolx + rtolx * (abs(dy) + abs(y))
            
            return (x, y, f, g, iter)
       end
    end
    return (x, y, f, g, -maxiter)
end


function rateandstate_vectorized(V_v, ψ, σn, τ_v, η, RSas, RSV0)
    # V and τ both stand for absolute value of slip rate and traction vecxtors. 
    Y_v = (1 ./ (2 .* RSV0)) .* exp.(ψ ./ RSas)
    f_v = RSas .* asinh.(V_v .* Y_v)
    dfdV_v = RSas .* (1 ./ sqrt.(1 .+ (V_v .* Y_v) .^ 2)) .* Y_v
  
    g_v = σn .* f_v .+ η .* V_v .- τ_v
    dgdV_v = σn .* dfdV_v .+ η
    return (g_v, dgdV_v)
end

function newtbndv_vectorized(rateandstate_vectorized, xL, xR, V_v, ψ, σn, τ_v, η, 
                        RSas, RSV0; ftol=1e-6, maxiter = 500, minchange = 0, atolx = 1e-4, rtolx=1e-4)
    fL_v = rateandstate_vectorized(xL, ψ, σn, τ_v, η, RSas, RSV0)[1]
    fR_v = rateandstate_vectorized(xR, ψ, σn, τ_v, η, RSas, RSV0)[1]

    if any(x -> x > 0, fL_v .* fR_v)
        # return (fill(typeof(V_v)(NaN), length(V_v)), fill(typeof(V_v)(NaN), length(V_v)), -maxiter)
        return (fill(NaN, length(V_v)), fill(NaN, length(V_v)), -maxiter)
        # typeof(V_v)(NaN) no longer supported in Julia
    end

    f_v = rateandstate_vectorized(V_v, ψ, σn, τ_v, η, RSas, RSV0)[1]
    df_v = rateandstate_vectorized(V_v, ψ, σn, τ_v, η, RSas, RSV0)[2]
    dxlr_v = xR .- xL

 

    for iter = 1:maxiter
        dV_v = -f_v ./ df_v
        V_v = V_v .+ dV_v
        
        mask = (V_v .< xL) .| (V_v .> xR) .| (abs.(dV_v) ./ dxlr_v .< minchange)
        V_v[mask] .= (xR[mask] .+ xL[mask]) ./ 2
        dV_v[mask] .= (xR[mask] .- xL[mask]) ./ 2

        f_v = rateandstate_vectorized(V_v, ψ, σn, τ_v, η, RSas, RSV0)[1]
        df_v = rateandstate_vectorized(V_v, ψ, σn, τ_v, η, RSas, RSV0)[2]
        
        mask_2 = f_v .* fL_v .> 0
        fL_v[mask_2] .= f_v[mask_2]
        xL[mask_2] .= V_v[mask_2]
        fR_v[.!mask_2] .= f_v[.!mask_2]
        xR[.!mask_2] .= V_v[.!mask_2]

        dxlr_v .= xR .- xL

        if all(abs.(f_v) .< ftol) && all(abs.(dV_v .< atolx .+ rtolx .* (abs.(dV_v) .+ abs.(V_v))))
            return (V_v, f_v, iter)
        end
    end
    return (V_v, f_v, -maxiter)

end


function map_jacobian(x, y, a, b)
    return [x y; a b]
end

function map_jacobian_inv(x, y, a, b)
    return inv([x y; a b])
end

function map_cat(x,y)
    return [x;y]
end

function get_first(a)
    return a[1]
end

function get_second(a)
    return a[2]
end

# Auxiliary functions for debugging
# Plut RS region
function plot_RS(vector)
    spy(reshape(vector, fN2, fN3)')
end
# Plot favorible region
function plot_f(vector)
    spy(reshape(vector, fN2_VW_favorable, fN3_VW)')
end


# Plot the slip in 2D from BP1 problem
function plot_slip(S, δNp, yf, stride_time)

    m = length(yf)
    no_time_steps = size(S.t)
    slip_final = S.u[end][end]

    for i = 1:stride_time:no_time_steps[1]

        slip_t = S.u[i][δNp+1:end] # slip at time t
        #pyplot()
        display(plot(slip_t, -yf, xtickfont=font(18),
            ytickfont=font(18),
            guidefont=font(18),
            legendfont=font(18), ylabel="Depth (km)", xlabel="Slip (m)", xlims=(0, slip_final)))
        sleep(0.1)
    end

    #nothing
end

# Plot slipt against time 
function plot_time_series(path, filenm, region, fieldname, col; headerlines=16, figuresize=(800,600))
    default(fontfamily = "Times", xtickfont = font(20), ytickfont = font(20), legendfont = font(20), guidefont = font(20))

    filename = string(path, filenm)
    A = readdlm(filename, skipstart=headerlines)

    if region == "fault"
        T = A[:, 1] ./ (60*60*24*30) # seconds to years
        slip = A[:, 2]
        V = A[:, 3]
        tau = A[:, 4]
        state = A[:, 5]

        if fieldname == "slip"
            plot(T, slip, color=col, linewidth=2, size=figuresize)
        elseif fieldname == "V"
            plot(T, V, color=col, linewidth=2, size=figuresize)
        elseif fieldname == "10toV"
            plot(T, 10 .^ V, color=col, linewidth=2, size=figuresize)
        elseif fieldname == "tau"
            plot(T, tau, color=col, linewidth=2, size=figuresize)
        else
            plot(T, state, color=col, linewidth=2, size=figuresize)
        end

    elseif region == "surface"
        T = A[:, 1] ./ (60*60*24*30) # seconds to years
        horiz = A[:, 2]
        vert = A[:, 3]
        Vhoriz = A[:, 4]
        Vvert = A[:, 5]

        if fieldname == "horiz"
            plot(T, horiz, color=col)
        elseif fieldname == "vert"
            plot(T, vert, color=col)
        elseif fieldname == "Vhoriz"
            plot(T, Vhoriz, color=col)
        else
            plot(T, Vvert, color=col)
        end
    end

    xlabel!("Time (months)")
    ylabel!(fieldname)
    output_figure_name = split(filename,'.')[1]
    savefig("$output_figure_name.png")
    return
end


########################## Test functions, subjected to changes ################################


function test_unpack(p)
    @unpack_namedtuple p
    @show Vp
    @show RHS
end


function create_text_files(path, station_strings, station_indices, δ, τb, θ, t)
    Vzero = 1e-20
    for n = 1:length(station_strings)
        XXX = path * "fltst_strk" * station_strings[n] * ".txt"
        ww = Array{Float64}(undef, 1, 8)
        ww[1] = t
        RS_index = RS_filter_2D_nzind[station_indices[n]]
        ww[2] = δ[2 * RS_index-1]
        ww[3] = δ[2 * RS_index]
        # ww[4] = log10(BP7_coeff.Vinit)
        ww[4] = log10(abs(V2_v[station_indices[n]]))
        ww[5] = log10(Vzero)
        ww[6] = τb[2 * RS_index - 1] # need to define this
        ww[7] = τb[2 * RS_index]
        ww[8] = log10(θ[station_indices[n]])  # 
        open(XXX, "w") do io
            write(io, "# This is the file header\n")
            write(io, "# problem=SEAS Benchmark BP7-QD\n")  # 
            write(io, "# code=Thrase\n")
            write(io, "# modeler=A. Chen, B. A. Erickson\n")
            write(io, "# date=$(Dates.format(now(), "yyyymmdd"))\n")
            write(io, "# element size=10 m\n")
            write(io, "#location=on fault, 0km along strike, 8km away from the fault, 0km depth")
            write(io, "# minimum_time_step=0.1\n")
            write(io, "# maximum_time_step=3.157e6\n")
            write(io, "# num_time_steps=2400\n")
            write(io, "# Column #1 = Time(s)\n")
            write(io, "# Column #2 = Slip_2(m)\n")
            write(io, "# Column #3 = Slip_3(m)\n")
            write(io, "# Column #4 = Slip_rate_2(log10 m/s)\n")
            write(io, "# Column #5 = Slip_rate_3(log10 m/s)\n")
            write(io, "# Column #6 = Shear_stress_2 (MPa)\n")
            write(io, "# Column #7 = Shear_stress_3 (MPa)\n")
            write(io, "# Column #8 = State (log10 s)\n")
            write(io, "t\t slip_2\t slip_3\t slip_rate_2\t slip_rate_3\t shear_stress_2\t shear_stress_3\t state\n")
            write(io, "# Here is the time-series data.\n")
            writedlm(io, ww)
        end
    end

    ww = Array{Float64}(undef, 1, 4)
    ww[1] = t
    # V2_A = V[2 * A_filter_2D_nzind - 1]
    # V3_A = V[2 * A_filter_2D]

    # V_A = hypot.(V2_A, V3_A)

    # V2_VW = V[2 * VW_filter_2D_nzind - 1]
    # V3_VW = V[2 * VW_filter_2D_nzind]

    # V_VW = hypot.(V2_VW, V3_VW)

    # ww[2] = maximum(V_A)
    # ww[3] = mean(V_A) * π * (1.5 * BP7_coeff.RVW)^2 * BP7_coeff.μ * 10^9
    # ww[3] = mean(V_VW) * π * (BP7_coeff.RVW) ^ 2 * BP7_coeff.μ * 10^9
    ww[2] = 0.0
    ww[3] = 0.0
    ww[4] = 0.0
        
    open(path * "global.dat","w") do io
        write(io, "# This is the file header: \n")
        write(io, "# problem=SEAS Benchmark BP7-QD-A\n")
        write(io, "# code=Thrase\n")
        write(io, "# modeler = A. Chen, B. A. Erickson\n")
        write(io, "# date=$(Dates.format(now(), "yyyymmdd"))\n")
        write(io, "# element_size = 10 m, DRS=$(BP7_coeff.DRS), Lx = $(Lx)\n")
        write(io, "# location = 1.5 times VW patch radius\n")
        write(io, "# minimum_time_step=1e-10\n")
        write(io, "# maximum_time_step=undef\n")
        write(io, "# num_time_steps=undef\n")
        write(io, "# Column #1 = Time (s)\n")
        write(io, "# Column #2 = Max_slip_rate (log10 m/s)\n")
        write(io, "# Column #3 = Moment_rate (N-m/s) \n")
        write(io, "# Column #4 = Moment_rate of VW patch (N-m/s)\n")
        write(io, "# The line below lists the names of the data fields\n")
        write(io, "t\t max_slip_rate\t moment_rate\t moment_rate_vw\n")
        write(io, "# Here is the time-seris data.\n")
        writedlm(io, ww)
    end
end



function write_to_file(path, ψδ, t, i, odeparam, station_strings, station_indices)
    @unpack_namedtuple odeparam;
    Vmax = 0.0
    # if t == (sim_years ./ 31556926)
   
    if isdefined(i, :fsallast)
        dψV = i.fsallast
        dψ, V, ψ, δ = create_view(dψV, ψδ)
        if mod(ctr[], odeparam.stride_time) == 0
            for n = 1:length(station_strings)
                XXX = path * "fltst_strk" * station_strings[n] * ".txt"
                ww = Array{Float64}(undef, 1, 8)
                ww[1] = t
                RS_index = RS_filter_2D_nzind[station_indices[n]] # Check if should be just station_indices[n]
                ww[2] = δ[2 * RS_index-1]
                ww[3] = δ[2 * RS_index]
                V2_real = abs((V[2 * RS_index-1]))
                V3_real = abs((V[2 * RS_index]))
                ww[4] = log10(V2_real)
                ww[5] = log10(V3_real)
                ww[6] = τfb[2 * RS_index - 1] # need to define this
                ww[7] = τfb[2 * RS_index]
                θ = BP7_coeff.DRS / BP7_coeff.V0 * exp((ψ[station_indices[n]] - BP7_coeff.f0)/BP7_coeff.b0)
                ww[8] = log10(θ)
                open(XXX, "a") do io
                    writedlm(io, ww)
                end
            end

            ww = Array{Float64}(undef, 1, 4)
            ww[1] = t
            V2_A = V[2 * A_filter_2D_nzind .- 1]
            V3_A = V[2 * A_filter_2D_nzind]
    
            V_A = hypot.(V2_A, V3_A)
    
            V2_VW = V[2 * VW_filter_2D_nzind .- 1]
            V3_VW = V[2 * VW_filter_2D_nzind]
    
            V_VW = hypot.(V2_VW, V3_VW)
    
            ww[2] = log(10, maximum(V_A))
            ww[3] = mean(V_A) * π * (1.5 * BP7_coeff.RVW)^2 * BP7_coeff.μ * 10^9
            ww[4] = mean(V_VW) * π * (BP7_coeff.RVW) ^ 2 * BP7_coeff.μ * 10^9
            open(path * "global.dat", "a") do io 
                writedlm(io, ww)
            end
        end
        if mod(ctr[], 1000) == 0
            restart_values = (t = t, 
            ψ = ψ,
            δ = δ,
            V = V)
            open(path * "restart_values.json", "w") do file
                JSON.print(file, restart_values)
            end
        end
        global ctr[] += 1
        @show ctr[]
    end
    Vmax
end


# function find_value_in_file(filename, value)
#     # Open to file
#     open(filename) do file
#         for line in each line(file)
#             columns = split(line)
#             if columns[i] == value
#         end
#     end
# end