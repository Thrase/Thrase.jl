# This file specifies the computational domain
using Dates
include("helper.jl")
const year_seconds = 31556926
sim_years = 1800

# calling parameterized constructor to set values for BP5
# if BP5_coeff is not defined here, the BP5-QD.jl will 
# call default constructor to construct BP5_coeff using
# values in helper.jl (line 51)
BP5_coeff = coefficients(
    2670,                   # ρ
    3.464,                  # cs
    0.25,                   # ν
    0.004,                  # a0 default value 0.004
    0.04,                   # amax default value 0.04
    0.03,                   # b0
    25,                     # σn
    0.14,                   # L or Dc 
    1E-9,                   # Vp
    1E-9,                   # Vinit initial value 1E-9, change it to 1E-3
    1E-6,                   # V0
    0.6,                    # f0
    2,                      # hs
    2,                      # ht
    12,                     # H
    60,                     # l
    40,                     # Wf
    100,                    # lf
    12,                     # w
    1000,                   # Δz
    1800                    # tf
)


include("Assembling_3D_matrices.jl")

# The entire domain is 128 km by 128 km by 128 km
Lx = Ly = Lz = 128
N_x = N_y = N_z = Int(128 / (BP5_coeff.Δz / 1000))
Nx = N_x + 1
Ny = N_y + 1
Nz = N_z + 1
n_levels = Int(log(2, N_x))

u1_filter_matrix = get_u1(Nx, Ny, Nz)
u2_filter_matrix = get_u2(Nx, Ny, Nz)
u3_filter_matrix = get_u3(Nx, Ny, Nz)


fN2 = BP5_coeff.lf / (BP5_coeff.Δz / 1000) + 1
fN2 = round(Int, fN2, RoundUp)
fN3 = BP5_coeff.Wf / (BP5_coeff.Δz / 1000) + 1
fN3 = round(Int, fN3, RoundUp)

# fNy and fNz represents the indices of the fault region
fNy_start = (Ly - BP5_coeff.lf) / (2 * BP5_coeff.Δz / 1000) + 1
fNy_start = round(Int, fNy_start, RoundUp)
fNy = (fNy_start, fNy_start + fN2 - 1)                    # y direction 
fNz = (1, fN3)                   # z direction

# VW region of the RS

fN2_VW = BP5_coeff.l /(BP5_coeff.Δz / 1000) + 1
fN2_VW = round(Int, fN2_VW, RoundUp)
fN2_VW_favorable = BP5_coeff.w / (BP5_coeff.Δz / 1000) + 1
fN2_VW_favorable = round(Int, fN2_VW_favorable, RoundUp)
fN3_VW = BP5_coeff.H / (BP5_coeff.Δz / 1000) + 1
fN3_VW = round(Int, fN3_VW, RoundUp)

# fNy_VW and fNz_VW represents the indices of the fault region
fNy_VW_start = (Ly - BP5_coeff.l) / (2 * BP5_coeff.Δz / 1000) + 1
fNy_VW_start = round(Int, fNy_VW_start, RoundUp)
fNy_VW = (fNy_VW_start, fNy_VW_start + fN2_VW - 1)
fNy_VW_favorable = (fNy_VW_start, fNy_VW_start + fN2_VW_favorable - 1)
fNz_VW_start = round(Int, (BP5_coeff.hs + BP5_coeff.ht) / (BP5_coeff.Δz / 1000), RoundUp) + 1
fNz_VW = (fNz_VW_start, fNz_VW_start + fN3_VW - 1)


# VW-VS transition region
fN2_VW_VS = (BP5_coeff.l + 2 * BP5_coeff.ht) / (BP5_coeff.Δz / 1000) + 1
fN2_VW_VS = round(Int, fN2_VW_VS, RoundUp)

fN3_VW_VS = (BP5_coeff.H + BP5_coeff.ht * 2) / (BP5_coeff.Δz / 1000) + 1
fN3_VW_VS = round(Int, fN3_VW_VS, RoundUp)

fNy_VW_VS_start = (Ly- BP5_coeff.l - BP5_coeff.ht) / (2 * BP5_coeff.Δz / 1000) + 1
fNy_VW_VS_start = round(Int, fNy_VW_VS_start, RoundUp)
fNy_VW_VS = (fNy_VW_VS_start, fNy_VW_VS_start + fN2_VW_VS - 1)
fNz_VW_VS_start = round(Int, BP5_coeff.hs / (BP5_coeff.Δz / 1000), RoundUp) + 1
fNz_VW_VS = (fNz_VW_VS_start, fNz_VW_VS_start + fN3_VW_VS - 1)

# fNy_VW_VS and fNz_VW represents the indices of the fault region

# Assembling matrices for 3D SBP-SAT
SBPp = 2                # SBPp order
(M, RHS, H_tilde, HI_tilde, analy_sol, source, traction_operators, 
    u_filters, Face_operators, sigmas, updators) = Assembling_3D_matrices(N_x, N_y, N_z;SBPp=SBPp,Lx=Lx,Ly=Ly,Lz=Lz);
M_GPU = CUDA.CUSPARSE.CuSparseMatrixCSR(M);

# set RHS to be zero at the beginning
RHS .= 0

sigma_11 = sigmas[1]
sigma_21 = sigmas[2]
sigma_31 = sigmas[3]
 

function create_ψVδ()
    dψV = zeros(fN2 * fN3 + 2 * (Ny) * (Nz))
    ψδ = zeros(fN2 * fN3 + 2 * (Ny) * (Nz))
    return dψV, ψδ
end

function create_view(dψV, ψδ)
    index_1 = fN2 * fN3

    dψ = @view dψV[1:index_1]
    V = @view dψV[index_1 + 1:end]
    ψ = @view ψδ[1:index_1]
    δ = @view ψδ[index_1 + 1:end]

    return dψ, V, ψ, δ
end

# getting rate-and-state 
function get_RS_indices(Nx, Ny, Nz, fNy, fNz)
    # fNy is a tuple of the start index and end index of fault in y(x2) direction 
    # fNz is a tuple of the start index and end index of fault in z(x3) direction 
    x_idx = spzeros(Nx)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    x_idx[1] = 1 # the fault is on the first face
    y_idx[fNy[1]:fNy[2]] .= 1
    z_idx[fNz[1]:fNz[2]] .= 1
    return kron(z_idx, y_idx, x_idx) # x->y->z is the index order, hence the kron order is z_idx<-y_idx<-x_idx
end

function get_RS_indices_2D(Ny, Nz, fNy, fNz)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    y_idx[fNy[1]:fNy[2]] .= 1
    z_idx[fNz[1]:fNz[2]] .= 1
    return kron(z_idx, y_idx)
end

function get_uniform_indices(Nx, Ny, Nz, fNy_VW, fNz_VW)
    x_idx = spzeros(Nx)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    x_idx[1] = 1 # the fault is on the first face
    y_idx[fNy_VW[1]:fNy_VW[2]] .= 1
    z_idx[fNz_VW[1]:fNz_VW[2]] .= 1
    return kron(z_idx, y_idx, x_idx)
end

function get_uniform_indices_2D(Ny, Nz, fNy_VW, fNz_VW)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    y_idx[fNy_VW[1]:fNy_VW[2]] .= 1
    z_idx[fNz_VW[1]:fNz_VW[2]] .= 1
    return kron(z_idx, y_idx)
end


function get_transition_indices(Nx, Ny, Nz, fNy_VW, fNz_VW, fNy_VW_VS, fNz_VW_VS)
    x_idx = spzeros(Nx)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    x_idx[1] = 1 # the fault is on the first face
    y_idx[fNy_VW_VS[1]:fNy_VW_VS[2]] .= 1
    z_idx[fNz_VW_VS[1]:fNz_VW_VS[2]] .= 1
    # return kron(z_idx, y_idx, x_idx)
    return kron(z_idx, y_idx, x_idx) - get_uniform_indices(Nx, Ny, Nz, fNy_VW, fNz_VW)
end

function get_favorable_indices(Nx, Ny, Nz, fNy_VW_favorable, fNz_VW)
    x_idx = spzeros(Nx)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    x_idx[1] = 1 # the fault is on the first face
    y_idx[fNy_VW_favorable[1]:fNy_VW_favorable[2]] .= 1
    z_idx[fNz_VW[1]:fNz_VW[2]] .= 1
    return kron(z_idx, y_idx, x_idx)
end

function get_favorable_indices_2D(Ny, Nz, fNy_VW_favorable, fNz_VW)
    y_idx = spzeros(Ny)
    z_idx = spzeros(Nz)
    y_idx[fNy_VW_favorable[1]:fNy_VW_favorable[2]] .= 1
    z_idx[fNz_VW[1]:fNz_VW[2]] .= 1
    return kron(z_idx, y_idx)
end

function get_favorable_indices_RS(fNy, fNz, fNy_VW_favorable, fNz_VW)
    y_idx = spzeros(fNy[2] - fNy[1] + 1)
    z_idx = spzeros(fNz[2] - fNz[1] + 1)
    y_idx[fNy_VW_favorable[1] - fNy[1] + 1: fNy_VW_favorable[2] - fNy[1] + 1] .= 1
    z_idx[fNz_VW[1]: fNz_VW[2]] .= 1
    return kron(z_idx, y_idx)
end



RS_filter = get_RS_indices(Nx, Ny, Nz, fNy, fNz)
RS_filter_nzind = RS_filter.nzind
RS_filter_2D = get_RS_indices_2D(Ny, Nz, fNy, fNz)
RS_filter_2D_nzind = RS_filter_2D.nzind
VW_filter = get_uniform_indices(Nx, Ny, Nz, fNy_VW, fNz_VW)
VW_filter_2D = get_uniform_indices_2D(Ny, Nz, fNy_VW, fNz_VW)
VW_filter_2D_nzind = VW_filter_2D.nzind
VW_favorable_filter = get_favorable_indices(Nx, Ny, Nz, fNy_VW_favorable, fNz_VW)
VW_favorable_filter_nzind = VW_favorable_filter.nzind
VW_favorable_filter_2D = get_favorable_indices_2D(Ny, Nz, fNy_VW_favorable, fNz_VW)
VW_favorable_filter_2D_nzind = VW_favorable_filter_2D.nzind
VW_favorable_filter_RS = get_favorable_indices_RS(fNy, fNz, fNy_VW_favorable, fNz_VW)
VW_favorable_filter_RS_nzind = VW_favorable_filter_RS.nzind
VW_VS_transition_filter = get_transition_indices(Nx, Ny, Nz, fNy_VW, fNz_VW, fNy_VW_VS, fNz_VW_VS)
VS_filter = RS_filter - VW_filter - VW_VS_transition_filter


@assert length(VW_filter.nzind) + length(VW_VS_transition_filter.nzind) + length(VS_filter.nzind) == fN2 * fN3




# Time series
# On-Fault series
fltst = [
    [0, -36, 0],
    [0, -16, 0],
    [0, 0, 0],
    [0, 16, 0],
    [0, 36, 0],
    [0, -24, 10],
    [0, -16, 10],
    [0, 0, 10],
    [0, 16, 10],
    [0, 0, 22]
]

function find_flt_indices(indices, lf, fN2)
    x2 = indices[2]
    x3 = indices[3]
    j = Int(round((x2 - (-lf / 2)) / (BP5_coeff.Δz / 1000))) + 1 # starting with 1
    k = Int(round((x3 - 0) / (BP5_coeff.Δz / 1000))) # starting with 0 (multiplied by fN2) no +1
    return j + k * fN2
end

# time_string = Dates.format(now(),"yyyymmddHHMM")
# path="./output/$time_string/"
path="./output"
station_indices = find_flt_indices.(fltst,BP5_coeff.lf, fN2)
station_strings = ["-36dp+00", "-16dp+00", "+00dp+00", "+16dp+00", "+36dp+00",
                    "-24dp+10", "-16dp+10", "+00dp+10","+16dp+10",
                    "+00dp+22"]

nothing # avoid printing out results 