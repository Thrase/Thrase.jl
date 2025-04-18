# This file specifies the computational domain
using Dates
include("helper.jl")
const year_seconds = 31556926
sim_years = 10

# calling parameterized constructor to set values for BP7
# if BP7_coeff is not defined here, the BP7-QD.jl will 
# call default constructor to construct BP7_coeff using
# values in helper.jl (line 51)
BP7_coeff = coefficients(
    2670,                   # ρ
    3.464,                  # cs
    32.04,                   # μ
    0.25,                   # ν
    0.004,                  # a0
    0.016,                   # amax
    0.01,                   # b0            value for b in this problem
    25,                     # σn
    0.00053,                   # DRS
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
    1.75,                      # Δτ0
    150,                    # Rnuc
    1,                      # T 
    -50,                    # y2
    -50                     # y3
)



include("Assembling_3D_matrices.jl")

# The entire domain is 1 km by 1 km by 1 km
Lx = Ly = Lz = 1000
N_x = N_y = N_z = Int(Lx / (BP7_coeff.Δz))
Nx = N_x + 1
Ny = N_y + 1
Nz = N_z + 1
# n_levels = Int(log(2, N_x))

u1_filter_matrix = get_u1(Nx, Ny, Nz)
u2_filter_matrix = get_u2(Nx, Ny, Nz)
u3_filter_matrix = get_u3(Nx, Ny, Nz)


fN2 = 2 * BP7_coeff.lf / (BP7_coeff.Δz ) + 1
fN2 = round(Int, fN2, RoundUp)
fN3 = 2 * BP7_coeff.Wf / (BP7_coeff.Δz ) + 1
fN3 = round(Int, fN3, RoundUp)

# fNy and fNz represents the indices of the fault region
fNy_start = (Ly - 2 * BP7_coeff.lf) / (2 * BP7_coeff.Δz) + 1
fNy_start = round(Int, fNy_start, RoundUp)
fNy = (fNy_start, fNy_start + fN2 - 1)                    # y direction 

# fNz = (1, fN3)                   # z direction
fNz_start = (Lz - 2 * BP7_coeff.Wf) / (2 * BP7_coeff.Δz) + 1
fNz_start = round(Int, fNz_start, RoundUp)
fNz = (fNz_start, fNz_start + fN3 - 1)

# VW region of the RS
function create_sparse_matrix(size, condition)
    """
    Creates a sparse matrix where elements satisfy the given condition.

    Args:
        size: A tuple representing the dimensions of the matrix (rows, columns).
        condition: A function that takes two arguments (row index, column index) 
                   and returns true if the element should be non-zero, false otherwise.

    Returns:
        A sparse matrix with non-zero elements determined by the condition.
    """

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for i in 1:size[1]
        for j in 1:size[2]
            if condition(i, j)
                push!(rows, i)
                push!(cols, j)
                # You can assign any value here, e.g., 1.0
                push!(vals, 1.0) 
            end
        end
    end

    sparse(rows, cols, vals, size[1], size[2])
end

function within_VW(i, j, BP7_coeff, Ly, Lz)
    y_coord = (i - 1) * BP7_coeff.Δz
    z_coord = (j - 1) * BP7_coeff.Δz
    y2 = Ly/2
    y3 = Lz/2
    return (y_coord - y2)^2 + (z_coord - y3)^2 <= BP7_coeff.RVW^2
end




function get_VW_indices_2D(Ny, Nz, BP7_coeff)
    matrix_size = (Ny, Nz)
    condition(i,j) = within_VW(i, j, BP7_coeff, Ly, Lz)
    sparse_matrix_2D = create_sparse_matrix(matrix_size, condition)
    return sparse_matrix_2D[:]
end

function get_VW_indices(Nx, Ny, Nz, BP7_coeff)
    idx = spzeros(Nx)
    idx[1] = 1
    matrix_size = (Ny, Nz)
    condition(i,j) = within_VW(i, j, BP7_coeff, Ly, Lz)
    sparse_matrix_2D = create_sparse_matrix(matrix_size, condition)
    sparse_matrix_3D = kron(sparse_matrix_2D,idx)[:]
    return sparse_matrix_3D
end


function within_A(i, j, BP7_coeff, Ly, Lz)
    y_coord = (i - 1) * BP7_coeff.Δz
    z_coord = (j - 1) * BP7_coeff.Δz
    y2 = Ly/2
    y3 = Lz/2
    return (y_coord - y2)^2 + (z_coord - y3)^2 <= (1.5 * BP7_coeff.RVW)^2 
end

function get_A_indices_2D(Ny, Nz, BP7_coeff)
    matrix_size = (Ny, Nz)
    condition(i,j) = within_A(i, j, BP7_coeff, Ly, Lz)
    sparse_matrix_2D = create_sparse_matrix(matrix_size, condition)
    return sparse_matrix_2D[:]
end


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



RS_filter = get_RS_indices(Nx, Ny, Nz, fNy, fNz)
RS_filter_nzind = RS_filter.nzind
RS_filter_2D = get_RS_indices_2D(Ny, Nz, fNy, fNz)
RS_filter_2D_nzind = RS_filter_2D.nzind


VW_filter = get_VW_indices(Nx, Ny, Nz, BP7_coeff)
VW_filter_nzind = VW_filter.nzind
VW_filter_2D = get_VW_indices_2D(Ny, Nz, BP7_coeff)
VW_filter_2D_nzind = VW_filter_2D.nzind

VS_filter_2D = RS_filter_2D - VW_filter_2D
VS_filter = RS_filter - VW_filter

# Get patched VW for moment output
A_filter_2D = get_A_indices_2D(Ny, Nz, BP7_coeff)
A_filter_2D_nzind = A_filter_2D.nzind

function find_indices_in_array(A, B)
    """
    Finds the indices of all elements in array A within array B.

    Args:
        A: The array of elements to find.
        B: The array to search within.

    Returns:
        An array of the same length as A, where each element is 
        the index of the corresponding element in B, or 0 if 
        the element is not found in B.
    """

    indices = zeros(Int, length(A))
    for (i, a) in enumerate(A)
        indices[i] = findfirst(==(a), B) 
    end
    return indices
end

VW_on_RS_filter = find_indices_in_array(VW_filter_2D_nzind, RS_filter_2D_nzind)
A_on_RS_filter = find_indices_in_array(A_filter_2D_nzind, RS_filter_2D_nzind)

@assert length(VW_filter.nzind) + length(VS_filter.nzind) == fN2 * fN3
sparse(reshape(VS_filter_2D, Ny, Nz))


# Time series
# On-Fault series
fltst = [
    [0, 0, 0],
    [0, -100, 0],
    [0, 0, 100],
    [0, 100, 0],
    [0, 0, -100],
    [0, -100, -100],
    [0, -100, 100],
    [0, 100, -100],
    [0, 100, 100],
    [0, -300, 0],
    [0, 0, 300],
    [0, 300, 0],
    [0, 0, -300]
]

function find_flt_indices(indices, lf, Wf, fN2)
    x2 = indices[2]
    x3 = indices[3]
    j = Int(round((x2 - (-lf)) / (BP7_coeff.Δz ))) + 1 # starting with 1
    k = Int(round((x3 - (-Wf)) / (BP7_coeff.Δz ))) # starting with 0 (multiplied by fN2) no +1
    return j + k * fN2
end

# time_string = Dates.format(now(),"yyyymmddHHMM")
# path="./output/$time_string/"
path="./output"
station_indices = find_flt_indices.(fltst,BP7_coeff.lf, BP7_coeff.Wf, fN2)
# station_strings = ["-36dp+00", "-16dp+00", "+00dp+00", "+16dp+00", "+36dp+00",
#                     "-24dp+10", "-16dp+10", "+00dp+10","+16dp+10",
#                     "+00dp+22"]
station_strings = ["+000dp+000",
                "-100dp+000",
                "+000dp+100",
                "+100dp+000",
                "+000dp-100",
                "-100dp-100",
                "-100dp+100",
                "+100dp-100",
                "+100dp+100",
                "-300dp+000",
                "+000dp+300",
                "+300dp+000",
                "+000dp-300"]

nothing # avoid printing out results 