# 3D array indexing in Julia
#
#    3
#    |
#    |
#    |___________  2
#   /
#  /
# 1
#
# Face 1: End, Face 2: Front
# Face 3: Left, Face 4: Right
# Face 5: Bottom, Face 6: Top 
#
# Example: julia> test_3D
# 3×3×3 Array{Float64, 3}:
# [:, :, 1] =
# 1.99217   -0.419911   1.27037
# -0.616689  -0.615541  -0.452607
# -0.144164  -1.93718    0.418154

# [:, :, 2] =
# 3.03034   0.37151    0.964773
# 1.40377  -0.338357  -0.2258
# 1.21607   0.57358   -0.364077

# [:, :, 3] =
# -0.417334   -0.675461   -1.81201
# -0.0343112   0.0751002   0.562684
# -0.326663    0.949942    0.587156
#
#

#   julia> test_3D[:]
#   27-element Vector{Float64}:
#   1.9921749921668979
#  -0.6166891690372321
#  -0.1441644800574074
#  -0.4199112178717434
#  -0.6155411265727337
#  -1.9371793724089905
#   1.270365705663953
#  -0.45260724271602304
#   0.41815434826401493
#   3.0303354534340596
#   1.4037652552774593
#   1.2160687591726536
#   0.3715100914804111
#  -0.3383569289529389
#   0.5735799032430334
#   0.9647725600192497
#  -0.22579997325200207
#  -0.3640769675163103
#  -0.4173335304749411
#  -0.03431124808637245
#  -0.326663027771571
#  -0.6754607583657409
#   0.07510016190239746
#   0.9499416235735255
#  -1.8120098781175735
#   0.5626840107886724
#   0.5871558963442842

using LinearAlgebra
using SparseArrays

function e(i,n)
    # A = Matrix{Float64}(I,n,n)
    A = spzeros(n)
    A[i] = 1
    return A
end

function e_t(i,n)
    A = spzeros(1,n)
    A[i] = 1
    return A
end

function eyes(n)
    # return Matrix{Float64}(I,n,n)
    return sparse(I,n,n)
end


function Diag(A)
    # Self defined function that is similar to Matlab Diag
    return Diagonal(A[:])
end


function get_bottom_face(Nx,Ny,Nz)
    # mat = kron(e(1,Nz)',sparse(eyes(Nx)),sparse(eyes(Ny)))
    # mat = kron(sparse(e(1,Nz)'),eyes(Nx),eyes(Ny))
    mat = kron(e_t(1,Nz),eyes(Nx),eyes(Ny))
    return mat
end

function get_bottom_face(A)
    Nx,Ny,Nz = size(A)
    return get_bottom_face(Nx,Ny,Nz)*A[:]
end

function get_top_face(Nx,Ny,Nz)
    # mat = kron(e(Nz,Nz)',sparse(eyes(Nx)),sparse(eyes(Ny)))
    mat = kron(e_t(Nz,Nz),eyes(Nx),eyes(Ny))
    return mat
end

function get_top_face(A)
    Nx,Ny,Nz = size(A)
    return get_top_face(Nx,Ny,Nz)*A[:]
end

function get_left_face(Nx,Ny,Nz)
    # mat = kron(sparse(eyes(Nz)),e(1,Ny)',sparse(eyes(Nx)))
    mat = kron(eyes(Nz),e_t(1,Ny),eyes(Nx))
    return mat
end

function get_left_face(A)
    Nx,Ny,Nz = size(A)
    return get_left_face(Nx,Ny,Nz)*A[:]
end

function get_right_face(Nx,Ny,Nz)
    # mat = kron(sparse(eyes(Nz)),e(Ny,Ny)',sparse(eyes(Nx)))
    mat = kron(eyes(Nz),e_t(Ny,Ny),eyes(Nx))
    return mat
end

function get_right_face(A)
    Nx,Ny,Nz = size(A)
    return get_right_face(Nx,Ny,Nz)*A[:]
end

function get_front_face(Nx,Ny,Nz)
    # mat = kron(sparse(eyes(Ny)),sparse(eyes(Nz)),sparse(e(1,Nx)'))
    mat = kron(eyes(Ny),eyes(Nz),e_t(Nx,Nx))
    return mat
end

function get_front_face(A)
    Nx,Ny,Nz = size(A)
    return get_front_face(Nx,Ny,Nz)*A[:]
end

function get_end_face(Nx,Ny,Nz)
    # mat = kron(sparse(eyes(Ny)),sparse(eyes(Nz)),e(Nx,Nx)')
    mat = kron(eyes(Ny),eyes(Nz),e_t(1,Nx))
    return mat
end

function get_end_face(A)
    Nx,Ny,Nz = size(A)
    return get_end_face(Nx,Ny,Nz)*A[:]
end




test_matrix = randn(2,3,4)
get_bottom_face(test_matrix)
get_top_face(test_matrix)
get_front_face(test_matrix)
get_end_face(test_matrix)
get_left_face(test_matrix)
get_right_face(test_matrix)
