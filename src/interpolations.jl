using SparseArrays
using LinearAlgebra



"""
    restriction_matrix_v0(nxf,nyf,nxc,nyc)
    Generating standard restriction matrix from fine grid (nxf,nyf,nzf) to coarse grid (nxc, nyc, nzc)
    The matrix size generated is ((nxc+1) * (nyc+1) * (nzc+1) by (nxf+1) * (nyf+1) * (nzf+1))
    
    # Examples
    ```julia
    julia> restriction_matrix_v0(4,4,4,2,2,2)
    27×125 SparseMatrixCSC{Float64, Int64} with 343 stored entries:
    ⎡⠙⠲⢌⡓⠦⣀⠀⢀⡀⠀⠀⠀⠈⠓⠦⣙⠲⢄⡀⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
    ⎢⠀⠀⠀⠉⠓⠈⠙⠢⣍⡓⢤⣀⠀⠀⠀⠈⠙⠂⠉⠓⢬⣙⠢⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠈⠑⠦⣄⠲⢤⡀⠀⠀⠀⠈⠁⠀⠉⠲⢤⡐⠦⣄⠀⠀⠀⠀⠀⠀⠀⠐⠦⣄⠲⢤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠲⠌⠓⠦⣙⠲⢄⡀⠀⠀⠀⠈⠓⠦⠙⠲⢌⡓⠦⣀⠀⠀⠀⠀⠙⠲⠌⠓⠦⣙⠲⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠂⠉⠓⠀⠀⠀⠀⠀⠀⠀⠀⠉⠓⠈⠙⠢⣄⡀⢤⣀⠀⠀⠀⠈⠙⠂⠉⠓⢤⣀⠠⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠲⢬⡑⠦⣄⠲⢤⡀⠀⠀⠀⠈⠑⠦⣍⠲⢤⡐⠦⣄⠀⠀⠀⎥
    ⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠲⠌⠓⠦⠀⠀⠀⠀⠀⠀⠀⠈⠓⠦⠙⠲⠄⎦
    ```
"""
function restriction_matrix_v0(nxf,nyf,nzf,nxc,nyc,nzc)
    # x
    restriction_matrix_x = spzeros(nxc+1, nxf+1)
    restriction_matrix_x[1,1] = 1/2
    restriction_matrix_x[1,2] = 1/2
    for i in 2:nxc
        restriction_matrix_x[i,2*i-1] = 1/2
        restriction_matrix_x[i,2*i-2] = 1/4
        restriction_matrix_x[i,2*i] = 1/4
    end
    restriction_matrix_x[end,end] = 1/2
    restriction_matrix_x[end,end-1] = 1/2

    # y
    restriction_matrix_y = spzeros(nyc+1, nyf+1)
    restriction_matrix_y[1,1] = 1/2
    restriction_matrix_y[1,2] = 1/2
    for i in 2:nyc
        restriction_matrix_y[i,2*i-1] = 1/2
        restriction_matrix_y[i,2*i-2] = 1/4
        restriction_matrix_y[i,2*i] = 1/4
    end
    restriction_matrix_y[end,end] = 1/2
    restriction_matrix_y[end,end-1] = 1/2

    # z
    restriction_matrix_z = spzeros(nzc+1, nzf+1)
    restriction_matrix_z[1,1] = 1/2
    restriction_matrix_z[1,2] = 1/2
    for i in 2:nzc
        restriction_matrix_z[i,2*i-1] = 1/2
        restriction_matrix_z[i,2*i-2] = 1/4
        restriction_matrix_z[i,2*i] = 1/4
    end
    restriction_matrix_z[end,end] = 1/2
    restriction_matrix_z[end,end-1] = 1/2

    restriction_matrix_ = kron(restriction_matrix_x,restriction_matrix_y, restriction_matrix_z)
    return restriction_matrix_
end



"""
    prolongation_matrix_v0(nxf,nyf,nxc,nyc)
    Generating standard prolongation matrix from coarse grid (nxc, nyc) to fine grid (nxf,nyf)
    The matrix size generated is (nxf+1) * (nyf+1)) by ((nxc+1) * (nyc+1) 

    # Examples
    ```julia
    julia> prolongation_matrix_v0(4,4,4,2,2,2)
    125×27 SparseMatrixCSC{Float64, Int64} with 343 stored entries:
    ⎡⣳⡀⠀⠀⠀⠀⠀⠀⎤
    ⎢⠘⣗⠀⠀⠀⠀⠀⠀⎥
    ⎢⠀⢾⡆⠀⠀⠀⠀⠀⎥
    ⎢⣄⠈⢧⠀⠀⠀⠀⠀⎥
    ⎢⢾⣆⢸⣧⡀⠀⠀⠀⎥
    ⎢⠀⣧⡀⣳⣄⠀⠀⠀⎥
    ⎢⠀⠳⡇⠘⢾⠀⠀⠀⎥
    ⎢⠀⠀⢱⡄⠈⠀⠀⠀⎥
    ⎢⠀⠀⠘⣟⠆⠀⠀⠀⎥
    ⎢⠀⠀⠀⢮⣧⠀⠀⠀⎥
    ⎢⠀⠀⢀⠀⠻⡀⠀⠀⎥
    ⎢⠀⠀⢸⢦⠀⡷⡄⠀⎥
    ⎢⠀⠀⠈⠻⡅⠙⡇⠀⎥
    ⎢⠀⠀⠀⠙⣿⠈⢿⠆⎥
    ⎢⠀⠀⠀⠀⠘⣆⠈⠃⎥
    ⎢⠀⠀⠀⠀⠀⢿⠆⠀⎥
    ⎢⠀⠀⠀⠀⠀⠨⣧⠀⎥
    ⎣⠀⠀⠀⠀⠀⠀⠹⡅⎦
    ```
"""
function prolongation_matrix_v0(nxf,nyf,nzf,nxc,nyc,nzc) 
    # SBP preserving prolongation matrix
    prolongation_matrix_x = spzeros(nxf+1,nxc+1)
    for i in 1:nxc
        prolongation_matrix_x[2*i-1,i] = 1
        prolongation_matrix_x[2*i,i] = 0.5
        prolongation_matrix_x[2*i,i+1] = 0.5
    end
    prolongation_matrix_x[end,end] = 1

    prolongation_matrix_y = spzeros(nyf+1,nyc+1)
    for i in 1:nyc
        prolongation_matrix_y[2*i-1,i] = 1
        prolongation_matrix_y[2*i,i] = 0.5
        prolongation_matrix_y[2*i,i+1] = 0.5
    end
    prolongation_matrix_y[end,end] = 1

    prolongation_matrix_z = spzeros(nzf+1,nzc+1)
    for i in 1:nyc
        prolongation_matrix_z[2*i-1,i] = 1
        prolongation_matrix_z[2*i,i] = 0.5
        prolongation_matrix_z[2*i,i+1] = 0.5
    end
    prolongation_matrix_z[end,end] = 1

    prolongation_matrix_ = kron(prolongation_matrix_x,prolongation_matrix_y,prolongation_matrix_z)
    return prolongation_matrix_
end

