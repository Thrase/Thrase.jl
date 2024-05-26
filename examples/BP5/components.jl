using LinearAlgebra
using SparseArrays

include("3D_face.jl")

function get_u1(Nx,Ny,Nz)
    output = eyes(Nx*Ny*Nz)
    e1 = spzeros(1,3)
    e1[1] = 1
    return kron(output,e1)   
end


function get_u2(Nx,Ny,Nz)
    output = eyes(Nx*Ny*Nz)
    e1 = spzeros(1,3)
    e1[2] = 1
    return kron(output,e1)   
end

function get_u3(Nx,Ny,Nz)
    output = eyes(Nx*Ny*Nz)
    e1 = spzeros(1,3)
    e1[3] = 1
    return kron(output,e1)   
end




