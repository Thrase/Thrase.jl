# Stripped down version of Thrase code (for training purposes)
# Does not integrate a coordinate transformation. 
# On a structured mesh of quads.
# Solves the 2D antiplane problem 0 = μ(u_xx + u_zz) + f(x, z)
# with Dirichlet BC on left/right and traction-free on top/bottom
# Does convergence tests with MMS (method of manufactured solutions)

using Thrase
using LinearAlgebra

include("ops_stripped.jl")
include("../utils_2D.jl")


# Define shear modulus μ
μ = 2 

# Define domain
xc = (0, 2)  # need to be integers for manufactured solution to workout correctly
zc = (0, 2)  # need to be integers for manufactured solution to workout correctly

# manufactured solution
function exact(μ, x, z)

    
    uexact = sin.(pi * x) * cos.(pi * z')
    uexact_x = pi * cos.(pi * x) * cos.(pi * z')
    uexact_xx = -pi^2 * sin.(pi * x) * cos.(pi * z')
    uexact_z = -pi * sin.(pi * x) * sin.(pi * z')
    uexact_zz = -pi^2 * sin.(pi * x) * cos.(pi * z')
    f = -μ .* (uexact_xx .+ uexact_zz)  # source data

    # boundary data
    gL = sin(pi  * (xc[1])) * cos.(pi * z)
    gR = sin(pi * (xc[end])) * cos.(pi * z)
    gT = -pi * μ * sin.(pi * x) * sin(pi * zc[end])
    gB = pi * μ * sin.(pi * x) * sin(pi * zc[1])

    return uexact[:], f[:], gL, gR, gT, gB

end


function antiplane_solve(Nx, Nz, μ, SBPp)

    # Define 2D domain (x, z) in [0, Lx] x [0, Lz], constant grid spacing
    x = Array(LinRange(xc[1], xc[2], Nx+1))  
    z = Array(LinRange(zc[1], zc[2], Nz+1))


    # Domain lengths 
    Lx = xc[2]-xc[1]
    Lz = zc[2]-zc[1]

    # Get data from manufactured solution
    uexact, f, gL, gR, gT, gB = exact(μ, x, z)
    
    # create operators
    (A, B, H̃, T, e) = get_operators(SBPp, Nx, Nz, μ; xc = xc, zc = zc)

    A = lu(A)  # LU matrix factorization

    # initialize vector b that stores boundary data (linear system will be Au = b, where b = B*g)
    b = zeros((Nx+1) * (Nz+1))
    
    # fill in initial boundary data into b (note that manuf. solution is such that gT = gB = 0)
 
    bdry_vec_strip!(b, B, x, z, gL, gR, gT, Lx, Lz)
    
    u = A \ (b-f) # solve linear system with a backsolve to obtain initial displacement.

    return (u, uexact, H̃, x, z)
end

# Convergence tests:

# Define spatial order of accuracy
SBPp = 2

# Define number of grid grid_points
Nx = [2^5, 2^6, 2^7, 2^8, 2^9]
Nz = Nx

err = zeros(length(Nx))
rates = zeros(length(Nx)-1)

for n = 1:length(Nx)

    (u, uexact, H̃, x, z) = antiplane_solve(Nx[n], Nz[n], μ, SBPp)
    diff = u .- uexact
    err[n] = sqrt(diff' * H̃ * diff) ./ sqrt(uexact' * H̃ * uexact)   # relative error
    
end

for n = 2:length(Nx)
    rates[n-1] = log2(err[n-1]/err[n])
end

@show err
@show rates

# p = (2, 4, 6, 8) should yield rates of (2, 3, 4, 5), respectively.