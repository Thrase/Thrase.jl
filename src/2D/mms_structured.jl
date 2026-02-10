# MMS integrates a coordinate transformation and variable μ. 
# On a structured mesh of quads.
# Solves the 2D antiplane problem 0 = d_dx(mu(x,z)du_dx) + d_dz(mu(x,z)du_dz) + f(x, z)
# with Dirichlet BC on left/right and traction-free on top/bottom
# Does convergence tests with MMS (method of manufactured solutions)

using Thrase
using LinearAlgebra

include("ops_BP1-QD_structured.jl")
include("odefun_BP1-QD_structured.jl")

# manufactured solution
function exact_mu(x, z)
    # Define shear modulus μ
    μ = 2 .+ 1 .* sin.(x .+ z)
    return μ
end

function exact_mu_x(x, z)
    # Define shear modulus derivative μ_x
    μ_x = 1 .* cos.(x .+ z)
    return μ_x
end

function exact_mu_z(x, z)
    # Define shear modulus derivative μ_z
    μ_z = 1 .* cos.(x .+ z)
    return μ_z
end

function exact(x, z, facecoord, nx, ny)
  
  
    uexact = sin.(pi * x) .* cos.(pi * z)
    uexact_x = pi * cos.(pi * x) .* cos.(pi * z)
    uexact_xx = -pi^2 * sin.(pi * x) .* cos.(pi * z)
    uexact_z = -pi * sin.(pi * x) .* sin.(pi * z)
    uexact_zz = -pi^2 * sin.(pi * x) .* cos.(pi * z)
    f = -exact_mu_x(x, z) .* uexact_x .- exact_mu_z(x, z) .* uexact_z .- exact_mu(x, z) .* (uexact_xx .+ uexact_zz)  # source data


    # boundary data
    gL = sin.(pi  .* facecoord[1][1]) .* cos.(pi * facecoord[2][1])
    gR = sin.(pi .* facecoord[1][2]) .* cos.(pi * facecoord[2][2])

    # gB = sin.(pi  .* facecoord[1][3]) .* cos.(pi * facecoord[2][3])
    # gT = sin.(pi .* facecoord[1][4]) .* cos.(pi * facecoord[2][4])
    mu_uexact_x_3 = exact_mu(facecoord[1][3], facecoord[2][3]) .* pi .* cos.(pi * facecoord[1][3]) .* cos.(pi * facecoord[2][3])
    mu_uexact_z_3 = exact_mu(facecoord[1][3], facecoord[2][3]) .* (-pi * sin.(pi * facecoord[1][3]) .* sin.(pi * facecoord[2][3]))
    gB = nx[3] .* mu_uexact_x_3 .+ ny[3] .* mu_uexact_z_3
    
    mu_uexact_x_4 = exact_mu(facecoord[1][4], facecoord[2][4]) .* pi .* cos.(pi * facecoord[1][4]) .* cos.(pi * facecoord[2][4])
    mu_uexact_z_4 = exact_mu(facecoord[1][4], facecoord[2][4]) .* (-pi * sin.(pi * facecoord[1][4]) .* sin.(pi * facecoord[2][4]))
    gT = nx[4] .* mu_uexact_x_4 .+ ny[4] .* mu_uexact_z_4

    return uexact[:], f[:], gL, gR, gB, gT

end


function antiplane_solve(Nr, Ns, SBPp)

    ################################## COORDINATE TRANSFORM ###################################
    # DOMAIN: ordered pairs that define the physical domain
    (x1,x2,x3,x4) = (-.2,1,0,2)  
    (y1,y2,y3,y4) = (-.3,0,1,2)

    # Details of transfinite interpolation: straight edges:
    # x-parameterization on each face
    ex = [(α) -> x1 * (1 .- α) / 2 + x3 * (1 .+ α) / 2,
    (α) -> x2 * (1 .- α) / 2 + x4 * (1 .+ α) / 2,
    (α) -> x1 * (1 .- α) / 2 + x2 * (1 .+ α) / 2,
    (α) -> x3 * (1 .- α) / 2 + x4 * (1 .+ α) / 2]

    # derivative of x-parameterization
    exα = [(α) -> -x1 / 2 + x3 / 2,
    (α) -> -x2 / 2 + x4 / 2,
    (α) -> -x1 / 2 + x2 / 2,
    (α) -> -x3 / 2 + x4 / 2]
    
    # y-parameterization on each face
    ey = [(α) -> y1 * (1 .- α) / 2 + y3 * (1 .+ α) / 2,
    (α) -> y2 * (1 .- α) / 2 + y4 * (1 .+ α) / 2,
    (α) -> y1 * (1 .- α) / 2 + y2 * (1 .+ α) / 2,
    (α) -> y3 * (1 .- α) / 2 + y4 * (1 .+ α) / 2]

    # derivative of y-parameterization
    eyα = [(α) -> -y1 / 2 + y3 / 2,
    (α) -> -y2 / 2 + y4 / 2,
    (α) -> -y1 / 2 + y2 / 2,
    (α) -> -y3 / 2 + y4 / 2]

    # Coordinate transformation with transfinite interpolation
    xt(r,s) = transfinite_blend(ex[1], ex[2], ex[3], ex[4], exα[1], exα[2], exα[3], exα[4], r, s)
    zt(r,s) = transfinite_blend(ey[1], ey[2], ey[3], ey[4], eyα[1], eyα[2], eyα[3], eyα[4], r, s)

    # Uncomment if coordinate transform user-specified, for example:
    #xt(r,s) = (2 * r, 2 * ones(size(r)), 0 * ones(size(s)))
    #zt(r,s) = (2 * s, 0 * ones(size(r)), 2 * ones(size(s)))

    metrics = create_metrics(Nr, Ns, exact_mu, xt, zt) # create coordinate transform

    ###################################################################### 
    # create finite difference operators on computational domain:
    (A, B, H̃, T, e, J, coord, facecoord, sJ, nx, ny) = get_operators(SBPp, Nr, Ns, metrics) 
    A = lu(A)  # LU matrix factorization

    # Get data from manufactured solution
    uexact, f, gL, gR, gB, gT = exact(coord[1], coord[2], facecoord, nx, ny)

    # initialize vector b that stores boundary data (linear system will be Au = b, where b = B*g)
    b = zeros((Nr+1) * (Ns+1))
    
    # fill in initial boundary data into b
    bdry_vec!(b, B, gL, gR, gB, gT, sJ)
    
    u = A \ (b-J*f) # solve linear system with a backsolve to obtain initial displacement.

    return (u, uexact, H̃, coord)
end

# Convergence tests:

# Define spatial order of accuracy
SBPp = 6

# Number of grid points
Nr = [2^5, 2^6, 2^7, 2^8, 2^9, 2^10]
Ns = Nr

errors = zeros(length(Nr))
rates = zeros(length(Nr)-1)

for n = 1:length(Nr)

    (u, uexact, H̃, coord) = antiplane_solve(Nr[n], Ns[n], SBPp)
    
    diff = u .- uexact
    errors[n] = sqrt(diff' * H̃ * diff) ./ sqrt(uexact' * H̃ * uexact)   # relative error
    
end

for n = 2:length(Nr)
    rates[n-1] = log2(errors[n-1]/errors[n])
end

@show errors
@show rates

# p = (2, 4, 6) should yield rates of (2, 3, 4), respectively.

# Plot error on coarsest grid
using Plots
pyplot();
(u, uexact, H̃, coord) = antiplane_solve(Nr[1], Ns[1], SBPp)
u = reshape(u, Nr[1]+1, Ns[1]+1)
uexact = reshape(uexact, Nr[1]+1, Ns[1]+1)
p1 = surface(coord[1], coord[2], uexact - u, xlabel = "x", ylabel = "z", showaxis = true,camera=(30,30))
display(p1)