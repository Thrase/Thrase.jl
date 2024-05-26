using DifferentialEquations
using Printf
using DelimitedFiles
using IterativeSolvers

# include("coefficients.jl")
include("helper.jl")
# include("domain.jl")

global const ctr = Ref{Int64}(1)

initialize_mg_struct_CUDA(mg_struct_CUDA, N_x, N_y, N_z, n_levels)

odeparam = (
    reject_step = [false],                          # to reject a step or not
    Vp = BP5_coeff.Vp,                              # plate rate
    M = M,                                          # LHS of the linear system M*u=RHS
    M_GPU = M_GPU,                                  # GPU array of the LHS system
    u = zeros(size(RHS)),                           # solution for the linear system 
    u_old = zeros(size(RHS)),                       # solution from the previous step
    u_GPU = CuArray(zeros(size(RHS))),              # GPU array for solutions
    Δτb = spzeros(2 * (N_x + 1) * (N_y + 1)),       # store the traction computed
    τb = spzeros(2 * (N_x + 1) * (N_y + 1)),        # shear stress vector \boldsymbol{τ} = [τ; τ_z]
    τfb = spzeros(2 * (N_x + 1) * (N_y + 1)),
    V2_v = fill(1e-9, fN2 * fN3),                   # acutal velocity (not in logical domain)
    V3_v = fill(1e-20, fN2 * fN3),                  # actual velocity (not in logical comain)
    V_v = fill(1e-9, fN2 * fN3),                    # norm of the velocity
    τ2_v = fill(13.0, fN2 * fN3),                     # traction in the second direction for the RS region
    τ3_v = fill(0.0, fN2 * fN3),                      # traction in the third direction for the RS region
    τ_v = fill(13.0, fN2 * fN3),                      # norm of the traction
    counter = [],                                   # counter for slip with Vmax >= threshold
    RHS = RHS,                                      # RHS of the linear system
    μshear = BP5_coeff.cs^2 * BP5_coeff.ρ ,         # constant?
    RSa = BP5_coeff.a0,                             # rate-and-state distance a0
    RSb = BP5_coeff.b0,                             # rate-and-state distance b0
    σn = BP5_coeff.σn,                              # effective normal stress
    η = BP5_coeff.cs * BP5_coeff.ρ / (2 * 1000) ,   # bug? should be \mu /(2 * cs) 
    RSV0 = BP5_coeff.V0,                            # rate-and-state reference slip rate
    τ0 = zeros((N_x + 1) * (N_y + 1)),              # pre-stress                                     # 
    RSL = BP5_coeff.L,                              # rate-and-state critical slip distance L
    RSLs = fill(BP5_coeff.L, fN2 * fN3),            # rate-and-state critical slip distance Ls (0.13/0.14)
    RSf0 = BP5_coeff.f0,                            # rate-and-state reference friction coefficient 
    N = N_x,                                        # number of grids in each direction, assuming idential of grid in x,y,z directions
    δNp = N_x + 1,                                  # number of grid points in each direction, assuming idential of grid in x,y,z directions
    Face_operators,                                 # getting face values from 3D SparseArrays
    updators,                                       # updating RHS values using SBP-SAT operators for Dirichlet Operations
    u_filters,                                      # filtering u1, u2, u3 from stacked u
    stride_time = 5,                                 # 
    RSas = zeros(fN2 * fN3),                         # RSas
    mg_struct_CUDA = mg_struct_CUDA
);



# ODE function
function odefun(dψV, ψδ, odeparam, t)
    @show t
    @unpack_namedtuple odeparam;
    if reject_step[1]
        return
    end

    ## Setting up ratees of change for state and slip
    dψ, V, ψ, δ = create_view(dψV, ψδ) # creating "views" to get dψ, V, ψ, δ
    dψ .= 0;
    V .= 0;
    ## End setting up dψV and ψδ

    # Updating RHS using δ
    RHS .= 0 # Resetting RHS
    RHS .+= updators[1] * δ[1:2:end] ./ 2 # divide slip by half to get displacements
    RHS .+= updators[2] * δ[2:2:end] ./ 2

    # Updating RHS using remote loading for face 2 for V2
    RHS .+= updators[3] * (fill(t .* Vp / 2, div(length(δ),2)))

    # End updating RHS using δ
    abstol_ = norm(RHS) * sqrt(eps(Float64))

    # Solving linear system using iterative methods
    if t == 0
        u_GPU, history = cg!(CuArray(u),M_GPU, CuArray(RHS), abstol=abstol_, log=true);    # solving with non preconditioned cg
    else
    # @show t, history.iters

    # End solving with cg!

    # Test solving with MGCG
        mg_struct_CUDA.b_mg[1] .= CuArray(RHS)
        # @assert norm(mg_struct_CUDA.b_mg[1]) != 0
        mg_struct_CUDA.x_CUDA[1] .= CuArray(u)
        mgcg_CUDA(mg_struct_CUDA,nx=N_x,ny=N_y,nz=N_z,n_levels=n_levels,precond=true,max_mg_iterations=1, v1=5, v2=5, v3=5, max_cg_iter=20, scaling_factor=1, rel_tol=1e-6, print_results=false)
        u_GPU = mg_struct_CUDA.x_CUDA[1]
    # End testing with MGCG test
    end

    u[:] .= Array(u_GPU)
    # End of solving 

    # updating values for traction
    Δτ = @view Δτb[1:2:length(Δτb)] 
    Δτz = @view Δτb[2:2:length(Δτb)]

    τ0 = @view τb[1:2:length(τb)]
    τz0 =  @view τb[2:2:length(τb)]

    Δτ .= Face_operators[1] * sigma_21 * u
    Δτz .= Face_operators[1] * sigma_31 * u
    # finish updating values for traction
    
    odeparam.τfb .= τb .+ Δτb 

    # getting tractions for the RS region
    τ2_v .= (τ0 + Δτ)[RS_filter_2D_nzind]
    τ3_v .= (τz0 + Δτz)[RS_filter_2D_nzind]
    τ_v .= hypot.(τ2_v, τ3_v)
    # end of getting tractions for the RS region

    # bisection guarded newton's method
    xL = fill(0.0, length(τ_v))
    xR = τ_v ./ η
    (V_v_tmp, f_v, iter) = newtbndv_vectorized(rateandstate_vectorized, xL, xR, V_v, ψ, σn, τ_v, η,
                                    RSas, RSV0; ftol=1e-6, maxiter=500, minchange=0, atolx = 1e-4, rtolx=1e-4)

    # end of bisection guarded newton's method

    # calculating V2_v and V3_v from V_v
    V_v .= V_v_tmp
    V2_v = V_v .* τ2_v ./ τ_v
    V3_v = V_v .* τ3_v ./ τ_v
    # end of calculating V2_v and V3_v from V_v

    
    # rejecting if V2 or V3 has infinite entries
    if !all(isfinite.(V2_v)) || !all(isfinite.(V3_v))
        println("V reject")
        reject_step[1] = true
        return
    end

    # or newton's method does not converge
    if iter < 0
        println("iter reject")
        reject_step[1] = true
        return
    end
    # end of rejecting from V2, V3, or iter



    if iter > 1
        @show iter
    end

    # out side of RS, V2 = Vp, V3 = 0
    V[1:2:end] .= Vp 
    V[2:2:end] .= 0


    V[2 .* RS_filter_2D_nzind .- 1] .= V2_v
    V[2 .* RS_filter_2D_nzind] .= V3_v

    # Updating ψ based on iteration convergence
    # dψ[n] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0) # BP1
    # dψ .= (RSb * RSV0 / RSL) .* (exp.((RSf0 .- ψ) ./ RSb) .- sqrt.(V2_v.^2 .+ V3_v.^2) ./ RSV0)
    if iter > 0
        dψ .= (RSb * RSV0 ./ RSLs) .* (exp.((RSf0 .- ψ) ./ RSb) .- sqrt.(V2_v.^2 .+ V3_v.^2) ./ RSV0)
    else
        dψ .= 0
    end
    # End of updating ψ

    # Rejecting dψ if there are infinite values
    if !all(isfinite.(dψ))
        println("ψ reject")
        dψ .= 0
        reject_step[1] = true
        return
    end
    # End of rejecting
    nothing
end