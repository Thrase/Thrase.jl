include("../diagonal_sbp.jl")
include("../3D_face.jl")
include("../analy_sol.jl")
include("../components.jl")
include("../coefficients.jl")

using LinearAlgebra
using IterativeSolvers
using BenchmarkTools
using Plots
using SparseArrays
using CUDA


function Assembling_3D_matrices(N_x, N_y, N_z;SBPp=2, Lx=128, Ly=128, Lz=128)
    Lx = Lx
    Ly = Ly
    Lz = Lz

    hx = 1 / N_x;
    hy = 1 / N_y;
    hz = 1 / N_z;

    x = collect(range(0,step=hx,1));
    y = collect(range(0,step=hy,1));
    z = collect(range(0,step=hz,1));

    Nx = N_x + 1;
    Ny = N_y + 1;
    Nz = N_z + 1;

    (D1x, HIx, H1x, r1x) = diagonal_sbp_D1(SBPp,N_x,xc=(0,Lx));
    (D2x, S0x, SNx, HI2x, H2x, r2x) = diagonal_sbp_D2(SBPp,N_x,xc=(0,Lx));


    (D1y, HIy, H1y, r1y) = diagonal_sbp_D1(SBPp,N_y,xc=(0,Ly));
    (D2y, S0y, SNy, HI2y, H2y, r2y) = diagonal_sbp_D2(SBPp,N_y,xc=(0,Ly));


    (D1z, HIz, H1z, r1z) = diagonal_sbp_D1(SBPp,N_y,xc=(0,Lz));
    (D2z, S0z, SNz, HI2z, H2z, r2z) = diagonal_sbp_D2(SBPp,N_y,xc=(0,Lz));

    BSx = sparse(SNx - S0x);
    BSy = sparse(SNy - S0y);
    BSz = sparse(SNz - S0z);

    H_tilde = kron(H1x,H1y,H1z);
    HI_tilde = kron(HIx, HIy, HIz);

    I_Nx = eyes(N_x+1);
    I_Ny = eyes(N_y+1);
    I_Nz = eyes(N_z+1);

    D2_x = kron(I_Nz,I_Ny,D2x);
    D2_y = kron(I_Nz,D2y,I_Nx);
    D2_z = kron(D2z,I_Nx,I_Ny);

    e_1x = e(1,N_x + 1);
    e_1y = e(1,N_y + 1);
    e_1z = e(1,N_z + 1);

    e_Nx = e(N_x + 1,N_x + 1);
    e_Ny = e(N_y + 1,N_y + 1);
    e_Nz = e(N_z + 1,N_z + 1);


    # H and HI Operators for End and Front
    # Needs to redefine the shape to be N^3 by N^3 
    # Size: N^2 by N^2
    H_1 = kron(H1z,H1y);
    HI_1 = kron(I_Nz,I_Ny,HIx);

    H_2 = kron(H1z, H1y);
    HI_2 = kron(I_Nz,I_Ny,HIx);

    # H and HI operators for Left and Right
    # Size: N^2 by N^2
    H_3 = kron(H1z, H1x);
    HI_3 = kron(I_Ny,HIy,I_Nx);

    H_4 = kron(H1z, H1x);
    HI_4 = kron(I_Ny,HIy,I_Nx);

    # H and HI operators for Bottom and Top
    # Size: N^2 by N^2

    H_5 = kron(H1y, H1x);
    HI_5 = kron(HIz,I_Ny,I_Nx);

    H_6 = kron(H1y, H1x);
    HI_6 = kron(HIz,I_Ny,I_Nx);

    # BS operators for 6 faces
    # Size: N^3 by N^3



    BS_Front = kron(I_Ny,I_Nz,BSx);
    BS_End = kron(I_Ny,I_Nz,BSx);
    BS_Left = kron(I_Nz,BSy,I_Nx);
    BS_Right = kron(I_Nz,BSy,I_Nx);
    # BS_Bottom = kron(I_Nx,I_Ny,BSz)
    # BS_Top = kron(I_Nx,I_Ny,BSz)
    BS_Bottom = kron(BSz,I_Nx,I_Ny);
    BS_Top = kron(BSz,I_Nx,I_Ny);

    Front_operator = get_front_face(N_x+1,N_y+1,N_z+1);
    End_operator = get_end_face(N_x+1,N_y+1,N_z+1);
    Left_operator = get_left_face(N_x+1,N_y+1,N_z+1);
    Right_operator = get_right_face(N_x+1,N_y+1,N_z+1);
    Top_operator = get_top_face(N_x+1,N_y+1,N_z+1);
    Bottom_operator = get_bottom_face(N_x+1,N_y+1,N_z+1);

    # Penalty Parameters
    tau_x = 13/hx;
    tau_y = 13/hy;
    tau_z = 13/hz;


    beta = -1;



    # numerical_sol_3D = reshape(numerical_sol,N_x+1,N_y+1,N_z+1)

    # analy_sol_3D u1 = u2 = u3 = sin(πx + πy + πz)

    u1_filter = get_u1(Nx,Ny,Nz);
    u2_filter = get_u2(Nx,Ny,Nz);
    u3_filter = get_u3(Nx,Ny,Nz);




    # replacing D1x D1y D1z with BSx BSy BSz
    p_px_new = kron(I_Nz,I_Ny,D1x);
    p_py_new = kron(I_Nz,D1x,I_Nx);
    p_pz_new = kron(D1x,I_Ny,I_Nx);

    p_px_hat_new = kron(I_Nz, I_Ny, BSx);
    p_py_hat_new = kron(I_Nz, BSy, I_Nx);
    p_pz_hat_new = kron(BSz, I_Ny, I_Nx);


    # Second order derivatives
    p2_px2_new = kron(I_Nz,I_Ny,D2x);
    p2_py2_new = kron(I_Nz,D2y,I_Nx);
    p2_pz2_new = kron(D2z,I_Ny,I_Nx);

    # crossterms

    p2_pypx_new = kron(I_Nz,D1y,D1x); # equivalent to p_py * p_px ? actually true 
    p2_pxpy_new = kron(I_Nz,D1y,D1x);

    p2_pzpy_new = kron(D1z,D1y,I_Nx);
    p2_pypz_new = kron(D1z,D1y,I_Nx);

    p2_pxpz_new = kron(D1z,I_Ny,D1x);
    p2_pzpx_new = kron(D1z,I_Ny,D1x);

    # Express σ tensors as operators on u vector (u1, u2, u3 stacked)
    sigma_11_new = (K_v - 2/3*μ_v) * (p_px_new * u1_filter + p_py_new * u2_filter + p_pz_new * u3_filter) + 2 * μ_v * p_px_new * u1_filter;
    sigma_12_new = μ_v*(p_py_new*u1_filter + p_px_new*u2_filter);
    sigma_13_new = μ_v*(p_pz_new*u1_filter + p_px_new*u3_filter);


    sigma_21_new = μ_v*(p_px_new*u2_filter + p_py_new*u1_filter); 
    sigma_22_new = (K_v - 2/3*μ_v) * (p_px_new * u1_filter + p_py_new * u2_filter + p_pz_new * u3_filter) + 2 * μ_v * p_py_new * u2_filter;
    sigma_23_new = μ_v * (p_pz_new * u2_filter + p_py_new * u3_filter);

    sigma_31_new = μ_v * (p_px_new * u3_filter + p_pz_new * u1_filter);
    sigma_32_new = μ_v * (p_py_new * u3_filter + p_pz_new * u2_filter);
    sigma_33_new = (K_v - 2/3 * μ_v) * (p_px_new * u1_filter + p_py_new * u2_filter + p_pz_new * u3_filter) + 2 * μ_v * p_pz_new * u3_filter;



    # # TO DO: Need to rewrite this part 
    u1_operator_new =  ( (K_v - 2/3 * μ_v) * (p2_px2_new * u1_filter + p2_pxpy_new * u2_filter + p2_pxpz_new * u3_filter) 
                + 2 * μ_v * p2_px2_new * u1_filter 
                + μ_v * (p2_py2_new * u1_filter + p2_pxpy_new * u2_filter)
                + μ_v * (p2_pz2_new * u1_filter + p2_pxpz_new * u3_filter)
    );

    u2_operator_new = ( μ_v * (p2_px2_new * u2_filter + p2_pxpy_new * u1_filter)
                + (K_v - 2/3 * μ_v) * (p2_pxpy_new * u1_filter + p2_py2_new * u2_filter + p2_pypz_new * u3_filter)
                + 2 * μ_v * p2_py2_new * u2_filter
                + μ_v * (p2_pz2_new * u2_filter + p2_pypz_new * u3_filter)
    );

    u3_operator_new = ( μ_v * (p2_px2_new * u3_filter + p2_pxpz_new * u1_filter)
                + μ_v * (p2_py2_new * u3_filter + p2_pypz_new * u2_filter)
                + (K_v - 2/3 * μ_v) * (p2_pxpz_new * u1_filter + p2_pypz_new * u2_filter + p2_pz2_new * u3_filter)
                + 2 * μ_v * p2_pz2_new * u3_filter
    );






    ### Assembling SBP terms according to the note

    ### Assembling components of stress tensors and boundary operators
    ### Face 1
    e_1 = End_operator';
    e_1T = End_operator;

    # new formulation

    T_11_1_new = (K_v + 4/3 * μ_v) * p_px_hat_new; #* u1_filter
    T_12_1_new = - (K_v - 2/3 * μ_v) * p_py_new; #* u1_filter # Not quite sure 
    T_13_1_new = - (K_v - 2/3 * μ_v) * p_pz_new; #* u1_filter

    T_21_1_new = - μ_v * p_py_new; #* u2_filter
    T_22_1_new = μ_v * p_px_hat_new; #* u2_filter
    T_23_1_new = 0;  #* u2_filter

    T_31_1_new = - μ_v * p_pz_new; #* u3_filter
    T_32_1_new = 0; #* u3_filter
    T_33_1_new = μ_v * p_px_hat_new; #* u3_filter


    ## TO DO Fix Z values it should be of the same size as T values
    Z_11_1_new = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx);#* u1_filter
    Z_12_1_new = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u1_filter
    Z_13_1_new = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u1_filter

    Z_21_1_new =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u2_filter
    Z_22_1_new = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx);#* u2_filter
    Z_23_1_new =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u2_filter ## 0 ?

    Z_31_1_new = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx);#* u3_filter
    Z_32_1_new = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx);#* u3_filter
    Z_33_1_new = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx); #* u3_filter


    ### Face 2
    e_2 = Front_operator';
    e_2T = Front_operator;


    T_11_2_new = (K_v + 4/3 * μ_v)  * p_px_hat_new; #* u1_filter
    T_12_2_new = (K_v - 2/3 * μ_v) * p_py_new; #* u1_filter
    T_13_2_new = (K_v - 2/3 * μ_v) * p_pz_new; #* u1_filter

    T_21_2_new = μ_v * p_py_new; #* u2_filter
    T_22_2_new = μ_v * p_px_hat_new; #* u2_filter
    T_23_2_new = 0;# * u2_filter

    T_31_2_new = μ_v * p_pz_new; #* u3_filter
    T_32_2_new = 0; #* u3_filter
    T_33_2_new = μ_v * p_px_hat_new; #* u3_filter

    ## TO DO Fix Z values
    Z_11_2_new = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx); #* u1_filter
    Z_12_2_new = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u1_filter
    Z_13_2_new = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u1_filter

    Z_21_2_new =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u2_filter
    Z_22_2_new = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx); #* u2_filter
    Z_23_2_new =  (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u2_filter ## 0 ?

    Z_31_2_new = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u3_filter
    Z_32_2_new = (d * β/ H1x[1]) * (K_v - 2/3 * μ_v + μ_v) * 0 * kron(I_Nz, I_Ny, I_Nx); #* u3_filter
    Z_33_2_new = (d * β / H1x[1]) * (K_v + 4/3 * μ_v + 2 * μ_v) * kron(I_Nz, I_Ny, I_Nx);#* u3_filter



    ### Face 3
    e_3 = Left_operator';
    e_3T = Left_operator;

    T_11_3_new = μ_v * p_py_hat_new; #* u1_filter
    T_12_3_new = - μ_v * p_px_new; #* u2_filter
    T_13_3_new = 0; # Face 1

    T_21_3_new = - (K_v - 2/3 * μ_v) * p_px_new; #* u1_filter
    T_22_3_new = (K_v + 4/3 * μ_v) * p_py_hat_new; #* u2_filter
    T_23_3_new = - (K_v - 2/3 * μ_v) * p_pz_new; #* u3_filter

    T_31_3_new = 0; #* u1_filter
    T_32_3_new = - μ_v * p_pz_new; #* u2_filter
    T_33_3_new = μ_v * p_py_hat_new; #* u3_filter

    ### Face 4
    e_4 = Right_operator';
    e_4T = Right_operator;

    T_11_4_new = μ_v * p_py_hat_new; #* u1_filter
    T_12_4_new = μ_v * p_px_new; #* u2_filter
    T_13_4_new = 0; #* u3_filter

    T_21_4_new = (K_v - 2/3 * μ_v) * p_px_new; #* u1_filter
    T_22_4_new = (K_v + 4/3 * μ_v) * p_py_hat_new; #* u2_filter
    T_23_4_new = (K_v - 2/3 * μ_v) * p_pz_new; #* u3_filter

    T_31_4_new = 0; #* u1_filter
    T_32_4_new = μ_v * p_pz_new; #* u2_filter
    T_33_4_new = μ_v * p_py_hat_new; #* u3_filter

    ### Face 5
    e_5 = Bottom_operator';
    e_5T = Bottom_operator;

    T_11_5_new = μ_v * p_pz_hat_new; #* u1_filter
    T_12_5_new = 0;#* u2_filter
    T_13_5_new = - μ_v * p_px_new; #* u3_filter

    T_21_5_new = 0;#* u1_filter
    T_22_5_new = μ_v * p_pz_hat_new; #* u2_filter
    T_23_5_new = - μ_v * p_py_new; #* u3_filter

    T_31_5_new = - (K_v - 2/3 * μ_v) * p_px_new; #* u1_filter
    T_32_5_new = - (K_v - 2/3 * μ_v) * p_py_new; #* u2_filter
    T_33_5_new = (K_v + 4/3 * μ_v) * p_pz_hat_new; #* u3_filter

    ### Face 6
    e_6 = Top_operator';
    e_6T = Top_operator;

    T_11_6_new = μ_v * p_pz_hat_new; #* u1_filter
    T_12_6_new = 0;#* u2_filter
    T_13_6_new = μ_v * p_px_new; #* u3_filter

    T_21_6_new = 0;#* u1_filter
    T_22_6_new = μ_v * p_pz_hat_new; #* u2_filter
    T_23_6_new = μ_v * p_py_new; #* u3_filter

    T_31_6_new = (K_v - 2/3 * μ_v) * p_px_new; #* u1_filter
    T_32_6_new = (K_v - 2/3 * μ_v) * p_py_new; #* u2_filter
    T_33_6_new = (K_v + 4/3 * μ_v) * p_pz_hat_new; #* u3_filter

    ### Assembling SBP terms for left-hand-side (LHS) traction condition

    SAT_1_LHS_new = beta * HI_tilde * (
            e_3 * H_3 * e_3T * (T_11_3_new * u1_filter .+ T_12_3_new * u2_filter .+ T_13_3_new * u3_filter)
        +   e_4 * H_4 * e_4T * (T_11_4_new * u1_filter .+ T_12_4_new * u2_filter .+ T_13_4_new * u3_filter)
        +   e_5 * H_5 * e_5T * (T_11_5_new * u1_filter .+ T_12_5_new * u2_filter .+ T_13_5_new * u3_filter)
        +   e_6 * H_6 * e_6T * (T_11_6_new * u1_filter .+ T_12_6_new * u2_filter .+ T_13_6_new * u3_filter)

    ); 

    SAT_2_LHS_new = beta * HI_tilde * (
            e_3 * H_3 * e_3T * (T_21_3_new * u1_filter .+ T_22_3_new * u2_filter .+ T_23_3_new * u3_filter)
        +   e_4 * H_4 * e_4T * (T_21_4_new * u1_filter .+ T_22_4_new * u2_filter .+ T_23_4_new * u3_filter)
        +   e_5 * H_5 * e_5T * (T_21_5_new * u1_filter .+ T_22_5_new * u2_filter .+ T_23_5_new * u3_filter)
        +   e_6 * H_6 * e_6T * (T_21_6_new * u1_filter .+ T_22_6_new * u2_filter .+ T_23_6_new * u3_filter)

    ); 


    SAT_3_LHS_new = beta * HI_tilde * (
            e_3 * H_3 * e_3T * (T_31_3_new * u1_filter .+ T_32_3_new * u2_filter .+ T_33_3_new * u3_filter)
        +   e_4 * H_3 * e_4T * (T_31_4_new * u1_filter .+ T_32_4_new * u2_filter .+ T_33_4_new * u3_filter)
        +   e_5 * H_3 * e_5T * (T_31_5_new * u1_filter .+ T_32_5_new * u2_filter .+ T_33_5_new * u3_filter)
        +   e_6 * H_3 * e_6T * (T_31_6_new * u1_filter .+ T_32_6_new * u2_filter .+ T_33_6_new * u3_filter)
    );




    ### Assembling SBP terms for left-hand-side (LHS) Dirichlet condition


    SAT_tilde_1_LHS_new = HI_tilde * (
            (T_11_1_new' .- Z_11_1_new') * (e_1 * H_1 * (e_1T)) * u1_filter
        +   (T_21_1_new' .- Z_21_1_new') * (e_1 * H_1 * (e_1T)) * u2_filter
        +   (T_31_1_new' .- Z_31_1_new') * (e_1 * H_1 * (e_1T)) * u3_filter
        +   (T_11_2_new' .- Z_11_2_new') * (e_2 * H_2 * (e_2T)) * u1_filter
        +   (T_21_2_new' .- Z_21_2_new') * (e_2 * H_2 * (e_2T)) * u2_filter
        +   (T_31_2_new' .- Z_31_2_new') * (e_2 * H_2 * (e_2T)) * u3_filter
    );

    SAT_tilde_2_LHS_new = HI_tilde * (
            (T_12_1_new' .- Z_12_1_new') * (e_1 * H_1 * (e_1T)) * u1_filter
        +   (T_22_1_new' .- Z_22_1_new') * (e_1 * H_1 * (e_1T)) * u2_filter
        +   (T_32_1_new' .- Z_32_1_new') * (e_1 * H_1 * (e_1T)) * u3_filter
        +   (T_12_2_new' .- Z_12_2_new') * (e_2 * H_2 * (e_2T)) * u1_filter
        +   (T_22_2_new' .- Z_22_2_new') * (e_2 * H_2 * (e_2T)) * u2_filter
        +   (T_32_2_new' .- Z_32_2_new') * (e_2 * H_2 * (e_2T)) * u3_filter
    );

    SAT_tilde_3_LHS_new = HI_tilde * (
            (T_13_1_new' .- Z_13_1_new') * (e_1 * H_1 * (e_1T)) * u1_filter
        +   (T_23_1_new' .- Z_23_1_new') * (e_1 * H_1 * (e_1T)) * u2_filter
        +   (T_33_1_new' .- Z_33_1_new') * (e_1 * H_1 * (e_1T)) * u3_filter
        +   (T_13_2_new' .- Z_13_2_new') * (e_2 * H_2 * (e_2T)) * u1_filter
        +   (T_23_2_new' .- Z_23_2_new') * (e_2 * H_2 * (e_2T)) * u2_filter
        +   (T_33_2_new' .- Z_33_2_new') * (e_2 * H_2 * (e_2T)) * u3_filter
    );



    # Forming analytical solutions
    u1 = form_analy_sol(;N = N_x)[1][:]; # u1 is the only non-zero component
    u2 = form_analy_sol(;N = N_x)[2][:]; # u2 = 0 for the test case
    u3 = form_analy_sol(;N = N_x)[3][:]; # u3 = 0 for the test case

    analy_sol = u1_filter' * u1 + u2_filter' * u2 + u3_filter' * u3;



    # Assembling boundary conditions
    # Getting boundary values

    # u1
    u1_Front_value = Front_operator' * u1_Front(y,z)[:]; # Dirichlet Conditions
    u1_End_value = End_operator' * u1_End(y,z)[:]; # Dirichlet Conditions

    u1_Top_value = Top_operator' * u1_Top(x,y)[:]; # Dirichlet Conditions
    u1_Bottom_value = Bottom_operator' * u1_Top(x,y)[:]; # Dirichlet Conditions

    u1_Left_value = Left_operator' * u1_y_Left(x,z)[:]; # Neumann Conditions
    u1_Right_value = Right_operator' * u1_y_Right(x,z)[:];

    # u2
    u2_Front_value = Front_operator' * u1_Front(y,z)[:]; # Dirichlet Conditions
    u2_End_value = End_operator' * u1_End(y,z)[:]; # Dirichlet Conditions

    u2_Top_value = Top_operator' * u1_Top(x,y)[:]; # Dirichlet Conditions
    u2_Bottom_value = Bottom_operator' * u1_Top(x,y)[:]; # Dirichlet Conditions

    u2_Left_value = Left_operator' * u1_y_Left(x,z)[:]; # Neumann Conditions
    u2_Right_value = Right_operator' * u1_y_Right(x,z)[:];

    # u3
    u3_Front_value = Front_operator' * u1_Front(y,z)[:]; # Dirichlet Conditions
    u3_End_value = End_operator' * u1_End(y,z)[:]; # Dirichlet Conditions

    u3_Top_value = Top_operator' * u1_Top(x,y)[:]; # Dirichlet Conditions
    u3_Bottom_value = Bottom_operator' * u1_Top(x,y)[:]; # Dirichlet Conditions

    u3_Left_value = Left_operator' * u1_y_Left(x,z)[:]; # Neumann Conditions
    u3_Right_value = Right_operator' * u1_y_Right(x,z)[:];



    # # Assembling left hand side

    E1 = (u1_filter' * H_tilde * u1_operator_new);

    E2 = (u2_filter' * H_tilde * u2_operator_new);

    E3 = (u3_filter' * H_tilde * u3_operator_new);

    E = (E1 + E2 + E3);

    # Assembling right hand side
    # Assembling source source_terms
    source_u1 = u1_filter' * H_tilde * ((K_v - 2/3 * μ_v) * (-π^2 * u1_analy(x,y,z)[:] -π^2 * u2_analy(x,y,z)[:]) 
                + 2 * μ_v * (-π^2 * u1_analy(x,y,z)[:])
                + μ_v * (-π^2 * u1_analy(x,y,z)[:] -π^2 * u2_analy(x,y,z)[:]) 
                + μ_v * (-π^2 * u1_analy(x,y,z)[:])
                );

    source_u2 = u2_filter' * H_tilde * (μ_v * (-π^2 * u1_analy(x,y,z)[:] - π^2 * u2_analy(x,y,z)[:])
                + (K_v - 2/3 * μ_v) * (-π^2 * u1_analy(x,y,z)[:] -π^2 * u2_analy(x,y,z)[:] ) 
                + 2 * μ_v * (-π^2 * u2_analy(x,y,z)[:])
                + μ_v * (-π^2 * u2_analy(x,y,z)[:])
                );

    source_u3 = u3_filter' * H_tilde * (μ_v * (-π^2 * u1_analy(x,y,z)[:])
                + μ_v * (-π^2 * u2_analy(x,y,z)[:])
                + (K_v - 2/3 * μ_v) * (-π^2 * u1_analy(x,y,z)[:] + -π^2 * u2_analy(x,y,z)[:])
                # u3 is set to be zero 
                );
    source = ( source_u1 + source_u2 + source_u3);

    # Assembling boundary data
    # Face 1: Dirichlet
    g₁¹ = u1_End(y,z);
    g₂¹ = u2_End(y,z);
    g₃¹ = zeros(Ny,Nz);

    # Face 2: Dirichlet
    g₁² = u1_Front(y,z);
    g₂² = u2_Front(y,z);
    g₃² = zeros(Ny,Nz);

    # Face 3: Neumann
    g₁³ = (u1_y_Left(x,z) + u2_x_Left(x,z));
    g₂³ = ((K_v - 2/3 * μ_v) * u1_x_Left(x,z) + (K_v + 4/3 * μ_v) * u2_y_Left(x,z));
    g₃³ = zeros(Nx,Nz);

    # Face 4: Neumann
    g₁⁴ = u1_y_Right(x,z) + u2_x_Right(x,z);
    g₂⁴ = (K_v - 2/3 * μ_v) * u1_x_Right(x,z) + (K_v + 4/3 * μ_v) * u2_y_Right(x,z);
    g₃⁴ = zeros(Nx,Nz);

    # Face 5: Neumann
    g₁⁵ = u1_z_Bottom(x,y);
    g₂⁵ = u2_z_Bottom(x,y);
    g₃⁵ = (K_v - 2/3 * μ_v) * (u1_x_Bottom(x,y) + u2_y_Bottom(x,y));

    # Face 6: Neumann
    g₁⁶ = u1_z_Top(x,y);
    g₂⁶ = u2_z_Top(x,y);
    g₃⁶ = (K_v - 2/3 * μ_v) * (u1_x_Top(x,y) + u2_y_Top(x,y));


    ### Assembling SBP terms for right-hand-side (RHS) traction condition
    SAT_1_RHS_new = beta * HI_tilde * (
            e_3 * H_3 * g₁³[:]
        +   e_4 * H_4 * g₁⁴[:]
        +   e_5 * H_5 * g₁⁵[:]
        +   e_6 * H_6 * g₁⁶[:]
    );

    SAT_2_RHS_new = beta * HI_tilde * (
            e_3 * H_3 * g₂³[:]
        +   e_4 * H_4 * g₂⁴[:]
        +   e_5 * H_5 * g₂⁵[:]
        +   e_6 * H_6 * g₂⁶[:]
    );

    SAT_3_RHS_new = beta * HI_tilde * (
            e_3 * H_3 * g₃³[:]
        +   e_4 * H_4 * g₃⁴[:]
        +   e_5 * H_5 * g₃⁵[:]
        +   e_6 * H_6 * g₃⁶[:]
    );



    ### Assembling SBP terms for right-hand-side (RHS) Dirichlet condition

    SAT_tilde_1_RHS_new =  HI_tilde * (
            (T_11_1_new' .- Z_11_1_new') * (e_1 * H_1 * g₁¹[:])
        +   (T_21_1_new' .- Z_21_1_new') * (e_1 * H_1 * g₂¹[:])
        +   (T_31_1_new' .- Z_31_1_new') * (e_1 * H_1 * g₃¹[:])
        +   (T_11_2_new' .- Z_11_2_new') * (e_2 * H_2 * g₁²[:]) # face 1 Dirichlet
        +   (T_21_2_new' .- Z_21_2_new') * (e_2 * H_2 * g₂²[:])
        +   (T_31_2_new' .- Z_31_2_new') * (e_2 * H_2 * g₃²[:])
    );

    SAT_tilde_2_RHS_new =  HI_tilde * (
            (T_12_1_new' .- Z_12_1_new') * (e_1 * H_1 * g₁¹[:])
        +   (T_22_1_new' .- Z_22_1_new') * (e_1 * H_1 * g₂¹[:])
        +   (T_32_1_new' .- Z_32_1_new') * (e_1 * H_1 * g₃¹[:])
        +   (T_12_2_new' .- Z_12_2_new') * (e_2 * H_2 * g₁²[:]) # face 1 Dirichlet
        +   (T_22_2_new' .- Z_22_2_new') * (e_2 * H_2 * g₂²[:])
        +   (T_32_2_new' .- Z_31_2_new') * (e_2 * H_2 * g₃²[:])
    );

    SAT_tilde_3_RHS_new =  HI_tilde * (
            (T_13_1_new' .- Z_13_1_new') * (e_1 * H_1 * g₁¹[:])
        +   (T_23_1_new' .- Z_23_1_new') * (e_1 * H_1 * g₂¹[:])
        +   (T_33_1_new' .- Z_33_1_new') * (e_1 * H_1 * g₃¹[:])
        +   (T_13_2_new' .- Z_13_2_new') * (e_2 * H_2 * g₁²[:]) # face 1 Dirichlet
        +   (T_23_2_new' .- Z_23_2_new') * (e_2 * H_2 * g₂²[:])
        +   (T_33_2_new' .- Z_33_2_new') * (e_2 * H_2 * g₃²[:])
    );

    # Traction updating operators
    δ2_update = - (
            u1_filter' * H_tilde * HI_tilde * (T_21_1_new' .- Z_21_1_new') * (e_1 * H_1)
        +   u2_filter' * H_tilde * HI_tilde * (T_22_1_new' .- Z_22_1_new') * (e_1 * H_1)
        +   u3_filter' * H_tilde * HI_tilde * (T_23_1_new' .- Z_23_1_new') * (e_1 * H_1)
    )

    δ3_update = - (
            u1_filter' * H_tilde * HI_tilde * (T_31_1_new' .- Z_31_1_new') * (e_1 * H_1)
        +   u2_filter' * H_tilde * HI_tilde * (T_32_1_new' .- Z_32_1_new') * (e_1 * H_1)
        +   u3_filter' * H_tilde * HI_tilde * (T_33_1_new' .- Z_33_1_new') * (e_1 * H_1)
    )

    face_2_V2_update = - ( # something not quite right here
            u1_filter' * H_tilde * HI_tilde * (T_21_2_new' .- Z_21_2_new') * (e_2 * H_2)
        +   u2_filter' * H_tilde * HI_tilde * (T_22_2_new' .- Z_22_2_new') * (e_2 * H_2)
        +   u3_filter' * H_tilde * HI_tilde * (T_23_2_new' .- Z_23_2_new') * (e_2 * H_2)
    )

    face_2_V3_update = - (
            u1_filter' * H_tilde * HI_tilde * (T_31_2_new' .- Z_31_2_new') * (e_2 * H_2)
        +   u2_filter' * H_tilde * HI_tilde * (T_32_2_new' .- Z_31_2_new') * (e_2 * H_2)
        +   u3_filter' * H_tilde * HI_tilde * (T_33_2_new' .- Z_33_2_new') * (e_2 * H_2)
    )



    # Assembling LHS of the linear system

    M_new = - (E + u1_filter' * H_tilde * SAT_1_LHS_new
            + u2_filter' * H_tilde * SAT_2_LHS_new 
            + u3_filter' * H_tilde * SAT_3_LHS_new 
            + u1_filter' * H_tilde * SAT_tilde_1_LHS_new 
            + u2_filter' * H_tilde * SAT_tilde_2_LHS_new 
            + u3_filter' * H_tilde * SAT_tilde_3_LHS_new
        );

    RHS_new = - (source + u1_filter' * H_tilde * SAT_1_RHS_new 
            + u2_filter' * H_tilde * SAT_2_RHS_new 
            + u3_filter' * H_tilde * SAT_3_RHS_new
            + u1_filter' * H_tilde * SAT_tilde_1_RHS_new 
            + u2_filter' * H_tilde * SAT_tilde_2_RHS_new 
            + u3_filter' * H_tilde * SAT_tilde_3_RHS_new);

    return M_new, RHS_new, H_tilde, HI_tilde, analy_sol, source, 
            [T_11_1_new, T_12_1_new, T_13_1_new, T_21_1_new, T_22_1_new, T_23_1_new, T_31_1_new, T_32_1_new, T_33_1_new], 
            [u1_filter, u2_filter, u3_filter], 
            [End_operator, Front_operator], 
            [sigma_11_new, sigma_21_new, sigma_31_new],
            [δ2_update, δ3_update, face_2_V2_update, face_2_V3_update];
end
