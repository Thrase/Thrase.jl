
const ϵ = 0
const ϵ1 = 0
const γ = 0
const γ1 = 0
const α = 0
const β = 1
const β1 = 1
const β2 = 1
const ν = 0
const c = 1
const c1 = 1
const b1 = 1
const b2 = 1
const L = 8


function exactu(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        uexact = (x .+ 1) .* (x .- 1) .* (y .+ 1) .* (y .- 1) .+ β .* sin.(pi * x ./ L) .* cos.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        uexact = 0 .* x
    end 
    return uexact[:]
end

function exactu_x(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return (2 .* x) .* (y .+ 1) .* (y .- 1) .+ β .* (pi/L) * cos.(pi * x ./ L) .* cos.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end
  
function exactu_xx(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 2 .* (y .+ 1) .* (y .- 1) .+ -β .* (pi/L)^2 * sin.(pi * x ./ L)  .* cos.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end

function exactu_y(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return (x .+ 1) .* (x .- 1) .* (2 .* y) .+ -β .* (pi/L) * sin.(pi * x ./ L)  .* sin.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end

function exactu_yy(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 2 .* (x .+ 1) .* (x .- 1) .+ -β .* (pi/L)^2 * sin.(pi * x ./ L) .* cos.(pi * y ./ L) .* cos.(pi * z ./ L)
    else
        return 0 .* x
    end
end
function exactu_z(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .- β .* (pi/L) .* sin.(pi * x ./ L) .* cos.(pi * y ./ L) .* sin.(pi * z ./ L)
    else
        return 0 .* x
    end
end
  
function exactu_zz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .- β .* (pi/L)^2 .* sin.(pi * x ./ L) .* cos.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end


function exactu_xy(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 2 .* x .* 2 .* y .+ β .* (pi/L) .* (-pi/L) .* cos.(pi * x ./ L) .* sin.(pi * y ./ L) .* cos.(pi * z ./ L) 

    else
        return 0 .* x
    end
end

function exactu_xz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ β .* (pi/L) * (-pi/L) * cos.(pi * x ./ L) .* cos.(pi * y ./ L) .* sin.(pi * z ./ L) 

    else
        return 0 .* x
    end
end

function exactu_yz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ β .* (pi/L)^2 .* sin.(pi * x ./ L) .* sin.(pi * y ./ L) .* sin.(pi * z ./ L) 
    else
        return 0 .* x
    end
end

#############################################

function exactv(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        vexact = 0 .* x .^ 3 .+ β1 .* sin.(pi * x ./ L) .* sin.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        vexact = 0 .* x
    end
    return vexact[:]
end

function exactv_x(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* 3 .* x .^ 2 .+ β1 .* (pi/L) * cos.(pi * x ./ L) .* sin.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end
  
function exactv_xx(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* 6 .* x .+ β1 .* (pi/L) .* (-pi/L) * sin.(pi * x ./ L) .* sin.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end

function exactv_y(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ β1 .* (pi/L) * sin.(pi * x ./ L) .* cos.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end
  
function exactv_yy(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ β1 .* (pi/L)*(-pi/L) .* sin.(pi * x ./ L) .* sin.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end

function exactv_z(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ β1 .* (-pi/L) * sin.(pi * x ./ L) .* sin.(pi * y ./ L) .* sin.(pi * z ./ L) 
    else
        return 0 .* x
    end
end
  
function exactv_zz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ β1 .* (-pi/L) * (pi/L) .* sin.(pi * x ./ L) .* sin.(pi * y ./ L) .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end
end

function exactv_xz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ β1 .* (pi/L) .* (-pi/L) .* cos.(pi * x ./ L) .* sin.(pi * y ./ L) .* sin.(pi * z ./ L) 
    else
        return 0 .* x
    end
end

function exactv_yz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return  0 .* x .+ β1 .* (pi/L) .* (-pi/L) .* sin.(pi * x ./ L) .* cos.(pi * y ./ L) .* sin.(pi * z ./ L) 

    else
        return 0 .* x
    end
end

function exactv_xy(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ β1 .* (pi/L)^2 .* cos.(pi * x ./ L) .* cos.(pi * y ./ L) .* cos.(pi * z ./ L) 

    else
        return 0 .* x
    end
end


#############################################


function exactw(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        wexact = z .* x .^ 3 .+ β2 .* sin.(pi * x ./ L) .+ β2 .* cos.(pi * y ./ L) .+ β2 .* cos.(pi * z ./ L) 
    else
        wexact = 0 .* x
    end
    return wexact[:]
end

function exactw_x(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return z .* 3 .* x .^ 2 .+ β2 .* (pi/L) * cos.(pi * x ./ L) 
    else
        return 0 .* x
    end
end
  
function exactw_xx(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return z .* 6 .* x .+ -β2 .* (pi/L)^2 * sin.(pi * x ./ L)
    else
        return 0 .* x
    end
end

function exactw_y(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ -β2 .* (pi/L) * sin.(pi * y ./ L) 
    else
        return 0 .* x
    end
end
  
function exactw_yy(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .+ -β2 .* (pi/L)^2 * cos.(pi * y ./ L) 
    else
        return 0 .* x
    end
end

function exactw_z(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return x .^ 3 .- β2 .* (pi/L) .* sin.(pi * z ./ L) 
    else
        return 0 .* x
    end
end
  
function exactw_zz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x .- β2 .* (pi/L)^2 .* cos.(pi * z ./ L) 
    else
        return 0 .* x
    end       
end

function exactw_xz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 3 .* x .^ 2
    else
        return 0 .* x
    end
end

function exactw_yz(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x
    else
        return 0 .* x
    end
end

function exactw_xy(x, y, z, e, EToDomain)
    if EToDomain[e] == 1
        return 0 .* x
    else
        return 0 .* x
    end
end



function exact_mu(x, y, z)
    # Define shear modulus μ
    μ = 4 .+ c .* sin.(b2 .* x ./ L .+ b1 .* y ./ L .+ b1 .* z ./ L)
    return μ
end

function exact_mu_x(x, y, z)
    # Define shear modulus derivative μ_x
    μ_x = b2 .* c ./ L .* cos.(b2 .* x ./ L .+ b1 .* y ./ L .+ b1 .* z ./ L)
    return μ_x
end

function exact_mu_y(x, y, z)
    # Define shear modulus derivative μ_z
    μ_z = b1 .* c ./ L .* cos.(b2 .* x ./ L .+ b1 .* y ./ L .+ b1 .* z ./ L)
    return μ_z
end


function exact_mu_z(x, y, z)
    # Define shear modulus derivative μ_z
    μ_z = b1 .* c ./ L .* cos.(b2 .* x ./ L .+ b1 .* y ./ L .+ b1 .* z ./ L)
    return μ_z
end



function exact_lambda(x, y, z)
    lam = 1 .+ c1 .* cos.(b2 .* x ./ L .+ b1 .* y ./ L .+ b1 .* z ./ L)
    return lam
end

function exact_lambda_x(x, y, z)
    lam_x = -c1 .* (b2 ./ L) .* sin.(b2 .* x ./ L .+ b1 .* y ./ L .+ b1 .* z ./ L)
    return lam_x
end

function exact_lambda_y(x, y, z)
    lam_z = -c1 .* (b1 ./ L) .* sin.(b2 .* x ./ L  .+ b1 .* y ./ L .+ b1 .* z ./ L)

    return lam_z
end

function exact_lambda_z(x, y, z)
    lam_z = -c1 .* (b1 ./ L) .* sin.(b2 .* x ./ L  .+ b1 .* y ./ L .+ b1 .* z ./ L)

    return lam_z
end

function exact_nu(x, y, z) # nu should range strictly between 0 and 1/2 
    # Define Poisson ratio
    nu = exact_lambda(x, y, z) ./ (2 .* exact_lambda(x, y, z) .+ 2 .* exact_mu(x, y, z))
    return nu
end



function source(x, y, z, e, EToDomain) # this should operate in local orientation
    c = exact_lambda(x, y, z) .+ 2 .* exact_mu(x, y, z)
    cx = exact_lambda_x(x, y, z) .+ 2 .* exact_mu_x(x, y, z)
    cy = exact_lambda_y(x, y, z) .+ 2 .* exact_mu_y(x, y, z)
    cz = exact_lambda_z(x, y, z) .+ 2 .* exact_mu_z(x, y, z)

    lam = exact_lambda(x, y, z)
    lamx = exact_lambda_x(x, y, z)
    lamy = exact_lambda_y(x, y, z)
    lamz = exact_lambda_z(x, y, z)

    mu = exact_mu(x, y, z)
    mux = exact_mu_x(x, y, z)
    muy = exact_mu_y(x, y, z)
    muz = exact_mu_z(x, y, z)
   
    ux = exactu_x(x, y, z, e, EToDomain)
    uy = exactu_y(x, y, z, e, EToDomain)
    uz = exactu_z(x, y, z, e, EToDomain)
    uxx = exactu_xx(x, y, z, e, EToDomain)
    uxy = exactu_xy(x, y, z, e, EToDomain)
    uxz = exactu_xz(x, y, z, e, EToDomain)
    uyy = exactu_yy(x, y, z, e, EToDomain)
    uyz = exactu_yz(x, y, z, e, EToDomain)
    uzz = exactu_zz(x, y, z, e, EToDomain)

    vx = exactv_x(x, y, z, e, EToDomain)
    vy = exactv_y(x, y, z, e, EToDomain)
    vz = exactv_z(x, y, z, e, EToDomain)
    vxx = exactv_xx(x, y, z, e, EToDomain)
    vxy = exactv_xy(x, y, z, e, EToDomain)
    vxz = exactv_xz(x, y, z, e, EToDomain)
    vyy = exactv_yy(x, y, z, e, EToDomain)
    vyz = exactv_yz(x, y, z, e, EToDomain)
    vzz = exactv_zz(x, y, z, e, EToDomain)

    wx = exactw_x(x, y, z, e, EToDomain)
    wy = exactw_y(x, y, z, e, EToDomain)
    wz = exactw_z(x, y, z, e, EToDomain)
    wxx = exactw_xx(x, y, z, e, EToDomain)
    wxy = exactw_xy(x, y, z, e, EToDomain)
    wxz = exactw_xz(x, y, z, e, EToDomain)
    wyy = exactw_yy(x, y, z, e, EToDomain)
    wyz = exactw_yz(x, y, z, e, EToDomain)
    wzz = exactw_zz(x, y, z, e, EToDomain)

    exx = ux 
    exy = 0.5*(uy + vx)
    exz = 0.5*(uz + wx)
    eyy = vy 
    eyz = 0.5(vz + wy)
    ezz = wz 
    ekk = exx + eyy + ezz

    exx_x = uxx 
    exx_y = uxy
    exx_z = uxz 

    eyy_x = vxy 
    eyy_y = vyy
    eyy_z = vyz
    
    ezz_x = wxz
    ezz_y = wyz 
    ezz_z = wzz 
 

    exy_x = 0.5*(uxy + vxx)
    exy_y = 0.5*(uyy + vxy)
    exz_x = 0.5*(uxz + wxx)
    exz_z = 0.5*(uzz + wxz)

    eyz_y = 0.5(vyz + wyy)
    eyz_z = 0.5(vzz + wyz)

   
    ekk_x = exx_x + eyy_x + ezz_x
    ekk_y = exx_y + eyy_y + ezz_y
    ekk_z = exx_z + eyy_z + ezz_z

    f1 = -(lamx .* ekk .+ lam .* ekk_x .+ 2 .* mux .* exx .+ 2 .* mu .* exx_x .+ 2 .* muy .* exy .+ 2 .* mu .* exy_y .+ 2 .* muz .* exz .+ 2 .* mu .* exz_z)
    f2 = -(2 .* mux .* exy .+ 2 .* mu .* exy_x .+ lamy .* ekk .+ lam .* ekk_y .+ 2 .* muy .* eyy .+ 2 .* mu .* eyy_y .+ 2 .* muz .* eyz .+ 2 .* mu .* eyz_z)
    f3 = -(2 .* mux .* exz .+ 2 .* mu .* exz_x .+ 2 .* muy .* eyz .+ 2 .* mu .* eyz_y .+ lamz .* ekk .+ lam .* ekk_z .+ 2 .* muz .* ezz .+ 2 .* mu .* ezz_z)

    # f1 = -((lam .+ 2 .* mu) .* uxx .+ (lam .+ mu) .* (vxy .+ wxz) .+ mu .* (uyy .+ uzz))
    # f2 = -((lam .+ 2 .* mu) .* vyy .+ (lam .+ mu) .* (uxy .+ wyz) .+ mu .* (vxx .+ vzz))
    # f3 = -((lam .+ 2 .* mu) .* wzz .+ (lam .+ mu) .* (uxz .+ vyz) .+ mu .* (wxx .+ wyy))
    
    f = [f1[:] f2[:] f3[:]]
    return f
end

function exact(x, y, z, e, EToDomain)
    e1 = exactu(x, y, z, e, EToDomain)
    e2 = exactv(x, y, z, e, EToDomain)
    e3 = exactw(x, y, z, e, EToDomain)

    return [e1 e2 e3]
end

function exact_traction(x, y, z, nx, ny, nz, e, EToDomain) #operates in local orientaion


    c = exact_lambda(x, y, z) .+ 2 * exact_mu(x, y, z)
    mu = exact_mu(x, y, z)
    lam = exact_lambda(x, y, z)

    sig_xx =   c .* exactu_x(x, y, z, e, EToDomain) .+ lam .* exactv_y(x, y, z, e, EToDomain) .+ lam .*  exactw_z(x, y, z, e, EToDomain) 
    sig_yy = lam .* exactu_x(x, y, z, e, EToDomain) .+   c .* exactv_y(x, y, z, e, EToDomain) .+ lam .*  exactw_z(x, y, z, e, EToDomain) 
    sig_zz = lam .* exactu_x(x, y, z, e, EToDomain) .+ lam .* exactv_y(x, y, z, e, EToDomain) .+   c .*  exactw_z(x, y, z, e, EToDomain) 


    sig_xy = mu .* exactu_y(x, y, z, e, EToDomain) .+ mu .*  exactv_x(x, y, z, e, EToDomain) 
    sig_xz = mu .* exactu_z(x, y, z, e, EToDomain) .+ mu .*  exactw_x(x, y, z, e, EToDomain) 
    sig_yz = mu .* exactv_z(x, y, z, e, EToDomain) .+ mu .*  exactw_y(x, y, z, e, EToDomain)  
        


    T1 = nx .* sig_xx[:] .+ ny .* sig_xy[:] .+ nz .* sig_xz[:]
    T2 = nx .* sig_xy[:] .+ ny .* sig_yy[:] .+ nz .* sig_yz[:]
    T3 = nx .* sig_xz[:] .+ ny .* sig_yz[:] .+ nz .* sig_zz[:]


    T = [T1 T2 T3]

    return T 

end

export exact_mu, exact_lambda, exact_nu, exact_mu_x, exact_mu_y, exact_mu_z, exact_lambda_x, exact_lambda_y, exact_lambda_z
export exactu, exactu_x, exactu_y, exactu_z, exactu_xx, exactu_xy, exactu_xz, exactu_yy, exactu_yz, exactu_zz
export exactv, exactv_x, exactv_y, exactv_z, exactv_xx, exactv_xy, exactv_xz, exactv_yy, exactv_yz, exactv_zz
export exactw, exactw_x, exactw_y, exactw_z, exactw_xx, exactw_xy, exactw_xz, exactw_yy, exactw_yz, exactw_zz
export source, exact, exact_traction