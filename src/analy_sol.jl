function u1_analy(x,y,z)
    a = [sin(π*i + π*j + π*k) for i ∈ x, j ∈ y, k ∈ z]
    return a
end


function u2_analy(x,y,z)
    a = [0 for i ∈ x, j ∈ y, k ∈ z]
    return a
end

function u3_analy(x,y,z)
    a = [0 for i ∈ x, j ∈ y, k ∈ z]
    return a
end


function u1_new(x,y,z)
    a = [sin(π*i + 2π*j + 3π*k) for i ∈ x, j ∈ y, k ∈ z]
    return a
end

# x_ex = 0:0.25:1
# y_ex = 0:0.25:1
# z_ex = 0:0.25:1

# ex_analy_sol = u(x_ex,y_ex,z_ex)


function form_analy_sol(;N = 2^3)
    x_ex = 0:1/N:1
    y_ex = 0:1/N:1
    z_ex = 0:1/N:1
    # ex_analy_sol = u(x_ex,y_ex,z_ex)
    # return (ex_analy_sol,x_ex,y_ex,z_ex)
    return (u1_analy(x_ex, y_ex, z_ex), u2_analy(x_ex, y_ex, z_ex), u3_analy(x_ex, y_ex, z_ex), x_ex,y_ex,z_ex)
end


# Dirichlet conditions
# u1
function u1_Front(y,z)
    a = [sin(π + π*j + π*k) for j ∈ y, k ∈ z]
    return a
end

function u1_End(y,z)
    a = [sin(π*j + π*k) for j ∈ y, k ∈ z]
    return a
end

function u1_Left(x,z)
    a = [sin(π*i + π*k) for i ∈ x, k ∈ z]
    return a
end

function u1_Right(x,z)
    a = [sin(π*i + π*k) for i ∈ x, k ∈ z]
    return a
end

function u1_Top(x,y)
    a = [sin(π*i + π*j + π) for i ∈ x, j ∈ y]
    return a
end


function u1_Bottom(x,y)
    a = [sin(π*i + π*j) for i ∈ x, j ∈ y]
    return a
end

# u2
function u2_Front(y,z)
    a = [0 for j ∈ y, k ∈ z]
    return a
end

function u2_End(y,z)
    a = [0 for j ∈ y, k ∈ z]
    return a
end

function u2_Left(x,z)
    a = [0 for i ∈ x, k ∈ z]
    return a
end

function u2_Right(x,z)
    a = [0 for i ∈ x, k ∈ z]
    return a
end

function u2_Top(x,y)
    a = [0 for i ∈ x, j ∈ y]
    return a
end


function u2_Bottom(x,y)
    a = [0 for i ∈ x, j ∈ y]
    return a
end


# Neumann conditions

# u1
function u1_x_Front(y,z)
    a = π * [cos(π + π*j + π*k) for j ∈ y, k ∈ z]
    return a 
end


function u1_x_End(y,z)
    a = -π * [cos(π*j + π*k) for j ∈ y, k ∈ z] # -1 normal direction
    return a 
end

function u1_y_End(y,z)
    a = -π * [cos(π*j + π*k) for j ∈ y, k ∈ z] # -1 normal direction
    return a 
end

function u1_y_Left(x,z)
    a = -π * [cos(π*i + π*k) for i ∈ x, k ∈ z] # -1 normal direction
    return a 
end

function u1_x_Left(x,z)
    a = -π * [cos(π*i + π*k) for i ∈ x, k ∈ z] # -1 normal direction
end

function u1_y_Right(x,z)
    a = π * [cos(π*i + π + π*k) for i ∈ x, k ∈ z] # +1 normal direction
    return a 
end

function u1_x_Right(x,z)
    a = π * [cos(π*i + π + π*k) for i ∈ x, k ∈ z] # +1 normal direction
    return a 
end

function u1_z_Top(x,y)
    a = π * [cos(π*i + π*j + π) for i ∈ x, j ∈ y] # +1 normal direction
    return a 
end

function u1_x_Top(x,y)
    a = π * [cos(π*i + π*j + π) for i ∈ x, j ∈ y] # +1 normal direction
    return a 
end


function u1_z_Bottom(x,y)
    a = -π * [cos(π*i + π*j) for i ∈ x, j ∈ y] # -1 normal direction
    return a 
end

function u1_x_Bottom(x,y)
    a = -π * [cos(π*i + π*j) for i ∈ x, j ∈ y] # -1 normal direction
    return a 
end

# u2

function u2_x_Front(y,z)
    a = π * [0 for j ∈ y, k ∈ z]
    return a 
end


function u2_x_End(y,z)
    a = -π * [0 for j ∈ y, k ∈ z] # -1 normal direction
    return a 
end

function u2_y_End(y,z)
    a = -π * [0 for j ∈ y, k ∈ z] # -1 normal direction
    return a 
end

function u2_y_Left(x,z)
    a = -π * [0 for i ∈ x, k ∈ z] # -1 normal direction
    return a 
end

function u2_x_Left(x,z)
    a = -π * [0 for i ∈ x, k ∈ z] # -1 normal direction
end

function u2_y_Right(x,z)
    a = π * [0 for i ∈ x, k ∈ z] # +1 normal direction
    return a 
end

function u2_x_Right(x,z)
    a = π * [0 for i ∈ x, k ∈ z] # +1 normal direction
    return a 
end

function u2_z_Top(x,y)
    a = π * [0 for i ∈ x, j ∈ y] # +1 normal direction
    return a 
end

function u2_y_Top(x,y)
    a = π * [0 for i ∈ x, j ∈ y] # +1 normal direction
    return a 
end


function u2_z_Bottom(x,y)
    a = -π * [0 for i ∈ x, j ∈ y] # -1 normal direction
    return a 
end

function u2_y_Bottom(x,y)
    a = -π * [0 for i ∈ x, j ∈ y] # -1 normal direction
    return a 
end