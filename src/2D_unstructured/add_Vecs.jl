A = zeros(2,2)
B = A

B[1, 1] = 1
B[1, 2] = 2
B[2, 1] = 3
B[2, 2] = 4

x = 3

# Julia funcitons return the last expression or statement after the word return
# recommend using docstring with type contract
"""
    add_vecs(x, y)

Add two 1D arrays x and y together.

# Examples
'''julia-repl
julia> add_vecs([1, 2], [3, 4])
2-element Array{Int64,1}:
4
5
'''
"""
function add_vecs(x, y)
    return x + y
end

function fuse_mult_add(a, x, y)
    return a * x + y
end

function add_vecs!(z, x, y)
    N = length(x)
    for i = 1:N
        z[i] = x[i] + y[i]
    end
end
