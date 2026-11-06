
function var_3D_D2q(p, Nqp, Nrp, Nsp, C, HIq; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    # Global index: (iq, ir, is) -> iq + (ir-1)*Nqp + (is-1)*Nqp*Nrp
    # For fixed (ir,is), D2 acts along iq: global = iq + offset, offset = (ir-1)*Nqp + (is-1)*Nqp*Nrp

    N = Nqp * Nrp * Nsp
    I_D2 = Int[]; J_D2 = Int[]; V_D2 = Float64[]
    I_S0 = Int[]; J_S0 = Int[]; V_S0 = Float64[]
    I_SN = Int[]; J_SN = Int[]; V_SN = Float64[]

    for is = 1:Nsp, ir = 1:Nrp
        B = C[:, ir, is]
        (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nqp-1, B; xc = (-1,1))
        offset = (ir-1)*Nqp + (is-1)*Nqp*Nrp

        rows, cols, vals = findnz(D2)
        append!(I_D2, rows .+ offset); append!(J_D2, cols .+ offset); append!(V_D2, vals)

        rows, cols, vals = findnz(S0)
        append!(I_S0, rows .+ offset); append!(J_S0, cols .+ offset); append!(V_S0, vals)

        rows, cols, vals = findnz(SN)
        append!(I_SN, rows .+ offset); append!(J_SN, cols .+ offset); append!(V_SN, vals)
    end

    D2q = sparse(I_D2, J_D2, V_D2, N, N)
    S0q = sparse(I_S0, J_S0, V_S0, N, N)
    SNq = sparse(I_SN, J_SN, V_SN, N, N)
    return D2q, S0q, SNq
end

function var_3D_D2r(p, Nqp, Nrp, Nsp, C, HIr; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    # Global index: (iq, ir, is) -> iq + (ir-1)*Nqp + (is-1)*Nqp*Nrp
    # For fixed (iq,is), D2 acts along ir: global = iq + (ir-1)*Nqp + (is-1)*Nqp*Nrp

    N = Nqp * Nrp * Nsp
    I_D2 = Int[]; J_D2 = Int[]; V_D2 = Float64[]
    I_S0 = Int[]; J_S0 = Int[]; V_S0 = Float64[]
    I_SN = Int[]; J_SN = Int[]; V_SN = Float64[]

    for is = 1:Nsp, iq = 1:Nqp
        B = C[iq, :, is]
        (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nrp-1, B; xc = (-1,1))

        rows, cols, vals = findnz(D2)
        grows = iq .+ (rows .- 1) .* Nqp .+ (is-1)*Nqp*Nrp
        gcols = iq .+ (cols .- 1) .* Nqp .+ (is-1)*Nqp*Nrp
        append!(I_D2, grows); append!(J_D2, gcols); append!(V_D2, vals)

        rows, cols, vals = findnz(S0)
        grows = iq .+ (rows .- 1) .* Nqp .+ (is-1)*Nqp*Nrp
        gcols = iq .+ (cols .- 1) .* Nqp .+ (is-1)*Nqp*Nrp
        append!(I_S0, grows); append!(J_S0, gcols); append!(V_S0, vals)

        rows, cols, vals = findnz(SN)
        grows = iq .+ (rows .- 1) .* Nqp .+ (is-1)*Nqp*Nrp
        gcols = iq .+ (cols .- 1) .* Nqp .+ (is-1)*Nqp*Nrp
        append!(I_SN, grows); append!(J_SN, gcols); append!(V_SN, vals)
    end

    D2r = sparse(I_D2, J_D2, V_D2, N, N)
    S0r = sparse(I_S0, J_S0, V_S0, N, N)
    SNr = sparse(I_SN, J_SN, V_SN, N, N)
    return D2r, S0r, SNr
end

function var_3D_D2s(p, Nqp, Nrp, Nsp, C, HIq; xc = (-1, 1))
    # C has not been diagonalized, e.g. send in C[1,1,1,1], which is size (Nqp, Nrp, Nsp)
    # Global index: (iq, ir, is) -> iq + (ir-1)*Nqp + (is-1)*Nqp*Nrp
    # For fixed (iq,ir), D2 acts along is: global = base + (is-1)*Nqp*Nrp, base = iq + (ir-1)*Nqp

    N = Nqp * Nrp * Nsp
    I_D2 = Int[]; J_D2 = Int[]; V_D2 = Float64[]
    I_S0 = Int[]; J_S0 = Int[]; V_S0 = Float64[]
    I_SN = Int[]; J_SN = Int[]; V_SN = Float64[]

    for ir = 1:Nrp, iq = 1:Nqp
        B = C[iq, ir, :]
        (D2, S0, SN, _, _, _, _) = variable_diagonal_sbp_D2(p, Nsp-1, B; xc = (-1,1))
        base = iq + (ir-1)*Nqp

        rows, cols, vals = findnz(D2)
        append!(I_D2, base .+ (rows .- 1) .* Nqp*Nrp)
        append!(J_D2, base .+ (cols .- 1) .* Nqp*Nrp)
        append!(V_D2, vals)

        rows, cols, vals = findnz(S0)
        append!(I_S0, base .+ (rows .- 1) .* Nqp*Nrp)
        append!(J_S0, base .+ (cols .- 1) .* Nqp*Nrp)
        append!(V_S0, vals)

        rows, cols, vals = findnz(SN)
        append!(I_SN, base .+ (rows .- 1) .* Nqp*Nrp)
        append!(J_SN, base .+ (cols .- 1) .* Nqp*Nrp)
        append!(V_SN, vals)
    end

    D2s = sparse(I_D2, J_D2, V_D2, N, N)
    S0s = sparse(I_S0, J_S0, V_S0, N, N)
    SNs = sparse(I_SN, J_SN, V_SN, N, N)
    return D2s, S0s, SNs
end

export var_3D_D2q, var_3D_D2r, var_3D_D2s
