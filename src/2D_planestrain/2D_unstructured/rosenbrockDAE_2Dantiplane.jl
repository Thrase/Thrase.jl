using Thrase
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqCallbacks
using DelimitedFiles
using Dates
using Plots
using SparseArrays


import NonlinearSolve as NLS

#using UnicodePlots
const year_seconds = 31556926

include("ops_BP3-QD_unstructured.jl")
include("odefun_BP9-QD_unstructured.jl")
include("../utils_2D.jl")

sim_years = 20

p = 2
N = 40

L = 40
h = L/N

(D, Hinv, H, r) = diagonal_sbp_D1(p, N; xc = (0,L))
(D2, S0, SN, _, _, _) = diagonal_sbp_D2(p, N; xc = (0,L))

myI = sparse(1.0 * I, N+1, N+1)
D2x = SparseArrays.kron(myI, D2)
D2y = SparseArrays.kron(D2, myI)
D2 = D2x + D2y

S0x = SparseArrays.kron(myI, S0)
SNx = SparseArrays.kron(myI, SN)
Sx = S0x + SNx
S0y = SparseArrays.kron(S0, myI)
SNy = SparseArrays.kron(SN, myI)
Sy = S0y + SNy

Hx = SparseArrays.kron(myI, H)
Hy = SparseArrays.kron(H, myI)
Hinvx = SparseArrays.kron(myI, Hinv)
Hinvy = SparseArrays.kron(Hinv, myI)


x = collect(0:h:L)
y = collect(-L:h:0)

mu = 32

α = -13*mu/h

E0 = spzeros(N+1, N+1)
E0[1, 1] = 1
EN = spzeros(N+1, N+1)
EN[end, end] = 1
e0 = spzeros(N+1)
e0[1] = 1
eN = spzeros(N+1)
eN[N+1] = 1

# restriction ops (lifts are their transpose)
EB = SparseArrays.kron(e0', myI)
ET = SparseArrays.kron(eN', myI)
EL = SparseArrays.kron(myI, e0')
ER = SparseArrays.kron(myI, eN')

b0 = mu*EL*S0x

A = mu*D2 + α*Hinvy*EL'*EL + α*Hinvy*ER'*ER + Hinvx*EB'*EB*Sy - Hinvx*ET'*ET*Sy
A = lu(A)

onebyN = spzeros(N+1, (N+1)^2)
NbyN = spzeros((N+1)^2, (N+1)^2)
Nbyone = spzeros((N+1)^2, N+1)
spz0 = spzeros(N+1, N+1)

c = 3.4
η = mu/c
σn_min = 50
σn_max = 50
σn = σn_min .-(σn_max - σn_min) .* y ./ L
b = .015

RSamin = 0.01
RSamax = 0.025
RSH1 = 15
RSH2 = 18
a = zeros(N+1)
for i = 1:N+1
    a[i] = RSamin - (RSamin - RSamax) *
                min(1, max(0, (RSH1 + y[i])/(RSH1 - RSH2)))
end





Dc = .032 
f0 = 0.6
V0 = 1e-6
Vp = 1e-9
n = 1

M = [NbyN Nbyone Nbyone Nbyone;
    onebyN I spz0 spz0;
    onebyN spz0 spz0 spz0;
    onebyN spz0 spz0 I]
 
# initial conditions:
u_0 = zeros(N+1, N+1)
δ_0 = 0 .* y
τ0 = σn .* RSamax .* asinh.((0.5) .* exp((f0 .+ b*log(1)) ./ RSamax)) .+ η .* V0
psi_0 = a .* log.(2 * sinh.((τ0 .- η*V0) ./ (σn .* a)))
V_0 = V0 .+ 0 .* y
Y_0 = u_0[:]
Y_0 = vcat(Y_0, δ_0)
Y_0 = vcat(Y_0, V_0[:])
Y_0 = vcat(Y_0, psi_0[:])


params = (A = A, L = L, x = x, N = N, α = α, EL = EL, ER = ER, Hinvx = Hinvx, Hinvy = Hinvy, D2 = D2, f0 = f0, V0 = V0, Vp = Vp, c = c, mu = mu, η = η, σn = σn, b = b, a = a, Dc = Dc, n = n, b0 = b0)

function rhs_britt(dY, Y, p, t)

    N = p.N 

    m = (N+1)^2
    u = Y[1:m]
    δ = Y[m+1:m+N+1]
    V = Y[m+N+2:m+2N+2]
    psi = Y[m+2N+3:m+3N+3]


    dY[1:m] = p.mu*p.D2*u .+ p.α*p.Hinvy*p.EL'*(p.EL*u-δ/2) + p.α*p.Hinvy*p.ER'*(ER*u - p.Vp*ones(N+1)*t/2) + Hinvx*EB'*EB*Sy*u - Hinvx*ET'*ET*Sy*u
    dY[m+1:m+N+1] = V
    dY[m+N+2:m+2N+2] = p.n .* (p.b0*u) .+ τ0 .- p.η .* V .- p.σn .* p.a .* asinh.((V ./ (2*p.V0)) .* exp.(psi ./ p.a))
    dY[m+2N+3:m+3N+3] = (p.b*p.V0/p.Dc) .* (exp.((p.f0 .- psi) ./ p.b) .- abs.(V) ./ p.V0)

 
    nothing
end


function my_ode!(dw, W, p, t)

    N = p.N 
    δ = W[1:N+1]
    psi = W[N+2:2N+2]
 

    b = p.α*p.Hinvy*p.EL'*(δ/2) + p.α*p.Hinvy*p.ER'*(p.Vp*ones(N+1)*t/2)
    u = p.A \ b[:]
    tau = p.n*p.b0*u

    V = zeros(N+1)
    for i = 1:N+1
        f(v, params) = tau[i] .+ τ0[i] .- p.η*v .- p.σn[i] .* p.a[i] .* asinh.((v/(2*p.V0)) * exp.((psi[i]) ./ p.a[i]))
        prob = NLS.NonlinearProblem(f, 0)
        Vsol = NLS.solve(prob, NLS.NewtonRaphson())
        V[i] = Vsol.u
    end

    dw[1:N+1] = V
    dw[N+2:2N+2] = (p.b .* p.V0/p.Dc) .* (exp.((p.f0 .- psi)/p.b) .- abs.(V) ./ p.V0)

  
    nothing
end

# Bad? Normal way:
tspan = (0, sim_years * year_seconds)
w0 = [δ_0;psi_0]
prob = ODEProblem(my_ode!, w0, tspan, params)
#@time sol = solve(prob,  Rodas5P(autodiff = AutoFiniteDiff()), reltol = 1e-6, abstol = 1e-6)
@time sol = solve(prob,  Tsit5(), reltol = 1e-6, abstol = 1e-6)
Plots.plot(sol, vars = (N+1), tspan = tspan)

#  Better? Rosenbrock DAE:
# tspan = (0, sim_years * year_seconds)
# f = ODEFunction(rhs_britt, mass_matrix = M)
# prob_mm = ODEProblem(f, Y_0, tspan, params)
# @time sol = solve(prob_mm, Rodas5P(), reltol = 1e-6, abstol = 1e-6)
# Plots.plot(sol, vars = (N+1)^2, tspan = tspan)
