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

sim_years = 1000

p = 2
N = 200

L = 20
h = L/N

(D, Hinv, H, r) = diagonal_sbp_D1(p, N; xc = (0,L))
(D2, S0, SN, _, _, _) = diagonal_sbp_D2(p, N; xc = (0,L))
x = collect(0:h:L)



mu = 32
b0 = mu*S0[1,:]

E0 = zeros(N+1, N+1)
E0[1, 1] = 1
EN = zeros(N+1, N+1)
EN[end, end] = 1


e0 = zeros(N+1)
e0[1] = 1

eN = zeros(N+1)
eN[N+1] = 1

onebyN = zeros(1, N+1)
NbyN = zeros(N+1, N+1)
Nbyone = zeros(N+1, 1)

c = 3.4
η = mu/c
σn = 50
b = .015
a = .01
α = -13*mu/h

#V = 1.123e-6
Dc = .008
f0 = 0.6
V0 = 1e-6
Vp = 1e-9
n = 1

M = [NbyN Nbyone Nbyone Nbyone;
    onebyN 1 0 0;
    onebyN 0 0 0;
    onebyN 0 0 1]
 
# initial conditions:
u_0 = ((η*V0+σn * a * asinh(0.5*exp(f0/a)))/mu)*x
u_L = ((η*V0+σn * a * asinh(0.5*exp(f0/a)))/mu)*L
δ_0 = 0
#theta_0 = Dc / V0
psi_0 = f0
V_0 = V0
Y_0 = u_0[:]
append!(Y_0, δ_0)
append!(Y_0, V_0)
append!(Y_0, psi_0)

A = mu*D2 + α*Hinv*E0 + α*Hinv*EN 
A = lu(A)



params = (A = A, L = L, x = x, u_L = u_L, N = N, α = α, e0 = e0, eN = eN, Hinv = Hinv, D2 = D2, f0 = f0, V0 = V0, Vp = Vp, c = c, mu = mu, η = η, σn = σn, b = b, a = a, Dc = Dc, n = n, b0 = b0)

function rhs_britt(dY, Y, p, t)

    N = p.N 
    u = Y[1:N+1]
    δ = Y[N+2]
    V = Y[N+3]
    psi = Y[N+4]
    # @show u
    # @show (δ, V, theta, t)

    dY[1:N+1] = p.mu*p.D2*u .+ p.α*p.Hinv*(u[1]-δ/2)*p.e0 + p.α*p.Hinv*(u[end] - p.u_L - p.Vp*t/2)*p.eN
    dY[N+2] = V
    # dY[N+3] = p.n*p.b0'*u - p.η*V - p.σn * (p.f0 + p.a * log(abs(V)/p.V0) + p.b * log(p.V0 * theta ./ p.Dc))
    dY[N+3] = p.n*p.b0'*u - p.η*V - p.σn * p.a * asinh((V/(2*p.V0)) * exp((psi)/p.a))
    dY[N+4] = (p.b*p.V0/p.Dc) * (exp((p.f0 - psi)/p.b) - abs(V)/p.V0)#1 - V*theta/p.Dc 

  
    nothing
end


function my_ode!(dw, W, p, t)

    N = p.N 
    δ = W[1]
    psi = W[2]
 
   #p.mu*p.D2*u .+ p.α*p.Hinv*(u[1]-δ/2)*p.e0 + p.α*p.Hinv*(u[end] - p.u_L - p.Vp*t/2)*p.eN
    
    b = p.α*p.Hinv*(δ/2)p.e0 + p.α*p.Hinv*(p.u_L + p.Vp*t/2)*p.eN 
    u = p.A \ b[:]
    #u = (p.u_L + p.Vp*t/2 - δ/2) * p.x ./ p.L
    tau = p.n*p.b0'*u

    f(V, params) = tau - p.η*V - p.σn * p.a * asinh((V/(2*p.V0)) * exp((psi)/p.a))
    prob = NLS.NonlinearProblem(f, 0)
    Vsol = NLS.solve(prob, NLS.NewtonRaphson())
    V = Vsol.u

    dw[1] = V
    dw[2] = (p.b*p.V0/p.Dc) * (exp((p.f0 - psi)/p.b) - abs(V)/p.V0)

  
    nothing
end

# Bad? Normal way:
# tspan = (0, sim_years * year_seconds)
# w0 = [δ_0;psi_0]
# prob = ODEProblem(my_ode!, w0, tspan, params)
# #@time sol = solve(prob,  Rodas5P(autodiff = AutoFiniteDiff()), reltol = 1e-6, abstol = 1e-6)
#  @time sol = solve(prob,  Tsit5(), reltol = 1e-6, abstol = 1e-6)

#  Plots.plot(sol, vars = (1), tspan = tspan)

#  Better? Rosenbrock DAE:
tspan = (0, sim_years * year_seconds)
f = ODEFunction(rhs_britt, mass_matrix = M)
prob_mm = ODEProblem(f, Y_0, tspan, params)
@time sol = solve(prob_mm, Rodas5P(), reltol = 1e-6, abstol = 1e-6)
Plots.plot(sol, vars = (N+2), tspan = tspan)
