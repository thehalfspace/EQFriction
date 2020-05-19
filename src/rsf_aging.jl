#######################################################
#   RATE- AND STATE- DEPENDENT FRICTION WITH AGING LAW
#######################################################

# Regularized rate-state-dependent friction
#   τ = σ_n * a * asinh[(V/2*Vo)exp(ψ/a)]
#   dψ_dt = (b*Vo/Lc) * [exp((fo - ψ)/b) - (V/Vo)]

using NLsolve
using LinearAlgebra

# Global constant friction parameters
const a = 0.019
const b = 0.015
const Vo = 1.0e-6
const fo = 0.6
const Lc = 0.008

function aging_law(V, θ_o, Δt)
    Lc*(θ_o + Δt)/(Lc + V*Δt)
end

function fault_strength(θ, τ_o, σ_n)
    # returns V
    2*Vo*sinh(τ_o/(σ_n*a)) * exp(-(fo + b*log(Vo*θ/Lc))/a)
end

function fric_solve()
    # Initial Conditions
    τ_o = repeat([22.5e6], 1000)
    σ_n = repeat([50e6], 1000)
    θ_o = repeat([1.0], 1000)
    V_o = repeat([1.0e-3], 1000)
    Δt = 1.0

    # Solution
    θ = zeros(1000)
    V = zeros(1000)
    result = 0.

    # at each point on fault
    for i in 1:length(τ_o)
        function f!(F,x)
            F[1] = aging_law(x[2], θ_o[i], Δt)              # θ
            F[2] = fault_strength(x[1], τ_o[i], σ_n[i])     # V
        end
        θ[i], V[i] = nlsolve(f!, [θ_o[i], V_o[i]]).zero
    end
    
    θ, V
end








