#######################################################
#   RATE- AND STATE- DEPENDENT FRICTION WITH AGING LAW
#######################################################

# Regularized rate-state-dependent friction
#   τ = σ_n * a * asinh[(V/2*Vo)exp(ψ/a)]
#   ψ = fo + b*log(V*θ/Lc)
   
# Aging law
#   dθ_dt = 1 - (Vθ/Lc)

using NLsolve
using LinearAlgebra

struct rsf_params{T<:AbstractVector}
    a::T
    b::T
    Vo::T
    fo::T
    Lc::T
end

struct rsf_initial_conditions{T<:AbstractVector}
    τ_o::T
    σ_n::T
    θ_o::T
    V_o::T
end

function aging_law(V, θ_o, Lc, Δt)
    Lc*(θ_o + Δt)/(Lc + V*Δt)
end

function fault_strength(θ, τ_o, σ_n, Vo, a, b, fo, Lc)
    # returns V
    2.0*Vo*sinh(τ_o/(σ_n*a)) * exp(-(fo + b*log(Vo*θ/Lc))/a)
end

function fric_solve(params::rsf_params, ics::rsf_initial_conditions, Δt)
    # Unpack friction parameters
    a = params.a
    b = params.b
    Vo = params.Vo
    fo = params.fo
    Lc = params.Lc

    # Unpack initial conditions
    τ_o = ics.τ_o
    σ_n = ics.σ_n
    θ_o = ics.θ_o
    V_o = ics.V_o

    # Number of points on the fault
    nFault = length(τ_o)

    # Solution
    θ::AbstractVector = zeros(nFault)
    V::AbstractVector = zeros(nFault)

    # at each point on fault
    for i in 1:nFault
        function f!(F,x)
            F[1] = aging_law(x[2], θ_o[i], Lc[i], Δt)           # θ
            F[2] = fault_strength(x[1], τ_o[i], σ_n[i], Vo[i],
                                  a[i], b[i], fo[i], Lc[i])     # V
        end
        θ[i], V[i] = nlsolve(f!, [θ_o[i], V_o[i]]).zero
    end
    
    θ, V
end

# Testing
function testing()
    nF = 100

    # Initialize friction parameters
    a = repeat([0.019], nF)
    b = repeat([0.015], nF)
    Vo = repeat([1.0e-6], nF)
    fo = repeat([0.6], nF)
    Lc = repeat([0.008], nF)
    params = rsf_params(a, b, Vo, fo, Lc)

    # Initial conditions
    τ_o = repeat([22.5e6], nF) 
    σ_n = repeat([50e6], nF) 
    θ_o = repeat([10.0], nF) 
    V_o = repeat([1.0e-3], nF) 
    ics = rsf_initial_conditions(τ_o, σ_n, θ_o, V_o)

    Δt = 0.01
    fric_solve(params, ics, Δt)
end
