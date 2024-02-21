# α
#######################################################
#   RATE- AND STATE- DEPENDENT FRICTION WITH AGING LAW
#######################################################

# Regularized rate-state-dependent friction
#   τ = σn * a * asinh[(V/2*Vo)exp(ψ/a)]
#   ψ = fo + b*log(V*θ/Lc)
   
# Aging law
#   dθ_dt = 1 - (Vθ/Lc)

using DifferentialEquations

struct rsf_params{T<:AbstractVector}
    a::T
    b::T
    Vo::T
    fo::T
    Lc::T
end

struct rsf_initial_conditions{T<:AbstractVector}
    τo::T
    σn::T
    θo::T
    Vo::T
end


function aging_law(V, θo, Lc, Δt)
    Lc*(θo + Δt)/(Lc + V*Δt)
end

function fault_strength(θ, τo, σn, Vo, a, b, fo, Lc)
    # returns V
    2.0*Vo*sinh(τo/(σn*a)) * exp(-(fo + b*log(Vo*θ/Lc))/a)
end

function friction_solver(params::rsf_params, ics::rsf_initial_conditions, Δt)
    p
end
