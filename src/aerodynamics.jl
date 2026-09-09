CL(u, p::Parameters, t) = CL(alpha(u, p, t), p.aerodynamics)
CD(u, p::Parameters, t) = CD(alpha(u, p, t), p.aerodynamics)

∂CL_∂α(u, p::Parameters, t) = ∂CL_∂α(alpha(u, p, t), p.aerodynamics)
∂CD_∂α(u, p::Parameters, t) = ∂CD_∂α(alpha(u, p, t), p.aerodynamics)
max_CL(p::Parameters) = max_CL(p.aerodynamics)
max_CD(p::Parameters) = max_CD(p.aerodynamics)
min_CD(p::Parameters) = min_CD(p.aerodynamics)
alpha_CL_max(p::Parameters) = alpha_CL_max(p.aerodynamics)
alpha_CD_max(p::Parameters) = alpha_CD_max(p.aerodynamics)
alpha_CL_CD_max(p::Parameters) = alpha_CL_CD_max(p.aerodynamics)
inverse_CL(p::Parameters) = inverse_CL(p.aerodynamics)
inverse_CD(p::Parameters) = inverse_CD(p.aerodynamics)


"""
    TrigAerodynamics <: Aerodynamics

CL and CD are trigonometric functions of α:

    CL = CLmax sin(2α)
    CD = CD₀ + CDmax sin(α)²
"""
struct TrigAerodynamics <: Aerodynamics
    CLmax::Float64
    CD₀::Float64
    CDmax::Float64
end
TrigAerodynamics() = TrigAerodynamics(1.5, 0.03, 2.0)

CL(α, p::TrigAerodynamics) = p.CLmax * sin(2α)
CD(α, p::TrigAerodynamics) = p.CD₀ + p.CDmax * sin(α)^2

∂CL_∂α(α, p::TrigAerodynamics) = 2p.CLmax * cos(2α)
∂CD_∂α(α, p::TrigAerodynamics) = p.CDmax * sin(2α)
max_CL(p::TrigAerodynamics) = p.CLmax
min_CD(p::TrigAerodynamics) = p.CD₀
max_CD(p::TrigAerodynamics) = p.CD₀ + p.CDmax
alpha_CL_max(::TrigAerodynamics) = π / 4
alpha_CD_max(::TrigAerodynamics) = π / 2
alpha_CL_CD_max(p::TrigAerodynamics) = atan(sqrt(p.CD₀ / (p.CD₀ + p.CDmax)))
inverse_CL(cl, p::TrigAerodynamics) = cl ≤ p.CLmax ? 0.5asin(cl / p.CLmax) : π / 2 - 0.5asin(cl / p.CLmax)
inverse_CD(cd, p::TrigAerodynamics) = asin(sqrt((cd - p.CD₀) / p.CDmax))

"""
    LinearAerodynamics <: Aerodynamics

Combines linear pre-stall aerodynamics with post-stall model.

Pre-stall:

    CL = CLα α              {α ≤ αstall}
    CD = CD₀ + k CLα² α²    {α ≤ αstall}

Post-stall model can be any other <:Aerodynamics.
"""
struct LinearAerodynamics <: Aerodynamics
    slope::Float64
    CD₀::Float64
    k::Float64
    αstall::Float64
    poststall::TrigAerodynamics

    LinearAerodynamics(slope, CD₀, k, αstall, CLmax, CDmax) = new(slope, CD₀, k, αstall, TrigAerodynamics(CLmax, CD₀, CDmax))
end
LinearAerodynamics(slope, k, αstall, post::TrigAerodynamics) = LinearAerodynamics(slope, post.CD₀, k, αstall, post.CLmax, post.CDmax)
LinearAerodynamics() = LinearAerodynamics(2π, 0.03, deg2rad(15), TrigAerodynamics(1.5, 0.03, 2.0))

CL(α, p::LinearAerodynamics) = α ≤ p.αstall ? p.slope * α : CL(α, p.poststall)
CD(α, p::LinearAerodynamics) = α ≤ p.αstall ? p.CD₀ + p.k * p.slope^2 * α^2 : CD(α, p.poststall)

∂CL_∂α(α, p::LinearAerodynamics) = α < p.αstall ? p.slope : ∂CL_∂α(α, p.poststall)
∂CD_∂α(α, p::LinearAerodynamics) = α < p.αstall ? 2p.k * p.slope^2 * α : ∂CD_∂α(α, p.poststall)
max_CL(p::LinearAerodynamics) = p.poststall.CLmax
max_CD(p::LinearAerodynamics) = p.CD₀ + p.poststall.CDmax
min_CD(p::LinearAerodynamics) = p.CD₀
alpha_CL_max(p::LinearAerodynamics) = p.slope * p.αstall > max_CL(p.poststall) ? p.αstall : alpha_CL_max(p.poststall)
alpha_CD_max(p::LinearAerodynamics) = p.k * p.slope^2 * p.αstall^2 > max_CD(p.poststall) ? p.αstall : alpha_CD_max(p.poststall)
# alpha_CL_CD_max(p::LinearAerodynamics) = p.k * p.slope^2 > p.CDmax ? sqrt(p.CD₀ / p.k * p.a^2) : atan(sqrt(p.CD₀ / (p.CD₀ + p.CDmax)))
inverse_CL(cl, p::LinearAerodynamics) = cl ≤ CL(p.αstall, p) ? cl / p.slope : inverse_CL(cl, p.poststall)
inverse_CD(cd, p::LinearAerodynamics) = cd ≤ CD(p.αstall, p) ? sqrt((cd - p.CD₀) / (p.k * p.slope^2)) : inverse_CD(cd, p.poststall)

"""
    SmoothPreStallAerodynamics <: Aerodynamics

Smooth representation incorporating pre-stall aerodynamics, adapted from Li et al., 2022 J. Fluid. Mech.

Ignores the rotational contribution to lift.

Only valid for -π/2 ≤ α ≤ π/2.
"""
struct SmoothPreStallAerodynamics <: Aerodynamics
    CL₁::Float64 # approx slope of pre-stall CL
    CL₂::Float64 # max post-stall CL
    CD₀::Float64 # skin friction drag
    CD₁::Float64 # approx slope of pre-stall CD
    CD₂::Float64 # max post stall CD
    αstall::Float64 # α₀ in Li et al.
    δ::Float64 # smoothness of transition between pre-stall and stall
end
SmoothPreStallAerodynamics() = SmoothPreStallAerodynamics(5.2, 0.95, 0.1, 5.0, 1.9, deg2rad(14), deg2rad(6)) # values from Li et al
SmoothPreStallAerodynamics(p::TrigAerodynamics) = SmoothPreStallAerodynamics(2π, p.CLmax, p.CD₀, 0.1, p.CDmax, deg2rad(14), deg2rad(6))
TrigAerodynamics(p::SmoothPreStallAerodynamics) = TrigAerodynamics(p.CL₂, p.CD₀, p.CD₂-p.CD₀)

f̃(α, αstall, δ) = 0.5 * (1 - tanh((sqrt(α^2 + 1e-16) - αstall) / δ)) # symmetric about α = 0, so α ∈ [-π/2, π/2]
C̃L(α, CL₁, CL₂, αstall, δ) = f̃(α, αstall, δ) * CL₁ * sin(α) + (1 - f̃(α, αstall, δ)) * CL₂ * sin(2α)
C̃D(α, CD₀, CD₁, CD₂, αstall, δ) = f̃(α, αstall, δ) * (CD₀ + CD₁ * sin(α)^2) + (1 - f̃(α, αstall, δ)) * (CD₂ * sin(α)^2)
∂f̃_∂α(α, αstall, δ) = -0.5α*sech((αstall - sqrt(α^2 + 1e-16)) / δ)^2 / (δ * sqrt(α^2 + 1e-16))

CL(α, p::SmoothPreStallAerodynamics) = C̃L(α, p.CL₁, p.CL₂, p.αstall, p.δ)
CD(α, p::SmoothPreStallAerodynamics) = C̃D(α, p.CD₀, p.CD₁, p.CD₂, p.αstall, p.δ)

∂CL_∂α(α, p::SmoothPreStallAerodynamics) =
    ∂f̃_∂α(α, p.αstall, p.δ) * (p.CL₁ * sin(α) - p.CL₂ * sin(2α)) + f̃(α, p.αstall, p.δ) * (p.CL₁ * cos(α) - 2p.CL₂ * cos(2α)) + 2p.CL₂ * cos(2α)
∂CD_∂α(α, p::SmoothPreStallAerodynamics) =
    ∂f̃_∂α(α, p.αstall, p.δ) * (p.CD₀ + (p.CD₁ - p.CD₂) * sin(α)^2) + f̃(α, p.αstall, p.δ) * sin(2α) * (p.CD₁ - p.CD₂) + p.CD₂ * sin(2α)
max_CL(p::SmoothPreStallAerodynamics) = max(CL(p.αstall, p), CL(π / 4, p))
max_CD(p::SmoothPreStallAerodynamics) = p.CD₂
min_CD(p::SmoothPreStallAerodynamics) = p.CD₀

