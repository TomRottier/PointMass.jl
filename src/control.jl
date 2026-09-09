alpha(u, p::Parameters, t) = alpha(u, p.control, t)
thrust(u, p::Parameters, t) = thrust(u, p.control, t)
thrust_angle(u, p::Parameters, t) = thrust_angle(u, p.control, t)


## constant alpha, no thrust
struct ConstantAlpha <: Control
    α::Float64
end

alpha(u, p::ConstantAlpha, t) = p.α
thrust(u, p::ConstantAlpha, t) = 0.0
thrust_angle(u, p::ConstantAlpha, t) = 0.0

## linear alpha, no thurst
struct LinearAlpha <: Control
    slope::Float64
    intercept::Float64
end

alpha(u, p::LinearAlpha, t) = p.slope * t + p.intercept
thrust(u, p::LinearAlpha, t) = 0.0
thrust_angle(u, p::LinearAlpha, t) = 0.0

## pendulum alpha, no thrust
struct PendulumAlpha <: Control
    r::Float64
end
alpha(u, p::Parameters{<:Aerodynamics,PendulumAlpha}, t) = inverse_CL(1 / p.control.r + cos(u[2]) / u[1]^2, p)
thrust(u, p::PendulumAlpha, t) = 0.0
thrust_angle(u, p::PendulumAlpha, t) = 0.0

## interpolated control - expects some interpolation object that can be called at time t to return corresponding value
struct InterpolatedControl{T1,T2,T3} <: Control
    α_itp::T1
    CT_itp::T2
    β_itp::T3
end
InterpolatedControl(α_itp) = InterpolatedControl(α_itp, CT -> 0.0, β -> 0.0)
PointMass.alpha(u, p::InterpolatedControl, t) = p.α_itp(t)
PointMass.thrust(u, p::InterpolatedControl, t) = p.CT_itp(t)
PointMass.thrust_angle(u, p::InterpolatedControl, t) = p.β_itp(t)
