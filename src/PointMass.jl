module PointMass

export Aerodynamics, Control, Parameters
export ν̇, γ̇, χ̇, ζ̇, dynamics!, dynamics, full_dynamics!, full_dynamics
export kinetic_energy, potential_energy, total_energy, aerodynamic_power, gravitational_power, thrust_power
export TrigAerodynamics, LinearAerodynamics, SmoothPreStallAerodynamics
export CL, CD, inverse_CL, inverse_CD, alpha_CL_max, alpha_CD_max, alpha_CL_CD_max, ∂CL_∂α, ∂CD_∂α, min_CD, max_CD, max_CL
export alpha, thrust, thrust_angle, ConstantAlpha, LinearAlpha, InterpolatedControl
export DimensionalModel, nominal_time, nominal_speed, nominal_length


# abstract type for aerodynamics models
abstract type Aerodynamics end

# abstract type for all control models
abstract type Control end

# struct for parameters for simulation
struct Parameters{A<:Aerodynamics,C<:Control}
    aerodynamics::A
    control::C
end


include("aerodynamics.jl")
include("control.jl")
include("dynamics.jl")
include("dimensional.jl")
include("analysis.jl")

end # module PointMass
