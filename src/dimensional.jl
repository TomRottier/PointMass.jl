# struct to hold parameters for dimensional model
struct DimensionalModel{T}
    m::T
    S::T
    ρ::T
    g::T
end

# default parameters (a pigeon)
DimensionalModel(m, S, ρ=1.225, g=9.81) = DimensionalModel(m, S, ρ, g)
DimensionalModel() = DimensionalModel(0.3, 0.04, 1.225, 9.81)

nominal_speed(m, S, ρ, g) = √(2m * g / (ρ * S))
nominal_time(m, S, ρ, g) = nominal_speed(m, S, ρ, g) / g
nominal_length(m, S, ρ) = 2m / (ρ * S)
nominal_speed(mdl::DimensionalModel) = nominal_speed(mdl.m, mdl.S, mdl.ρ, mdl.g)
nominal_time(mdl::DimensionalModel) = nominal_time(mdl.m, mdl.S, mdl.ρ, mdl.g)
nominal_length(mdl::DimensionalModel) = nominal_length(mdl.m, mdl.S, mdl.ρ)
