## non dimensional dynamics in polar coordinates
ν̇(u, p, t) = thrust(u, p, t) * cos(thrust_angle(u, p, t)) - CD(u, p, t) * u[1]^2 - sin(u[2])
γ̇(u, p, t) = thrust(u, p, t) * sin(thrust_angle(u, p, t)) / u[1] + CL(u, p, t) * u[1] - cos(u[2]) / u[1]
χ̇(u, p, t) = u[1] * cos(u[2])
ζ̇(u, p, t) = u[1] * sin(u[2])

dynamics(u, p, t) = [ν̇(u, p, t), γ̇(u, p, t)]
full_dynamics(u, p, t) = [ν̇(u, p, t), γ̇(u, p, t), χ̇(u, p, t), ζ̇(u, p, t)]
dynamics!(du, u, p, t) = begin
    du[1] = ν̇(u, p, t)
    du[2] = γ̇(u, p, t)
    return nothing
end
full_dynamics!(du, u, p, t) = begin
    du[1] = ν̇(u, p, t)
    du[2] = γ̇(u, p, t)
    du[3] = χ̇(u, p, t)
    du[4] = ζ̇(u, p, t)
    return nothing
end
projectile_dynamics(u, p, t) = [-sin(u[2]), -cos(u[2]) / u[1]]
projectile_dynamics!(du, u, p, t) = [-sin(u[2]), -cos(u[2]) / u[1]]

## cartesian components of non-dimensional dynamics
ν̇x(u, p, t) = thrust(u, p, t) * cos(u[2] + thrust_angle(u, p, t)) - CD(u, p, t) * u[1]^2 * cos(u[2]) - CL(u, p, t) * u[1]^2 * sin(u[2])
ν̇y(u, p, t) = thrust(u, p, t) * sin(u[2] + thrust_angle(u, p, t)) - CD(u, p, t) * u[1]^2 * sin(u[2]) + CL(u, p, t) * u[1]^2 * cos(u[2]) - 1


## calculated, non-dimensional quantities
kinetic_energy(u, p, t) = 0.5 * u[1]^2
potential_energy(u, p, t) = u[4]
total_energy(u, p, t) = kinetic_energy(u, p, t) + potential_energy(u, p, t)

aerodynamic_power(u, p, t) = -u[1]^3 * CD(u, p, t)
gravitational_power(u, p, t) = -sin(u[2]) * u[1]
thrust_power(u, p, t) = u[1] * thrust(u, p, t) * cos(thrust_angle(u, p, t))