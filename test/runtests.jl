using Test, PointMass, OrdinaryDiffEq

function simulate(f, u0, p, tspan; kwargs...)
    prob = ODEProblem(f, u0, tspan, p; kwargs...)
    return solve(prob, Tsit5(); saveat=0.01, abstol=1e-8, reltol=1e-8)
end

@testset "all" verbose = true begin
    # test parameters
    p = Parameters(
        TrigAerodynamics(1.5, 0.03, 2.0),
        ConstantAlpha(0.2)
    )
    tspan = (0.0, 5.0)

    # energy conservation
    @testset "energy conservation" verbose = true begin
        p = Parameters(
            TrigAerodynamics(1.5, 0.0, 0.0),
            ConstantAlpha(0.2)
        )
        for _ in 1:10
            u₀ = rand(4)
            sol = simulate(full_dynamics!, u₀, p, tspan)
            TE₀ = total_energy(sol(0.0), p, 0.0)
            TE₁ = total_energy(sol(sol.t[end]), p, sol.t[end])

            @test TE₀ ≈ TE₁
        end
    end

    # power conservation
    @testset "power conservation" verbose = true begin
        for _ in 1:10
            u₀ = rand(4)
            sol = simulate(dynamics!, u₀, p, tspan)
            net_power = map(sol.t) do t
                u = sol(t)
                dke_dt = u[1] * ν̇(u, p, t)
                aero_power = aerodynamic_power(u, p, t)
                grav_power = gravitational_power(u, p, t)

                return aero_power + grav_power - dke_dt
            end
            @test net_power ≈ fill(0.0, size(net_power)) atol = 1e-7
        end
    end

    # test negative velocity
    # @testset "negative velocity" verbose = true begin
    #     @test dynamics([1.0, 0.0], p, 0.0) ≈ dynamics([-1.0, π], p, 0.0)
    # end
end
