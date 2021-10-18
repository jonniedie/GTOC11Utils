function solve_with(vars, p; kwargs...)
	@unpack λ, t = vars
	u0 = ComponentArray([p.chaser; p.target; λ], getaxes(p.prob.u0))
	prob = remake(p.prob; u0, tspan=(zero(t), t), p)
	return solve(prob; kwargs...)
end

function loss(vars, p)
	sol = solve_with(vars, p)
	@unpack x, λ = sol[end]
	@unpack chaser, target = x
	return sum(abs2, SVector{6}(chaser) - SVector{6}(target))
end

"""
    low_thrust_transfer(chaser, target; kwargs...)

Calculate low thrust transfer between a chaser with initial conditions `chaser` and a target with
initial conditions `target`. Initial conditions should be given as `[rx, ry, rz, vx, vy, vz]`

## Keyword arguments
| Name | Default | Description |
|:---- |:-------|:----------- |
| `λ` | `-1e12rand(6).-5e11` | Initial guess for costate initial conditions
| `λ_lb` | `fill(-Inf,6)` | Lower bound for costate initial condition
| `λ_ub` | `fill(Inf,6)` | Upper bound for costate initial condition
| `t` | `8.0` | Initial guess for stop time
| `t_lb` | `0.0` | Lower bound for stop time
| `t_ub` | `15.0` | Upper bound for stop time
| `Γ` | `ustrip(Γ)` | Acceleration magnitude (in AU/yr)
| `μ` | `ustrip(μ)` | Sun's gravitational constant (in AU³/yr²)
| `alg`| `Fminbox(NelderMead())` | Optimization algorithm (from Optim.jl or similar)
| `autodiff` | `GalacticOptim.AutoFiniteDiff()` | Autodiff method (from GalacticOptim.jl)
"""
function low_thrust_transfer(chaser, target;
                            λ=-1e12rand(6).-5e11, λ_lb=fill(-Inf,6), λ_ub=fill(Inf,6),
                            t=8.0, t_lb=0.0, t_ub=15.0,
                            Γ=ustrip(Γ), μ=ustrip(μ),
                            alg=Fminbox(NelderMead()), autodiff=GalacticOptim.AutoFiniteDiff(),
                            opt_kwargs...
                            )

	x = ComponentArray(OptInput(; λ, t))
	lb = ComponentArray(OptInput(; λ=λ_lb, t=t_lb))
	ub = ComponentArray(OptInput(; λ=λ_ub, t=t_ub))

    p = (; Γ, μ, chaser, target, prob=sim_prob)
    f = OptimizationFunction(loss, autodiff)
    prob = OptimizationProblem(f, x, p; lb, ub)

    return solve(prob, alg; opt_kwargs...)
end
