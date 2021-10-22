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
    dr = SVector{3}(chaser.r) - SVector{3}(target.r)
    dṙ = SVector{3}(chaser.ṙ) - SVector{3}(target.ṙ)
    return 5e6*1/2*dr'dr + 1e6*1/2*dṙ'dṙ
	# return sum(abs2, SVector{6}(chaser) - SVector{6}(target))
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

    p = (; Γ, μ, chaser, target, prob=opt_rel_prob)
    f = OptimizationFunction(loss, autodiff)
    prob = OptimizationProblem(f, x, p; lb, ub)

    return solve(prob, alg; opt_kwargs...)
end


## Everything above this is garbage

function nl_fun(u, p)
	@unpack λ, t = u
	@unpack u0, xf, prob, α, out = p

	u0 = copy(u0)
	u0.λ = λ
	prob = remake(prob; u0, tspan=(prob.tspan[1], t))
	sol = solve(prob)
	uf = sol[end]

    out.λ.r .= (xf.r .- uf.x.r).*AU .|> km .|> ustrip
    out.λ.ṙ .= (xf.ṙ .- uf.x.ṙ).*(AU/yr) .|> km/s .|> ustrip
    out.t = (-1 + uf.λ'uf.x) * α
    return out
end

get_candidate_solutions(station, asteroids::AbstractMatrix, args...; kwargs...) = get_candidate_solutions(station, collect(eachrow(asteroids)), args...; kwargs...)
function get_candidate_solutions(station, asteroids, back_time, args...; n_candidates=1, trans_scale=1e-8, autodiff=:forward, saveat=ustrip(yr(1d)), kwargs...)
    @assert back_time>0 "Second argument should be a positive number representing the time before current time. Got $back_time"

    ## Solve reverse problem
    # Choose the first five final costates at random and calculate last from Hamiltonian
    λf = @SVector(rand(5)) .- 0.5
    λf = [λf; -(1 + station[1:end-1]'λf)/station[end]]

    # Set up problem
    uf = ComponentArray(OneVehicleSimState(x=collect(station), λ=λf))
    t0 = 0.0
    back_prob = remake(opt_prob; u0=uf, tspan=(back_time, t0))

    # Solve
    back_sol = solve(back_prob)

    # Get initial state and costate
    u0 = back_sol[end]
    back_station = u0.x
    λ0 = u0.λ


    ## Get closest asteroids at that point
    sorted = sort(asteroids; by=asteroid->sum(abs2, asteroid-back_station))
    besties = sorted[1:n_candidates]


    ## Loop over chosen asteroids
    return map(besties) do bestie
        u0.x = bestie

        ## Solve for the optimal trajectory
        # Set up forward ODE problem
        tspan = (t0, back_time)
        forward_prob = remake(back_prob; u0=u0, tspan=tspan)

        # Set up nonlinear problem
        nl_u = ComponentArray(; λ=λ0, t=back_time)
        nl_p = (
            u0 = u0,
            xf = uf.x,
            prob = forward_prob,
            α = trans_scale,
            out = copy(nl_u),
        )
        nl_prob = NonlinearProblem(nl_fun, nl_u, nl_p)

        # Solve
        nl_sol = solve(nl_prob; autodiff=autodiff, kwargs...)
        u0.λ = nl_sol.u.λ

        solve(remake(forward_prob; u0=u0, tspan=(t0, nl_sol.u.t)), args...; saveat=saveat)
    end
end