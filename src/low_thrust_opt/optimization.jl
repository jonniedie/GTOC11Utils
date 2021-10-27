function solve_with(vars, p; alg=DEFAULT_ALG, kwargs...)
	@unpack λ, t = vars
	u0 = ComponentArray([p.chaser; p.target; λ], getaxes(p.prob.u0))
	prob = remake(p.prob; u0, tspan=(zero(t), t), p)
	return solve(prob, alg; DEFUALT_SIM_ARGS..., kwargs...)
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
                            opt_alg=Fminbox(NelderMead()), autodiff=GalacticOptim.AutoFiniteDiff(),
                            opt_kwargs...
                            )

	x = ComponentArray(OptInput(; λ, t))
	lb = ComponentArray(OptInput(; λ=λ_lb, t=t_lb))
	ub = ComponentArray(OptInput(; λ=λ_ub, t=t_ub))

    p = (; Γ, μ, chaser, target, prob=opt_rel_prob)
    f = OptimizationFunction(loss, autodiff)
    prob = OptimizationProblem(f, x, p; lb, ub)

    return solve(prob, opt_alg; opt_kwargs...)
end


## Everything above this is garbage
function nl_fun(u, p)
	@unpack λ, t = u
	@unpack station_initial, prob, α, out, alg = p

    station_final = propagate(t, station_initial)

	u0 = copy(prob.u0)
	u0.λ = λ
	prob = remake(prob; u0, tspan=(prob.tspan[1], t))
	sol = solve(prob, alg)
	uf = sol[end]

    @. out.λ.r = 1e-10(1e10station_final[1:3] - 1e10uf.x.r)*DEFAULT_DISTANCE_UNIT |> km |> ustrip
    @. out.λ.ṙ = (station_final[4:6] - uf.x.ṙ)*(DEFAULT_DISTANCE_UNIT/DEFAULT_TIME_UNIT) |> km/s |> ustrip
    out.t = (1 + uf.λ'uf.x) * α
    return out
end

get_candidate_solutions(station, asteroids::AbstractMatrix, time_guess; kwargs...) = get_candidate_solutions(station, collect(eachrow(asteroids)), time_guess; kwargs...)
function get_candidate_solutions(station, asteroids, t0, time_guess;
                                 n_candidates=1,
                                 trans_scale=1e-8,
                                 alg=DEFAULT_ALG,
                                 autodiff=:forward,
                                 saveat=ustrip(DEFAULT_TIME_UNIT(1d)),
                                 kwargs...)
    ## Solve reverse problem
    #t0 = 0.0 # station[1]    
    #station_state_initial = station[2:end]
    tf = t0+time_guess
    station_state_final = propagate(tf, station)

    # Choose the first five final costates at random and calculate last from Hamiltonian
    λf = @SVector(rand(6)) .- 0.5
    max_i = argmax(station_state_final)
    # max_i = 6
    @set! λf[max_i] = -(1 + station_state_final[Not(max_i)]'*λf[Not(max_i)]) / station_state_final[max_i]

    # Set up problem
    uf = ComponentArray(OneVehicleSimState(x=collect(station_state_final), λ=λf))
    back_prob = remake(opt_prob; u0=uf, tspan=(tf, t0))

    # # Solve
    back_sol = solve(back_prob, Tsit5())		

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
        tspan = (t0, tf)
        forward_prob = remake(back_prob; u0=u0, tspan=tspan)

        # Set up nonlinear problem
        nl_u = ComponentArray(; λ=λ0, t=tf)
        nl_p = (
            station_initial = station,
            prob = forward_prob,
            α = trans_scale,
            out = copy(nl_u),
            alg = alg,
        )
        nl_prob = NonlinearProblem(nl_fun, nl_u, nl_p)

        # Solve
        nl_sol = solve(nl_prob; autodiff=autodiff, kwargs...)
        u0.λ = nl_sol.u.λ

				tff = nl_sol.u.t
				print("$t0,$tf\n")
        solve(remake(forward_prob; u0=u0, tspan=(t0, nl_sol.u.t)), alg; saveat=saveat)
    end
end
