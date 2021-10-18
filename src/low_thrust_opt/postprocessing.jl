"""
    run_solution(out; kwargs...)

Get the simulated results using optimal values `out` returned from `low_thrust_transfer`.
"""
run_solution(out; kwargs...) = solve_with(out.u, out.prob.p; kwargs...)


"""
    get_accelerations(out::AbstractOptimizationSolution, args...)
    get_accelerations(sol::AbstractODESolution, t=sol.t)

Get accelerations from either a simulated solution `sol` or directly from an optimization result `out`.
"""
get_accelerations(out::AbstractOptimizationSolution, args...) = get_accelerations(run_solution(out), args...)
function get_accelerations(sol::AbstractODESolution, t=sol.t)
    sol_log = get_log(sol, t)
    return reduce(hcat, sol_log.a)'
end