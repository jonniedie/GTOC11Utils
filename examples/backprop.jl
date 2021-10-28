using ComponentArrays
using DifferentialEquations
using GTOC11Utils
using GTOC11Utils: ustrip, km, AU, m, s, yr, d, NLSolveJL, NewtonRaphson, @SVector
using LinearAlgebra
using Plots
plotlyjs()


##
unitize_state(state; pos_units=km, vel_units=m/s) = ComponentArray(r=pos_units.((state.r)AU), ṙ=vel_units.((state.ṙ)AU/yr))

## Data ingestion
data = read_asteroids_file("data/Candidate_Asteroids.txt")
asteroids = GTOC11Utils.strip_ephemeris.(eachrow(data))

# Chose a random asteroid's state as the station
station = @SVector [
    ustrip(yr, data[1,:epoch]),
    2.5,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
]


## Solve the back problem for asteroids
t0_guess = station[1]
tf = t0_guess + 1.75
@time sols = get_candidate_solutions(station, asteroids, tf, t0_guess;
    n_candidates = 20,
    trans_scale=1,
    # alg=Tsit5(),
    nl_alg=NLSolveJL(
        autoscale=false,
        # autodiff=:forward,
        method=:trust_region,
        # ftol=1e-10,
        # xtol=1e-8,
    ),
    # nl_alg=NewtonRaphson(
    #     autodiff=false,
    # ),
    maxiters=100,
    ftol=1e-10,
    xtol=1e-10,
    show_trace=true,
);

# final_stations = [collect(propagate(sol.t[end], station)) for sol in sols]
station_initial_state = state_vec(propagate(t0_guess, station))
station_final_state = state_vec(propagate(tf, station))


## Evaluate convenvergence
final_states = [sol[end].x for sol in sols]
state_errors = final_states .- Ref(station_final_state)
er = [(e.r)AU .|> km for e in state_errors]
ev = [(e.ṙ)AU/yr .|> m/s for e in state_errors]
count([all(abs.(x).<2m/s) for x in ev] .& [all(abs.(x).<1km) for x in er])




## Solve the forward problem for the tf
prob = remake(GTOC11Utils.spacecraft_prob, u0=station_initial_state, tspan=(t0_guess, tf))
tols = ComponentArray(r=fill(ustrip(AU(10km)), 3), ṙ=fill(ustrip((AU/yr)(0.01m/s)), 3)) * 1e-5
@time station_sol = solve(prob; abstol=tols);

# sundials_prob = remake(prob; u0=station_sol[end], tspan=(0.0, back_time))
# @time station_sol = solve(sundials_prob, GTOC11Utils.DEFAULT_ALG; GTOC11Utils.DEFAULT_SIM_ARGS...)

# error = station_sol[1] - sundials_sol[end] |> unitize_state



## Plot
asteroid_states = propagate.(tf, asteroids)
defaults = (linewidth=1.5, xlims=(-5,5), ylims=(-5,5), zlims=(-5,5), aspect_ratio=1, size=(900,900))
plt = plot(; defaults...)
scatter!(eachcol(reduce(hcat,asteroid_states)'[:,1:3])...; color=:gray, markersize=0.2, opacity=0.5, label="asteroids")
plot!(plt, station_sol; vars=(1,2,3), color=:red, lw=3, label="station", defaults...)
primary = true
for sol in sols
    plot!(plt, sol; vars=(1,2,3), primary=primary, color=:goldenrod, label="candidates", defaults...)
    primary = false
end
display(plt)