using ComponentArrays
using DifferentialEquations
using GTOC11Utils
using GTOC11Utils: ustrip, km, AU, m, s, yr, d
using Plots
plotlyjs()


##
unitize_state(state; pos_units=km, vel_units=m/s) = ComponentArray(r=pos_units.((state.r)AU), ṙ=vel_units.((state.ṙ)AU/yr))

## Data ingestion
data = read_asteroids_file("data/Candidate_Asteroids.txt")
# asteroids = ustrip.(reduce(hcat, propagate.(0yr, eachrow(data))))'
asteroids = propagate.(0yr, eachrow(data))

# Chose a random asteroid's state as the station
station = NamedTuple(rand(eachrow(data)))
station_state = state_vec(propagate(0yr, station))


## Solve the back problem for asteroids
back_time = 1.5
@time sols = get_candidate_solutions(station_state, asteroids, back_time;
    n_candidates = 20,
    trans_scale=1e-8,
    autoscale=false,
);


## Solve the forward problem for the station
prob = remake(GTOC11Utils.spacecraft_prob, u0=station_state, tspan=(0.0, -back_time))
tols = ComponentArray(r=fill(ustrip(AU(10km)), 3), ṙ=fill(ustrip((AU/yr)(0.01m/s)), 3)) * 1e-5
@time station_sol = solve(prob; abstol=tols)

sundials_prob = remake(prob; u0=station_sol[end], tspan=(0.0, back_time))
@time sundials_sol = solve(sundials_prob, GTOC11Utils.DEFAULT_ALG; GTOC11Utils.DEFAULT_SIM_ARGS...)

error = station_sol[1] - sundials_sol[end] |> unitize_state


## Evaluate convenvergence
final_states = [sol[end].x for sol in sols]
state_errors = final_states .- Ref(station_state)
er = [(e.r)AU .|> km for e in state_errors]
ev = [(e.ṙ)AU/yr .|> m/s for e in state_errors]


## Plot
defaults = (linewidth=1.5, xlims=(-5,5), ylims=(-5,5), zlims=(-5,5), aspect_ratio=1)
plt = scatter(eachcol(reduce(hcat,asteroids)'[:,1:3])...; color=:gray, markersize=0.2, opacity=0.5, label="asteroids")
# plt = scatter(eachcol(asteroids[:,1:3])...; color=:gray, markersize=0.2, opacity=0.5, label="asteroids")
plot!(station_sol; vars=(1,2,3), color=:red, lw=3, label="station")
primary = true
for sol in sols
    plot!(sol; vars=(1,2,3), primary=primary, color=:goldenrod, label="candidates", defaults...)
    primary = false
end
display(plt)