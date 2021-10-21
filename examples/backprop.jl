using DifferentialEquations
using GTOC11Utils
using GTOC11Utils: yr
using Plots
plotlyjs()


## Data ingestion
data = read_asteroids_file("data/asteroids.txt")
asteroids = ustrip.(reduce(hcat, propagate.(0yr, eachrow(data))))'

# Chose a random asteroid's state as the station
station = rand(eachrow(data))
station_state = state_vec(ustrip.(propagate(0yr, station)))


## Solve
# Solve the back problem for asteroids
back_time = 1.5
sols = get_candidate_solutions(station_state, ast_mat, back_time; n_candidates=20);

# Solve the forward problem for the station
prob = remake(GTOC11Utils.spacecraft_prob, u0=station_state, tspan=(0.0, -back_time))
station_sol = solve(prob)


## Plot
defaults = (linewidth=1.5, xlims=(-5,5), ylims=(-5,5), zlims=(-5,5), aspect_ratio=1)
plt = scatter(eachcol(asteroids[:,1:3])...; color=:gray, markersize=0.2, opacity=0.5, widen=true)
plot!(station_sol; vars=(1,2,3), color=:red, lw=3, label="station")
for sol in sols
    plot!(sol; vars=(1,2,3), label=false, primary=false, color=:lightblue, defaults...)
end
display(plt)