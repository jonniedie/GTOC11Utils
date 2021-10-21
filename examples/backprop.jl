using DifferentialEquations
using GTOC11Utils
using GTOC11Utils: ustrip, km, AU, m, s, yr
using Plots
plotlyjs()


## Data ingestion
data = read_asteroids_file("data/Candidate_Asteroids.txt")
asteroids = ustrip.(reduce(hcat, propagate.(0yr, eachrow(data))))'

# Chose a random asteroid's state as the station
station = rand(eachrow(data))
station_state = state_vec(ustrip.(propagate(0yr, station)))


## Solve
# Solve the back problem for asteroids
back_time = 2.0
@time sols = get_candidate_solutions(station_state, asteroids, back_time;
    n_candidates = 100,
    autodiff = :forward,
);

# Solve the forward problem for the station
prob = remake(GTOC11Utils.spacecraft_prob, u0=station_state, tspan=(0.0, -back_time))
station_sol = solve(prob)


## Evaluate convenvergence
final_states = [sol[end].x for sol in sols]
state_errors = final_states .- Ref(station_state)
er = [(e.r)AU .|> km for e in state_errors]
ev = [(e.ṙ)AU/yr .|> km/s for e in state_errors]


## Plot
defaults = (linewidth=1.5, xlims=(-5,5), ylims=(-5,5), zlims=(-5,5), aspect_ratio=1, size=(900,900))
plt = scatter(eachcol(asteroids[:,1:3])...; color=:gray, markersize=0.2, opacity=0.5, label="asteroids")
plot!(station_sol; vars=(1,2,3), color=:red, lw=3, label="station")
primary = true
for sol in sols
    plot!(sol; vars=(1,2,3), primary=primary, color=:goldenrod, label="candidates", defaults...)
    primary = false
end
display(plt)