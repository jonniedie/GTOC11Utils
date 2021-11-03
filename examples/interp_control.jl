using ComponentArrays
using DifferentialEquations
using GTOC11Utils
using GTOC11Utils: ustrip, AU, km, m, yr, d, s, °, rad
using StaticArrays
using Plots
plotlyjs()


## Utility Functions
unitize_state(state; pos_units=km, vel_units=m/s) = ComponentArray(r=pos_units.((state.r)AU), ṙ=vel_units.((state.ṙ)AU/yr))

my_scatter!(v::AbstractVector{<:Real}; kwargs...) = scatter!(v[1:1], v[2:2], v[3:3]; kwargs...)
my_scatter!(v::AbstractVector{<:AbstractVector}; kwargs...) = my_scatter!(reduce(hcat, v)'; kwargs...)
my_scatter!(m::AbstractMatrix; kwargs...) = scatter!(m[:,1], m[:,2], m[:,3]; kwargs...)


## Data Ingestion
# Read asteroids in from file, put in canonical units, and strip the units
data = read_asteroids_file("data/Candidate_Asteroids.txt")
asteroids = GTOC11Utils.strip_ephemeris.(eachrow(data))


## Inputs
# Station
radius = 2.5AU
mean_anomoly = 0°

# Asteroid
asteroid = rand(asteroids)

# Time
Δt_guess = 3.0yr


## Solve
# Get final time and initial time guess
tf = asteroid[1]
t0_guess = tf - ustrip(yr, Δt_guess)

# Set station ephemeris
station = @SVector [tf, ustrip(AU, radius), 0, 0, 0, 0, ustrip(rad, mean_anomoly)]

# Solve optimization problem
@time sol = interp_control_transfer(station, asteroid, tf, t0_guess;
    time_divisions = 10,
)
display(sol.opt.original)

# Get error
asteroid_final = unitize_state(sol.sim[end])
station_state = unitize_state(state_vec(propagate(tf, station)))
e = asteroid_final - station_state


## Plotting
plot()
plot!(sol.sim; vars=(1,2,3), label="Asteroid trajectory", lw=2)
my_scatter!(ustrip.(AU.(station_state.r)); label="Station final state", xlims=(-3,3), ylims=(-3,3), zlims=(-3,3))
