using ComponentArrays
using DifferentialEquations
using GTOC11Utils
using GTOC11Utils: ustrip, AU, km, m, yr, d, s, °,rad, @SVector
using LinearAlgebra
using Plots
plotlyjs()


## Utility Functions
unitize_state(state; pos_units=km, vel_units=m/s) = ComponentArray(r=pos_units.((state.r)AU), ṙ=vel_units.((state.ṙ)AU/yr))


## Data Ingestion
data = read_asteroids_file("data/Candidate_Asteroids.txt")
asteroids = GTOC11Utils.strip_ephemeris.(eachrow(data))


## Inputs
# Station
radius = 2.5AU
mean_anomoly = 0°

# Asteroid (just pick a random one)
asteroid = rand(asteroids)

# Time
Δt_guess = 3.0yr


## Solve
# Get final time and initial time guess
tf = asteroid[1]
t0_guess = tf - ustrip(yr, Δt_guess)

# Set station ephemeris
station = @SVector [tf, ustrip(AU, radius), 0, 0, 0, 0, ustrip(rad, mean_anomoly)]

# Solve the optimization problem
@time sol = GTOC11Utils.low_thrust_transfer2(station, asteroid, tf, t0_guess)
