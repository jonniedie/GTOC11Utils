module GTOC11Utils

using Base: @kwdef
using ComponentArrays: ComponentArrays, ComponentArray, getdata, getaxes
using ConcreteStructs: @concrete
using DataFrames: DataFrame, DataFrameRow
using DataInterpolations: LagrangeInterpolation
using DelimitedFiles: readdlm
using DiffEqSensitivity: InterpolatingAdjoint
using DifferentialEquations: ODEProblem, Tsit5, solve, remake
using DifferentialEquations.SciMLBase: AbstractODESolution, AbstractOptimizationSolution
import ForwardDiff
using GalacticOptim: GalacticOptim, OptimizationFunction, OptimizationProblem
using InvertedIndices: Not
using LinearAlgebra: norm, normalize, ×, ⋅, ldiv!, lu
using NLopt: Opt, equality_constraint!, optimize
using NonlinearSolve: NonlinearProblem, NewtonRaphson, DEFAULT_LINSOLVE
using Optim: Fminbox, NelderMead
using Plots: plot, plot!, scatter!, quiver!
using Rotations: RotZXZ, UnitQuaternion
using SciMLNLSolve: NLSolveJL, CMINPACK
using Setfield: @set!
using SimulationLogs: @log, get_log, scope, scope!
using StaticArrays: @SVector, @SMatrix, SVector
using Sundials: ARKODE, Explicit, FEHLBERG_13_7_8
using UnPack: @unpack
using Unitful: @unit, Quantity, NoDims, FreeUnits, ustrip, unit, m, km, s, d, °, kg, rad
using UnitfulAstro: yr


export get_log, scope, scope!

include("units_and_constants.jl")

include("utils.jl")
export read_asteroids_file, state_vec, delta_v_with_radius

include("orbit_conversions.jl")
export kepler_to_cartesian, cartesian_to_kepler, propagate, strip_ephemeris

include("low_thrust_opt/types.jl")
export OneVehicleSimState, TwoVehicleSimState, OptInput

include("low_thrust_opt/equations_of_motion.jl")

include("low_thrust_opt/optimization.jl")
export low_thrust_transfer, get_candidate_solutions

include("low_thrust_opt/postprocessing.jl")
export run_solution, get_accelerations

include("low_thrust_opt/plotting.jl")
export chase_plot, acceleration_plot


end # module
