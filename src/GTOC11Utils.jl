module GTOC11Utils

using Base: @kwdef
using ComponentArrays: ComponentArrays, ComponentArray, getdata, getaxes
using ConcreteStructs: @concrete
using DataFrames: DataFrame, DataFrameRow
using DelimitedFiles: readdlm
using DifferentialEquations: ODEProblem, solve, remake
using DifferentialEquations.SciMLBase: AbstractODESolution, AbstractOptimizationSolution
using GalacticOptim: GalacticOptim, OptimizationFunction, OptimizationProblem, NonlinearProblem
using LinearAlgebra: norm, ×, ⋅
using Optim: Fminbox, NelderMead
using Plots: plot, plot!, scatter!, quiver!
using Rotations: RotZXZ
using SimulationLogs: @log, get_log, scope, scope!
using StaticArrays: @SVector, @SMatrix, SVector
using UnPack: @unpack
using Unitful: @unit, ustrip, m, km, s, d, °, kg
using UnitfulAstro: yr


export get_log, scope, scope!

include("utils.jl")
export read_asteroids_file, kepler_to_cartesian, cartesian_to_kepler, propagate, state_vec

include("constants.jl")

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
