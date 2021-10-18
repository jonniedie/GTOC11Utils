module GTOC11Utils

using Base: @kwdef
using ComponentArrays: ComponentArrays, ComponentArray, getdata, getaxes
using ConcreteStructs: @concrete
using DifferentialEquations: ODEProblem, solve, remake
using DifferentialEquations.SciMLBase: AbstractODESolution, AbstractOptimizationSolution
using GalacticOptim: GalacticOptim, OptimizationFunction, OptimizationProblem
using LinearAlgebra: norm, ×
using Optim: Fminbox, NelderMead
using Plots: plot, plot!, scatter!, quiver!
using SimulationLogs: @log, get_log, scope, scope!
using StaticArrays: @SVector, @SMatrix, SVector
using UnPack: @unpack
using Unitful: @unit, ustrip, m, km, s, d
using UnitfulAstro: yr


export get_log, scope, scope!

include("utils.jl")

include("constants.jl")

include("low_thrust_opt/types.jl")

include("low_thrust_opt/equations_of_motion.jl")

include("low_thrust_opt/optimization.jl")
export low_thrust_transfer

include("low_thrust_opt/postprocessing.jl")
export run_solution, get_accelerations

include("low_thrust_opt/plotting.jl")
export chase_plot, acceleration_plot


end # module
