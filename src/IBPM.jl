module IBPM

export IBProblem, MultiGrid, StateData, solve, solve!

using LinearAlgebra
using SparseArrays
using FFTW
using LinearMaps
using IterativeSolvers
using InplaceOps  # @! macro
using Plots


"""
    gridstep(problem::AbstractIBProblem)
    gridstep(grid::Grid)

The minimum grid step size of a problem or grid.
"""
function gridstep end

"""
    timestep(problem::AbstractIBProblem)

The time step size of a problem.
"""
function timestep end

#Caution, include order matters!
include("fluid-domain/fluid-domain-include.jl")
include("structure-domain/structure-domain-include.jl")
include("interface-coupling/interface-coupling-include.jl")
include("pre-processing/pre-processing-include.jl")
include("fluid-operators/fluid-operators-include.jl")
include("timestepping/timestepping-include.jl")
include("plotting/plotting-include.jl")
include("experimental/sfd.jl")

end
