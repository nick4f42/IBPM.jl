module IBPM

export IBProblem, MultiGrid
export InterpNone, InterpQuadratic
export StateCallback, StateData, solve, solve!

using LinearAlgebra
using SparseArrays
using FFTW
using LinearMaps
using IterativeSolvers
using InplaceOps  # @! macro
using Plots

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
