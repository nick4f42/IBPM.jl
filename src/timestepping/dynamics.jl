module Dynamics2D

using ...IBPM: IBProblem, motions_of, Static, MotionFunction, MovingGrid

export AbstractBasis, GridBasis, RelativeBasis
export AbstractPoint, GridOrigin, OffsetPoint
export AbstractFrame, GridFrame, Frame
export grid_basis, grid_origin, grid_frame
export lab_basis, lab_origin, lab_frame

abstract type AbstractBasis end
struct GridBasis <: AbstractBasis end

abstract type AbstractPoint end
struct GridOrigin <: AbstractPoint end

struct RelativeBasis <: AbstractBasis
    # TODO: Implement
end

struct OffsetPoint <: AbstractPoint
    # TODO: Implement
end

# The frame structs are separated to avoid type instabilities in quantities.jl
abstract type AbstractFrame end
struct GridFrame <: AbstractFrame end
struct Frame{P<:AbstractPoint,B<:AbstractBasis} <: AbstractFrame
    # TODO: Implement
end

framedef(::GridOrigin, ::GridBasis) = GridFrame()
framedef(origin::AbstractPoint, basis::AbstractBasis) = Frame(origin, basis)

grid_basis(::IBProblem) = GridBasis()
grid_origin(::IBProblem) = GridOrigin()
grid_frame(::IBProblem) = GridFrame()

lab_basis(prob::IBProblem) = prob |> motions_of |> first |> lab_basis
lab_basis(::Union{Static,MotionFunction}) = GridBasis()
lab_basis(m::MovingGrid) = error("not implemented")

lab_origin(prob::IBProblem) = prob |> motions_of |> first |> lab_origin
lab_origin(::Union{Static,MotionFunction}) = GridOrigin()
lab_origin(m::MovingGrid) = error("not implemented")

lab_frame(x) = framedef(lab_origin(x), lab_basis(x))

end # module
