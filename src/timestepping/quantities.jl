"""
    IBPM.Quantities

Types for various quantities and methods to retrieve those quantities from a state.

The [`Quantity`](@ref) type is the supertype of all quantities.
"""
module Quantities

using Interpolations
using StaticArrays
using LinearAlgebra

include("dynamics.jl")
using .Dynamics2D
using ..IBPM: IBState, IBProblem, AbstractIBProblem, MultiGrid, grid_of

export DynamicsScalar, DynamicsVector
export grid_basis, grid_origin, grid_frame, lab_basis, lab_origin, lab_frame

export Quantity, BodyQuantity, DomainVector, DomainScalar, gridranges
export Dense, Discrete
export StreamFunction, Vorticity, Velocity, VelocityX, VelocityY, VelocityNorm,
    CFL, Slip, SurfaceStress, BodyImpulse, BodyPoints, DragCoef, LiftCoef

cross2d((x1, y1), (x2, y2)) = x1 * y2 - y1 * x2
rot90((x, y)) = SVector(-y, x)
function rot((x, y), θ)
    c = cos(θ)
    s = sin(θ)
    return SVector(c * x - s * y, s * x + c * y)
end

abstract type DomainKind end
struct Discrete <: DomainKind end
struct Dense{T} <: DomainKind
    interp::T
end

abstract type Quantity end
abstract type BodyQuantity <: Quantity end
abstract type DomainVector{N} <: Quantity end
const DomainScalar = DomainVector{1}

struct StreamFunction <: DomainScalar end
struct Vorticity <: DomainScalar end
struct Velocity <: DomainScalar end
struct VelocityX <: DomainScalar end
struct VelocityY <: DomainScalar end
struct VelocityNorm <: DomainVector{2} end

struct CFL <: Quantity end
struct Slip <: Quantity end

struct SurfaceStress <: BodyQuantity end
struct BodyImpulse <: BodyQuantity end
struct BodyPoints <: BodyQuantity end
struct DragCoef <: BodyQuantity end
struct LiftCoef <: BodyQuantity end

(T::Type{<:Quantity})(args...; kw...) = T()(args...; kw...)

for (Qty, field) in (
    (:CFL, :cfl),
    (:Slip, :slip),
)
    @eval (::$Qty)(::IBProblem) = state::IBState -> state.$field
end

for (Qty, field) in (
    (:SurfaceStress, :fb),
    (:BodyImpulse, :F̃b),
    (:BodyPoints, :xb),
    (:DragCoef, :CD),
    (:LiftCoef, :CL),
)
    @eval (::$Qty)(::IBProblem, bodies=(:)) = state::IBState -> state.$field[bodies]
end

default_dense_kind(::IBProblem) = (Dense ∘ BSpline ∘ Linear)()

function (qty::DomainVector)(
    prob::AbstractIBProblem;
    kind::DomainKind=default_dense_kind(prob),
    frame::AbstractFrame=GridFrame(),
    kw...
)
    return qty(prob, kind, frame; kw...)
end

function (::StreamFunction)(
    prob::IBProblem, ::Discrete, ::GridFrame;
    subgrids=1:grid_of(prob).mg
)
    grid = grid_of(prob)
    return state::IBState -> reshape(state.ψ, grid.nx - 1, grid.ny - 1, :)[:, :, subgrids]
end

function (q::StreamFunction)(prob::IBProblem, d::Dense, frame::GridFrame; kw...)
    return _interp_quantity(q, prob, d; kw...)
end

function gridranges(grid::MultiGrid, ::StreamFunction; subgrids=1:grid.mg)
    h = grid.h
    len = grid.len

    return map(subgrids) do i
        fac = 2.0^(i - 1)
        δ = h * fac
        xlen = len * fac

        ylen = xlen * (grid.ny / grid.nx)
        offx = fac * len / 2.0 - len / 2.0 + grid.offx
        offy = fac * (grid.ny * h) / 2.0 - (grid.ny * h) / 2.0 + grid.offy

        xs = range(-offx + δ, xlen - offx - δ, length=grid.nx - 1)
        ys = range(-offy + δ, ylen - offy - δ, length=grid.ny - 1)

        (xs, ys)
    end
end

function (::Vorticity)(
    prob::IBProblem, ::Discrete, ::GridFrame;
    subgrids=1:grid_of(prob).mg
)
    grid = grid_of(prob)

    return function (state::IBState)
        ω = reshape(state.Γ, grid.nx - 1, grid.ny - 1, :)[:, :, subgrids]
        for (i, subgrid) in zip(axes(ω, 3), subgrids)
            δ = grid.h * 2.0^(subgrid - 1)
            ω[:, :, i] ./= δ^2
        end

        return ω
    end
end

function (q::Vorticity)(prob::IBProblem, d::Dense, frame::GridFrame; kw...)
    return _interp_quantity(q, prob, d; kw...)
end

function gridranges(grid::MultiGrid, ::Vorticity; subgrids=1:grid.mg)
    h = grid.h
    len = grid.len

    return map(subgrids) do i
        fac = 2.0^(i - 1)
        δ = h * fac
        xlen = len * fac

        ylen = xlen * (grid.ny / grid.nx)
        offx = fac * len / 2.0 - len / 2.0 + grid.offx
        offy = fac * (grid.ny * h) / 2.0 - (grid.ny * h) / 2.0 + grid.offy

        xs = range(-offx + δ, xlen - offx - δ, length=grid.nx - 1)
        ys = range(-offy + δ, ylen - offy - δ, length=grid.ny - 1)

        (xs, ys)
    end
end

function (::VelocityX)(
    prob::IBProblem, ::Discrete, ::GridFrame;
    subgrids=1:grid_of(prob).mg
)
    grid = grid_of(prob)
    nu = grid.ny * (grid.nx + 1)

    return function (state::IBState)
        u = sum((state.q, state.q0)) do x
            @views reshape(x[1:nu, :], grid.nx + 1, grid.ny, :)[:, :, subgrids]
        end

        for (i, subgrid) in zip(axes(u, 3), subgrids)
            δ = grid.h * 2.0^(subgrid - 1)
            u[:, :, i] ./= δ
        end

        # TODO: For some reason the top left velocity on grid 1 is erroneous.

        # (Nick) The erroneous velocity was cropped out in plotting-utils.jl, see:
        # https://github.com/NUFgroup/IBPM.jl/blob/aca48b50c9ca92d24dab5659edaf7196d1328f63/src/plotting/plotting-utils.jl#L146-L152
        # It won't be cropped out here since it messes with the `gridranges` function
        #return @view u[2:end-1, 2:end-1, axes(u)[3:end]...]

        return u
    end
end

function (q::VelocityX)(prob::IBProblem, d::Dense, frame::GridFrame; kw...)
    return _interp_quantity(q, prob, d; kw...)
end

function gridranges(grid::MultiGrid, ::VelocityX; subgrids=1:grid.mg)
    h = grid.h
    len = grid.len

    return map(subgrids) do i
        fac = 2.0^(i - 1)
        δ = h * fac
        xlen = len * fac

        ylen = xlen * (grid.ny / grid.nx)
        offx = fac * len / 2.0 - len / 2.0 + grid.offx
        offy = fac * (grid.ny * h) / 2.0 - (grid.ny * h) / 2.0 + grid.offy

        x = range(-offx, xlen - offx, length=grid.nx + 1)
        y = range(-offy + δ / 2.0, ylen - offy - δ / 2.0, length=grid.ny)

        (x, y)
    end
end


function (::VelocityY)(
    prob::IBProblem, ::Discrete, ::GridFrame;
    subgrids=1:grid_of(prob).mg
)
    grid = grid_of(prob)
    nu = grid.ny * (grid.nx + 1)

    return function (state::IBState)
        v = sum((state.q, state.q0)) do x
            @views reshape(x[nu+1:end, :], grid.nx, grid.ny + 1, :)[:, :, subgrids]
        end

        for (i, subgrid) in zip(axes(v, 3), subgrids)
            δ = grid.h * 2.0^(subgrid - 1)
            v[:, :, i] ./= δ
        end

        return v
    end
end

function (q::VelocityY)(prob::IBProblem, d::Dense, frame::GridFrame; kw...)
    return _interp_quantity(q, prob, d; kw...)
end

function gridranges(grid::MultiGrid, ::VelocityY; subgrids=1:grid.mg)
    h = grid.h
    len = grid.len

    return map(subgrids) do i
        fac = 2.0^(i - 1)
        δ = h * fac
        xlen = len * fac

        ylen = xlen * (grid.ny / grid.nx)
        offx = fac * len / 2.0 - len / 2.0 + grid.offx
        offy = fac * (grid.ny * h) / 2.0 - (grid.ny * h) / 2.0 + grid.offy

        x = range(-offx + δ / 2.0, xlen - offx - δ / 2.0, length=grid.nx)
        y = range(-offy, ylen - offy, length=grid.ny + 1)

        (x, y)
    end
end

function (::StreamFunction)(
    prob::IBProblem, d::Dense, frame::Frame;
    subgrids=1:grid_of(prob).mg
)
    ψ_grid = StreamFunction()(prob, d, GridFrame(); subgrids)

    error("not implemented")
    # rf =
    # vf =
    # θf =
    # ωf =

    return function (state::IBState)
        t = state.t
        ψ = ψ_grid(state)
        return function (x, y)
            r = SVector(x, y)
            rg = rf(t) + rot(r, θf(t))
            return ψ(rg...) - cross2d(vf(t) / 2 + ωf(t) * rot90(r) / 3, r)
        end
    end
end

function (::Vorticity)(
    prob::IBProblem, d::Dense, frame::Frame;
    subgrids=1:grid_of(prob).mg
)
    vort_grid = Vorticity()(prob, d, GridFrame(); subgrids)

    error("not implemented")
    # rf =
    # vf =
    # θf =
    # ωf =

    return function (state::IBState)
        t = state.t
        vort = vort_grid(state)
        return function (x, y)
            r = SVector(x, y)
            rg = rf(t) + rot(r, θf(t))
            vort(rg...) - 2 * ωf(t)
        end
    end
end

function (q::Velocity)(prob::IBProblem, d::Dense, frame::GridFrame; kw...)
    u = VelocityX()(prob, d, GridFrame(); kw...)
    v = VelocityY()(prob, d, GridFrame(); kw...)
    return function(state::IBState)
        us = u(state)
        vs = v(state)
        return (x, y) -> SVector(us(x, y), vs(x, y))
    end
end

function (::Velocity)(
    prob::IBProblem, d::Dense, frame::Frame;
    subgrids=1:grid_of(prob).mg
)
    # fluid velocity in the grid frame
    u_grid = Velocity()(prob, d, GridFrame(); subgrids)

    error("not implemented")
    # rf =
    # vf =
    # θf =
    # ωf =

    return function (state::IBState)
        t = state.t
        u = u_grid(state)
        return function (x, y)
            r = SVector(x, y)
            rg = rf(t) + rot(r, θf(t))
            ug = u(rg...)
            return rot(ug, -θf(t)) - vf(t) - ωf(t) * rot90((x, y))
        end
    end
end

function (q::VelocityNorm)(prob::IBProblem, d::Dense, frame::AbstractFrame; kw...)
    u = Velocity()(prob, d, GridFrame(); kw...)
    return function(state::IBState)
        us = u(state)
        (x, y) -> norm(us(x, y))
    end
end

# Interpolates a quantity in the grid frame.
function _interp_quantity(
    qty::DomainVector{N}, prob::IBProblem, d::Dense; subgrids=1:grid_of(prob).mg
) where {N}
    subgrids = sort(subgrids)
    ranges = gridranges(grid_of(prob), qty; subgrids)
    f = qty(prob, Discrete(), GridFrame(); subgrids)

    return function (state::IBState)
        return _subdomain_interp!(f(state), ranges, d.interp, Val(N))
    end
end

# Interpolation function for a scalar quantity
# Overwrites `a`
function _subdomain_interp!(a::AbstractArray, ranges, interp, ::Val{1})
    # Interp function for each subdomain
    funcs = map(eachslice(a; dims=3), ranges) do a_sub, (xs, ys)
        itp = extrapolate(interpolate!(a_sub, interp), Flat())
        scale(itp, xs, ys)
    end
    bounds = map(r -> extrema.(r), ranges)

    return function (x, y)
        i = _smallest_grid_index(bounds, (x, y))
        return funcs[i](x, y)
    end
end

# Interpolation function for a vector quantity
# Overwrites `a`
function _subdomain_interp!(a::AbstractArray, ranges, interp, ::Val{N}) where {N}
    # Interp function for each subdomain
    funcs = map(eachslice(a; dims=3), ranges) do slice, (xs, ys)
        # Each component of the vector
        fs = eachslice(slice; dims=3) do a_sub
            itp = extrapolate(interpolate!(a_sub, interp), Flat())
            scale(itp, xs, ys)
        end
        SVector{N}(fs)
    end
    bounds = map(r -> extrema.(r), ranges)

    return function (x, y)
        i = _smallest_grid_index(bounds, (x, y))
        return map(f -> f(x, y), funcs[i])
    end
end

# The smallest index of the grid that might contain position pos
# Assumes bounds are sorted from smallest to largest
function _smallest_grid_index(bounds, pos::NTuple{N}) where {N}
    for (i, b::NTuple{N,NTuple{2}}) in pairs(bounds)
        if all(c1 < c < c2 for ((c1, c2), c) in zip(b, pos))
            return i
        end
    end
    return lastindex(bounds)
end

end # module
