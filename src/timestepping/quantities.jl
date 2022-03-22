module Quantities

using IdentityRanges
using ..IBPM: IBState, IBProblem, MultiGrid

export PreferView, AlwaysCopy
export gridranges


"""
    CopyPreference

Whether to return a view or copy of a quantity. Can be one of:
- [`PreferView`](@ref)
- [`AlwaysCopy`](@ref)
"""
abstract type CopyPreference end

"""
    PreferView <: CopyPreference

Return a view to a state's quantity if possible.
"""
struct PreferView <: CopyPreference end

"""
    AlwaysCopy <: CopyPreference

Never return a view to a state's quantity.
"""
struct AlwaysCopy <: CopyPreference end


"""
    surface_stress(state, [copy])

See [`Quantity`](@ref).
"""
function surface_stress end

"""
    body_impulse(state, [copy])

See [`Quantity`](@ref).
"""
function body_impulse end

"""
    drag_coef(state, [copy])

See [`Quantity`](@ref).
"""
function drag_coef end

"""
    lift_coef(state, [copy])

See [`Quantity`](@ref).
"""
function lift_coef end

"""
    cfl(state, [copy])

See [`Quantity`](@ref).
"""
function cfl end

"""
    slip(state, [copy])

See [`Quantity`](@ref).
"""
function slip end

"""
    body_points(state, [copy])

See [`Quantity`](@ref).
"""
function body_points end

"""
    streamfunction(state, [copy])

See [`Quantity`](@ref).
"""
function streamfunction end

"""
    vorticity(state, [copy]; [subgrids])

See [`Quantity`](@ref).
"""
function vorticity end

"""
    x_velocity(state, [copy]; [subgrids])

See [`Quantity`](@ref).
"""
function x_velocity end

"""
    y_velocity(state, [copy]; [subgrids])

See [`Quantity`](@ref).
"""
function y_velocity end


"""
    quantities

A tuple of all supported quantity functions.
"""
const quantities = (
    :surface_stress,
    :body_impulse,
    :drag_coef,
    :lift_coef,
    :cfl,
    :slip,
    :body_points,
    :streamfunction,
    :vorticity,
    :x_velocity,
    :y_velocity,
)

"""
    Quantity

Type of all quantity functions. Every `quantity` has the following signature:

    quantity(state::IBState, copy=AlwaysCopy(); ...)

`copy` is a [`CopyPreference`](@ref) that specifies whether to prefer views or always
copy. By default, a copy is always made.

---

Quantities that apply to each subdomain have the following signature:

    quantity(state, [copy]; [subgrids], ...) -> q::OffsetArray

`subgrids` specifies the subdomains to return the quantities over. The last dimension of the
returned array specifies the subdomain such that `q[..., i]` is the quantity on the `i`'th
grid.

---

The following quantities are supported:

$(join(("[`$q`](@ref)" for q in quantities), ", ")).
"""
const Quantity = Union{map(typeof ∘ eval, quantities)...}


"""
    gridranges(quantity::Quantity, problem::IBProblem; [subgrids]) -> r

The grid ranges `(x, y)` that correspond to the return value of `quantity` for each
subdomain in `subgrids`.

For each subdomain `k` in `subgrids`, the value at `quantity(...)[i,j,k]` has the coordinate
`(x, y) = (r[k][i], r[k][j])`.
"""
function gridranges(quantity, problem::IBProblem; kwargs...)
    return gridranges(quantity, problem.model.grid; kwargs...)
end

const IndexRange = AbstractUnitRange{<:Integer}


for (name, field) in (
    (:surface_stress, :fb),
    (:body_impulse, :F̃b),
    (:drag_coef, :CD),
    (:lift_coef, :CL),
    (:cfl, :cfl),
    (:slip, :slip),
    (:body_points, :xb),
)
    @eval $name(state::IBState, ::PreferView) = state.$field
end


function vorticity(
    state::IBState, ::AlwaysCopy; subgrids::IndexRange=1:state.grid[].mg
)
    grid = state.grid[]

    Γ = reshape(state.Γ, grid.nx-1, grid.ny-1, :)

    ω = Γ[:, :, IdentityRange(subgrids)]
    for i in subgrids
        δ = grid.h * 2.0 ^ (i - 1)
        ω[:, :, i] ./= δ^2
    end

    return ω
end

function gridranges(
    ::typeof(vorticity), grid::MultiGrid; subgrids::IndexRange=1:grid.mg
)
    h = grid.h
    len = grid.len

    return map(IdentityRange(subgrids)) do i
        fac = 2.0 ^ (i - 1)
        δ = h * fac
        xlen = len * fac

        ylen = xlen * (grid.ny / grid.nx)
        offx = fac * len/2.0 - len/2.0 + grid.offx
        offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

        x = range(-offx+δ, xlen-offx-δ, length=grid.nx-1)
        y = range(-offy+δ, ylen-offy-δ, length=grid.ny-1)

        (x, y)
    end
end


function streamfunction(
    state::IBState, ::PreferView; subgrids::IndexRange=1:state.grid[].mg
)
    grid = state.grid[]
    ψ = reshape(state.ψ, grid.nx-1, grid.ny-1, :)
    return @view ψ[:, :, IdentityRange(subgrids)]
end

function gridranges(
    ::typeof(streamfunction), grid::MultiGrid; subgrids::IndexRange=1:grid.mg
)
    h = grid.h
    len = grid.len

    return map(IdentityRange(subgrids)) do i
        fac = 2.0 ^ (i - 1)
        δ = h * fac
        xlen = len * fac

        ylen = xlen * (grid.ny / grid.nx)
        offx = fac * len/2.0 - len/2.0 + grid.offx
        offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

        x = range(-offx+δ, xlen-offx-δ, length=grid.nx-1)
        y = range(-offy+δ, ylen-offy-δ, length=grid.ny-1)

        (x, y)
    end
end

function x_velocity(
    state::IBState, ::AlwaysCopy; subgrids::IndexRange=1:state.grid[].mg
)
    grid = state.grid[]
    nu = grid.ny * (grid.nx + 1)

    u = sum((state.q, state.q0)) do x
        @views reshape(
            x[1:nu, :], grid.nx+1, grid.ny, :
        )[:, :, IdentityRange(subgrids)]
    end

    for i in subgrids
        δ = grid.h * 2.0 ^ (i - 1)
        u[:, :, i] ./= δ
    end

    # TODO: For some reason the top left velocity on grid 1 is erroneous.
    # The erroneous velocity was cropped out in plotting-utils.jl, but I'm keeping it in
    # here for now. Cropping it out would require changing the corresponding `gridranges`.

    return u
end

function gridranges(
    ::typeof(x_velocity), grid::MultiGrid; subgrids::IndexRange=1:grid.mg
)
    nx, ny = grid.nx, grid.ny

    h = grid.h
    len = grid.len

    return map(IdentityRange(subgrids)) do i
        fac = 2.0 ^ (i - 1)
        δ = h * fac
        xlen = len * fac

        ylen = xlen * (grid.ny / grid.nx)
        offx = fac * len/2.0 - len/2.0 + grid.offx
        offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

        x = range(-offx, xlen-offx, length=grid.nx+1)
        y = range(-offy+δ/2.0, ylen-offy-δ/2.0, length=grid.ny)

        (x, y)
    end
end


function y_velocity(
    state::IBState, ::AlwaysCopy; subgrids::IndexRange=1:state.grid[].mg
)
    grid = state.grid[]
    nu = grid.ny * (grid.nx + 1)
    nv = grid.nx * (grid.ny + 1)

    v = sum((state.q, state.q0)) do x
        @views reshape(
            x[nu.+(1:nv), :], grid.nx, grid.ny+1, :
        )[:, :, IdentityRange(subgrids)]
    end

    for i in subgrids
        δ = grid.h * 2.0 ^ (i - 1)
        v[:, :, i] ./= δ
    end

    return v
end

function gridranges(
    ::typeof(y_velocity), grid::MultiGrid; subgrids::IndexRange=1:grid.mg
)
    nx, ny = grid.nx, grid.ny

    h = grid.h
    len = grid.len

    return map(IdentityRange(subgrids)) do i
        fac = 2.0 ^ (i - 1)
        δ = h * fac
        xlen = len * fac

        ylen = xlen * (grid.ny / grid.nx)
        offx = fac * len/2.0 - len/2.0 + grid.offx
        offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

        x = range(-offx+δ/2.0, xlen-offx-δ/2.0, length=grid.nx)
        y = range(-offy, ylen-offy, length=grid.ny+1)

        (x, y)
    end
end


for func in quantities
    @eval begin
        export $func

        # By default, make a copy of the quantity.
        $func(state::IBState; kw...) = $func(state, AlwaysCopy(); kw...)
    end

    # Only one of the following needs to be defined for each quantity.
    # The other has this default definition.
    f = eval(func)
    if !hasmethod(f, (IBState, PreferView))
        @eval function $func(state::IBState, ::PreferView; kw...)
            return $func(state, AlwaysCopy(); kw...)
        end
    elseif !hasmethod(f, (IBState, AlwaysCopy))
        @eval function $func(state::IBState, ::AlwaysCopy; kw...)
            return $func(state, PreferView(); kw...) |> deepcopy
        end
    end
end

end # module
