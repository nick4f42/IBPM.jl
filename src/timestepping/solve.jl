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


"""
    InterpKind

Used to control how [`StateData`](@ref) does or doesn't interpolate data.

Currently, the following kinds are supported:
- [`InterpNone`](@ref)
- [`InterpQuadratic`](@ref)
"""
abstract type InterpKind end

"""
    InterpNone <: InterpKind

Do not interpolate. Use the nearest time and state.
"""
struct InterpNone <: InterpKind end

"""
    InterpQuadratic <: InterpKind

Quadratically interpolate between states.
"""
struct InterpQuadratic <: InterpKind end


"""
    StateCallback(f::Function, problem::IBProblem; at::AbstractVector)

In [`solve`](@ref), call `f(t, state)` each time step in `at` where `t` and `state` are the
time and [`IBState`](@ref).

If `at` is a vector of floats, `f` is called at the nearest time steps. If it is a vector
of integers, `f` is called at each i'th time step.

---

    StateCallback(
        f::Function,
        problem::IBProblem,
        quantities::Function...;
        at::AbstractVector,
        interp::InterpKind = InterpNone()
    )

In [`solve`](@ref), call `f(t, q...)` with the time and return values from `quantities`.

Each function in `quantities` must have only one method with a signature `(state,)` or `(t,
state)`. Each function must not return a view of `state` unless `interp` is `InterpNone`.

If `interp` is not `InterpNone`, the return value of each `quantity` is interpolated to
approximate the values between time steps. The return value must either be a number or
define broadcasted arithmetic operators (`.+`, `.*`, etc).

The [`IBPM.Quantities`](@ref) module defines helpful functions that may be passed as
`quantity`.
"""
struct StateCallback
    callback::Function
    function StateCallback(
        f::Function,
        problem::IBProblem,
        quantities::Function...;
        at::AbstractVector,
        interp::InterpKind = InterpNone()
    )
        qty_funcs = map(quantities) do qty
            nargs = _infer_nargs(qty)
            if nargs == 1
                (_, state) -> qty(state)
            elseif nargs == 2
                qty
            else
                throw(ArgumentError(
                    "quantity function must have signature (state,) or (t, state)"
                ))
            end
        end

        times = _interp_times(problem, at, interp)
        return new(_state_callback(f, qty_funcs, times, problem, interp))
    end
end

function StateCallback(f::Function, problem::IBProblem; at::AbstractVector)
    return StateCallback(f, problem, identity; at, interp=InterpNone())
end

function _infer_nargs(quantity::Function)
    method_iter = Iterators.flatten((
        methods(quantity, (IBState,)),
        methods(quantity, (AbstractFloat, IBState))
    ))
    nargs = Iterators.map(m -> m.nargs-1, method_iter)

    isempty(nargs) && return nothing
    n = first(nargs)
    return all(==(n), nargs) ? n : nothing
end


function _state_callback(
    f::Function, quantities, times::AbstractVector{Float64},
    ::IBProblem, ::InterpNone
)
    t_interp = Iterators.Stateful(times)

    return function (t, state)
        if peek(t_interp) == t
            popfirst!(t_interp)
            vals = map(qty -> qty(t, state), quantities)
            f(t, vals...)
        end
        nothing
    end
end

function _state_callback(
    f::Function, quantities, times::AbstractVector{Float64},
    problem::IBProblem, ::InterpQuadratic
)
    t_interp = Iterators.Stateful(times)
    t_sim = problem.t
    dt = step(t_sim)

    cache = []

    i = 0
    return function (t0, state)
        isempty(t_interp) && return

        i += 1
        t = peek(t_interp)
        t2 = t_sim[min(i+2, length(t_sim))]

        t > t2 && return

        current = map(qty -> qty(t0, state), quantities)

        if t > t0 || i < 3
            push!(cache, current)
            return
        end

        while true
            x = (t - t0) / dt
            vals = map(cache..., current) do y...
                _quadratic_interp(y, x)
            end
            f(t, vals...)

            popfirst!(t_interp)
            isempty(t_interp) && return
            t = peek(t_interp)
            t > t0 && break
        end

        popfirst!(cache)
        t ≤ t2 && push!(cache, current)
    end
end

function _quadratic_interp end

let interp_expr = :(
    x * (x + 1) / 2 * a
    - x * (x + 2) * b
    + (x + 1) * (x + 2) / 2 * c
)
    function _factory(T, expr)
        # Quadratically interpolate values a, b, c at x.
        # a, b, c are at x values -2, -1, and 0 respectively.
        @eval function _quadratic_interp((a, b, c)::NTuple{3,$T}, x)
            return $expr
        end
    end

    _factory(Number, interp_expr)
    _factory(Any, :(@. $interp_expr))
end

function _interp_times(
        problem::IBProblem,
        times::AbstractVector{<:Integer},
        ::InterpKind
    )
    return problem.t[times]
end

function _interp_times(
        problem::IBProblem,
        times::AbstractVector{<:AbstractFloat},
        ::InterpNone
    )
    return _nearest_values(problem.t, times)
end

function _interp_times(
        problem::IBProblem,
        times::AbstractVector{<:AbstractFloat},
        ::InterpKind
    )
    if first(times) < problem.t[1] || last(times) > problem.t[end]
        throw(DomainError(times, "times must be within the problem's time span"))
    end
    return times
end

function _nearest_values(a::AbstractVector{T}, b::AbstractVector{T}) where T
    # Returns a vector of unique closest values in `a` to values in `b`.
    # Assumes `a` and `b` are sorted.

    vals = T[]
    isempty(a) && return vals

    b_iter = Iterators.Stateful(b)

    for (x1, x2) in zip(a, Iterators.drop(a, 1))
        isempty(b_iter) && break

        y = peek(b_iter)
        y - x1 < x2 - y || continue

        push!(vals, x1)
        popfirst!(b_iter)
        while !isempty(b_iter) && (y = peek(b_iter); y - x1 < x2 - y)
            popfirst!(b_iter)
        end
    end

    isempty(b_iter) || push!(vals, last(a))

    return vals
end

struct StateData{T, V<:AbstractVector} <: AbstractVector{T}
    quantity::Function
    problem::IBProblem
    interp::InterpKind
    data::Vector{T}
    t::V

    function StateData{T}(
        quantity::Function, problem::IBProblem;
        saveat::AbstractVector, interp::InterpKind=InterpNone()
    ) where T
        times = _interp_times(problem, saveat, interp)
        data = sizehint!(T[], length(times))

        return new{T, typeof(times)}(quantity, problem, interp, data, times)
    end
end

Base.size(s::StateData) = size(s.data)
Base.getindex(s::StateData, i::Int) = s.data[i]
Base.IndexStyle(::StateData) = IndexLinear()

function saveto!(s::StateData)
    return StateCallback(s.problem, s.quantity; at=s.t, interp=s.interp) do _, qty
        push!(s.data, qty)
    end
end


function solve!(f::Function, state::IBState, problem::IBProblem)
    for t in problem.t
        advance!(state, problem, t)

        if !all(isfinite, state.CL)
            @error "Lift coefficient is isfinite, terminating loop."
            break
        end

        f(t, state) === false && break
    end

    return state
end

function solve!(f, state, problem, callbacks::AbstractVector{StateCallback})
    return solve!(state, problem) do t, state
        for cb in callbacks
            cb.callback(t, state)
        end

        f(t, state)
    end
end

solve!(args...) = solve!((_,_)->nothing, args...)

function solve(f::Function, problem::IBProblem, args...)
    return solve!(f, IBState(problem), problem, args...)
end

function solve(problem::IBProblem, args...)
    return solve!(IBState(problem), problem, args...)
end


"""
    compute_cfl(state, prob)

Compute the CFL number (uΔt/Δx) based on the fine-grid flux

Note that this uses working memory that is also used in `nonlinear!`
"""
function compute_cfl(state, prob)
    Δt, Δx = prob.scheme.dt, prob.model.grid.h
    qwork = prob.model.work.q5
    @views @. qwork = abs( state.q[:, 1] )
    return maximum(qwork)*Δt/Δx
end
