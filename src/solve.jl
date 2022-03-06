struct StateData{D, S<:AbstractVector} <: AbstractVector{D}
    func::Function
    data::Vector{D}
    t::Vector{Float64}  # actual times
    saveat::S  # requested times to save at

    function StateData(
            func::Function, ::Type{D}, nargs::Integer;
            saveat::S=[]
        ) where {D, S<:AbstractVector}

        f = if nargs == 1
            (_, state) -> func(state)
        elseif nargs == 2
            func
        else
            throw(ArgumentError("nargs must be 1 or 2"))
        end

        new{D, S}(f, D[], Float64[], saveat)
    end
end

function StateData(func::Function, ::Type{D}; kwargs...) where D
    nargs, _ = _infer_signature(func)
    StateData(func, D, nargs; kwargs...)
end

function StateData(func::Function; kwargs...)
    nargs, return_type = _infer_signature(func)
    StateData(func, return_type, nargs; kwargs...)
end

StateData(; kwargs...) = StateData(identity, Any, 1; kwargs...)

function _infer_signature(func::Function)
    # Infer the number of arguments and return type of func

    # func can either be a function of time and state (t, state)
    # or a sole function of state (state,)

    types = Base.return_types(func, (Float64, State))
    nargs = 2
    if isempty(types)
        types = Base.return_types(func, (State,))
        nargs = 1
        if isempty(types)
            throw(ArgumentError("Function does not match signtaure (t, state) or (state,)"))
        end
    end

    return_type = Union{types...}
    if !isconcretetype(return_type)
        return_type = Any
    end

    return (nargs, return_type)
end

Base.size(s::StateData) = size(s.data)
Base.getindex(s::StateData, i::Int) = s.data[i]
Base.IndexStyle(::StateData) = IndexLinear()

function Base.empty!(s::StateData)
    empty!(s.data)
    empty!(s.t)
    s
end

function Base.sizehint!(s::StateData, n)
    sizehint!(s.data, n)
    sizehint!(s.t, n)
    s
end


function solve!(
        datalist::AbstractVector{<:StateData}, problem::IBProblem;
        callback = (_,_)->nothing
    )
    t = problem.t
    Δt = step(t)

    next_save_times = map(datalist) do data
        itr = Iterators.Stateful(isempty(data.saveat) ? t : data.saveat)
        empty!(data)
        sizehint!(data, length(itr))
        itr
    end

    state = IBState(problem)

    for t_k in t
        advance!(state, problem, t_k)
        all(isfinite, state.CL) || break

        callback(t_k, state)

        for (data, t_saves) in zip(datalist, next_save_times)
            while !isempty(t_saves) && peek(t_saves) < t_k + Δt / 2
                popfirst!(t_saves)
                push!(data.t, t_k)
                push!(data.data, data.func(t_k, state))
            end
        end
    end

    state
end

function solve!(data::StateData, problem; kwargs...)
    solve!([data], problem; kwargs...)
end

solve(problem; kwargs...) = solve!(StateData[], problem; kwargs...)


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
