module Plotting

using Plots

using ..IBPM: IBProblem, IBState
using ..IBPM.Quantities


@recipe function f(problem::IBProblem, state::IBState, quantity; subgrids=missing)
    a = if ismissing(subgrids)
        quantity(state, PreferView())
    else
        quantity(state, PreferView(); subgrids)
    end
    (problem, a, quantity)
end

@recipe function f(
    problem::IBProblem,
    a::AbstractArray,
    quantity::Union{typeof.((vorticity, streamfunction))...}
)
    aspect_ratio --> :equal
    legend --> false
    color --> :seaborn_icefire_gradient

    seriestype --> :contourf
    levels --> 30
    linewidth --> 0

    subgrids = axes(a, 3)
    ranges = gridranges(quantity, problem; subgrids)

    x, y = last(ranges)
    xlims --> extrema(x)
    ylims --> extrema(y)

    # Iterate over each subdomain in reverse order so that the smallest is displayed on top
    for i in Iterators.reverse(subgrids)
        @series @views (ranges[i]..., transpose(a[:, :, i]))
    end
end

@recipe function f(
    problem::IBProblem, pt_vec::AbstractVector{<:AbstractMatrix}, ::typeof(body_points)
)
    for pts in pt_vec
        @series Shape(pts[:, 1], pts[:, 2])
    end
end

end # module
