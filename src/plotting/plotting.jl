using Plots


"""
    fluidplot(x, y, f; ...)

Plot quantity `f` at discrete coordinates `x` and `y`.

Keywords are passed to the `Plots` package.

# Arguments
- `x`: Vector of x coordinates.
- `y`: Vector of y coordinates.
- `f`: Callable `f(x, y)` or array of the quantity to display.
"""
@userplot FluidPlot

@recipe function f(p::FluidPlot)
    if length(p.args) != 3
        throw(ArgumentError("fluidplot expects 3 arguments (x, y, f)"))
    end
    x, y = p.args

    tick_direction --> :out
    aspect_ratio --> :equal
    grid --> false

    seriestype --> :heatmap
    seriescolor --> :seaborn_icefire_gradient

    xlims --> extrema(x)
    ylims --> extrema(y)

    p.args
end

"""
    bodyplot(points)

Plot bodies given the return value of [`BodyPoints`](@ref).

Keywords are passed to the `Plots` package.
"""
@userplot BodyPlot

@recipe function f(p::BodyPlot)
    if length(p.args) != 1
        throw(ArgumentError("bodyplot expects 1 argument: points"))
    end
    points = first(p.args)

    seriescolor --> :gray
    linecolor --> :gray
    legend --> false

    for xy in (ndims(points) == 1 ? points : (points,))
        @series Shape(xy[:, 1], xy[:, 2])
    end
end
