include("../src/IBPM.jl")

import Plots
using Plots: plot, @animate, gif

using .IBPM
using .IBPM.Quantities: lift_coef

Δx = 0.02  # Grid step size
boundary = ((-1.0, 3.0), (-2.0, 2.0))  # x and y boundaries of domain
gridcount = 5  # Number of sub-domains
grid = MultiGrid(Δx, boundary; mg=gridcount)

r = 0.5
bodies = [ IBPM.make_cylinder(r, Δx, 0.0, 0.0) ]

Δt = 0.004  # Time step size
T = 10.0  # Final time
times = range(0, T, step=Δt)

# Δt can be autodetermined by supplying a tuple of initial and final times.
# The default aims for a CFL of 0.1 with a 5x safety factor on max velocity.
# times = (0.0, 1.0)

Re = 100.0  # Reynolds number
freestream = t -> (1.0, 0.0)  # Freestream velocity as a function of time

# lift_coef returns the lift coefficient for each body, so take the first
lift = StateData(state -> lift_coef(state)[1])
# By default, the entire state is saved
states = StateData(saveat=LinRange(0, T, 31))

problem = IBProblem(grid, bodies, times; Re, freestream)
solve!([lift, states], problem)

anim = @animate for (t, state) in zip(states.t, states)
    print("t = ", t, '\r')
    IBPM.plot_state(
        problem, state, t, var=:omega,
        xlims=(-4.0, 10.0), ylims=(-3.0, 3.0), clims=(-5.0, 5.0), clevs=40
    )
end

gif(anim, "cyl_100.gif", fps=10)
