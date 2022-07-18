using IBPM
using IBPM: Bodies
using IBPM.Quantities

using Printf
using Plots

xlim = (-1.0, 2.0) # x bounds
ylim = (-2.0, 2.0) # y bounds
dx = 0.01          # Grid step size
mg = 5             # Number of sub-domains
mgrid = MultiGrid(dx, (xlim, ylim); mg=mg)

r = 0.5 # cylinder radius
spacing = 1.0 # spacing between the cylinders
y = spacing / 2 + r # y coordinate of cylinder center
bodies = [Bodies.cylinder((0.0, y), r, dx) for y in (y, -y)]

# Freestream velocity
freestream(t) = (1.0, 0.0)

Re = 200.0 # Reynolds number
dt = 0.002 # Time step size
T = 24.0   # Final time

# Specify the problem using the grid, bodies, freestream, etc
prob = IBProblem(mgrid, bodies, freestream, Re=Re, dt=dt)

# Create functions vorticity(state) and bodypoints(state)
vorticity = Vorticity(prob)
bodypoints = BodyPoints(prob)

# Initialize animation
anim = Animation()

# Create callback that adds to animation
save_anim = at_times(range(0, T, length=240)) do state
    # Plot the vorticity and bodies
    xs = range(-3, 10, step=dx)
    ys = range(-3, 3, step=dx)
    fluidplot(xs, ys, vorticity(state); clims=(-5, 5))
    bodyplot!(bodypoints(state))

    # Save an animation frame
    frame(anim)
end

# Create callback to show progress
show_progress = each_timestep() do state
    @printf "\r%.2f%%" (100 * state.t / T)
end

# Solve the problem and update the animation
solve(prob, (0.0, T); call=[save_anim, show_progress])

# Save the animation to disk
gif(anim, "$(@__DIR__)/double_cyl.gif")
