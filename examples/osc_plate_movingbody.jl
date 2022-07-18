using IBPM
using IBPM: Bodies
using IBPM.Quantities

using Printf
using Plots
using ProgressMeter

xlim = (-4.0, 4.0) # x bounds
ylim = (-2.0, 2.0) # y bounds
dx = 0.02          # Grid step size
mg = 3             # Number of subdomains
mgrid = MultiGrid(dx, (xlim, ylim); mg=mg)

# Specify motion
xc(t) = cos(t)  # x position
yc(t) = 0.0     # y potision
θ(t)  = 0.0     # Angular position
uc(t) = -sin(t) # x velocity
vc(t) = 0.0     # y velocity
ω(t)  = 0.0     # Angular velocity
motion = IBPM.MotionFunction([xc, yc, θ], [uc, vc, ω])

# Create plate body
L = 1.0          # Length
a = deg2rad(90)  # Angle with x-axis
x, y = 0.0, -L/2 # Coordinates of plate edge
bodies = [Bodies.plate((x, y), L, a, dx; motion=motion)]

# Other parameters
Re = 200.0   # Reynolds number
dt = 1e-3    # Time step size
T = 1.0 * 2π # Final time

# Specify the problem using the grid, bodies, etc
prob = IBProblem(mgrid, bodies; Re=Re, dt=dt)

# Create functions vorticity(state) and bodypoints(state)
vorticity = Vorticity(prob)
bodypoints = BodyPoints(prob)

# Initialize animation
anim = Animation()

# Create callback that adds to animation
save_anim = at_times(range(0, T, length=240)) do state
    # Plot the vorticity and bodies
    xs = range(-4.0, 4.0, step=dx)
    ys = range(-2.0, 2.0, step=dx)
    fluidplot(xs, ys, vorticity(state); clims=(-5, 5))
    bodyplot!(bodypoints(state); linewidth=3)

    # Save an animation frame
    frame(anim)
end

# Create callback to show progress using the ProgressMeter package
progress = Progress(round(Int, T/dt) + 1)
show_progress = each_timestep() do state
    next!(progress)
end

# Solve the problem and update the animation
solve(prob, (0.0, T); call=[save_anim, show_progress])

# Save the animation to disk
gif(anim, "$(@__DIR__)/osc_plate_movingbody.gif")
