using IBPM
using IBPM: Bodies
using IBPM.Quantities

using Printf
using Plots
using ProgressMeter

xlim = (-1.5, 6.5) # x bounds
ylim = (-2.0, 2.0) # y bounds
dx = 0.02          # Grid step size
mg = 3             # Number of subdomains
mgrid = MultiGrid(dx, (xlim, ylim), mg=mg)

# Specify motion
ω = 2π * 0.1               # Pitch frequency
A = deg2rad(40.0)          # Pitch amplitude
θ(t) = -A * sin(ω * t)     # Angular position
Ω(t) = -ω * A * cos(ω * t) # Angular velocity
U(t) = 1.0                 # x velocity
V(t) = 0.0                 # y velocity
motion = IBPM.MovingGrid(U, V, θ, Ω)

# Create airfoil body
x0 = 0.25
nb = 48 # Number of body points
bodies = [Bodies.naca_airfoil(x0, nb, "0012", motion=motion)]

# Other parameters
Re = 200.0         # Reynolds number
dt = 1e-3          # Time step size
T = 2.0 * (2π / ω) # Final time

# Specify the problem using the grid, bodies, etc
prob = IBProblem(mgrid, bodies; Re=Re, dt=dt)

# Create functions vorticity(state) and bodypoints(state)
vorticity = Vorticity(prob; frame=lab_frame(prob))
bodypoints = BodyPoints(prob)

# Initialize animation
anim = Animation()

# Create callback that adds to animation
save_anim = at_times(range(0, T, length=240)) do state
    # Plot the vorticity and bodies
    xs = range(-1.5, 6.5, step=dx)
    ys = range(-2.0, 2.0, step=dx)
    fluidplot(xs, ys, vorticity(state); clims=(-5, 5))
    bodyplot!(bodypoints(state))

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
gif(anim, "$(@__DIR__)/pitching_naca.gif")
