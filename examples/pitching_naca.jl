include("../src/IBPM.jl")

# using .IBPM
# using .IBPM.Quantities: vorticity, body_points

using Plots

xlims = (-1.5, 6.5)
ylims = (-2.0, 2.0)
mg = 3   # num domains
Δx = 0.02
mgrid = IBPM.MultiGrid(Δx, (xlims, ylims), mg=mg)

# Other parameters
Re = 200.0
Δt = 1e-3 * 3

# Initialize motion
ω = 2π*0.1          # Pitch frequency
A = 40.0 * π/180.0  # Pitch amplitude, degrees
θ(t) = -A*sin(ω*t)
θ̇(t) = -ω*A*cos(ω*t);
U(t) = 1.0
V(t) = 0.0
motion = IBPM.MovingGrid(U, V, θ, θ̇)

# Create airfoil
x0 = 0.25
nb = 48;  # Number of body points
bodies = [IBPM.make_naca(x0, nb, "0012", motion=motion)]

# T = 2.0*(2π/ω)
T = 1.5*(2π/ω)

prob = IBPM.IBProblem(mgrid, bodies, (0, T), Δt; Re)

anim = Animation()

boundary = (xlims, ylims)
discretization = (20, 20)
particlePlot = IBPM.ParticlePlot(discretization, mgrid)
save_anim = IBPM.StateCallback(prob; at=range(0, T, length=200)) do t, state
    plot(prob, state, IBPM.Quantities.vorticity; subgrids=1:1, clims=(-5, 5))
    IBPM.plotParticles(T/(200-1), particlePlot, state, mgrid)
    plot!(prob, state, IBPM.Quantities.body_points; color=:black)
    frame(anim)
end

IBPM.solve(prob, [save_anim]) do t, state
    percent = round(100*t/T, digits=1)
    print("solving... ", percent, "%\r")
end

gif(anim, "IBPM.jl/examples/pitching_naca_pp_12.gif", fps=30)