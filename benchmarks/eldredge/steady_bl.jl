### Steady flat plate boundary layer
include("config.jl")  # Set up grid and other variables
#using BSON: @save, @load

Δt = 1e-2
T = 100.0
t = 0:Δt:T

# Create plate
bodies = [IBPM.make_plate( L, α, grid.h, x0, y0; n=nb)]

prob = IBPM.IBProblem(grid, bodies, Δt, Re, freestream=freestream);

#@load "benchmarks/eldredge/steady.bson" state
state = IBPM.IBState(prob);

run_sim!(t, state, prob, output=20) #advance to final time

#@save "benchmarks/eldredge/results/steady.bson" state
