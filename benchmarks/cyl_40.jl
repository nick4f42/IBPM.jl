include("../src/IBPM.jl")

# Define grid
xlims = (-1.0, 3.0)
ylims = (-2.0, 2.0)

mg = 1   # num domains
Δx = 0.02
grid =  IBPM.MultiGrid(Δx, (xlims, ylims), mg=mg)

# Other parameters
Re = 40.0
Δt = 1e-2

Uinf = 1.0;   # Free-stream flow
r = 0.5; # Cylinder radius

cyls = [IBPM.make_cylinder( r, grid.h, 0.0, 0.0 )]

 #freestream conditions
freestream = (Ux=t->1.0, Uy=t->0.0, inclination=t->0.0)

function run_sim!(t, state, prob; output=1, callback=(state, prob)->nothing)
	for i=1:length(t)
		IBPM.advance!(state, prob, t[i])
        if mod(i,output) == 0
			callback(state, prob);  # Primitive callback, can be used for plotting or other output
            @show (t[i], state.CD, state.CL, state.cfl)
        end
	end
end

prob = IBPM.IBProblem(grid, cyls, Δt, Re, freestream=freestream);

state = IBPM.IBState(prob);
T = 10
t = 0:Δt:T

run_sim!(t[1:2], state, prob) # Pre-compile

# Advance to final time
runtime = @elapsed run_sim!(t, state, prob; output=20)
println(runtime)
