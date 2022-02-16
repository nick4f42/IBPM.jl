module ibpm

using LinearAlgebra
using SparseArrays
using FFTW
using LinearMaps
using IterativeSolvers
using InplaceOps  # @! macro
using Plots

export IBPM_advance

#Caution, include order matters!
include("fluid-domain/fluid-domain-include.jl")
include("structure-domain/structure-domain-include.jl")
include("interface-coupling/interface-coupling-include.jl")
include("pre-processing/pre-processing-include.jl")
include("fluid-operators/fluid-operators-include.jl")
include("timestepping/timestepping-include.jl")
include("plotting/plotting-include.jl")
include("experimental/sfd.jl")

"""
    IBPM_advance(Re, boundary::Tuple, body, freestream, Δx, mg, Δt, T, save_info)
 
Convenience function to solve the full problem and plot final solution
   for flow over a single cylinder
 
For more control, just use this as a template - see benchmarks and examples
 
# Arguments
- `Re::Any`: Reynolds Number.
- `boundary::Tuple`: Describes the boundaries of the interesed simulation domain. Defined by a tuple of format: (left, right, bottom, top).
- `body::Vector{NamedTuple{}}`: Defines simulated bodies. See [`make_body`](@ref).
- `freestream`: Optional. Describes the incomming freestream. Defined by a named typle of format: (Ux=f(t), Uy=f(t), inclination=f(t))
- `Δx`: Optional. Grid sizing. Any number type.
- `mg`: Optional. Number of subdomains in each grid block. Default is 5.
- `Δt`: Optional. Timestep interval. Choose Δt such that CFL<1 for stable/accurate simulation solution.
- `T`: Optional. Total simulation time. Default is 20*dt. <<<dt is not defined anywhere???>>>
- `save_info`: Optional. See [`init_save_data`](@ref).
 
# Returns
- `prob`: Stucture of type AbstractIBProblem. Defines what type of problem is being solved. Can be structs IBProblem, LinearizedIBProblem, SFDProblem. 
        See [`IBProblem`](@ref).
- `data`: Solution data to simulation at defined time intervals. See [`IBState`](@ref).
- `runtime`: Total runtime of simulation.
"""
function IBPM_advance(Re, boundary, body, freestream=(Ux=t->0.0,);
    Δx=missing, mg=5, Δt=missing, T=20.0*dt, save_info=missing)

    #--extract user params to sim variables
        Δx, Δt, T, freestream = read_user_vars(Δt, Δx, freestream, Re, T)
    #--

    #--Build flow grid
        grid =  make_grid(Δx, boundary, mg=mg)
    #--

    #-build body
        #@aditya: made a template function make_body
        #(in structure-domain/bodies.jl) that makes a body using some flags
        #specified by the user. The function is a wrapper that calls routines
        #within sample-bodies.jl. make_body does not support calling make_naca
        #because
        #    (1) that function defines the airfoil using the number of body
        #       points as the user prescribed variable, rather than \Delta x
        #    (2) the function doesn't allow for AoA to be specified.
        #It would be nice to address (1) and (2) and let the user have a naca
        #airfoil be made by specifying the type key as :naca, and perhaps adding
        #new key, params, with subkeys spec (e.g. params.spec="0012" for a
        #NACA0012) and \alpha.
        #Would be nice to also modify the file to allow the user to
        #specify the body as a function of the arc-length s.
        body = make_body( body, Δx )
    #--

	#--Initialize problem types and models
	    prob = IBProblem(grid, body, Δt, Re, freestream=freestream);
	    state = IBState(prob);
	#--

	#Time over which simulation will be run
	t = 0:Δt:T

	#initialize user desired save variables as a NamedTuple called data
	data = init_save_data( t, save_info, state )

	#run simulation over desired time window
    # run_sim(t[1:2], state, prob) # Pre-compilation for benchmarking
    runtime = @elapsed run_sim!(t, state, prob, data=data) #advance to final time

    return prob, data, runtime
end

"""
    compute_cfl(state, prob)

Compute the CFL number (uΔt/Δx) based on the fine-grid flux

Note that this uses working memory that is also used in `nonlinear!`

# Arguments
- `state`: Stores values of properties of flow at single time instance. Of Struct IBState
- `prob`: Stucture of type AbstractIBProblem. Defines what type of problem is being solved. Can be structs IBProblem, LinearizedIBProblem, SFDProblem.

# Returns
- `cfl`: The Courant-Friedrichs-Lewy (CFL) number, used to verify the CFL condition.
"""
function compute_cfl(state, prob)
	Δt, Δx = prob.scheme.dt, prob.model.grid.h
	qwork = prob.model.work.q5
	@views @. qwork = abs( state.q[:, 1] )
	return maximum(qwork)*Δt/Δx
end

"""
    run_sim!(t, state, prob)

Computes the solution to the prescribed problem for a range of times. 

# Arguments
- `t`: An array of desired simulation times. 
- `state`: Stores values of properties of flow at single time instance. Of Struct IBState
- `prob`: Stucture of type AbstractIBProblem. Defines what type of problem is being solved. Can be structs IBProblem, LinearizedIBProblem, SFDProblem.

# Returns
- `cfl`: The Courant-Friedrichs-Lewy (CFL) number, used to verify the CFL condition.
"""
function run_sim!(t, state, prob;
	display_freq=25,
	data::Array{user_var, 1})
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
		if ~all(isfinite.(state.CL))
			break
		end
        if mod(i,display_freq) == 0
			state.cfl = compute_cfl(state, prob)
            @show (t[i], state.CD, state.CL, state.cfl)
        end

		#save stuff if desired
	    save_data!( t[i], t, prob, state, data  )
	end
end

end
