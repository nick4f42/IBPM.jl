"""
IBProblem has the info of IBModel as well as the problem structure (e.g., the
explicit time stepping scheme and information about the implicit treatment via
the A and B matrices and their inverses).

Maybe some opportunity for restructuring...

Looking towards possible compatibility with DifferentialEquations.jl, this
        would be similar to the ODEProblem
"""
abstract type AbstractIBProblem end

"""
    IBProblem <: AbstractIBProblem

Initialize the problem structure (matrices used, bodies and simulation
parameters, time steppping scheme, ...)

Note: the scheme actually speaks to the terms that are explicitly treated. This
is a projection method the directly enforces the no-slip condition, so some terms
are implicitly treated. This information is not contained in scheme, but in the
`A`, `Ainv`, `B`, and `Binv` matrices.

# Constructor
    IBProblem(
        grid::Grid,
        bodies::Array{<:Body, 1},
        t::Union{AbstractRange, Tuple},
        Re::Float64;
        freestream::Function = t -> (0.0, 0.0)
    )

# Arguments
- `grid::Grid`: Grid struct of type Grid which defines and discretizes the domain.
- `bodies::Array{<:Body, 1}`: 1D array of bodies created for simulation.
- `t::Union{AbstractRange, Tuple}`: Time stepping interval. Can be passed in as an AbstractRange or a Tuple of the start and end time.
- `Re::Float64`: Reynolds number.
- `freestream::Function`: Optional. Freestream velocity. Default: `t -> (0.0, 0.0)`.

# Fields
- `model::IBModel`: Contains information about simulation parameters and stores static matrices. 
- `scheme::ExplicitScheme`: Explicit time stepping scheme.
- `t::AbstractRange`: Time steps for simulation.
- `A`: A matrix in form of LinearMap.
- `Ainv`: A inverse matrix in form of LinearMap.
- `Binv`: B inverse matrix in form of LinearMap.
"""
mutable struct IBProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    t::typeof(range(0.0,1.0))
    A
    Ainv
    Binv
    function IBProblem(
            grid::Grid, bodies::Vector{<:Body}, t::AbstractRange{Float64};
            Re::Float64, freestream::Function
        )
        model = IBModel(grid, bodies, Re; freestream=freestream)
        scheme = AB2(step(t))   # Explicit time-stepping for nonlinear terms
        A, Ainv, Binv = get_AB(model, step(t))

        new(model, scheme, t, A, Ainv, Binv)
    end
end

function IBProblem(
        grid::Grid, bodies::Vector{<:Body}, t_range::Tuple{Float64, Float64};
        Re::Float64, freestream::Function
    )
    # TODO: 5000 time steps to find max is sort of arbitrary. Consider changing?
    U_max = maximum(t_k -> hypot(freestream(t_k)...), LinRange(t_range..., 5000))

    # satisfy CFL 0.2 constraint, with safety factor on max velocity
    Δx = gridstep(grid)
    Δt = 0.1 * Δx / (5.0 * U_max)
    t = range(t_range..., step=Δt)

    IBProblem(grid, bodies, t; Re, freestream)
end

gridstep(problem::IBProblem) = gridstep(problem.model.grid)
timestep(problem::IBProblem) = step(problem.t)

"""
    LinearizedIBProblem <: AbstractIBProblem

Create a linearized problem from the base_state and associated IBProblem

NOTE: The `freestream` value in the nonlinear IBProblem will be incorrect
for the linearized case, but this field is not used in
base_flux!(..., prob::LinearizedIBProblem, ...)

Modified IBProblem to include base state.  Only modification to the code
is the direct product called by the `nonlinear!` function

# Constructor
    LinearizedIBProblem(
        base_state::IBState,
        base_prob::IBProblem,
        dt::Float64
    )

# Arugments
- `base_state::IBState`: 
- `base_prob::IBProblem`: 
- `dt::Float64`: Time stepping interval.

# Fields
- `model::IBModel`: Contains information about simulation parameters and stores static matrices. 
- `scheme::ExplicitScheme`: Explicit time stepping scheme.
- `base_state::IBState`: 
- `QB::Array{Float64, 2}`: 
- `ΓB::Array{Float64, 2}`: 
- `A`: A matrix in form of LinearMap.
- `Ainv`: A inverse matrix in form of LinearMap.
- `Binv`: B inverse matrix in form of LinearMap.
"""
mutable struct LinearizedIBProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    base_state::IBState
    QB::Array{Float64, 2}
    ΓB::Array{Float64, 2}
    A
    Ainv
    Binv
    function LinearizedIBProblem(
                    base_state::IBState,
                    base_prob::IBProblem,
                    dt::Float64)
        prob = new()

        prob.model = deepcopy(base_prob.model)
        prob.scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
        prob.A, prob.Ainv, prob.Binv = get_AB(prob.model, dt)
        prob.base_state = deepcopy(base_state)

        # Averaged base flux used in mean flow advection
        #  Note this calls the IBProblem version of avg_flux, not Linearized
        #  so that the "background" or free-stream flux is accounted for
        prob.QB = zeros(prob.model.grid.nq, prob.model.grid.mg)
        for lev=1:prob.model.grid.mg
            prob.QB[:, lev] = copy(avg_flux(base_state, base_prob; lev=lev))
        end
        prob.ΓB = copy(base_state.Γ)

        return prob
    end
end
