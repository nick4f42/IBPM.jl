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
    grid_of(x)

The [`Grid`](@ref) cooresponding to a problem or model.
"""
function grid_of end

"""
    model_of(problem)

The [`SolnModel`](@ref) corresponding to `problem`.
"""
function model_of end

"""
    scheme_of(problem)

The [`ExplicitScheme`](@ref) corresponding to `problem`.
"""
function scheme_of end

"""
    motions_of(x)

The [`Motion`](@ref) subtype corresponding to a problem or model.
"""
function motions_of end

"""
    gridstep(x)

The minimum grid step size of a problem, model, or grid.
"""
function gridstep end

"""
    timestep(x)

The time step size of a problem or scheme.
"""
function timestep end

grid_of(prob::AbstractIBProblem) = grid_of(model_of(prob))
timestep(prob::AbstractIBProblem) = timestep(scheme_of(prob))

"""
    IBProblem <: AbstractIBProblem

Initialize the problem structure (matrices used, bodies and simulation
parameters, time steppping scheme, ...)

Note: the scheme actually speaks to the terms that are explicitly treated. This
is a projection method the directly enforces the no-slip condition, so some terms
are implicitly treated. This information is not contained in scheme, but in the
`A`, `Ainv`, `B`, and `Binv` matrices.

---

    IBProblem(
        grid::Grid,
        bodies::Vector{<:Body{M}} where M,
        freestream = t -> (0.0, 0.0);
        Re::Float64,
        dt::Float64,
    )

# Arguments

- `grid`: Discretization of the fluid domain.
- `bodies`: The bodies to simulate.
- `freestream`: The freestream velocity as a function of time. `freestream(t)` returns a
        tuple of each velocity component at `t`. Defaults to zero.
- `Re`: Reynolds number.
- `dt`: The time stepsize.
"""
mutable struct IBProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    A
    Ainv
    Binv
    function IBProblem(
        grid::Grid,
        bodies::Vector{<:Body{M}} where M,
        freestream = t -> (0.0, 0.0);
        Re::Float64,
        dt::Float64,
    )
        model = IBModel(grid, bodies, Re; freestream=freestream)
        scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
        A, Ainv, Binv = get_AB(model, dt)
        new(model, scheme, A, Ainv, Binv)
    end
end

model_of(prob::IBProblem) = prob.model
scheme_of(prob::IBProblem) = prob.scheme
motions_of(prob::IBProblem) = motions_of(prob.model)

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
