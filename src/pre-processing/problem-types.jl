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

---

    IBProblem(
        grid::Grid,
        bodies::Vector{<:Body},
        t_span::NTuple{2,Float64},
        [dt::Float64];
        Re::Float64,
        freestream::Function = t -> (0.0, 0.0)
    )

# Arguments

- `grid`: Discretization of the fluid domain.
- `bodies`: The bodies to simulate.
- `t_span`: The starting and ending time.
- `dt`: The time step. Auto-determined by default.
- `Re`: Reynolds number.
- `freestream`: The freestream velocity as a function of time. `freestream(t)` returns a
        tuple of each velocity component at `t`. Defaults to zero.

By default, `dt` is chosen to aim for a CFL of 0.1 with a saftey factor of 5 on the max
velocity.
"""
mutable struct IBProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    t::typeof(range(0.0,1.0))
    A
    Ainv
    Binv
    function IBProblem(
            grid::Grid,
            bodies::Vector{<:Body},
            t_span::NTuple{2,Real},
            dt::Float64;
            Re::Float64,
            freestream::Function = t -> (0.0, 0.0)
        )
        model = IBModel(grid, bodies, Re; freestream=freestream)
        scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
        A, Ainv, Binv = get_AB(model, dt)

        t = range(t_span..., step=dt)
        new(model, scheme, t, A, Ainv, Binv)
    end
end

function IBProblem(
        grid::Grid, bodies::Vector{<:Body}, t_span::NTuple{2,Real};
        Re::Float64, freestream::Function
    )
    # TODO: 5000 time steps to find max is sort of arbitrary. Consider changing?
    U_max = maximum(t_k -> hypot(freestream(t_k)...), LinRange(t_span..., 5000))

    # satisfy CFL 0.2 constraint, with safety factor on max velocity
    Δx = gridstep(grid)
    dt = 0.1 * Δx / (5.0 * U_max)

    IBProblem(grid, bodies, t_span, dt; Re, freestream)
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
