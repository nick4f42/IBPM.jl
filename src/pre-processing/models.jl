"""
From a method-of-lines perspective, this should contain everything related
    to the continuous-time problem... analogous to defining the RHS of an ODE.
    Then everything related to time-stepping goes in the IBProblem

At the moment, the distinction between this and IBProblem is not useful,
    but could potentially be useful for interfacing with DifferentialEquations
"""

"""
Matrices that can be precomputed

Correspondence with Taira & Colonius (2007)
    Note that not all matrices defined here are explicitly constructed
C  - Basic curl operator for single-grid
        Call curl! function for multigrid to take into account boundary conditions
C' - Transforms velocity flux to circulation: Γ = C'*q
E  - Maps fluxes to body motion, i.e. u_B = E*q
        Note that H = -E' is the regularization operator
A  - Implicit time-stepping operator for velocity flux
        A = I - (dt/2/h^2)*Lap

NOTE: Most of the functionality for this is in fluid-operators/lin.jl

# Constructor
    IBMatrices(grid::Grid, bodies::Vector{<:Body})

# Arguments
- `grid::Grid`: Grid struct of type Grid which defines and discretizes the domain.
- `bodies::Array{Body, 1}`: 1D array of bodies created for simulation.
"""
mutable struct IBMatrices
    C::LinearMap  # Basic curl operator for single-grid in form of linear map.
    Λ::Array{Float64, 2}  # Laplcacian eigenvalues.
    Δinv::LinearMap  # Inverse of Laplacian in form of linear map.
    E::LinearMap  # Interpolation/regularization matrix in form of linear map.
    dst_plan::Any  # Plans needed for optimized discrete sine transform (DST) inversion.

    function IBMatrices(grid::T, bodies::Array{V, 1}) where T <: Grid where V <: Body
        mats = new()
        Γwork = zeros(grid.nx-1, grid.ny-1)
        mats.C = LinearMap( (q, ψ) -> IBPM.curl!(q, ψ, grid),       # Forward
                            (Γ, q) -> IBPM.rot!(Γ, q, grid, Γwork), # Transpose
                            grid.nq, grid.nΓ; ismutating=true)
        #mats.Lap = mats.C'*mats.C/Re    # Laplacian
        mats.Λ = lap_eigs(grid)

        # Plan DST for inverse Laplacian
        #  Have to do this separately for every grid level because of FFT alignment issue
        Γtemp = ones(Float64, grid.nΓ, grid.mg)
        mats.dst_plan = get_dst_plan(ones(grid.nx-1, grid.ny-1))
        mats.Δinv = get_lap_inv(grid, mats.Λ, mats.dst_plan)

        # Interpolation/regularization matrix
        mats.E = IBPM.setup_reg(grid, bodies)   # interface-coupling/interface-oupling.jl
        return mats
    end
end

"""
A collection of different explicit time stepping schemes.
Types of explicit schemes include:
- AdamsBashforth
- TODO: RungeKuttaChebyshev

Is passed as an arugment into IBProblem, see [`Main.IBPM.IBProblem`](@ref).
"""
abstract type ExplicitScheme end

"""
    AdamsBashforth <: ExplicitScheme

The group of multistep methods known as Adams-Bashforth methods. Currently supports:
- Two step Adams-Bashforth Method
"""
struct AdamsBashforth <: ExplicitScheme
    dt::Float64
    β::Array{Float64, 1}
end

timestep(ab::AdamsBashforth) = ab.dt

# TODO: Make AdamsBashforth constructor to generate β automatically

"""
    AB2(dt::Float64)

Special case constructor for type::AdamsBashforth for second-order Adams-Bashforth scheme.

# Arguments

- `dt::Float64`: Time stepping interval.
"""
function AB2(dt::Float64)
    return AdamsBashforth(dt, [1.5, -0.5])
end


"""
Pre-allocate memory to certain vectors that can be re-used throughout the
computation process.

# Constructor
    WorkingMemory(grid::Grid)

# Arguments
- `grid::Grid`: Grid struct of type Grid which defines and discretizes the domain.
"""
mutable struct WorkingMemory
    q1::Array{Float64, 2}  # Trial flux qs.
    q2::Array{Float64, 1}  # boundary_forces.
    q3::Array{Float64, 1}  # Nonlinear.
    q4::Array{Float64, 1}  # Nonlinear.
    q5::Array{Float64, 1}  # Nonlinear.
    q6::Array{Float64, 1}  # Linear stability analysis.
    Γ1::Array{Float64, 2}  # Trial circulation Γs.
    Γ2::Array{Float64, 1}  # trial_state, project_circ.
    Γ3::Array{Float64, 1}  # vort2flux.
    Γbc::Array{Float64, 1}  # Poisson boundary conditions for multigrid
    rhsbc::Array{Float64, 1}  # Time-stepping boundary conditions for multigrid
    function WorkingMemory(grid::Grid)
        work = new()
        work.q1 = zeros(grid.nq, grid.mg)  # Trial flux qs
        work.Γ1 = zeros(grid.nΓ, grid.mg)  # Trial circulation Γs
        work.q2 = zeros(grid.nq) # boundary_forces
        work.Γ2 = zeros(grid.nΓ) # trial_state, project_circ
        work.q3 = zeros(grid.nq) # nonlinear
        work.q4 = zeros(grid.nq) # nonlinear
        work.q5 = zeros(grid.nq) # nonlinear
        work.q6 = zeros(grid.nq) # Used for linear stability analysis
        work.Γ3 = zeros(grid.nΓ) # vort2flux
        work.Γbc = zeros(2*(grid.nx+1)+2*(grid.ny+1))
        work.rhsbc = zeros(grid.nΓ)
        return work
    end
end


"""
SolnModel contains information about simulation parameters and stores
all static (non-time varying) matrices.
"""
abstract type SolnModel end

"""
    IBModel{T <: Grid, V <: Body} <: SolnModel

Contains information about simulation parameters and stores static matrices. 
Is created using a type of grid, simulation bodies, Reynolds number, and freestream velocity.
Is passed as an arugment into IBProblem, see [`Main.IBPM.IBProblem`](@ref).

# Constructor

    IBModel(
        grid::Grid,
        bodies::Array{Body, 1},
        Re::Number;
        freestream::Function,
        xc=0.0,
        yc=0.0
    )

# Arguments
- `grid::Grid`: Grid struct of type Grid which defines and discretizes the domain.
- `bodies::Array{Body, 1}`: Array of bodies.
- `Re::Number`: Reynolds number.
- `freestream::Function`: Free-stream velocity.
- `xc=0.0`: x coordinate offset for calculating rotational fluxes.
- `yc=0.0`: y coordinate offset for calculating rotational fluxes.
"""
struct IBModel{G<:Grid, B<:Body} <: SolnModel
    grid::G
    bodies::Vector{B}           # Array of bodies
    Re::Float64                 # Reynolds number
    freestream::Function        # Free-stream velocity
    mats::IBMatrices            # Various precomputed sparse matrices
    work::WorkingMemory
    XX::Union{Array{Float64, 2}, Nothing}  #  x-locations for computing rotational fluxes
    YY::Union{Array{Float64, 2}, Nothing}  #  y-locations for computing rotational fluxes
    function IBModel(grid::G,
                     bodies::Vector{B},
                     Re::Number;
                     freestream::Function,
                     xc=0.0,
                     yc=0.0) where {G<:Grid, M<:Motion, B<:Body{M}}
        mats = IBMatrices(grid, bodies)
        work = WorkingMemory(grid)

        # TODO: Put in different function??
        "Precompute grid locations for cross products (used for rotational flux)"
        if M == MovingGrid
            # Define x-xc, y-yc on all grids
            @assert length(bodies) == 1
            motion = bodies[1].motion
            nx, ny, h = grid.nx, grid.ny, grid.h
            XX = zeros(nx*(ny+1), grid.mg)  # Number of y-fluxes
            YY = zeros(ny*(nx+1), grid.mg)  # Number of x-fluxes

            for lev=1:grid.mg
                hc = h*2^(lev-1);  # Coarse grid spacing

                ### y-coordinates for calculating x-fluxes
                y = @. ((1:ny)-0.5-ny/2)*hc + ny/2*h - grid.offy
                YY[:, lev] = (ones(nx+1)*(y.-yc)')[:]

                ### x-coordinates for calculating y-fluxes
                x = @. ((1:nx)-0.5-nx/2)*hc + nx/2*h - grid.offx
                XX[:, lev] = ((x.-xc)*ones(ny+1)')[:]
            end
        else
            XX, YY = nothing, nothing
        end
        return new{G,B}(grid, bodies, Float64(Re), freestream, mats, work, XX, YY)
    end
end

grid_of(model::IBModel) = model.grid
motions_of(model::IBModel) = (body.motion for body in model.bodies)
