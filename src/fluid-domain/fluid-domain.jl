"""
A collection of types of grids used to discretize and solve the linear system.
Types of flow grids include:
- MultiGrid
- TODO: Quadtree
"""
abstract type Grid end

"""
    MultiGrid <: Grid

A type of Grid which uses a hierarchy of spacial discretizations to recursively solve until the solution converges. 
While faster than brute force solving large linear systems, MultiGrid typically uses relatively high memory allocation.

# Fields
- `nx`: Number of x points in discretized domain.
- `ny`: Number of y points in discretized domain.
- `nΓ`: Number of ciruclation points in discretized domain.
- `nq`: Number of vel flux points in discretized domain.
- `mg`: Number of domains in MultiGrid.
- `offx::Float64`: Offset of left border of domain.
- `offy::Float64`: Offset of bottom border of domain.
- `len::Float64`: Length of domain in x direction.
- `h::Float64`: Grid cell size.
- `split_flux::Any`: Return views to 2D array of fluxes. #THIS MAY BE UNCELAR, WHAT DOES THIS DO?#
- `LEFT::Int`: Constant offset for indexing left boundary condition.
- `RIGHT::Int`: Constant offset for indexing right boundary condition.
- `BOT::Int`: Constant offset for indexing bottom boundary condition.
- `TOP::Int`: Constant offset for indexing top boundary condition.
"""
struct MultiGrid <: Grid
    nx::Int
    ny::Int
    nΓ::Int
    nq::Int
    mg::Int
    offx::Float64
    offy::Float64
    len::Float64
    h::Float64
    split_flux::Any
    LEFT::Int
    RIGHT::Int
    BOT::Int
    TOP::Int
end

"""
    MultiGrid(
        h::Float64,
        boundary::NTuple{2, Tuple{Float64, Float64}};
        mg::Int = 1
    )

Constructor for struct MultiGrid. Returns an instance of struct MultiGrid, see [`Main.IBPM.MultiGrid`](@ref).

# Arugments
- `h::Float64`: Grid spacing.
- `boundary::NTuple{2, Tuple{Float64, Float64}}`: Tuple defining the boundaries of the domain.
- `mg::Int`: Optional. Number of domains in MultiGrid. Default: `1`.
"""
function MultiGrid(
        h::Float64,
        boundary::NTuple{2, Tuple{Float64, Float64}};
        mg::Int = 1
    )
    #back out variables that the software needs from user defined vars
    xlims, ylims = boundary
    offx = -xlims[1]
    offy = -ylims[1]
    xlen = xlims[2] - xlims[1]
    ylen = ylims[2] - ylims[1]

    nx = round(Int, xlen / h)
    ny = round(Int, ylen * nx / xlen)

    nΓ  = (nx-1)*(ny-1)  # Number of circulation points

    nu = ny * (nx+1); nv = nx * (ny+1);  # num of (flux) points
    # Total num of vel (flux) points
    nq = nu + nv;

    h = xlen / nx;  # Grid spacing

    "Return views to 2D arrays of fluxes"
    split_flux(q; lev=1) = reshape(@view(q[1:nu, lev]), nx+1, ny),
                           reshape(@view(q[nu+1:end, lev]), nx, ny+1)

    # Predefine constant offsets for indexing boundary conditions
    left = 0;  right = ny+1
    bot = 2*(ny+1); top = 2*(ny+1) + nx+1

    return MultiGrid(nx, ny, nΓ, nq, mg, offx, offy, xlen, h,
        split_flux, left, right, bot, top)
end

"""
    gridstep(grid::MultiGrid)

Returns the grid spacing for the passed in Grid, to be used in IBProblem (see [`Main.IBPM.IBProblem`](@ref)).
"""
gridstep(grid::MultiGrid) = grid.h