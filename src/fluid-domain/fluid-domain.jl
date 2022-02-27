"""
Types of flow grids

UniformGrid can probably be eliminated once the code is fairly stable

MOVE TO fluid-domain??
"""
abstract type Grid end

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

    nx = Int64(round(xlen/h))
    ny = Int64(round(ylen*nx/xlen))

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
