using LinearAlgebra: norm  # FOR DEBUGGING


"""
    advance!(y::IBState, x::IBState, prob::AbstractIBProblem, t::Float64)

Advance state x forward in time and save in y
"""
function advance!(y::IBState,
                  x::IBState,
                  prob::T where T<:AbstractIBProblem,
                  t::Float64)
    copy!(y, x)  # Copy all fields over, can now mutate y
    advance!(y, prob, t)
end

"""
    advance(x::IBState, prob::AbstractIBProblem, t::Float64)

Advance state x forward in time (not mutating)
"""
function advance(x::IBState,
                 prob::T where T<:AbstractIBProblem,
                 t::Float64)
    y = similar(x)
    advance!(y, x, prob, t)
    return y
end

"""
    advance!(state::IBState, prob::AbstractIBProblem, t::Float64)

Advance state forward in time.
"""
function advance!(state::IBState,
                  prob::T where T<:AbstractIBProblem,
                  t::Float64)
    grid = prob.model.grid
    update_bodies!(state, prob, t)

    if MotionType(prob.model.bodies) == MovingGrid
        base_flux!(state, prob, t)
    end

    # Alias working memory for notational clarity
    #   This leaves work.Γ2, Γ3 and work.q1, q2 available
    qs = prob.model.work.q1  # Trial flux
    Γs = prob.model.work.Γ1  # Trial circulation

    #Computes trial circulation Γs and associated strmfcn and vel flux that
    #don't satisfy no-slip (from explicitly treated terms)
    get_trial_state!(qs, Γs, state, prob)

    # Update surface quantities to be able to trim off part of circ
    # that doesn't satisfy no slip
    @views boundary_forces!(state.F̃b, qs[:, 1], state.q0[:, 1], prob)
    update_stress!(state, prob) #Compute integral quantities and store in state

    # --update circulation , vel-flux, and strmfcn on fine grid
    #   to satisfy no-slip updates state.Γ, state.ψ, state.q
    project_circ!(Γs, state, prob)

    # Interpolate values from finer grid to center region of coarse grid
    vort2flux!( state.ψ, state.q, state.Γ, prob.model, grid.mg );

    return nothing
end


"""
    get_trial_state!(qs, Γs, state, prob)

Compute trial circulation Γs that doesn't satisfy no-slip BCs

Combine explicit Laplacian and nonlinear terms into a rhs
   then invert implicit part to return trial circulation Γs

High-level version of AB2:
rhs = A*Γ .-
      3*dt/2 * nonlin .+
      dt/2 * nonlin_prev .+
      dt/2 * rhsbc

Then do Ainv of that to back out trial circ
"""
function get_trial_state!(qs::AbstractArray,
                          Γs::AbstractArray,
                          state::IBState,
                          prob::AbstractIBProblem)
    dt = prob.scheme.dt
    grid = prob.model.grid
    rhsbc = prob.model.work.rhsbc
    #rhs = @view(work.Γ2[:, 1])  # RHS of discretized equation
    rhs = prob.model.work.Γ2  # RHS of discretized equation
    bc = prob.model.work.Γbc

    for lev=grid.mg:-1:1
        bc .*= 0.0; rhsbc .*= 0.0
        hc = grid.h * 2^( lev - 1);

        if lev < grid.mg
            @views get_bc!(bc, state.Γ[:, lev+1], grid)

            fac = 0.25*dt/ ( prob.model.Re * hc^2 )
            apply_bc!( rhsbc, bc, fac, grid )
        end

        #compute the nonlinear term for the current time step
        @views nonlinear!( state.nonlin[1][:, lev], state, bc, lev, prob );

        @views mul!( rhs, prob.A[lev], state.Γ[:, lev] )

        for n=1:length(prob.scheme.β)
            rhs .+= (prob.scheme.β[n]*dt)*@view(state.nonlin[n][:, lev])
        end

        # Include boundary conditions
        rhs .+= rhsbc

        # Trial circulation  Γs = Ainv * rhs
        @views mul!(Γs[:, lev], prob.Ainv[lev], rhs)
    end

    # Store nonlinear solution for use in next time step
    state.nonlin[2] .= state.nonlin[1]

    vort2flux!( state.ψ, qs, Γs, prob.model, grid.mg )
    return nothing
end

"""
    boundary_forces!(F̃b, qs, q0, prob)

Solve the Poisson equation (25) in Colonius & Taira (2008).

Dispatch based on the type of motion in the problem - allows precomputing
    regularization and interpolation where possible.
"""
function boundary_forces!(F̃b, qs, q0, prob)
    boundary_forces!( MotionType(prob.model.bodies), F̃b, qs, q0, prob)

    return nothing
end

"""
    boundary_forces!(::Union{Type{Static}, Type{MovingGrid}},
                     F̃b, qs, q0, prob)

Solve modified Poisson problem for uB = 0 and bc2 = 0
```
 Bf̃ = Eq = ECψ
```
"""
function boundary_forces!(::Union{Type{Static}, Type{MovingGrid}},
                          F̃b, qs, q0, prob)
    Q = prob.model.work.q2  # Net flux

    broadcast!(+, Q, qs, q0)        # qs + q0
    mul!(F̃b, prob.model.mats.E, Q)  # E*(qs .+ state.q0)... using fb here as working array
    F̃b .= prob.Binv*F̃b              # Allocates a small amount of memory

    return nothing
end

"""
    boundary_forces!(::Type{V} where V <: Motion, F̃b, qs, q0, prob)

Solve the Poisson problem for bc2 = 0 (???) with nonzero boundary velocity ub
```
Bf̃ = Eq - ub
   = ECψ - ub
```
"""
function boundary_forces!(::Type{V} where V <: Motion,
                          F̃b, qs, q0, prob)
    # Working memory for in-place operations (small allocation)
    F̃work = similar(F̃b)
    Q = prob.model.work.q2   # Net flux

    broadcast!(+, Q, qs, q0)                # qs + q0
    mul!(F̃work, prob.model.mats.E, Q)       # E*(qs .+ state.q0)
    F̃work .-= get_ub(prob.model.bodies)*prob.model.grid.h   # Enforce no-slip conditions
    mul!(F̃b, prob.Binv, F̃work);

    return nothing
end

"""
project_circ!(Γs, state, prob)

Update circulation to satisfy no-slip condition.

Dispatch based on the type of motion in the problem.

This allows precomputing regularization and interpolation where possible.
"""
function project_circ!(Γs::AbstractArray,
                       state::IBState,
                       prob::AbstractIBProblem)
    project_circ!(MotionType(prob.model.bodies), Γs, state, prob)
    return nothing
end

function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState,
                       prob::AbstractIBProblem)
    """
    High-level version:
        Γ = Γs - Ainv * (E*C)'*F̃b
    """
    Γwork = prob.model.work.Γ2 # Working memory
    E, C = prob.model.mats.E, prob.model.mats.C

    # Low-level version:
    state.Γ .= Γs
    @views mul!( Γs[:, 1], (E*C)', state.F̃b[:, 1])  # Γ = ∇ x (E'*fb)
    @views mul!( Γwork, prob.Ainv[1], Γs[:, 1])

    state.Γ[:, 1] .-= Γwork

    return nothing
end


"""
    base_flux!(state::IBState, prob::AbstractIBProblem, t::Float64)
Set background flux based on `prob.model.bodies[].motion`
"""
function base_flux!(state::IBState,
                    prob::AbstractIBProblem,
                    t::Float64)
    base_flux!(MotionType(prob.model.bodies), state, prob, t)
end

"Freestream flux for linearized problem (uniform zero)"
function base_flux!(::Type{Static},
                    state::IBState,
                    prob::LinearizedIBProblem,
                    t::Float64)
    grid = prob.model.grid
    nu = grid.ny*(grid.nx+1);  # Number of x-flux points

    # Set all to zero since base flux is accounted for in base state
    for lev = 1 : grid.mg
        state.q0[ 1:nu, lev ] .= 0.0      # x-flux
        state.q0[ nu+1:end, lev ] .= 0.0  # y-flux
    end
end

"Initialize irrotational freestream flux when not time-varying"
function base_flux!(::Type{T} where T <: InertialMotion,
                    state::IBState,
                    prob::IBProblem,
                    t::Float64)
    grid = prob.model.grid
    Ux = prob.model.freestream.Ux(t)
    Uy = prob.model.freestream.Uy(t)
    α = prob.model.freestream.inclination(t)

    nu = grid.ny*(grid.nx+1);  # Number of x-flux points
    for lev = 1 : grid.mg
        # Coarse grid spacing
        hc = grid.h * 2^( lev - 1 );

        # write fluid velocity flux in body-fixed frame
        state.q0[ 1:nu, lev ] .= (Ux*cos(α) - Uy*sin(α))* hc      # x-flux
        state.q0[ nu+1:end, lev ] .= (Ux*sin(α) + Uy*cos(α))*hc  # y-flux
    end
end

"Update time-varying background flux for moving grid"
function base_flux!(::Type{MovingGrid},
                    state::IBState,
                    prob::IBProblem,
                    t::Float64)
    @assert length(prob.model.bodies) == 1 # Assumes only one body
    grid = prob.model.grid
    motion = prob.model.bodies[1].motion
    XX, YY = prob.model.XX, prob.model.YY;
    nu = grid.ny*(grid.nx+1);  # Number of x-flux points
    nq = grid.nq

    ### Rotational part
    Ω = -motion.θ̇(t)
    α = -motion.θ(t)

    ### Potential flow part (note θ = -α for angle of attack)
    Ux0 = motion.U(t)*cos(α) - motion.V(t)*sin(α)
    Uy0 = motion.U(t)*sin(α) + motion.V(t)*cos(α)

    ## Add in underlying freestream components
    Uxf = prob.model.freestream.Ux(t)
    Uyf = prob.model.freestream.Uy(t)
    αf = prob.model.freestream.inclination(t)

    Ux0 += Uxf*cos(αf)-Uyf*sin(αf)
    Uy0 += Uxf*sin(αf)+Uyf*cos(αf)

    state.q0 .*= 0.0
    for lev=1:grid.mg

        hc = grid.h*2.0^(Float64(lev)-1.0)  # Coarse grid spacing

        ### x-fluxes
        @views state.q0[1:nu, lev] .= YY[:, lev]
        @views state.q0[1:nu, lev] .*= -hc*Ω

        ### y-fluxes
        @views state.q0[(nu+1):nq, lev] .= XX[:, lev]
        @views state.q0[(nu+1):nq, lev] .*= hc*Ω

        ### Irrotational part
        @views state.q0[1:nu, lev] .+= hc*Ux0      # x-flux
        @views state.q0[(nu+1):nq, lev] .+= hc*Uy0  # y-velocity

    end
end




"""
    update_stress!(state, prob)

Store surface stresses and integrated forces.

Mutates "state"
"""
function update_stress!(state::IBState,
                        prob::AbstractIBProblem)
    nb, nf = get_body_info(prob.model.bodies)
    h = prob.model.grid.h
    dt = prob.scheme.dt
    state.F̃b = reshape(state.F̃b, sum(nb), 2)

    # Store surface stress and integrated forces
    nbod_tally = 0; #  Used to keep a tally of which body we're on
    for j = 1 : length(prob.model.bodies)
        ds = prob.model.bodies[j].ds
        # surface stresses
        state.fb[j] .= state.F̃b[nbod_tally .+ (1:nb[j]), :] *(h / dt) ./ [ds ds] ;

        # integrated forces
        state.CD[j] = 2 * sum( ds .* state.fb[j][:, 1] ) ;
        state.CL[j] = 2 * sum( ds .* state.fb[j][:, 2] ) ;

        # update body index
        nbod_tally += nb[j];
    end

    state.F̃b = reshape(state.F̃b, 2*sum(nb), 1)
    return nothing
end

"""
    update_bodies!(prob, t)

Update the immersed bodies and coupling matrices (if applicable).
"""
function update_bodies!(state::IBState, prob::AbstractIBProblem, t::Float64)
    update_bodies!(MotionType(prob.model.bodies), state, prob, t)
    return nothing
end

" No motion for static or moving grid cases "
function update_bodies!(::Union{Type{Static}, Type{MovingGrid}},
                        state::IBState,
                        prob::AbstractIBProblem,
                        t::Float64)
    return nothing
end

#DEPRECATED because of MotionFunction?
# " For a rotating cylinder, update the velocity, but not the operators "
# function update_bodies!(::Type{RotatingCyl}, prob::IBProblem, t::Float64)
#     model = prob.model
#     bodies, grid = prob.model.bodies, prob.model.grid
#     for j=1:length(bodies)
#         move_body!(bodies[j], t)
#     end
#     return nothing
# end

" For moving grid, update the bodies and the coupling operators"
function update_bodies!(::Type{MotionFunction}, state::IBState, prob::IBProblem, t::Float64)
    model = prob.model
    bodies, grid = prob.model.bodies, prob.model.grid
    for j=1:length(bodies)
        move_body!(bodies[j], t)
        state.xb[j] = bodies[j].xb
    end

    # For arbitrary motion in a fixed frame, have to update operators
    model.mats.E = setup_reg( grid, bodies )
    prob.Binv = get_Binv(model, prob.Ainv[1])

    return nothing
end
