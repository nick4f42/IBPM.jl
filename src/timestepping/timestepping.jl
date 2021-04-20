using LinearAlgebra: norm  # FOR DEBUGGING

"""
    advance!(state::IBState, prob::IBProblem, t::Float64)

Advance state forward in time.
"""
function advance!(state, prob, t) end

function advance!(state::IBState{UniformGrid},
                  prob::IBProblem,
                  t::Float64)
    grid = prob.model.grid

    # Move bodies and update coupling matrices (E)
    update_bodies!(prob, t)

    if MotionType(prob.model.bodies) == MovingGrid
        base_flux!(state, prob, t)
    end

    # Alias working memory for notational clarity
    qs = prob.work.q1  # Trial flux
    Γs = prob.work.Γ1  # Trial circulation

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
    vort2flux!( state.ψ, state.q, state.Γ, prob.model );

    #--A few simulation quantities of interest
    # get CFL (u * dt / dx) :
    state.cfl = maximum( abs.( (1/grid.h^2) * state.q * prob.scheme.dt ) ) ;
end



function advance!(state::IBState{MultiGrid},
                  prob::IBProblem,
                  t::Float64)
    grid = prob.model.grid
    update_bodies!(prob, t)

    if MotionType(prob.model.bodies) == MovingGrid
        base_flux!(state, prob, t)
    end

    # Alias working memory for notational clarity
    #   This leaves work.Γ2, Γ3 and work.q1, q2 available
    qs = prob.work.q1  # Trial flux
    Γs = prob.work.Γ1  # Trial circulation

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
    #println("Final circulation: ", sum(state.Γ.^2))
    vort2flux!( state.ψ, state.q, state.Γ, prob.model, grid.mg );

    #println("Final circulation: ", sum(state.Γ.^2))
    #println("Final flux: ", sum(state.q.^2))

    #--A few simulation quantities of interest
    # get CFL (u * dt / dx) :
    dt = prob.scheme.dt
    # TODO: DOES THIS ALLOCATE??
    state.cfl = maximum( @. abs( (1/(grid.h^2))*state.q[:, 1]*dt ) ) ;
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
function get_trial_state!(qs, Γs, state, prob) end

function get_trial_state!(qs::AbstractArray,
                          Γs::AbstractArray,
                          state::IBState{UniformGrid},
                          prob::IBProblem)
    grid = prob.model.grid  # SHOULDN'T NEED THIS HERE
    dt = prob.scheme.dt
    work = prob.work
    rhs = work.Γ2  # RHS of discretized equation

    #compute the nonlinear term for the current time step
    nonlinear!( state.nonlin[1], state, prob );

    # Explicit part of Laplacian
    mul!(rhs, prob.A, state.Γ)

    # Explicit nonlinear terms from multistep scheme
    for n=1:length(prob.scheme.β)
        work.Γ3 .= state.nonlin[n]
        rmul!(work.Γ3, prob.scheme.β[n]*dt)
        rhs .+= work.Γ3
    end

    # Trial circulation  Γs = Ainv * rhs
    mul!(Γs, prob.Ainv, rhs);

    # Trial velocity  (note ψ is used here as a dummy variable)
    vort2flux!( state.ψ, qs, Γs, prob.model )

    # Store current nonlinear term
    state.nonlin[2] .= state.nonlin[1];
end



function get_trial_state!(qs::AbstractArray,
                          Γs::AbstractArray,
                          state::IBState{MultiGrid},
                          prob::IBProblem)
    """
    I think Γbc might be used in Fortran the same as Γs in Matlab
    """
    dt = prob.scheme.dt
    work = prob.work
    grid = prob.model.grid
    rhsbc = work.bc
    rhs = work.Γ2  # RHS of discretized equation


    # TODO: PREALOCATE
    lastbc = zeros(2*(grid.nx+1)+2*(grid.ny+1), grid.mg)
    Γbc = zeros(2*(grid.nx+1)+2*(grid.ny+1))

    for lev=1:grid.mg-1
        @views get_bc!(lastbc[:, lev], state.Γ[:, lev+1], 0.25, grid)
    end
    #println(sum(lastbc.^2))

    for lev=grid.mg:-1:1
        #println("=== TRIAL STATE, LEVEL ", lev, " ===")

        if lev < grid.mg
            @views get_bc!(Γbc, state.Γ[:, lev+1], 0.25, grid)
        else
            Γbc .= 0.0
        end

        #compute the nonlinear term for the current time step
        #println("=== NONLIN ", lev, " ===")
        @views nonlinear!( state.nonlin[1][:, lev], state, lastbc[:, lev], prob, lev );

        rhsbc .= 0.0
        hc = grid.h * 2^( lev - 1);
        vfac = 0.5*dt/ ( prob.model.Re * hc^2 * 4*grid.nx*grid.ny )
        apply_bc!( rhsbc, @view(lastbc[:, lev]), vfac, grid )
        apply_bc!( rhsbc, Γbc, vfac, grid )

        #println(sum(rhsbc.^2))

        # Combine explicit Laplacian and nonlinear terms into a rhs
        #println("Trial flux, level ", lev)
        #println(sum(rhs[:, lev].^2))
        #mul!(rhs[:, lev], prob.A[lev], state.Γ[:, lev])
        #println(sum(rhs[:, lev].^2))
        rhs[:, lev] = prob.A[lev] * state.Γ[:, lev]

        for n=1:length(prob.scheme.β)
            work.Γ3[:, lev] .= state.nonlin[n][:, lev]/hc^2
            rmul!(work.Γ3, prob.scheme.β[n]*dt)
            rhs .-= work.Γ3
        end

        # Include boundary conditions
        # TODO: is this actually faster??
        #   High-level: rhs += 0.5*rhsbc
        work.Γ3[:, lev] .= rhsbc
        rmul!(work.Γ3, 0.5*dt)
        rhs .+= work.Γ3

        # Trial circulation  Γs = Ainv * rhs
        # TODO: use @view to do in-place multiplication
        # Doesn't work here because of view and FFT plan for even indices... WHY??
        Γs[:, lev] .= prob.Ainv[lev] * rhs[:, lev]

        #println(sum(Γs[:, lev].^2))
    end

    # Store nonlinear solution for use in next time step
    state.nonlin[2] .= state.nonlin[1]

    # THIS WAS JUST lev=2 IN THE MATLAB... DO WE NEED ALL?
    vort2flux!( state.ψ, qs,  Γs, prob.model, grid.mg )

    #println(sum(Γs.^2))
end

"""
    boundary_forces!(F̃b, qs, q0, prob)

Solve the Poisson equation (25) in Colonius & Taira (2008).

Dispatch based on the type of motion in the problem - allows precomputing
    regularization and interpolation where possible.
"""
function boundary_forces!(F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    boundary_forces!( MotionType(prob.model.bodies), F̃b, qs, q0, prob)
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
                          F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    #println("=== BOUNDARY FORCES ===")

    E = prob.model.mats.E
    h = prob.model.grid.h
    # Working memory for in-place operations
    qwork = @view(prob.work.q2[:, 1])
    broadcast!(+, qwork, qs, q0)                     # qs + q0
    mul!(F̃b, E, qwork)                               # E*(qs .+ state.q0)... using fb here as working array
    #println(sum(q0.^2))
    #println(sum(qs.^2))
    #println(sum(F̃b.^2))
    F̃b .= prob.Binv*F̃b                         # Allocates a small amount of memory

    #println(sum(F̃b.^2))
    #sleep(100)
end

"""
    boundary_forces!(::Type{RotatingCyl}, F̃b, qs, q0, prob)

Solve the Poisson problem for bc2 = 0 with special case of rotating cylinder.

In this case the points don't need to move, but they do have nonzero velocity
```
Bf̃ = Eq - ub
   = ECψ - ub
```
"""
function boundary_forces!(::Type{RotatingCyl},
                          F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    E = prob.model.mats.E
    h = prob.model.grid.h

    # Working memory for in-place operations
    F̃work = similar(F̃b)
    qwork = @view(prob.work.q2[:, 1])

    broadcast!(+, qwork, qs, q0)           # qs + q0
    mul!(F̃work, E, qwork)                  # E*(qs .+ state.q0)
    F̃work .-= get_ub(prob.model.bodies)*prob.model.grid.h   # Enforce no-slip conditions
    mul!(F̃b, prob.Binv, F̃work/h);
end

"""
    project_circ!(Γs, state, prob)

Update circulation to satisfy no-slip condition.

Dispatch based on the type of motion in the problem.

This allows precomputing regularization and interpolation where possible.
"""
function project_circ!(Γs::AbstractArray,
                     state::IBState,
                     prob::IBProblem)
    project_circ!(MotionType(prob.model.bodies), Γs, state, prob)
end

function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{UniformGrid},
                       prob::IBProblem)
    """
    High-level version:
        state.Γ[:, 1] .= Γs .- prob.Ainv[1] * (mats.RET*fb_til_dt)
    """
    Γwork = @view(prob.work.Γ3[:, 1]) # Working memory
    E, C = prob.model.mats.E, prob.model.mats.C
    fb_til_dt = state.F̃b

    # Low-level version:
    state.Γ .= Γs   # Now Γs is free for working memory
    @views mul!( Γs[:, 1], (E*C)', fb_til_dt)  # Γ = ∇ x (E'*fb)
    @views mul!( Γwork, prob.Ainv, Γs[:, 1])
    state.Γ[:, 1] .-= Γwork
end

function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{MultiGrid},
                       prob::IBProblem)
    """
    High-level version:
        state.Γ[:, 1] .= Γs .- prob.Ainv[1] * (mats.RET*fb_til_dt)
    """
    #println("=== PROJECT CIRC ===")
    Γwork = @view(prob.work.Γ3[:, 1]) # Working memory
    E, C = prob.model.mats.E, prob.model.mats.C
    fb_til_dt = state.F̃b

    # Low-level version:
    state.Γ .= Γs
    @views mul!( Γs[:, 1], (E*C)', fb_til_dt)  # Γ = ∇ x (E'*fb)
    #println(sum(Γs[:, 1].^2))
    @views mul!( Γwork, prob.Ainv[1], Γs[:, 1])  # This is the only difference with the UniformGrid version
    #println(sum(Γwork.^2))

    state.Γ[:, 1] .-= Γwork

end

"""
    update_stress!(state, prob)

Store surface stresses and integrated forces.

Mutates "state"
"""
function update_stress!(state::IBState,
                        prob::IBProblem)
    nb, nf = get_body_info(prob.model.bodies)
    h = prob.model.grid.h
    dt = prob.scheme.dt
    fb_til_dt = state.F̃b

    # Store surface stress and integrated forces
    nbod_tally = 0; #  Used to keep a tally of which body we're on
    for j = 1 : length(prob.model.bodies)
        ds = prob.model.bodies[j].ds
        # surface stresses
        state.fb[j] .= fb_til_dt[nbod_tally .+ (1:nf[j])] *(h / dt) ./ [ds; ds] ;

        # integrated forces
        state.CD[j] = 2 * sum( ds .* state.fb[j][1 : nb[j]] ) ;
        state.CL[j] = 2 * sum( ds .* state.fb[j][1 .+ nb[j] : nf[j] ] ) ;

        # update body index
        nbod_tally += nf[j];
    end
end

"""
    AB2(dt::Float64)

Initialize second-order Adams-Bashforth scheme.
"""
function AB2(dt::Float64)
    return AdamsBashforth(dt, [1.5, -0.5])
end


"""
    update_bodies!(prob, t)

Update the immersed bodies and coupling matrices (if applicable).

TODO: break out by multiple dispatch... but don't duplicate code
"""
function update_bodies!(prob::IBProblem, t::Float64)
    model = prob.model
    bodies, grid = prob.model.bodies, prob.model.grid
    motion = MotionType(bodies)

    if motion != Static
        for j=1:length(bodies)
            move_body!(bodies[j], t)
        end
    end

    # For arbitrary motion in an inertial frame, have to update operators
    if motion == MotionFunction
        model.mats.E = coupling_mat( grid, bodies )
        model.mats.RET = (model.mats.E*model.mats.C)'
        prob.Binv = get_B(model, prob.Ainv)
    end
end
