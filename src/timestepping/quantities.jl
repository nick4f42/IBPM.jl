module Quantities

using ..IBPM: IBState


"""
    vel_flux(state)
"""
function vel_flux end

"""
    vel_potential(state)
"""
function vel_potential end

"""
    circulation(state)
"""
function circulation end

"""
    streamfunction(state)
"""
function streamfunction end

"""
    surface_stress(state)
"""
function surface_stress end

"""
    body_impulse(state)
"""
function body_impulse end

"""
    drag_coef(state)
"""
function drag_coef end

"""
    lift_coef(state)
"""
function lift_coef end

"""
    cfl(state)
"""
function cfl end

"""
    slip(state)
"""
function slip end

"""
    body_points(state)
"""
function body_points end


for (name, field) in (
    (:vel_flux, :q),
    (:vel_potential, :q0),
    (:circulation, :Γ),
    (:streamfunction, :ψ),
    (:surface_stress, :fb),
    (:body_impulse, :F̃b),
    (:drag_coef, :CD),
    (:lift_coef, :CL),
    (:cfl, :cfl),
    (:slip, :slip),
    (:body_points, :xb),
)
    @eval $name(state::IBState) = deepcopy(state.$field)
end


end # module
