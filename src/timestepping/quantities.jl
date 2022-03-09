module Quantities

using ..IBPM: IBState


"""
    lift_coef(state)
"""
function lift_coef end

"""
    drag_coef(state)
"""
function drag_coef end

# TODO: Add functions for remaining quantities in IBState

for (name, field) in (
    (:lift_coef, :CL),
    (:drag_coef, :CD),
)
    @eval $name(state::IBState) = deepcopy(state.$field)
end


end # module
