"""
    compute_cfl(state, prob)

Compute the CFL number (uΔt/Δx) based on the fine-grid flux

Note that this uses working memory that is also used in `nonlinear!`
"""
function compute_cfl(state, prob)
    Δt, Δx = prob.scheme.dt, prob.model.grid.h
    qwork = prob.model.work.q5
    @views @. qwork = abs( state.q[:, 1] )
    return maximum(qwork)*Δt/Δx
end
