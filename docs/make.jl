using Documenter

include("../src/ibpm.jl")

using .ibpm

makedocs(
    sitename = "IBPM Documentation"
)

deploydocs(
    repo = "github.com/nick4f42/IBPM.jl.git"
)