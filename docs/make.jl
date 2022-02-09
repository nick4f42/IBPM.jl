using Documenter

include("../src/ibpm.jl")

using .ibpm

makedocs(
    sitename = "IBPM Documentation",
    pages = [
        "ibpm.jl" => "ibpm.md",
        "models.jl" => "models.md"
    ]
)

deploydocs(
    repo = "github.com/nick4f42/IBPM.jl.git"
)