using Documenter

include("../src/ibpm.jl")

using .ibpm

makedocs(
    sitename = "IBPM Documentation",
    pages = 
    [
        "IBPM.jl" => "index.md",
        "Files" => [
            "ibpm.jl" => "files/ibpm.md",
            "timestepping" => [
                "save_info.jl" => "files/save_info.md",
                "timestepping.jl" => "files/timestepping.md"
            ],
            "preprocessing" => [
                "models.jl" => "files/models.md",
                "problem-types.jl" => "files/problem-types.md",
                "read-user-vars.jl" => "files/read-user-vars.md",
                "state-types.jl" => "files/state-types.md",
            ],
            "structure-domain" => [
                "body-types.jl" => "files/body-types.md",
                "motion-types.jl" => "files/motion-types.md",
                "move-body-utils.jl" => "files/move-body-utils.md",
                "sample-bodies.jl" => "files/sample-bodies.md",
                "structure-domain.jl" => "files/structure-domain.md"
            ]
        ],
        "Guide" => [
            "Guide1" => "guide/guide.md"
            ]
    ]
)

deploydocs(
    repo = "github.com/nick4f42/IBPM.jl.git"
)