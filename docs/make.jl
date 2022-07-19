using Documenter
using IBPM

makedocs(
    sitename="IBPM",
    pages=
    [
        "Home" => "index.md",
        "Getting Started" => [
            "Tutorial" => "getting-started/tutorial.md",
        ],
        "Manual" => [
            "experimental" => [
                "sfd.jl" => "manual/experimental/sfd.md"
            ],
            "fluid-domain" => [
                "fluid-domain.jl" => "manual/fluid-domain/fluid-domain.md"
            ],
            "fluid-operators" => [
                "curl.jl" => "manual/fluid-operators/curl.md",
                "dst-inversion.jl" => "manual/fluid-operators/dst-inversion.md",
                "laplacian.jl" => "manual/fluid-operators/laplacian.md",
                "mg-utils.jl" => "manual/fluid-operators/mg-utils.md",
                "nonlinlear.jl" => "manual/fluid-operators/nonlinear.md"
            ],
            "interface-coupling" => [
                "interface-coupling.jl" => "manual/interface-coupling/interface-coupling.md"
            ],
            "plotting" => [
                "plotting-utils.jl" => "manual/plotting/plotting-utils.md"
            ],
            "pre-processing" => [
                "models.jl" => "manual/pre-processing/models.md",
                "problem-types.jl" => "manual/pre-processing/problem-types.md",
                "state-types.jl" => "manual/pre-processing/state-types.md"
            ],
            "structure-domain" => [
                "body-types.jl" => "manual/structure-domain/body-types.md",
                "motion-types.jl" => "manual/structure-domain/motion-types.md",
                "move-body-utils.jl" => "manual/structure-domain/move-body-utils.md",
                "sample-bodies.jl" => "manual/structure-domain/sample-bodies.md",
                "structure-domain.jl" => "manual/structure-domain/structure-domain.md"
            ],
            "timestepping" => [
                "quantities.jl" => "manual/timestepping/quantities.md"
                "solve.jl" => "manual/timestepping/solve.md"
                "timestepping.jl" => "manual/timestepping/timestepping.md"
            ]
        ]
    ]
)

deploydocs(
    repo="github.com/nick4f42/IBPM.jl.git",
    devbranch="UI-overhaul",
)
