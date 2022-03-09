using Documenter

include("../src/ibpm.jl")

using .ibpm

makedocs(
    sitename = "IBPM Documentation",
    pages = 
    [
        "Immersed Bounudary Projection Method: Compuational Fluid Dynamics Solver" => "index.md",
        "Manual" => [
            "Getting Started: Installation and More" => "manual/getting_started.md",
            "Simple Static Cylinder Tutorial (Before)" => "manual/examples/cyl_before.md",
            "Simple Static Cylinder Tutorial (After)" => "manual/examples/cyl_after.md",
            "Problem Options" => "manual/problem_options.md",
            "Save Data Options" => "manual/save_options.md",
            "Plotting" => "manual/plotting.md",
            "Examples" => "manual/examples/examples.md"
        ],
        "Files" => [
            "ibpm.jl" => "files/ibpm.md",
            "timestepping" => [
                "save_info.jl" => "files/timestepping/save_info.md",
                "timestepping.jl" => "files/timestepping/timestepping.md"
            ],
            "pre-processing" => [
                "models.jl" => "files/pre-processing/models.md",
                "problem-types.jl" => "files/pre-processing/problem-types.md",
                "read-user-vars.jl" => "files/pre-processing/read-user-vars.md",
                "state-types.jl" => "files/pre-processing/state-types.md",
            ],
            "structure-domain" => [
                "body-types.jl" => "files/structure-domain/body-types.md",
                "motion-types.jl" => "files/structure-domain/motion-types.md",
                "move-body-utils.jl" => "files/structure-domain/move-body-utils.md",
                "sample-bodies.jl" => "files/structure-domain/sample-bodies.md",
                "structure-domain.jl" => "files/structure-domain/structure-domain.md"
            ],
            "plotting" => [
                "plotting-utils.jl" => "files/plotting/plotting-utils.md"
            ],
            "interface-coupling" => [
                "interface-coupling.jl" => "files/interface-coupling/interface-coupling.md"
            ],
            "fluid-operators" => [
                "curl.jl" => "files/fluid-operators/curl.md",
                "dst-inversion.jl" => "files/fluid-operators/dst-inversion.md",
                "laplacian.jl" => "files/fluid-operators/laplacian.md",
                "mg-utils.jl" => "files/fluid-operators/mg-utils.md",
                "nonlinlear.jl" => "files/fluid-operators/nonlinear.md"
            ],
            "fluid-domain" => [
                "fluid-domain.jl" => "files/fluid-domain/fluid-domain.md"
            ],
            "experimental" => [
                "sfd.jl" => "files/experimental/sfd.md"
            ]
        ]
    ]
)

deploydocs(
    repo = "github.com/nick4f42/IBPM.jl.git"
)