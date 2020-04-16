using Documenter, eFEMpart

makedocs(;
    modules=[eFEMpart],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pseastham/eFEMpart.jl/blob/{commit}{path}#L{line}",
    sitename="eFEMpart.jl",
    authors="Patrick Eastham",
    assets=String[],
)

deploydocs(;
    repo="github.com/pseastham/eFEMpart.jl",
)
