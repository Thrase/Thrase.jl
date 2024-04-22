using Thrase
using Documenter

DocMeta.setdocmeta!(Thrase, :DocTestSetup, :(using Thrase); recursive=true)

makedocs(;
    modules=[Thrase],
    authors="Brittany A. Erickson",
    repo="https://github.com/Thrase/Thrase.jl/blob/{commit}{path}#{line}",
    sitename="Thrase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Thrase.github.io/Thrase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Numerical Method" => "method.md",
        "Install and Quick Start" => "quickstart.md",
        "Input Parameters" => "parameters.md",
    ],
)

deploydocs(;
    repo="github.com/Thrase/Thrase.jl",
    devbranch="main",
)
