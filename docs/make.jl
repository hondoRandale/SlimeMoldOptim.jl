using SlimeMoldOptim
using Documenter

DocMeta.setdocmeta!(SlimeMoldOptim, :DocTestSetup, :(using SlimeMoldOptim); recursive=true)

makedocs(;
    modules=[SlimeMoldOptim],
    authors="hondoRandale <Jules.Rasetaharison@tutanota.com>",
    repo="https://github.com/hondoRandale/SlimeMoldOptim.jl/blob/{commit}{path}#{line}",
    sitename="SlimeMoldOptim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hondoRandale.github.io/SlimeMoldOptim.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hondoRandale/SlimeMoldOptim.jl",
    devbranch="main",
)
