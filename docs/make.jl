using Documenter

include("setup.jl")

const doctestsetup = quote
    using LinearCombinations, SimplicialSets
end

DocMeta.setdocmeta!(SimplicialSets, :DocTestSetup, doctestsetup; recursive = true)

makedocs(sitename = "SimplicialSets.jl",
    modules = [SimplicialSets],
    pages = [
        "index.md",
        "simplices.md",
        "basic.md",
        "ez.md",
        "surjection.md",
        "helpers.md",
    ],
    format = Documenter.HTML(),
    warnonly = true)
