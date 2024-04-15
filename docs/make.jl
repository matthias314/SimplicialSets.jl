using Documenter

include("setup.jl")

const doctestsetup = quote
    using LinearCombinations, SimplicialSets
end

DocMeta.setdocmeta!(SimplicialSets, :DocTestSetup, doctestsetup; recursive = true)

makedocs(sitename = "SimplicialSets.jl",
    modules = [SimplicialSets],
    format = Documenter.HTML(),
    warnonly = true)
