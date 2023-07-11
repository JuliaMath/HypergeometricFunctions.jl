using ComplexPhasePortrait, Documenter, HypergeometricFunctions, Images

DocMeta.setdocmeta!(HypergeometricFunctions, :DocTestSetup, :(using HypergeometricFunctions); recursive=true)

makedocs(
    modules = [HypergeometricFunctions],
    sitename = "HypergeometricFunctions.jl",
    authors = "Richard Mikael Slevinsky",
    strict = false
)

deploydocs(repo = "github.com/JuliaMath/HypergeometricFunctions.jl.git")
