using Documenter, BioSemiBDF
makedocs(
  modules=[BioSemiBDF],
  format=:html,
  sitename="BioSemBDF.jl",
  doctest=false)

#= deploydocs = Deps.pip("mkdocs", "python-markdown.math"), =#
#=   repo = "github.com/igmmgi/BioSemiBDF.git", =#
#=   julia = "1.0.1", =#
#=   osname = "linux") =#

