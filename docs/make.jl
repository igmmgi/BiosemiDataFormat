using Documenter

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using BiosemiDataFormat

# Set up the documentation
makedocs(
  sitename="BiosemiDataFormat",
  format=Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    edit_link=:commit,
    assets=String[],
  ),
  modules=[BiosemiDataFormat],
  pages=[
    "Home" => "index.md",
    "API Reference" => "api.md",
  ],
  repo="https://github.com/igmmgi/BiosemiDataFormat.jl",
  doctest=true,
  checkdocs=:exports,
)
