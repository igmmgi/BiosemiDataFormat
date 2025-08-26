using Documenter

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using BiosemiDataFormat

# Set up the documentation
makedocs(
  sitename="BiosemiDataFormat",
  format=Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    assets=String[],
  ),
  modules=[BiosemiDataFormat],
  pages=[
    "Home" => "index.md",
    "API Reference" => "api.md",
  ],

  doctest=true,
  checkdocs=:exports,
)

deploydocs(;
    repo = "github.com/igmmgi/BiosemiDataFormat.jl.git",
    versions = ["stable" => "v^", "v#.#", "dev" => "master"],
    push_preview = true,
)
