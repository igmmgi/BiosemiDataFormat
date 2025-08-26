using Documenter

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using BiosemiDataFormat

# Set up the documentation
makedocs(
    sitename = "BiosemiDataFormat",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = :commit,
        repolink = "https://github.com/igmmgi/BiosemiDataFormat.jl",
        assets = String[],
    ),
    modules = [BiosemiDataFormat],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    repo = "github.com/igmmgi/BiosemiDataFormat.jl.git",
    doctest = true,
    checkdocs = :exports,
)
