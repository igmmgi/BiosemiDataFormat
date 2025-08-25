using Documenter

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using BiosemiDataFormat

# Set up the documentation
makedocs(
    sitename = "BiosemiDataFormat",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://igmmgi.github.io/BiosemiDataFormat.jl/",
        edit_link = :commit,
        assets = String[],
    ),
    modules = [BiosemiDataFormat],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md",
    ],
    repo = "github.com/igmmgi/BiosemiDataFormat.jl.git",
    doctest = true,
    checkdocs = :exports,
)
