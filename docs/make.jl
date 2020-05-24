using Documenter
using ClimateSatellite

format = Documenter.HTML(
    collapselevel = 1,
       prettyurls = get(ENV,"CI",nothing) == "true"
)

makedocs(
    modules  = [ClimateSatellite],
    doctest  = false,
    format   = Documenter.HTML(
        collapselevel = 1,
        prettyurls    = false
    ),
    authors  = "Nathanael Wong",
    sitename = "ClimateSatellite.jl",
    pages    = [
        "Home"      => "index.md",
        "Tutorials" => [
            "Basic Downloading"    => "tutorials/downloads.md",
            "Preliminary Analysis" => "tutorials/analysis.md",
            "Using GeoRegions"     => "tutorials/georegions.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/JuliaClimate/ClimateSatellite.jl.git",
)
