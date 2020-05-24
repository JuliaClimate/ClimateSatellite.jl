using Documenter
using ClimateSatellite

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
            "Data Downloads"       => "tutorials/downloads.md",
            "Using GeoRegions"     => "tutorials/georegions.md",
            # "Region Extraction"    => "tutorials/extract.md",
            # "Preliminary Analysis" => "tutorials/analysis.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/JuliaClimate/ClimateSatellite.jl.git",
)
