# ClimateSatellite.jl
*Download satellite data for climate applications using Julia!*

`ClimateSatellite.jl` aims to streamline the basic processing of satellite data for climate applications and research.  This includes:
* Downloading satellite datasets of interest (please see below for supported datasets)
* Preliminary analysis of downloaded satellite data, including
  * Yearly and monthly `mean`, `std`, `maximum` and `minimum` of daily averages
  * Yearly and monthly `mean`, `std`, `maximum` and `minimum` of sub-daily raw data (where available)
  * All of the above in the 2D-spatial domain, as well as meridional-averaged and zonal-averaged domains.
* Packing downloaded data into `Int16` format to save disk-space
* Extraction of data for specific geographic regions (GeoRegions)
  * See [`GeoRegions.jl`](https://github.com/JuliaClimate/GeoRegions.jl) for more information

## Supported Datasets
The following datasets are currently supported by `ClimateSatellite.jl`:

|    ID    | Mission | Product | Version |
|  :---:   |  :---:  |  :---   |  :---:  |
| gpmimerg |   PMM   | Global Precipitation Mission - IMERGv6 Research | pre-v1.0 |
| gpmlate  |   PMM   | Global Precipitation Mission - IMERGv6 NRT Late Run | pre-v1.0 |
| gpmearly |   PMM   | Global Precipitation Mission - IMERGv6 NRT Early Run | pre-v1.0 |
| mtpw2m   |  MIMIC  | Total Precipitable Water v2m | pre-v1.0 |

I plan to release support for the following datasets:

|    ID    | Mission | Product | Version |
|  :---:   |  :---:  |  :---   |  :---:  |
| trmm3b42 |  PMM  | Tropical Rainfall Measuring Mission - 3B42v7 Research | v1.0 |
| rsstmi |  RSS  | TRMM Microwave Imager | v1.1 |
| rssgmi |  RSS  | GPM Microwave Imager | v1.1 |

## Installation
`ClimateSatellite.jl` can be installed using Julia's built-in package manager as follows:

```
julia> ]
(@v1.4) pkg> add ClimateSatellite
```

You can update `ClimateSatellite.jl` to the latest version using
```
(@v1.4) pkg> update ClimateSatellite
```

And if you want to get the latest release without waiting for me to update the Julia Registry (although this generally isn't necessary since I make a point to release patch versions as soon as I find bugs or add new working features), you may fix the version to the `master` branch of the GitHub repository:
```
(@v1.4) pkg> add ClimateSatellite#master
```

## Documentation

The documentation for `ClimateSatellite.jl` is divided into three components:
1. Tutorials - meant as an introduction to the package
2. How-to Examples - geared towards those looking for specific examples of what can be done
3. API Reference - comprehensive summary of all exported functionalities

!!! tip "A Note on the Examples:"
    All the output for the coding examples were produced using my computer with key security information (such as login info) omitted.  The examples cannot be run online because the file size requirements are too big.  Copying and pasting the code examples (with relevant directory and login information changes) should produce the same results.

## Getting help
If you are interested in using `ClimateSatellite.jl` or are trying to figure out how to use it, please feel free to ask me questions and get in touch!  Please feel free to [open an issue](https://github.com/JuliaClimate/ClimateSatellite.jl/issues/new) if you have any questions, comments, suggestions, etc!
