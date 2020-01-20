"""
This file contains all the back-end scripts in ClimateSatellite.jl that are for the download
and extraction of data from the GPM satellite of the Precipitation Measurement Mission:
    * retrieval of file names
    * downloading data from arthurhou and jsimpson servers
    * extraction of chosen regions from global dataset

"""

function gpmftp(info::Dict)

    if info["product"] != "IMERG"
          return pmmftpopen("jsimpson")
    else; return pmmftpopen("arthurhou")
    end

end

function gpmfdwn(
    regions::Array{AbstractString,1}, year::Integer, info::Dict, email::AbstractString
)

    tdir = clisattmp(info); if !isdir(tdir) mkpath(tdir); end

    ftp = gpmftp(info);

    dvec = Date(year,1,1):Days(1):Date(year,12,31)
    for date in dvec

        gpmftpcd(date,ftp);

    end

    H5 = gpmfdt(date); gpmftpcd(date,ftp);

    @info "$(Dates.now()) - Downloading $(info["source"]) $(info["product"]) $(info["variable"]) data for $(Date(date))"
    for ii = 1 : length(fH5)
        fH5ii = fH5[ii];
        if !isfile("$(sroot)/tmp/$(fH5ii)"); gpmfget(ftp,fH5ii,sroot);
        else
            if overwrite
                @info "$(Dates.now()) - GPM Research (Final) precipitation data file $(fH5ii) already exists.  Overwriting."
                mimicget(ftp,fH5ii,sroot);
            else; @info "$(Dates.now()) - GPM Research (Final) precipitation data file $(fH5ii) already exists.  Not overwriting."
            end
        end
    end
    pmmftpclose(ftp);

end
