"""
This file contains all the back-end scripts in ClimateSatellite.jl that are for the download
and extraction of data from TRMM 3B42 product of the Precipitation Measurement Mission:
    * retrieval of file names
    * downloading data from arthurhou and jsimpson servers
    * extraction of chosen regions from global dataset

"""

function trmmlonlat()
    lon = convert(Array,-179.875:0.25:179.875); lat = convert(Array,49.875:-0.25:-49.875);
    return lon,lat
end

function trmmhdf(date::TimeType)

    yr = Dates.year(date); mo = Dates.month(date); ndy = daysinmonth(date);
    fname = Array{AbstractString,2}(undef,8,ndy);

    if date > Date(2010,9,30); suffix = "7A"; else; suffix = "7"; end

    @debug "$(Dates.now()) - Creating list of data files to download ..."
    for dy = 1 : ndy; dateii = Date(yr,mo,dy);
        for hh = 1 : 8

            hr = (hh-1)*3; hr = @sprintf("%02d",hr);
            fname[hh,dy] = "3B42.$(ymd2str(dateii)).$(hr).$(suffix).HDF"

        end
    end

    return fname

end
