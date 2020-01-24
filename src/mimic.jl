"""
This file contains all the back-end scripts in ClimateSatellite.jl that are for the download
and extraction of precipitable water data from the Morphed Integrated Microwave Imagery at
CIMSS (MIMIC)-TPW2m product:
    * retrieval of file names
    * downloading data from ftp server
    * extraction of chosen regions from global dataset

"""

function mimiclonlat()
    lon = convert(Array,-180:0.25:180); lat = convert(Array,-90:0.25:90); pop!(lon);
    return lon,lat
end

function mimicnclist(date::TimeType)

    yr = Dates.year(date); mo = Dates.month(date); ndy = daysinmonth(date);
    fname = Array{<:AbstractString,1}(undef,48*ndy);

    @debug "$(Dates.now()) - Creating list of data files to download ..."
    for dy = 1 : ndy; date = Date(yr,mo,ii);
        for hh = 1 : 24

            fname[hh,dy] = "comp$(ymd2str(date)).$(@sprintf("%02d",hh-1))0000.nc";

        end
    end

    furl = "ftp://ftp.ssec.wisc.edu/pub/mtpw2/data/$(yrmo2str(date))"

    return fname,furl

end

function mimicget(
    url::AbstractString,
    file::AbstractString,
    tdir::AbstractString,
    overwrite::Bool
)
    if overwrite && !isfile(joinpath(tdir,file))
        try download("$(url)/$(file)",joinpath(tdir,file));
            @debug "$(Dates.now()) - Downloaded MIMIC tropospheric precipitable water data file $(file)"
        catch;
            @warn "$(Dates.now()) - MIMIC tropospheric precipitable water data $(file) does not exist."
        end
    end

end

function mimicretrieve(
    fname::Array{<:AbstractString,2},
    furl::Array{<:AbstractString,2},
    tdir::AbstractString, info::Dict, overwrite::Bool
)

    @info "$(Dates.now()) - Starting data download of $(info["product"]) data for $(date) ..."

    for ii = 1 : length(fname); mimicget(furl,fname[ii],tdir,overwrite); end

    @info "$(Dates.now()) - Data download of $(info["product"]) data for $(date) has been completed."

end

function mimicextract(
    fname::Array{<:AbstractString,2}, fol::AbstractString, info::Dict,
    reg::AbstractString="GLB"
)

    lon,lat = mimiclonlat(); rlon,rlat,rinfo = regiongridvec(reg,lon,lat);
    nlon = length(lon); nrlon = length(rlon);
    nlat = length(lat); nrlat = length(rlat);

    data = zeros(Int16,nrlon,nrlat,length(fname));
    rawi = zeros(nlon,nlat); raw = zeros(nlon,nlat); tmp = zeros(nrlon,nrlat);

    @info "$(Dates.now()) - Extracted regional $(info["product"]) data for $(rinfo["fullname"])."

    for ii = 1 : length(fname)

        fii = joinpath(fol,fname[ii]);
        if isfile(fii)

            try
                ncread!(fii,"tpwGrid",raw);
                tmp .= regionextractgrid(raw,rinfo,lon,lat,rawi)
                real2int16!(data[:,:,ii],tmp,offset=info["offset"],scale=info["scale"]);
            catch
                @warn "$(Dates.now()) - Unable to extract/open $(fii).  Setting MIMIC tropospheric precipitable water data values set to NaN."
                data[:,:,ii] .= -32768;
            end

        else

            @info "$(Dates.now()) - $(fii) does not exist.  All values set to NaN."
            data[:,:,ii] .= -32768;

        end

    end

    @info "$(Dates.now()) - Extracted regional $(info["product"]) data for $(rinfo["fullname"])."
    return data,[rlon,rlat]

end

function mimicdwn(
    regions::Array{AbstractString,1}, date::TimeType, info::Dict; overwrite::Bool
)

    tdir = clisattmp(info); if !isdir(tdir) mkpath(tdir); end
    fname,furl = mimicnclist(date); mimicretrieve(fname,furl,tdir,info,overwrite);

    for reg in regions
        data,grid = mimicextract(fname,tdir,reg); clisatsave(data,grid,reg,info,date)
    end

end
