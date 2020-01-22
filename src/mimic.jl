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

function mimicnclist(yr::Integer, mo::Integer, ndy::Integer)

    fname = Array{AbstractString,2}(undef,48,ndy)
    for dy = 1 : ndy; date = Date(yr,mo,ii);
        for hh = 1 : 24

            fname[hh,dy] = "comp$(ymd2str(date)).$(@sprintf("%02d",hh-1))0000.nc";

        end
    end

    furl = "ftp://ftp.ssec.wisc.edu/pub/mtpw2/data/$(yrmo2str(date))"

    return fname,furl

end

function mimicget(url,file::AbstractString,tdir::AbstractString)

    try download("$(url)/$(file)",joinpath(tdir,file));
        @debug "$(Dates.now()) - Downloaded MIMIC tropospheric precipitable water data file $(file)"
    catch;
        @info "$(Dates.now()) - MIMIC tropospheric precipitable water data $(file) does not exist."
    end

end

function mimicretrieve(
    fname::Array{AbstractString,2}, furl::Array{AbstractString,2}, tdir::AbstractString
)

    for ii = 1 : length(fname); mimicget(furl,fname[ii],tdir); end

end

function mimicextract(
    fname::Array{AbstractString,2}, fol::AbstractString, reg::AbstractString="GLB"
)

    lon,lat = mimiclonlat(); nlon = length(lon); nlat = length(lat)
    rlon,rlat = regiongridvec(reg,lon,lat); nrlon = length(rlon); nrlat = length(rlat);

    data = Array{Int16,3}(undef,nrlon,nrlat,length(fH5));
    raw  = zeros(1440,721);
    tmp  = Array{Real,2}(undef,nrlon,nrlat);

    for ii = 1 : length(fname)

        fii = joinpath(fol,fname[ii]);
        if isfile(fH5ii)

            try
                ncread!(fncii,"tpwGrid",raw);
                tmp .= regionextractgrid(raw,reg,lon,lat,rawi)
                real2int16!(data[:,:,ii],tmp,offset=60,scale=60/32767);
            catch
                @warn "$(Dates.now()) - Unable to extract/open $(fncii).  Setting MIMIC tropospheric precipitable water data values set to NaN."
                data[:,:,ii] .= -32768;
            end

        else

            @info "$(Dates.now()) - $(fii) does not exist.  All values set to NaN."
            data[:,:,ii] .= -32768;

        end

    end

    return data,[rlon,rlat]

end

function gpmdwn(
    regions::Array{AbstractString,1}, yr::Integer, info::Dict
)

    tdir = clisattmp(info); if !isdir(tdir) mkpath(tdir); end

    for mo = 1 : 12;

        dateii = Date(yr,mo); ndy = daysinmonth(yr,mo);
        fname,furl = mimicnclist(yr,mo,ndy); mimicretrieve(fname,furl,tdir);

        for reg in regions
            data,grid = mimicextract(fname,tdir,reg); clisatsave(data,grid,reg,info,dateii)
        end

    end

end
