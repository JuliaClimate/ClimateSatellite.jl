"""
This file contains all the back-end scripts in ClimateSatellite.jl that are for the download
and extraction of data from the GPM satellite of the Precipitation Measurement Mission:
    * retrieval of file names
    * downloading data from arthurhou and jsimpson servers
    * extraction of chosen regions from global dataset

"""

function gpmftp(info::Dict)

    if info["short"] != "gpmimerg"
          return pmmftpopen("jsimpson",info["email"])
    else; return pmmftpopen("arthurhou",info["email"])
    end

end

function gpmlonlat()
    lon = convert(Array,-179.95:0.1:179.95); lat = convert(Array,-89.95:0.1:89.95);
    return lon,lat
end

function gpmh5!(info::Dict)

    if info["short"] == "gpmimerg"
        info["prefix"] = "3B-HHR.MS.MRG.3IMERG";   info["suffix"] = "V06B.HDF5";
    elseif info["short"] == "gpmlate"
        info["prefix"] = "3B-HHR-L.MS.MRG.3IMERG"; info["suffix"] = "V06B.RT-H5";
    elseif info["short"] == "gpmearly"
        info["prefix"] = "3B-HHR-E.MS.MRG.3IMERG"; info["suffix"] = "V06B.RT-H5";
    end

end

function gpmh5list(date::TimeType, ndy::Integer, info::Dict)

    fH5 = Array{<:AbstractString,2}(undef,48,ndy);
    yr = Dates.year(date); mo = Dates.month(date);

    @debug "$(Dates.now()) - Creating list of data files to download ..."
    for dy = 1 : ndy; dateii = Date(yr,mo,dy);
        for hh = 1 : 48

            hr = (hh-1)/2; mi = mod(hr,1); hr = hr - mi;
            hr = @sprintf("%02d",hr);
            id = @sprintf("%04d",(hh-1)*30);

            if mi == 0; dtstr = "$(ymd2str(dateii))-S$(hr)0000-E$(hr)2959.$(id)"
            else;       dtstr = "$(ymd2str(dateii))-S$(hr)3000-E$(hr)5959.$(id)"
            end

            fH5[hh,dy] = "$(info["prefix"]).$(dtstr).$(info["suffix"])"

        end
    end

    return fH5

end

function gpmftpcd(date::TimeType,ftp,info::Dict)

    @info "$(Dates.now()) - Entering IMERG directory for $(date) ..."

    if     info["short"] == "gpmimerg"
        cd(ftp,"gpmdata/$(yrmo2dir(date))")
    elseif info["short"] == "gpmlate"
        cd(ftp,"NRTPUB/imerg/late/$(yrmo2str(date))/")
    elseif info["short"] == "gpmearly"
        cd(ftp,"NRTPUB/imerg/early/$(yrmo2str(date))/")
    end

end

function gpmget(ftp,file::AbstractString,tdir::AbstractString,overwrite::Bool)

    if overwrite && !isfile(joinpath(tdir,file))
        try download(ftp,file,joinpath(tdir,file));
            @debug "$(Dates.now()) - Downloaded data file $(file)"
        catch; @warn "$(Dates.now()) - Data file $(file) does not exist."
        end
    end

end

function gpmretrieve(
    fH5::Array{String,2},
    date::TimeType, ndy::Integer,
    tdir::AbstractString, info::Dict,
    overwrite::Bool
)

    @info "$(Dates.now()) - Starting data download of $(info["product"]) data for $(date) ..."

    ftp = gpmftp(info); gpmftpcd(date,ftp,info);

    if info["short"] == "gpmimerg"

        for ii = 1 : ndy
            cd(ftp,"$(@sprintf("%02d",ii))/imerg");
            for jj = 1 : 48; gpmget(ftp,fH5[jj,ii],tdir,overwrite); end
            cd(ftp,"../../")
        end

    else

        for ii = 1 : length(fH5); gpmget(ftp,fH5[ii],tdir,overwrite); end

    end

    pmmftpclose(ftp);
    @info "$(Dates.now()) - Data download of $(info["product"]) data for $(date) has been completed."

end

function gpmextract(
    fH5::Array{<:AbstractString,2}, fol::AbstractString, info::Dict,
    reg::AbstractString="GLB"
)

    lon,lat = gpmlonlat(); rlon,rlat,rinfo = regiongridvec(reg,lon,lat);
    nlon = length(lon); nrlon = length(rlon);
    nlat = length(lat); nrlat = length(rlat);

    data = zeros(Int16,nrlon,nrlat,length(fH5));
    rawi = zeros(nlon,nlat,1); raw = zeros(nlon,nlat,1); tmp = zeros(nrlon,nrlat,1);

    @info "$(Dates.now()) - Extracting regional $(info["product"]) data for $(rinfo["fullname"]) ..."

    for ii = 1 : length(fH5)

        fH5ii = joinpath(fol,fH5[ii]);
        if isfile(fH5ii)

            rawii = h5read("$(fH5ii)","/Grid/precipitationCal");
            @debug "$(Dates.now()) - raw GPM precipitation data is given in (lat,lon) instead of (lon,lat).  Permuting to (lon,lat)"
            permutedims!(raw,rawii,[2,1,3]);
            tmp .= regionextractgrid(raw,rinfo,lon,lat,rawi)
            real2int16!(data[:,:,ii],tmp,offset=163.835,scale=1/200);

        else

            @info "$(Dates.now()) - $(fH5ii) does not exist.  All values set to NaN."
            data[:,:,ii] .= -32768;

        end

    end

    @info "$(Dates.now()) - Extracted regional $(info["product"]) data for $(rinfo["fullname"])."
    return data,[rlon,rlat]

end

function gpmdwn(
    regions::Array{<:AbstractString,1}, date::TimeType, info::Dict; overwrite::Bool
)

    tdir = clisattmp(info); if !isdir(tdir) mkpath(tdir); end; gpmh5!(info)
    ndy = daysinmonth(date); fH5 = gpmh5list(date,ndy,info);
    #gpmretrieve(fH5,date,ndy,tdir,info,overwrite);

    for reg in regions
        data,grid = gpmextract(fH5,tdir,info,reg); #clisatsave(data,grid,reg,info,date)
    end

end
