"""
This file contains all the back-end scripts in ClimateSatellite.jl that are for the download
and extraction of data from the GPM satellite of the Precipitation Measurement Mission:
    * retrieval of file names
    * downloading data from arthurhou and jsimpson servers
    * extraction of chosen regions from global dataset

"""

function gpmftp(info::Dict)

    if info["product"] != "IMERG"
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

function gpmh5list(yr::Integer, mo::Integer, ndy::Integer, info::Dict)

    fH5 = Array{AbstractString,2}(undef,48,ndy)
    for dy = 1 : ndy; date = Date(yr,mo,ii);
        for hh = 1 : 48

            hr = (hh-1)/2; mi = mod(hr,1); hr = hr - mi;
            hr = @sprintf("%02d",hr);
            id = @sprintf("%04d",(hh-1)*30);

            if mi == 0; dtstr = "$(ymd2str(date))-S$(hr)0000-E$(hr)2959.$(id)."
            else;       dtstr = "$(ymd2str(date))-S$(hr)3000-E$(hr)5959.$(id)."
            end

            fH5[hh,dy] = "$(info["prefix"]).$(dtstr).$(info["suffix"])"

        end
    end
    return fH5

end

function gpmftpcd(date::TimeType,ftp)

    @info "$(Dates.now()) - Entering IMERG directory for $(yrmostr(date))."

    if     info["short"] == "gpmimerg"
        cd(ftp,"gpmdata/$(yrmo2dir(date))")
    elseif info["short"] == "gpmlate"
        cd(ftp,"NRTPUB/imerg/late/$(yrmo2str(date))/")
    elseif info["short"] == "gpmearly"
        cd(ftp,"NRTPUB/imerg/early/$(yrmo2str(date))/")
    end

end

function gpmget(ftp,file::AbstractString,tdir::AbstractString)

    try download(ftp,file,joinpath(tdir,file));
        @debug "$(Dates.now()) - Downloaded GPM Near-RealTime (Late) precipitation data file $(file)"
    catch; @info "$(Dates.now()) - GPM Near-RealTime (Late) precipitation data $(file) does not exist."
    end

end

function gpmretrieve(
    fH5::Array{AbstractString,2},
    date::TimeType, ndy::Integer,
    tdir::AbstractString, info::Dict
)

    ftp = gpmftp(info); gpmftpcd(yr,mo,ftp);

    if info["short"] == "gpmimerg"

        for ii = 1 : ndy
            cd(ftp,"$(@sprintf("%02d",ii))/imerg");
            for jj = 1 : 48; gpmget(ftp,fH5[jj,ii],tdir); end
            cd(ftp,"../../")
        end

    else

        for ii = 1 : length(fH5); gpmget(ftp,fH5[ii],tdir); end

    end

    pmmftpclose(ftp);

end

function gpmextract(
    fH5::Array{AbstractString,2}, fol::AbstractString, reg::AbstractString="GLB"
)

    lon,lat = gpmlonlat(); nlon = length(lon); nlat = length(lat)
    rlon,rlat = regiongridvec(reg,lon,lat); nrlon = length(rlon); nrlat = length(rlat);

    data = Array{Int16,3}(undef,nrlon,nrlat,length(fH5));
    raw  = Array{Real,2}(undef,nlon,nlat);
    rawi = Array{Real,2}(undef,nlon,nlat);

    for ii = 1 : length(fH5)

        fH5ii = joinpath(fol,fH5[ii]);
        if isfile(fH5ii)

            rawii = h5read("$(fH5ii)","/Grid/precipitationCal");
            @debug "$(Dates.now()) - raw GPM precipitation data is given in (lat,lon) instead of (lon,lat).  Permuting to (lon,lat)"
            permutedims!(raw,rawii,[2,1,3]);
            tmp .= regionextractgrid(raw,reg,lon,lat,rawi)
            real2int16!(data[:,:,ii],tmp,offset=163.835,scale=1/200);

        else

            @info "$(Dates.now()) - $(fH5ii) does not exist.  All values set to NaN."
            data[:,:,ii] .= -32768;

        end

    end

    return data,[rlon,rlat]

end

function gpmdwn(
    regions::Array{AbstractString,1}, yr::Integer, info::Dict
)

    tdir = clisattmp(info); if !isdir(tdir) mkpath(tdir); end; gpmh5!(info)

    for mo = 1 : 12;

        dateii = Date(yr,mo); ndy = daysinmonth(yr,mo);
        fH5 = gpmh5list(yr,mo,ndy,info); gpmretrieve(fH5,dateii,ndy,tdir,info);

        for reg in regions
            data,grid = gpmextract(fH5,tdir,reg); clisatsave(data,grid,reg,info,dateii)
        end

    end

end
