"""
This file contains all the functions that are used in the retrieval of GPM
precipitable water datasets.  This includes the downloading and and retrieval
data from specific areas, as well as plotting of the data.

-- gpmdwn(date,region="GLB",root)
    downloads hour gridded precipitable water for a given date (and if
    specified, region)
-- gpmdt(date)
    extracts the year, month and day and generates:
    -- a list of files to download
    -- the url where the files are to be downloaded from

"""

# Setup Functions
function gpmroot(root::AbstractString)
    info(logger,"Making root folder for GPM datasets.")
    sroot = "$(root)/GPM/"; mkpath(sroot)
    if !isdir(sroot)
        notice(logger,"GPM directory does not exist.  Creating now.");
        mkpath(sroot);
    end
    info(logger,"The root folder for the GPM precipitation data is in $(sroot).")
    return sroot
end

function gpmlonlat()
    lon = convert(Array,-179.95:0.1:179.95); lat = convert(Array,-89.95:0.1:89.95);
    return lon,lat
end

function gpmncfile(date::Date,reg::AbstractString)
    return "gpm_$(reg)_prcp_$(Date(date)).nc"
end

function gpmhdf5(date)

    fname = Array{String}(undef,48)
    for ii = 1 : 48
        hr = (ii-1)/2; mi = mod(hr,1); hr = hr - mi;
        hr = @sprintf("%02d",hr);
        id = @sprintf("%04d",(ii-1)*30);

        if mi == 0; fname[ii] = "3B-HHR.MS.MRG.3IMERG.$(ymd2str(date))-S$(hr)0000-E$(hr)2959.$(id).V06B.HDF5"
        else;       fname[ii] = "3B-HHR.MS.MRG.3IMERG.$(ymd2str(date))-S$(hr)3000-E$(hr)5959.$(id).V06B.HDF5"
        end

    end
    return fname

end

function gpmfol(reg::AbstractString="GLB",date::Date,sroot::AbstractString)
    fol = "$(sroot)/$(reg)/$(yr2str(date))/$(mo2str(date))/"
    if !isdir(fol)
        notice(logger,"GPM data directory for the $(regionname(reg)) region, year $(yr2str(date)) and month $(mo2str(date)) does not exist.");
        info(logger,"Creating data directory $(fol)."); mkpath(fol);
    end
    return fol
end

# FTP Functions
function gpmftpopen()
    email = "natgeo.wong%40outlook.com"
    info(logger,"Opening FTP request to arthurhou.pps.eosdis.nasa.gov.")
    return FTP("ftp://$(email):$(email)@arthurhou.pps.eosdis.nasa.gov")
end

function gpmftpcd(date::Date,ftp)
    info(logger,"Entering IMERG directory for $(ymd2str(date)).")
    cd(ftp,"gpmdata/$(yr2str(date))/$(mo2str(date))/$(dy2str(date))/imerg")
end

function gpmftpclose(ftp)
    info(logger,"Closing FTP request to arthurhou.pps.eosdis.nasa.gov.")
    close(ftp)
end

# GPM Processing Functions
function gpmdt(date::Date)

    info(logger,"Extracting year, month and day and time from $(date).")
    fname = gpmhdf5(date)

    info(logger,"Extracted the list of GPM files to be downloaded for $(Date(date)).")
    return fname

end

function gpmget(ftp,file,sroot::AbstractString)
    try download(ftp,"$(file)","$(sroot)/tmp/$(file)")
        debug(logger,"Downloaded GPM precipitation data file $(file)")
    catch; notice(logger,"GPM precipitation data $(file) does not exist.")
    end
end

function gpmdwn(date,sroot::AbstractString,overwrite=false)

    tdir = "$(sroot)/tmp/";
    if !isdir(tdir) mkpath(tdir); end

    fH5 = gpmdt(date);
    ftp = gpmftpopen(); gpmftpcd(date,ftp);
    info(logger,"Downloading GPM precipitation data for $(Date(date))")
    for ii = 1 : length(fH5)
        fH5ii = fH5[ii];
        if !isfile("$(sroot)/tmp/$(fH5ii)"); gpmget(ftp,fH5ii,sroot);
        else
            if overwrite
                notice(logger,"GPM precipitation data file $(fH5ii) already exists.  Overwriting.")
                mimicget(ftp,fH5ii,sroot);
            else; notice(logger,"GPM precipitation data file $(fH5ii) already exists.  Not overwriting.")
            end
        end
    end
    gpmftpclose(ftp);

end

function gpmextract(date::Date,sroot::AbstractString,
    reg::AbstractString="GLB")

    data = zeros(1800,3600,48) #default grid step size is 0.1x0.1 in GPM
    fH5  = gpmdt(date);

    info(logger,"Extracting and compiling GPM precipitation data from HDF5 files.")
    for ii = 1 : length(fH5)
        fH5ii = "$(sroot)/tmp/$(fH5[ii])";
        if isfile(fH5ii)
              data[:,:,ii] = h5read("$(fH5ii)","/Grid/precipitationCal");
        else; notice(logger,"$(fH5ii) does not exist.  GPM precipitation data values set to NaN.")
              data[:,:,ii] .= NaN;
        end
    end

    if reg != "GLB"
        info(logger,"We do not wish to extract GPM precipitation data for the entire globe.")
        info(logger,"Finding grid-point boundaries ...")
        lon,lat = gpmlonlat();
        bounds = regionbounds(reg); igrid = regiongrid(bounds,lon,lat);

        info(logger,"Extracting GPM precipitation data for the region.")
        rdata,rpnts = regionextractgrid(reg,lon,lat,data)
    else; rdata = data; lon,lat = gpmlonlat(); rpnts = [lat[end],lat[1],lon[end],lon[1]];
    end

    return rdata,rpnts

end

function gpmsave(data,rpnts,date::Date,sroot::AbstractString,reg::AbstractString="GLB")

    fnc = gpmncfile(date,reg); data = permutedims(data,[2,1,3])
    nlon = size(data,1); lon = convert(Array,rpnts[4]:0.1:rpnts[3]);
    nlat = size(data,2); lat = convert(Array,rpnts[2]:0.1:rpnts[1]);
    if nlon != length(lon); error(logger,"nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error(logger,"nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_prcp = "prcp"; att_tpw = Dict("units" => "mm/hr");
    var_lon  = "lon";  att_lon = Dict("units" => "degree");
    var_lat  = "lat";  att_lat = Dict("units" => "degree");

    if isfile(fnc)
        notice(logger,"Unfinished netCDF file $(fnc) detected.  Deleting.");
        rm(fnc);
    end

    info(logger,"Creating GPM precipitation netCDF file $(fnc) ...")
    nccreate(fnc,var_prcp,"nlon",nlon,"nlat",nlat,"t",48,atts=att_prcp,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    info(logger,"Saving GPM precipitation data to netCDF file $(fnc) ...")
    ncwrite(data,fnc,var_prcp);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat);

    fol = gpmfol(date,sroot);
    info(logger,"Moving $(fnc) to data directory $(fol)")

    if isfile("$(fol)/$(fnc)"); notice(logger,"An older version of $(fnc) exists in the $(fol) directory.  Overwriting.") end

    mv(fnc,"$(fol)/$(fnc)",force=true);

end

function gpmrmtmp(date::Date,sroot::AbstractString)

    fH5 = gpmdt(date);
    info(logger,"Deleting raw GPM precipitation data files.")
    for ii = 1 : length(fH5)
        fH5ii = "$(sroot)/tmp/$(fH5[ii])";
        if isfile(fH5ii); rm(fH5ii) end
    end

end

# Compiled Function
function gpmrun(date::Date,sroot::AbstractString,reg::AbstractString="GLB")
    gpmdwn(date,sroot); data,pnts = gpmextract(date,sroot,reg);
    gpmsave(data,pnts,date,sroot,reg); gpmrmtmp(date,sroot);
end
