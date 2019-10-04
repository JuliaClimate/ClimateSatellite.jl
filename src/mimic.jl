"""
This file contains all the functions that are used in the retrieval of MIMIC
precipitable water datasets.  This includes the downloading and and retrieval
data from specific areas, as well as plotting of the data.



-- mimicdwn(date,region="GLB",root)
    downloads hour gridded precipitable water for a given date (and if
    specified, region)
-- mimicdt(date)
    extracts the year, month and day and generates:
    -- a list of files to download
    -- the url where the files are to be downloaded from

"""

# Setup Functions
function mimicroot(root::AbstractString)
    info(logger,"Finding root folder for MIMIC datasets.")
    sroot = "$(root)/MIMIC/";
    if !isdir(sroot)
        notice(logger,"MIMIC directory does not exist.  Creating now.");
        mkpath(sroot);
    end
    info(logger,"The root folder for MIMIC is in $(sroot).")
    return sroot
end

function mimiclonlat()
    lon = convert(Array,-180:0.25:180); lat = convert(Array,-90:0.25:90); pop!(lon);
    return lon,lat
end

function mimicfile(date::Date,reg::AbstractString="GLB")
    return "mimic_$(reg)_tpw_$(ymd2str(date)).nc"
end

function mimicfol(date::Date,sroot::AbstractString,reg::AbstractString="GLB")
    fol = "$(sroot)/$(reg)/$(yrmo2dir(date))/"
    if !isdir(fol)
        notice(logger,"MIMIC data directory for the $(regionname(reg)) region, year $(yr2str(date)) and month $(mo2str(date)) does not exist.");
        info(logger,"Creating data directory $(fol)."); mkpath(fol);
    end
    return fol
end

# MIMIC Processing Functions
function mimicdt(date::Date)

    info(logger,"Extracting year, month and day and time from $(date).")
    fname = Array{String}(undef,24)
    for hr = 1 : 24
        fname[hr] = "comp$(ymd2str(date)).$(@sprintf("%02d",hr-1))0000.nc";
    end

    furl = "ftp://ftp.ssec.wisc.edu/pub/mtpw2/data/$(yrmo2str(date))/"

    info(logger,"Extracted the list of MIMIC files for $(Date(date)).")
    return fname,furl

end

function mimicget(url,file,sroot::AbstractString)
    try download("$(url)$(file)","$(sroot)/tmp/$(file)")
        debug(logger,"Downloaded MIMIC tropospheric precipitable water data file $(file)")
    catch; notice(logger,"MIMIC tropospheric precipitable water data $(file) does not exist.")
    end
end

function mimicdwn(date,sroot::AbstractString,overwrite=false)

    tdir = "$(sroot)/tmp/"
    if !isdir(tdir) mkpath(tdir); end

    fnc,url = mimicdt(date);
    info(logger,"Downloading MIMIC tropospheric precipitable water data for $(Date(date))")
    for ii = 1 : length(fnc)
        fncii = fnc[ii];
        if !isfile("$(sroot)/tmp/$(fncii)"); mimicget(url,fncii,sroot);
        else
            if overwrite
                notice(logger,"MIMIC tropospheric precipitable water data file $(fncii) already exists.  Overwriting.")
                mimicget(url,fncii,sroot);
            else; notice(logger,"MIMIC tropospheric precipitable water data file $(fncii) already exists.  Not overwriting.")
            end
        end
    end

end

function mimicextract(date::Date,sroot::AbstractString,
    reg::AbstractString="GLB")

    data = zeros(1440,721,24) #default grid step size is 0.25x0.25 in MIMIC
    fnc,url = mimicdt(date);

    info(logger,"Extracting and compiling MIMIC tropospheric precipitable water data from raw netCDF files.")
    for ii = 1 : 24
        fncii = "$(sroot)/tmp/$(fnc[ii])";
        if isfile(fncii)
              data[:,:,ii] = ncread(fncii,"tpwGrid");
        else; notice(logger,"$(fncii) does not exists.  MIMIC tropospheric precipitable water data values set to NaN.")
              data[:,:,ii] .= NaN;
        end
    end

    if reg != "GLB"
        info(logger,"We do not wish to extract MIMIC tropospheric precipitable water data for the entire globe.")
        info(logger,"Finding grid-point boundaries ...")
        lon,lat = mimiclonlat();
        bounds = regionbounds(reg); igrid = regiongrid(bounds,lon,lat);

        info(logger,"Extracting MIMIC tropospheric precipitable water data for the region.")
        rdata,rpnts = regionextractgrid(reg,lon,lat,data)
    else; rdata = data; lon,lat = mimiclonlat(); rpnts = [lat[end],lat[1],lon[end],lon[1]];
    end

    return rdata,rpnts

end

function mimicsave(data,rpnts,date::Date,sroot::AbstractString,reg::AbstractString="GLB")

    fnc = mimicfile(date,reg);
    nlon = size(data,1); lon = convert(Array,rpnts[4]:0.25:rpnts[3]);
    nlat = size(data,2); lat = convert(Array,rpnts[2]:0.25:rpnts[1]);
    if nlon != length(lon); error(logger,"nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error(logger,"nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_tpw = "tpw"; att_tpw = Dict("units" => "mm");
    var_lon = "lon"; att_lon = Dict("units" => "degree");
    var_lat = "lat"; att_lat = Dict("units" => "degree");

    if isfile(fnc)
        notice(logger,"Unfinished netCDF file $(fnc) detected.  Deleting.");
        rm(fnc);
    end

    info(logger,"Creating MIMIC tropospheric precipitable water netCDF file $(fnc) ...")
    nccreate(fnc,var_tpw,"nlon",nlon,"nlat",nlat,"t",24,atts=att_tpw,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    info(logger,"Saving MIMIC tropospheric water vapour data to netCDF file $(fnc) ...")
    ncwrite(data,fnc,var_tpw);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat)

    fol = mimicfol(date,sroot);
    info(logger,"Moving $(fnc) to data directory $(fol)")

    if isfile("$(fol)/$(fnc)"); notice(logger,"An older version of $(fnc) exists in the $(fol) directory.  Overwriting.") end

    mv(fnc,"$(fol)/$(fnc)",force=true);

end

function mimicrmtmp(date::Date,sroot::AbstractString)

    fnc,url = mimicdt(date);
    info(logger,"Deleting raw MIMIC tropospheric water vapour files.")
    for ii = 1 : 24
        fncii = "$(sroot)/tmp/$(fnc[ii])";
        if isfile(fncii); rm(fncii) end
    end

end

# Compiled Function
function mimicrun(date::Date,sroot::AbstractString,reg::AbstractString="GLB")
    mimicdwn(date,sroot); data,pnts = mimicextract(date,sroot,reg);
    mimicsave(data,pnts,date,sroot,reg); mimicrmtmp(date,sroot);
end
