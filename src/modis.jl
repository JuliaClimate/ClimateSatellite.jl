"""
This file contains all the functions that are used in the retrieval of MODIS
precipitable water datasets.  This includes the downloading and and retrieval
data from specific areas, as well as plotting of the data.

-- modisdwn(date,region="GLB",root)
    downloads hour gridded precipitable water for a given date (and if
    specified, region)
-- modisdt(date)
    extracts the year, month and day and generates:
    -- a list of files to download
    -- the url where the files are to be downloaded from

"""

# Setup Functions
function modisroot(root::AbstractString)
    infl("Finding root folder for MODIS datasets.")
    sroot = "$(root)/MODIS/";
    if !isdir(sroot)
        infl("MODIS directory does not exist.  Creating now.");
        mkpath(sroot);
    end
    infl("The root folder for MODIS is in $(sroot).")
    return sroot
end

function modislonlat()
    lon = convert(Array,-180:0.25:180); lat = convert(Array,-90:0.25:90); pop!(lon);
    return lon,lat
end

function modisfile(date::Date,reg::AbstractString="GLB")
    return "modis_$(reg)_tpw_$(ymd2str(date)).nc"
end

function modisfol(date::Date,sroot::AbstractString)
    fol = "$(sroot)/$(reg)/$(yr2str(date))/$(mo2str(date))/"
    if !isdir(fol)
        infl("MODIS data directory for the $(regionname(reg)) region, year $(yr2str(date)) and month $(mo2str(date)) does not exist.");
        infl("Creating data directory $(fol)."); mkpath(fol);
    end
    return fol
end

# MODIS Processing Functions
function modisdt(date::Date)

    infl("Extracting year, month and day and time from $(date).")
    fname = Array{String}(undef,24)
    for hr = 1 : 24
        fname[hr] = "comp$(ymd2str(date)).$(@sprintf("%02d",hr-1))0000.nc";
    end

    furl = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MOD08_D3/$(yrdy2dir(date))/"

    infl("Extracted the list of MODIS files for $(Date(date)).")
    return fname,furl

end

function modisget(url,file,sroot::AbstractString)
    try download("$(url)$(file)","$(sroot)/tmp/$(file)")
        debl("Downloaded MODIS tropospheric precipitable water data file $(file)")
    catch; infl("MODIS tropospheric precipitable water data $(file) does not exist.")
    end
end

function modisdwn(date,sroot::AbstractString,overwrite=false)

    tdir = "$(sroot)/tmp/"
    if !isdir(tdir) mkpath(tdir); end

    fnc,url = modisdt(date);
    infl("Downloading MODIS tropospheric precipitable water data for $(Date(date))")
    for ii = 1 : length(fnc)
        fncii = fnc[ii];
        if !isfile("$(sroot)/tmp/$(fncii)"); modisget(url,fncii,sroot);
        else
            if overwrite
                infl("MODIS tropospheric precipitable water data file $(fncii) already exists.  Overwriting.")
                modisget(url,fncii,sroot);
            else; infl("MODIS tropospheric precipitable water data file $(fncii) already exists.  Not overwriting.")
            end
        end
    end

end

function modisextract(date::Date,sroot::AbstractString,
    reg::AbstractString="GLB")

    data = zeros(1440,721,24) #default grid step size is 0.25x0.25 in MODIS
    fnc,url = modisdt(date);

    infl("Extracting and compiling MODIS tropospheric precipitable water data from raw netCDF files.")
    for ii = 1 : 24
        fncii = "$(sroot)/tmp/$(fnc[ii])";
        if isfile(fncii)
              data[:,:,ii] = ncread(fncii,"tpwGrid");
        else; infl("$(fncii) does not exists.  MODIS tropospheric precipitable water data values set to NaN.")
              data[:,:,ii] .= NaN;
        end
    end

    if reg != "GLB"
        infl("We do not wish to extract MODIS tropospheric precipitable water data for the entire globe.")
        infl("Finding grid-point boundaries ...")
        lon,lat = modislonlat();
        bounds = regionbounds(reg); igrid = regiongrid(bounds,lon,lat);

        infl("Extracting MODIS tropospheric precipitable water data for the region.")
        rdata,rpnts = regionextractgrid(reg,lon,lat,data)
    else; rdata = data; lon,lat = modislonlat(); rpnts = [lat[end],lat[1],lon[end],lon[1]];
    end

    return rdata,rpnts

end

function modissave(data,rpnts,date::Date,sroot::AbstractString,reg::AbstractString="GLB")

    fnc = modisfile(date,reg);
    nlon = size(data,1); lon = convert(Array,rpnts[4]:0.25:rpnts[3]);
    nlat = size(data,2); lat = convert(Array,rpnts[2]:0.25:rpnts[1]);
    if nlon != length(lon); errl("nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); errl("nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_tpw = "tpw"; att_tpw = Dict("units" => "mm");
    var_lon = "lon"; att_lon = Dict("units" => "degree");
    var_lat = "lat"; att_lat = Dict("units" => "degree");

    if isfile(fnc)
        infl("Unfinished netCDF file $(fnc) detected.  Deleting.");
        rm(fnc);
    end

    infl("Creating MODIS tropospheric precipitable water netCDF file $(fnc) ...")
    nccreate(fnc,var_tpw,"nlon",nlon,"nlat",nlat,"t",24,atts=att_tpw,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    infl("Saving MODIS tropospheric water vapour data to netCDF file $(fnc) ...")
    ncwrite(data,fnc,var_tpw);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat)

    fol = modisfol(date,sroot);
    infl("Moving $(fnc) to data directory $(fol)")

    if isfile("$(fol)/$(fnc)"); infl("An older version of $(fnc) exists in the $(fol) directory.  Overwriting.") end

    mv(fnc,"$(fol)/$(fnc)",force=true);

end

function modisrmtmp(date::Date,sroot::AbstractString)

    fnc,url = modisdt(date);
    infl("Deleting raw MODIS tropospheric water vapour files.")
    for ii = 1 : 24
        fncii = "$(sroot)/tmp/$(fnc[ii])";
        if isfile(fncii); rm(fncii) end
    end

end

# Compiled Function
function modisrun(date::Date,sroot::AbstractString,reg::AbstractString="GLB")
    modisdwn(date,sroot); data,pnts = modisextract(date,sroot,reg);
    modissave(data,pnts,date,sroot,reg); modisrmtmp(date,sroot);
end
