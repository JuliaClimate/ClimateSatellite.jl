"""
This file contains all the functions that are used in the retrieval of MIMIC
precipitable water datasets.  This includes the downloading and and retrieval
data from specific areas, as well as plotting of the data.

"""

# Setup Functions
function mimicroot(root::AbstractString)
    @info "$(Dates.now()) - Finding root folder for MIMIC datasets."
    sroot = "$(root)/MIMIC/";
    if !isdir(sroot)
        @info "$(Dates.now()) - MIMIC directory does not exist.  Creating now."
        mkpath(sroot);
    end
    @info "$(Dates.now()) - The root folder for MIMIC is in $(sroot)."
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
        @info "$(Dates.now()) - MIMIC data directory for the $(regionfullname(reg)) region, year $(yr2str(date)) and month $(mo2str(date)) does not exist."
        @info "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end
    return fol
end

# MIMIC Processing Functions
function mimicdt(date::Date)

    @info "$(Dates.now()) - Extracting year, month and day and time from $(date)."
    fname = Array{String}(undef,24)
    for hr = 1 : 24
        fname[hr] = "comp$(ymd2str(date)).$(@sprintf("%02d",hr-1))0000.nc";
    end

    furl = "ftp://ftp.ssec.wisc.edu/pub/mtpw2/data/$(yrmo2str(date))/"

    @info "$(Dates.now()) - Extracted the list of MIMIC files for $(Date(date))."
    return fname,furl

end

function mimicget(url,file,sroot::AbstractString)
    try download("$(url)$(file)","$(sroot)/tmp/$(file)");
        @debug "$(Dates.now()) - Downloaded MIMIC tropospheric precipitable water data file $(file)"
    catch; @info "$(Dates.now()) - MIMIC tropospheric precipitable water data $(file) does not exist."
    end
end

function mimicdwn(date,sroot::AbstractString,overwrite=false)

    tdir = "$(sroot)/tmp/"
    if !isdir(tdir) mkpath(tdir); end

    fnc,url = mimicdt(date);
    @info "$(Dates.now()) - Downloading MIMIC tropospheric precipitable water data for $(Date(date))"
    for ii = 1 : length(fnc)
        fncii = fnc[ii];
        if !isfile("$(sroot)/tmp/$(fncii)"); mimicget(url,fncii,sroot);
        else
            if overwrite
                @info "$(Dates.now()) - MIMIC tropospheric precipitable water data file $(fncii) already exists.  Overwriting."
                mimicget(url,fncii,sroot);
            else; @info "$(Dates.now()) - MIMIC tropospheric precipitable water data file $(fncii) already exists.  Not overwriting."
            end
        end
    end

end

function mimicextract(date::Date,sroot::AbstractString,
    reg::AbstractString="GLB")

    data = zeros(1440,721,24) #default grid step size is 0.25x0.25 in MIMIC
    fnc,url = mimicdt(date);

    @info "$(Dates.now()) - Extracting and compiling MIMIC tropospheric precipitable water data from raw netCDF files."
    for ii = 1 : 24
        fncii = "$(sroot)/tmp/$(fnc[ii])";
        if isfile(fncii)
            try
                @info "$(Dates.now()) - Extracting MIMIC tropospheric precipitable water data from $(fncii) ..."
                data[:,:,ii] = ncread(fncii,"tpwGrid");
            catch
                @warn "$(Dates.now()) - Unable to extract/open $(fncii).  Setting MIMIC tropospheric precipitable water data values set to NaN."
                data[:,:,ii] .= NaN;
            end
        else
            @info "$(Dates.now()) - $(fncii) does not exists.  MIMIC tropospheric precipitable water data values set to NaN."
            data[:,:,ii] .= NaN;
        end
    end

    @debug "$(Dates.now()) - NetCDF.jl's ncread causes memory leakage.  Using ncclose() as a workaround."
    ncclose()

    if reg != "GLB"
        @info "$(Dates.now()) - We do not wish to extract MIMIC tropospheric precipitable water data for the entire globe."
        @info "$(Dates.now()) - Finding grid-point boundaries ..."
        lon,lat = mimiclonlat();
        bounds = regionbounds(reg); igrid = regiongrid(bounds,lon,lat);

        @info "$(Dates.now()) - Extracting MIMIC tropospheric precipitable water data for the region."
        rdata,rgrid = regionextractgrid(data,reg,lon,lat)
    else; rdata = data; rgrid = [mimiclonlat()];
    end

    return rdata,rgrid

end

function mimicsave(data,rgrid,date::Date,sroot::AbstractString,reg::AbstractString="GLB")

    fnc = mimicfile(date,reg);
    nlon = size(data,1); lon = rgrid[1];
    nlat = size(data,2); lat = rgrid[2];
    if nlon != length(lon); error("$(Dates.now()) - nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error("$(Dates.now()) - nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_tpw = "tpw"; att_tpw = Dict("units" => "mm");
    var_lon = "lon"; att_lon = Dict("units" => "degree");
    var_lat = "lat"; att_lat = Dict("units" => "degree");

    if isfile(fnc)
        @info "$(Dates.now()) - Unfinished netCDF file $(fnc) detected.  Deleting."
        rm(fnc);
    end

    @info "$(Dates.now()) - Creating MIMIC tropospheric precipitable water netCDF file $(fnc) ..."
    nccreate(fnc,var_tpw,"nlon",nlon,"nlat",nlat,"t",24,atts=att_tpw,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    @info "$(Dates.now()) - Saving MIMIC tropospheric water vapour data to netCDF file $(fnc) ..."
    ncwrite(data,fnc,var_tpw);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat)

    @debug "$(Dates.now()) - NetCDF.jl's ncread causes memory leakage.  Using ncclose() as a workaround."
    ncclose()

    fol = mimicfol(date,sroot,reg);
    @info "$(Dates.now()) - Moving $(fnc) to data directory $(fol)"

    if isfile("$(fol)/$(fnc)"); @info "$(Dates.now()) - An older version of $(fnc) exists in the $(fol) directory.  Overwriting." end

    mv(fnc,"$(fol)/$(fnc)",force=true);

end

function mimicrmtmp(date::Date,sroot::AbstractString)

    fnc,url = mimicdt(date);
    @info "$(Dates.now()) - Deleting raw MIMIC tropospheric water vapour files."
    for ii = 1 : 24
        fncii = "$(sroot)/tmp/$(fnc[ii])";
        if isfile(fncii); rm(fncii) end
    end

end

# Compiled Function
function mimicrun(date::Date,sroot::AbstractString,reg::AbstractArray=["GLB"])
    sroot = mimicroot(sroot); cd(sroot); mimicdwn(date,sroot);
    for regii in reg
        data,grid = mimicextract(date,sroot,regii);
        mimicsave(data,grid,date,sroot,regii);
    end
    mimicrmtmp(date,sroot);
end
