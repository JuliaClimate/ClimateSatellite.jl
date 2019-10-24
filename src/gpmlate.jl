"""
    This file contains all the functions that are used in the retrieval of GPM
    precipitable water datasets.  Current functionalities include:
        - downloading of raw data files
        - extraction of regional data (where applicable)
        - saving of data (global/regional)

    The following functionalities are forthcoming:
        - preliminary mapping of data
"""

# Setup Functions

"""
    gpmlroot(root) -> AbstractString

    Returns the path where GPM data will be downloaded and stored based on the
    given input "root", which must be a string
"""
function gpmlroot(root::AbstractString)
    @info "$(Dates.now()) - Making root folder for GPM Near-RealTime (Late) datasets."
    sroot = "$(root)/GPM-LATE/"; mkpath(sroot)
    if !isdir(sroot)
        @info "$(Dates.now()) - GPM Near-RealTime (Late) directory does not exist.  Creating now."
        mkpath(sroot);
    end
    @info "$(Dates.now()) - The root folder for the GPM Near-RealTime (Late) precipitation data is in $(sroot)."
    return sroot
end

"""
    gpmllonlat() -> Array,Array

    Returns two vector arrays for the longitude and latitude respectively.
    Longitude defaults are within [-180,180].  Grid spacing 0.1 x 0.1.
    Measurements are taken in the middle of grid points, therefore at .05º
    No inputs are required.
"""
function gpmllonlat()
    lon = convert(Array,-179.95:0.1:179.95); lat = convert(Array,-89.95:0.1:89.95);
    return lon,lat
end

"""
    gpmlncfile(date,reg) -> AbstractString

    Returns the ncfile name that the data extracted will be saved to, based on
    given date and region variable inputs.
"""
function gpmlncfile(date::Date,reg::AbstractString)
    return "gpml_$(reg)_prcp_$(ymd2str(date)).nc"
end

"""
    gpmlhdf5(date) -> Array{String}

    Returns an array of strings that contain the names of the raw HDF5 files for
    a given date input.  This array contains 48 different names, because GPM
    saves data every 30 minutes.
"""
function gpmlhdf5(date)

    fname = Array{String}(undef,48)
    for ii = 1 : 48
        hr = (ii-1)/2; mi = mod(hr,1); hr = hr - mi;
        hr = @sprintf("%02d",hr);
        id = @sprintf("%04d",(ii-1)*30);

        if mi == 0; fname[ii] = "3B-HHR-L.MS.MRG.3IMERG.$(ymd2str(date))-S$(hr)0000-E$(hr)2959.$(id).V06B.RT-H5"
        else;       fname[ii] = "3B-HHR-L.MS.MRG.3IMERG.$(ymd2str(date))-S$(hr)3000-E$(hr)5959.$(id).V06B.RT-H5"
        end

    end
    return fname

end

"""
    gpmlfol(date,root,reg) -> AbstractString

    Returns a string that is the path to which the extracted data will be saved.
    If folder does not exist, then the data directory will be created.
    Default value for "reg" is GLB (i.e. global).
"""
function gpmlfol(date::Date,sroot::AbstractString,reg::AbstractString="GLB")
    fol = "$(sroot)/$(reg)/$(yrmo2dir(date))/"
    if !isdir(fol)
        @info "$(Dates.now()) - GPM Near-RealTime (Late) data directory for the $(regionname(reg)) region, year $(yr2str(date)) and month $(mo2str(date)) does not exist."
        @info "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end
    return fol
end

# FTP Functions
"""
    gpmlftpcd(date,ftpID)

    Moves from the home FTP directory of the Precipitation Measurement Mission
    into the relevant GPM directory for the given date.
"""
function gpmlftpcd(date::Date,ftp)
    @info "$(Dates.now()) - Entering IMERG directory for $(ymd2str(date))."
    cd(ftp,"NRTPUB/imerg/late/$(yrmo2str(date))/")
end

# GPM Processing Functions
function gpmldt(date::Date)

    @info "$(Dates.now()) - Extracting year, month and day and time from $(date)."
    fname = gpmlhdf5(date)

    @info "$(Dates.now()) - Extracted the list of GPM files to be downloaded for $(Date(date))."
    return fname

end

function gpmlget(ftp,file,sroot::AbstractString)
    try download(ftp,"$(file)","$(sroot)/tmp/$(file)")
        @debug "$(Dates.now()) - Downloaded GPM Near-RealTime (Late) precipitation data file $(file)"
    catch; @info "$(Dates.now()) - GPM Near-RealTime (Late) precipitation data $(file) does not exist."
    end
end

function gpmldwn(date,sroot::AbstractString,overwrite=false)

    tdir = "$(sroot)/tmp/";
    if !isdir(tdir) mkpath(tdir); end

    fH5 = gpmldt(date);
    ftp = pmmnrtftpopen(); gpmlftpcd(date,ftp);
    @info "$(Dates.now()) - Downloading GPM Near-RealTime (Late) precipitation data for $(Date(date))"
    for ii = 1 : length(fH5)
        fH5ii = fH5[ii];
        if !isfile("$(sroot)/tmp/$(fH5ii)"); gpmlget(ftp,fH5ii,sroot);
        else
            if overwrite
                @info "$(Dates.now()) - GPM Near-RealTime (Late) precipitation data file $(fH5ii) already exists.  Overwriting."
                mimicget(ftp,fH5ii,sroot);
            else; @info "$(Dates.now()) - GPM Near-RealTime (Late) precipitation data file $(fH5ii) already exists.  Not overwriting."
            end
        end
    end
    pmmftpclose(ftp);

end

function gpmlextract(date::Date,sroot::AbstractString,
    reg::AbstractString="GLB")

    data = zeros(1800,3600,48) #default grid step size is 0.1x0.1 in GPM
    fH5  = gpmldt(date);

    @info "$(Dates.now()) - Extracting and compiling GPM Near-RealTime (Late) precipitation data from HDF5 files."
    for ii = 1 : length(fH5)
        fH5ii = "$(sroot)/tmp/$(fH5[ii])";
        if isfile(fH5ii)
              data[:,:,ii] = h5read("$(fH5ii)","/Grid/precipitationCal");
        else; @info "$(Dates.now()) - $(fH5ii) does not exist.  GPM Near-RealTime (Late) precipitation data values set to NaN."
              data[:,:,ii] .= NaN;
        end
    end

    @debug "$(Dates.now()) - NetCDF.jl's ncread causes memory leakage.  Using ncclose() as a workaround."
    ncclose()

    @info "$(Dates.now()) - raw GPM precipitation data is given in (lat,lon) instead of (lon,lat).  Permuting to (lon,lat)"
    data = permutedims(data,[2,1,3]);

    if reg != "GLB"
        @info "$(Dates.now()) - We do not wish to extract GPM Near-RealTime (Late) precipitation data for the entire globe."
        @info "$(Dates.now()) - Finding grid-point boundaries ..."
        lon,lat = gpmllonlat();
        bounds = regionbounds(reg); igrid = regiongrid(bounds,lon,lat);

        @info "$(Dates.now()) - Extracting GPM Near-RealTime (Late) precipitation data for the region."
        rdata,rgrid = regionextractgrid(reg,lon,lat,data)
    else; rdata = data; rgrid = [gpmllonlat()];
    end

    return rdata,rgrid

end

function gpmlsave(data,rgrid,date::Date,sroot::AbstractString,reg::AbstractString="GLB")

    fnc = gpmlncfile(date,reg);
    nlon = size(data,1); lon = rgrid[1];
    nlat = size(data,2); lat = rgrid[2];
    if nlon != length(lon); error("$(Dates.now()) - nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error("$(Dates.now()) - nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_prcp = "prcp"; att_prcp = Dict("units" => "mm/hr");
    var_lon  = "lon";  att_lon  = Dict("units" => "degree");
    var_lat  = "lat";  att_lat  = Dict("units" => "degree");

    if isfile(fnc)
        @info "$(Dates.now()) - Unfinished netCDF file $(fnc) detected.  Deleting."
        rm(fnc);
    end

    @info "$(Dates.now()) - Creating GPM Near-RealTime (Late) precipitation netCDF file $(fnc) ..."
    nccreate(fnc,var_prcp,"nlon",nlon,"nlat",nlat,"t",48,atts=att_prcp,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    @info "$(Dates.now()) - Saving GPM Near-RealTime (Late) precipitation data to netCDF file $(fnc) ..."
    ncwrite(data,fnc,var_prcp);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat);

    @debug "$(Dates.now()) - NetCDF.jl's ncread causes memory leakage.  Using ncclose() as a workaround."
    ncclose()

    fol = gpmlfol(date,sroot,reg);
    @info "$(Dates.now()) - Moving $(fnc) to data directory $(fol)"

    if isfile("$(fol)/$(fnc)"); @info "$(Dates.now()) - An older version of $(fnc) exists in the $(fol) directory.  Overwriting." end

    mv(fnc,"$(fol)/$(fnc)",force=true);

end

function gpmlrmtmp(date::Date,sroot::AbstractString)

    fH5 = gpmldt(date);
    @info "$(Dates.now()) - Deleting raw GPM Near-RealTime (Late) precipitation data files."
    for ii = 1 : length(fH5)
        fH5ii = "$(sroot)/tmp/$(fH5[ii])";
        if isfile(fH5ii); rm(fH5ii) end
    end

end

# Compiled Function
function gpmlrun(date::Date,sroot::AbstractString,reg::AbstractArray=["GLB"])
    sroot = gpmlroot(sroot); cd(sroot); gpmldwn(date,sroot);
    for regii in reg
        data,grid = gpmlextract(date,sroot,regii);
        gpmlsave(data,grid,date,sroot,regii);
    end
    gpmlrmtmp(date,sroot);
end
