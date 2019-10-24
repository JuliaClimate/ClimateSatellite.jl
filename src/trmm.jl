"""
    This file contains all the functions that are used in the retrieval of TRMM
    precipitable water datasets.  Current functionalities include:
        - downloading of raw data files
        - extraction of regional data (where applicable)
        - saving of data (global/regional)

    The following functionalities are forthcoming:
        - preliminary mapping of data
"""

# Setup Functions

"""
    trmmroot(root) -> AbstractString

    Returns the path where TRMM data will be downloaded and stored based on the
    given input "root", which must be a string
"""
function trmmroot(root::AbstractString)
    @info "$(Dates.now()) - Making root folder for TRMM datasets."
    sroot = "$(root)/TRMM/"; mkpath(sroot)
    if !isdir(sroot)
        @info "$(Dates.now()) - TRMM directory does not exist.  Creating now."
        mkpath(sroot);
    end
    @info "$(Dates.now()) - The root folder for the TRMM precipitation data is in $(sroot)."
    return sroot
end

"""
    trmmlonlat() -> Array,Array

    Returns two vector arrays for the longitude and latitude respectively.
    Longitude defaults are within [-180,180].  Grid spacing 0.25 x 0.25.
    Measurements are taken in the middle of grid points.
    No inputs are required.
"""
function trmmlonlat()
    lon = convert(Array,-179.875:0.25:179.875);
    lat = convert(Array,49.875:-0.25:-49.875);
    return lon,lat
end

"""
    trmmncfile(date,reg) -> AbstractString

    Returns the ncfile name that the data extracted will be saved to, based on
    given date and region variable inputs.
"""
function trmmncfile(date::Date,reg::AbstractString)
    return "trmm_$(reg)_prcp_$(ymd2str(date)).nc"
end

"""
    trmmhdf5(date) -> Array{String}

    Returns an array of strings that contain the names of the raw HDF5 files for
    a given date input.  This array contains 48 different names, because TRMM
    saves data every 30 minutes.
"""
function trmmhdf(date)

    fname = Array{String}(undef,8)
    if date > Date(2010,9,30)
        for ii = 1 : 8
            hr = (ii-1)*3; hr = @sprintf("%02d",hr);
            fname[ii] = "3B42.$(ymd2str(date)).$(hr).7A.HDF"
        end
    else
        for ii = 1 : 8
            hr = (ii-1)*3; hr = @sprintf("%02d",hr);
            fname[ii] = "3B42.$(ymd2str(date)).$(hr).7.HDF"
        end
    end
    return fname

end

"""
    trmmfol(date,root,reg) -> AbstractString

    Returns a string that is the path to which the extracted data will be saved.
    If folder does not exist, then the data directory will be created.
    Default value for "reg" is GLB (i.e. global).
"""
function trmmfol(date::Date,sroot::AbstractString,reg::AbstractString="GLB")
    fol = "$(sroot)/$(reg)/$(yrmo2dir(date))/"
    if !isdir(fol)
        @info "$(Dates.now()) - TRMM data directory for the $(regionname(reg)) region, year $(yr2str(date)) and month $(mo2str(date)) does not exist."
        @info "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end
    return fol
end

# FTP Functions
"""
    trmmftpcd(date,ftpID)

    Moves from the home FTP directory of the Precipitation Measurement Mission
    into the relevant TRMM directory for the given date.
"""
function trmmftpcd(date::Date,ftp)
    ftpdir = "trmmdata/ByDate/V07/$(ymd2dir(date))/";
    @info "$(Dates.now()) - Entering TRMM directory $(ftpdir) for $(ymd2str(date))."
    cd(ftp,"trmmdata/ByDate/V07/$(ymd2dir(date))/")
end

# TRMM Processing Functions
function trmmdt(date::Date)

    @info "$(Dates.now()) - Extracting year, month and day and time from $(date)."
    fname = trmmhdf(date)

    @info "$(Dates.now()) - Extracted the list of TRMM files to be downloaded for $(Date(date))."
    return fname

end

function trmmget(ftp,file,sroot::AbstractString)
    try download(ftp,"$(file).gz","$(sroot)/tmp/$(file).gz")
        @debug "$(Dates.now()) - Downloaded TRMM precipitation data file $(file).gz"
        trmmgunzip(file,sroot);
    catch; @info "$(Dates.now()) - TRMM precipitation data $(file).gz does not exist."
    end
end

function trmmgunzip(file,sroot::AbstractString)
    if isfile("$(sroot)/tmp/$(file)");
        @debug "$(Dates.now()) - TRMM precipitation HDF4 file $(sroot)/tmp/$(file) exists.  Overwriting."
        rm("$(sroot)/tmp/$(file)")
    end
    @debug "$(Dates.now()) - Unzipping tar.gz file $(file).gz into $(file)"
    run(pipeline(`gunzip $(sroot)/tmp/$(file).gz`));
end

function trmmh4read(file,variable)
    @debug "$(Dates.now()) - Extracting data from HDF4 file $(file) into a python object."
    dobj = h4read(file,variable); data = zeros(1440,400)
    @debug "$(Dates.now()) - Extracting data from a python object into numeric array."
    for ii = 0 : 1439
        for jj = 0 : 399
            data[ii+1,jj+1] = get(dobj,(ii,jj));
        end
    end
    return data
end

function trmmdwn(date,sroot::AbstractString,overwrite=false)

    tdir = "$(sroot)/tmp/";
    if !isdir(tdir) mkpath(tdir); end

    fHDF = trmmdt(date);
    ftp  = pmmfinftpopen(); trmmftpcd(date,ftp);
    @info "$(Dates.now()) - Downloading TRMM precipitation data for $(Date(date))"
    for ii = 1 : length(fHDF)
        fHii = fHDF[ii];
        if !isfile("$(sroot)/tmp/$(fHii).gz"); trmmget(ftp,fHii,sroot);
        else
            if overwrite
                @info "$(Dates.now()) - TRMM precipitation data .gz file $(fHii).gz already exists.  Overwriting."
                trmmget(ftp,fHii,sroot);
            else; @info "$(Dates.now()) - TRMM precipitation data .gz file $(fHii).gz already exists.  Not overwriting."
                trmmgunzip(fHii,sroot);
            end
        end
    end
    pmmftpclose(ftp);

end

function trmmextract(date::Date,sroot::AbstractString,
    reg::AbstractString="GLB")

    data = zeros(1440,400,8) #default grid step size is 0.1x0.1 in TRMM
    fHDF = trmmdt(date);

    @info "$(Dates.now()) - Extracting and compiling TRMM precipitation data from HDF5 files."
    for ii = 1 : length(fHDF)
        fHii = "$(sroot)/tmp/$(fHDF[ii])";
        if isfile(fHii)
              data[:,:,ii] = trmmh4read(fHii,"precipitation");
        else; @info "$(Dates.now()) - $(fHii) does not exist.  TRMM precipitation data values set to NaN."
              data[:,:,ii] .= NaN;
        end
    end

    @debug "$(Dates.now()) - NetCDF.jl's ncread causes memory leakage.  Using ncclose() as a workaround."
    ncclose()

    if reg != "GLB"
        @info "$(Dates.now()) - We do not wish to extract TRMM precipitation data for the entire globe."
        @info "$(Dates.now()) - Finding grid-point boundaries ..."
        lon,lat = trmmlonlat();
        bounds = regionbounds(reg); igrid = regiongrid(bounds,lon,lat);

        @info "$(Dates.now()) - Extracting TRMM precipitation data for the region."
        rdata,rgrid = regionextractgrid(data,reg,lon,lat)
    else; rdata = data; rgrid = [trmmlonlat()];
    end

    return rdata,rgrid

end

function trmmsave(data,rgrid,date::Date,sroot::AbstractString,reg::AbstractString="GLB")

    fnc = trmmncfile(date,reg);
    nlon = size(data,1); lon = rgrid[1];
    nlat = size(data,2); lat = rgrid[2];
    if nlon != length(lon); error("$(Dates.now()) - nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error("$(Dates.now()) - nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_prcp = "prcp"; att_prcp = Dict("units" => "mm/hr");
    var_lon  = "lon";  att_lon  = Dict("units" => "degree");
    var_lat  = "lat";  att_lat  = Dict("units" => "degree");

    if isfile(fnc)
        @info "$(Dates.now()) - Unfinished netCDF file $(fnc) detected.  Deleting.";
        rm(fnc);
    end

    @info "$(Dates.now()) - Creating TRMM precipitation netCDF file $(fnc) ..."
    nccreate(fnc,var_prcp,"nlon",nlon,"nlat",nlat,"t",8,atts=att_prcp,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    @info "$(Dates.now()) - Saving TRMM precipitation data to netCDF file $(fnc) ..."
    ncwrite(data,fnc,var_prcp);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat);

    @debug "$(Dates.now()) - NetCDF.jl's ncread causes memory leakage.  Using ncclose() as a workaround."
    ncclose()

    fol = trmmfol(date,sroot,reg);
    @info "$(Dates.now()) - Moving $(fnc) to data directory $(fol)"

    if isfile("$(fol)/$(fnc)"); @info "$(Dates.now()) - An older version of $(fnc) exists in the $(fol) directory.  Overwriting." end

    mv(fnc,"$(fol)/$(fnc)",force=true);

end

function trmmrmtmp(date::Date,sroot::AbstractString)

    fHDF = trmmdt(date);
    @info "$(Dates.now()) - Deleting raw TRMM precipitation data files."
    for ii = 1 : length(fHDF)
        fHii = "$(sroot)/tmp/$(fHDF[ii])";    if isfile(fHii); rm(fHii) end
        fHii = "$(sroot)/tmp/$(fHDF[ii]).gz"; if isfile(fHii); rm(fHii) end
    end

end

# Compiled Function
function trmmrun(date::Date,sroot::AbstractString,reg::AbstractArray=["GLB"])
    sroot = trmmroot(sroot); cd(sroot); trmmdwn(date,sroot);
    for regii in reg
        data,grid = trmmextract(date,sroot,regii);
        trmmsave(data,grid,date,sroot,regii);
    end
    trmmrmtmp(date,sroot);
end
