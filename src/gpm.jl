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
    gpmroot(root) -> AbstractString

    Returns the path where GPM data will be downloaded and stored based on the
    given input "root", which must be a string
"""
function gpmroot(root::AbstractString)
    infl("Making root folder for GPM datasets.")
    sroot = "$(root)/GPM/"; mkpath(sroot)
    if !isdir(sroot)
        infl("GPM directory does not exist.  Creating now.");
        mkpath(sroot);
    end
    infl("The root folder for the GPM precipitation data is in $(sroot).")
    return sroot
end

"""
    gpmlonlat() -> Array,Array

    Returns two vector arrays for the longitude and latitude respectively.
    Longitude defaults are within [-180,180].  Grid spacing 0.1 x 0.1.
    Measurements are taken in the middle of grid points, therefore at .05º
    No inputs are required.
"""
function gpmlonlat()
    lon = convert(Array,-179.95:0.1:179.95); lat = convert(Array,-89.95:0.1:89.95);
    return lon,lat
end

"""
    gpmncfile(date,reg) -> AbstractString

    Returns the ncfile name that the data extracted will be saved to, based on
    given date and region variable inputs.
"""
function gpmncfile(date::Date,reg::AbstractString)
    return "gpm_$(reg)_prcp_$(Date(date)).nc"
end

"""
    gpmhdf5(date) -> Array{String}

    Returns an array of strings that contain the names of the raw HDF5 files for
    a given date input.  This array contains 48 different names, because GPM
    saves data every 30 minutes.
"""
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

"""
    gpmfol(date,root,reg) -> AbstractString

    Returns a string that is the path to which the extracted data will be saved.
    If folder does not exist, then the data directory will be created.
    Default value for "reg" is GLB (i.e. global).
"""
function gpmfol(date::Date,sroot::AbstractString,reg::AbstractString="GLB")
    fol = "$(sroot)/$(reg)/$(yrmo2dir(date))/"
    if !isdir(fol)
        infl("GPM data directory for the $(regionname(reg)) region, year $(yr2str(date)) and month $(mo2str(date)) does not exist.");
        infl("Creating data directory $(fol)."); mkpath(fol);
    end
    return fol
end

# FTP Functions
"""
    gpmftpcd(date,ftpID)

    Moves from the home FTP directory of the Precipitation Measurement Mission
    into the relevant GPM directory for the given date.
"""
function gpmftpcd(date::Date,ftp)
    infl("Entering IMERG directory for $(ymd2str(date)).")
    cd(ftp,"gpmdata/$(ymd2dir(date))/imerg")
end

# GPM Processing Functions
function gpmdt(date::Date)

    infl("Extracting year, month and day and time from $(date).")
    fname = gpmhdf5(date)

    infl("Extracted the list of GPM files to be downloaded for $(Date(date)).")
    return fname

end

function gpmget(ftp,file,sroot::AbstractString)
    try download(ftp,"$(file)","$(sroot)/tmp/$(file)")
        debl("Downloaded GPM precipitation data file $(file)")
    catch; infl("GPM precipitation data $(file) does not exist.")
    end
end

function gpmdwn(date,sroot::AbstractString,overwrite=false)

    tdir = "$(sroot)/tmp/";
    if !isdir(tdir) mkpath(tdir); end

    fH5 = gpmdt(date);
    ftp = ppmftpopen(); gpmftpcd(date,ftp);
    infl("Downloading GPM precipitation data for $(Date(date))")
    for ii = 1 : length(fH5)
        fH5ii = fH5[ii];
        if !isfile("$(sroot)/tmp/$(fH5ii)"); gpmget(ftp,fH5ii,sroot);
        else
            if overwrite
                infl("GPM precipitation data file $(fH5ii) already exists.  Overwriting.")
                mimicget(ftp,fH5ii,sroot);
            else; infl("GPM precipitation data file $(fH5ii) already exists.  Not overwriting.")
            end
        end
    end
    ppmftpclose(ftp);

end

function gpmextract(date::Date,sroot::AbstractString,
    reg::AbstractString="GLB")

    data = zeros(1800,3600,48) #default grid step size is 0.1x0.1 in GPM
    fH5  = gpmdt(date);

    infl("Extracting and compiling GPM precipitation data from HDF5 files.")
    for ii = 1 : length(fH5)
        fH5ii = "$(sroot)/tmp/$(fH5[ii])";
        if isfile(fH5ii)
              data[:,:,ii] = h5read("$(fH5ii)","/Grid/precipitationCal");
        else; infl("$(fH5ii) does not exist.  GPM precipitation data values set to NaN.")
              data[:,:,ii] .= NaN;
        end
    end

    infl("raw GPM precipitation data is given in (lat,lon) instead of (lon,lat).  Permuting to (lon,lat)")
    data = permutedims(data,[2,1,3]);

    if reg != "GLB"
        infl("We do not wish to extract GPM precipitation data for the entire globe.")
        infl("Finding grid-point boundaries ...")
        lon,lat = gpmlonlat();
        bounds = regionbounds(reg); igrid = regiongrid(bounds,lon,lat);

        infl("Extracting GPM precipitation data for the region.")
        rdata,rpnts = regionextractgrid(reg,lon,lat,data)
    else; rdata = data; lon,lat = gpmlonlat(); rpnts = [lat[end],lat[1],lon[end],lon[1]];
    end

    return rdata,rpnts

end

function gpmsave(data,rpnts,date::Date,sroot::AbstractString,reg::AbstractString="GLB")

    fnc = gpmncfile(date,reg);
    nlon = size(data,1); lon = convert(Array,rpnts[4]:0.1:rpnts[3]);
    nlat = size(data,2); lat = convert(Array,rpnts[2]:0.1:rpnts[1]);
    if nlon != length(lon); errl("nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); errl("nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_prcp = "prcp"; att_tpw = Dict("units" => "mm/hr");
    var_lon  = "lon";  att_lon = Dict("units" => "degree");
    var_lat  = "lat";  att_lat = Dict("units" => "degree");

    if isfile(fnc)
        infl("Unfinished netCDF file $(fnc) detected.  Deleting.");
        rm(fnc);
    end

    infl("Creating GPM precipitation netCDF file $(fnc) ...")
    nccreate(fnc,var_prcp,"nlon",nlon,"nlat",nlat,"t",48,atts=att_prcp,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    infl("Saving GPM precipitation data to netCDF file $(fnc) ...")
    ncwrite(data,fnc,var_prcp);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat);

    fol = gpmfol(date,sroot);
    infl("Moving $(fnc) to data directory $(fol)")

    if isfile("$(fol)/$(fnc)"); infl("An older version of $(fnc) exists in the $(fol) directory.  Overwriting.") end

    mv(fnc,"$(fol)/$(fnc)",force=true);

end

function gpmrmtmp(date::Date,sroot::AbstractString)

    fH5 = gpmdt(date);
    infl("Deleting raw GPM precipitation data files.")
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