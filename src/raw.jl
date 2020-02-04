"""
This file stores the major front-end functions of ClimateSatellite.jl that are important for
raw data downloading and extraction.

Current functionalities stored in this file include:
- Downloading of satellite data (general-purpose function that calls the specific backend)
- Saving of satellite data into NetCDF files
- Extraction of downloaded satellite data for the following
    - Domain-name specific (previously defined in ClimateEasy)
    - Point location (nearest-neighbour)
    - Within a given boundary found in a specific domain

"""

function clisatdownload(
    productID::AbstractString, date::TimeType;
    email::AbstractString, path::AbstractString="",
    regions::Array{<:AbstractString,1}=["GLB"],
    overwrite::Bool=false
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    info = Dict{Any,Any}("root"=>dataroot,"email"=>replace(email,"@"=>"%40"));
    clisatinfo!(info,productID);

    if info["source"] == "PMM"
        if     isprod(info,"gpm");  gpmdwn(regions,date,info,overwrite=overwrite);
        elseif isprod(info,"3b42"); trmmdwn(regions,date,info,overwrite=overwrite);
        end
    elseif info["source"] == "MIMIC"; mtpwdwn(regions,date,info);
    elseif info["source"] == "RSS"
        if     isprod(info,"trmm"); rtmidwn(regions,date,info,email,overwrite=overwrite);
        elseif isprod(info,"gpm");  rgmidwn(regions,date,info,email,overwrite=overwrite);
        elseif isprod(info,"smif"); rsmidwn(regions,date,info,email,overwrite=overwrite);
        elseif isprod(info,"wind"); rwnddwn(regions,date,info,email,overwrite=overwrite);
        elseif isprod(info,"amsr"); rmsrdwn(regions,date,info,email,overwrite=overwrite);
        end
    end

end

function clisatrawsave(
    data::Array{<:Real,3}, grid::Vector{Any},
    region::AbstractString, info::Dict, date::TimeType
)

    @info "$(Dates.now()) - Saving $(info["source"]) $(info["product"]) data for the $(regionfullname(region)) region..."

    fol = clisatrawfol(info,date,region);
    fnc = joinpath(fol,clisatrawname(info,date,region)); lon,lat,t,tunit = grid;
    nlon,nlat,nt = size(data); nvar = length(info["variable"]);
    yrstr = @sprintf("%04d",Dates.year(date)); mostr = Dates.month(date);

    if nlon != length(lon)
        error("$(Dates.now()) - nlon is $(nlon) but lon contains $(length(lon)) elements")
    end
    if nlat != length(lat)
        error("$(Dates.now()) - nlat is $(nlat) but lat contains $(length(lat)) elements")
    end

    var_var = Vector{AbstractString}(undef,nvar)
    att_var = map(x->Dict(),1:nvar)

    for ii = 1 : nvar
        var_var[ii]                  = info["varID"][ii];
        att_var[ii]["units"]         = info["units"][ii];
        att_var[ii]["standard_name"] = info["standard"][ii];
        att_var[ii]["long_name"]     = info["variable"][ii];
        att_var[ii]["scale_factor"]  = info["scale"][ii];
        att_var[ii]["add_offset"]    = info["offset"][ii];
        att_var[ii]["missing_value"] = -32768;
    end

    att_lon = Dict("units"=>"degrees_east","long_name"=>"longitude");
    att_lat = Dict("units"=>"degrees_north","long_name"=>"latitude");
    att_t   = Dict("calendar"=>"gregorian","long_name"=>"time",
                   "units"=>"$(tunit) since $(yrstr)-$(mostr)-1 0:0:0");

    if isfile(fnc)
        @info "$(Dates.now()) - Unfinished netCDF file $(fnc) detected.  Deleting."
        rm(fnc);
    end

    ## Write data here

    @debug "$(Dates.now()) - Creating $(info["source"]) $(info["product"]) netCDF file $(fnc) ..."

    ds = Dataset(fnc,"c");
    ds.dim["longitude"] = nlon; ds.dim["latitude"] = nlat; ds.dim["time"] = nt;

    @debug "$(Dates.now()) - Saving $(info["source"]) $(info["product"]) data to netCDF file $(fnc) ..."

    if nvar != 1
        for ii = 1 : nvar
            defVar(ds,var_var[ii],data[:,:,:,ii],("longitude","latitude","time"),
                   atts=att_var[ii]);
        end
    else;
        defVar(ds,var_var[1],data,("longitude","latitude","time"),atts=att_var[1]);
    end

    defVar(ds,"longitude",lon,("longitude",),attrib=att_lon)
    defVar(ds,"latitude",lat,("latitude",),attrib=att_lat)
    defVar(ds,"time",t,("time",),attrib=att_t)

    close(ds);

    ## Write data end

    @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data for the $(regionfullname(region)) region has been saved into file $(fnc) and moved to the data directory $(fol)."

end

# Root Functions
function clisatrawall(
    productID::AbstractString, varname::AbstractString,
    start::TimeType, finish::TimeType;
    path::AbstractString="",
    region::AbstractString="GLB",
    unpack::Bool=true
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data for the entire $(regionfullname(region)) region ..."

    if !isdir(clisatregfol(info,region))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the $(regionfullname(region))")
    end

    clisatvarinfo!(info,productID,varname=varname);
    dvec,dys,dyf,ndy = clisatextractdate(start,finish);
    nt = info["dayfreq"]; ndates = length(dvec);

    lon,lat = clisatlonlat(info); rlon,rlat,rinfo = regiongridvec(region,lon,lat);
    nlon = length(rlon); nlat = length(rlat);
    data = zeros(Int16,nlon,nlat,ndy*nt);

    for ii = 1 : ndates; dateii = dvec[ii];

        fol = clisatrawfol(info,dateii,region);
        if !isdir(fol)
            error("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)).")
        end

        fnc = joinpath(fol,clisatrawname(info,dateii,region));
        ds  = Dataset(fnc,"r"); vds = ds[varname];

        if     ii == 1 && ii != ndates
            moday = daysinmonth(dateii);
            ibeg = (dys-1)*nt + 1; iend = (moday+1-dys)*nt;
            data[:,:,1:iend] = vds.var[:,:,ibeg:end];
        elseif ii != 1 && ii == ndates
            ibeg = iend+1; iend = dyf*nt;
            data[:,:,ibeg:end] = vds.var[:,:,1:iend];
        elseif ii == 1 && ii == ndates
            ibeg = (dys-1)*nt + 1; iend = dyf*nt;
            data = vds.var[:,:,ibeg:iend];
        else
            moday = daysinmonth(dateii);
            ibeg = iend+1; iend = ibeg-1 + moday*nt;
            data[:,:,ibeg:iend] = vds.var[:];
        end

    end

    @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data for the entire $(regionfullname(region)) region has been extracted."

    if !unpack; return datavec,info,[rlon,rlat]
    else; offset = info["offset"]; scale = info["scale"];
          return datavec.*scale.+offset,info,[rlon,rlat]
    end

end

function clisatrawpoint(
    productID::AbstractString, varname::AbstractString,
    start::TimeType, finish::TimeType;
    coord::Array{<:Real,1},
    path::AbstractString="",
    region::AbstractString="GLB",
    unpack::Bool=true
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    if length(coord) != 2
        error("$(Dates.now()) - Coordinate vector must be in the form [lon,lat]")
    end

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data at coordinates $(coord) ..."

    if !isdir(clisatregfol(info,region))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the $(regionfullname(region))")
    end

    clisatvarinfo!(info,productID,varname=varname);
    dvec,dys,dyf,ndy = clisatextractdate(start,finish);
    nt = info["dayfreq"]; ndates = length(dvec);

    lon,lat = clisatlonlat(info); rlon,rlat,rinfo = regiongridvec(region,lon,lat);
    plon,plat = coord; ilon,ilat = regionpoint(plon,plat,rlon,rlat);
    data = zeros(Int16,ndy*nt);

    for ii = 1 : ndates; dateii = dvec[ii];

        fol = clisatrawfol(info,dateii,region);
        if !isdir(fol)
            error("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)).")
        end

        fnc = joinpath(fol,clisatrawname(info,dateii,region));
        ds  = Dataset(fnc,"r"); vds = ds[varname];

        if     ii == 1 && ii != ndates
            moday = daysinmonth(dateii);
            ibeg = (dys-1)*nt + 1; iend = (moday+1-dys)*nt;
            data[1:iend] = vds.var[ilon,ilat,ibeg:end];
        elseif ii != 1 && ii == ndates
            ibeg = iend+1; iend = dyf*nt;
            data[ibeg:end] = vds.var[ilon,ilat,1:iend];
        elseif ii == 1 && ii == ndates
            ibeg = (dys-1)*nt + 1; iend = dyf*nt;
            data = vds.var[ilon,ilat,ibeg:iend];
        else
            moday = daysinmonth(dateii);
            ibeg = iend+1; iend = ibeg-1 + moday*nt;
            data[ibeg:iend] = vds.var[ilon,ilat,:];
        end

    end

    @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data for the coordinates $(coord) has been extracted."

    if !unpack; return datavec,info
    else; offset = info["offset"]; scale = info["scale"];
          return datavec.*scale.+offset,info
    end

end

function clisatrawgrid(
    productID::AbstractString, varname::AbstractString,
    start::TimeType, finish::TimeType;
    grid::Array{<:Real,1},
    path::AbstractString="",
    region::AbstractString="GLB",
    unpack::Bool=true
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    if length(grid) != 4
        error("$(Dates.now()) - Grid vector must be in the form [N,S,E,W]")
    end

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data for the [N,S,E,W] bounds $(grid) ..."

    if !isdir(clisatregfol(info,region))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the $(regionfullname(region))")
    end

    clisatvarinfo!(info,productID,varname=varname);
    dvec,dys,dyf,ndy = clisatextractdate(start,finish);
    nt = info["dayfreq"]; ndates = length(dvec);

    lon,lat = clisatlonlat(info); isgridinregion(grid,region);
    rlon,rlat,rinfo = regiongridvec(region,lon,lat);
    glon,glat,ginfo = regiongridvec(grid,rlon,rlat); iWE,iNS = ginfo["IDvec"];
    nlon = length(glon); nlat = length(glat);
    data = zeros(Int16,nlon,nlat,ndy*nt);

    for ii = 1 : ndates; dateii = dvec[ii];

        fol = clisatrawfol(info,dateii,region);
        if !isdir(fol)
            error("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)).")
        end

        fnc = joinpath(fol,clisatrawname(info,dateii,region));
        ds  = Dataset(fnc,"r"); vds = ds[varname];

        if     ii == 1 && ii != ndates
            moday = daysinmonth(dateii);
            ibeg = (dys-1)*nt + 1; iend = (moday+1-dys)*nt;
            data[:,:,1:iend] = vds.var[iWE,iNS,ibeg:end];
        elseif ii != 1 && ii == ndates
            ibeg = iend+1; iend = dyf*nt;
            data[:,:,ibeg:end] = vds.var[iWE,iNS,1:iend];
        elseif ii == 1 && ii == ndates
            ibeg = (dys-1)*nt + 1; iend = dyf*nt;
            data = vds.var[iWE,iNS,ibeg:iend];
        else
            moday = daysinmonth(dateii);
            ibeg = iend+1; iend = ibeg-1 + moday*nt;
            data[:,:,ibeg:iend] = vds.var[iWE,iNS,:];
        end

    end

    @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data within the [N,S,E,W] bounds $(grid) has been extracted."

    if !unpack; return datavec,info,[glon,glat]
    else; offset = info["offset"]; scale = info["scale"];
          return datavec.*scale.+offset,info,[glon,glat]
    end

end
