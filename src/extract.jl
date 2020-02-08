function clisatsubregion(
    productID::AbstractString, date::TimeType;
    path::AbstractString="",
    region::AbstractString="GLB"
    overwrite::Bool=false
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);

    parent = regionparent(region);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data for the entire $(regionfullname(region)) region ..."

    if !isdir(clisatregfol(info,parent))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the parent region: $(regionfullname(parent))")
    end

    folp = clisatrawfol(info,date,parent);
    fncp = joinpath(fol,clisatrawname(info,date,parent));
    if !isdir(fol)
        error("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)) in the parent region: $(regionfullname(parent))")
    end

    pds  = Dataset(fncp,"r");
    plon = ds["longitude"].var[:]; plat = ds["latitude"].var[:]; time = ds["time"].var[:];

    rlon,rlat,rinfo = regiongridvec(region,plon,plat); iWE,iNS = rinfo["IDvec"];
    nlon = length(rlon); nlat = length(rlat); nt = length(time);
    nvar = length(info["variable"]);

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

    att_lon = ds["longitude"].attrib;
    att_lat = ds["latitude"].attrib;
    att_t   = ds["time"].attrib;

    folr = clisatrawfol(info,date,region);
    fncr = joinpath(fol,clisatrawname(info,date,region));

    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fncr) detected.  Overwriting ..."
        rm(fnc);
    end

    ## Write data here

    @debug "$(Dates.now()) - Creating $(info["source"]) $(info["product"]) netCDF file $(fnc) ..."

    ds = Dataset(fncr,"c");
    ds.dim["longitude"] = nlon; ds.dim["latitude"] = nlat; ds.dim["time"] = nt;

    @debug "$(Dates.now()) - Saving $(info["source"]) $(info["product"]) data to netCDF file $(fnc) ..."

    if nvar != 1
        for ii = 1 : nvar
            v = defVar(ds,var_var[ii],Int16,("longitude","latitude","time"),
                       attrib=att_var[ii]);
            v.var[:] = ds[var_var[ii]].var[iWE,iNS,:];
        end
    else;
        v = defVar(ds,var_var[1],Int16,("longitude","latitude","time"),attrib=att_var[1]);
        v.var[:] = ds[var_var[1]].var[iWE,iNS,:];
    end

    defVar(ds,"longitude",lon,("longitude",),attrib=att_lon)
    defVar(ds,"latitude",lat,("latitude",),attrib=att_lat)
    defVar(ds,"time",t,("time",),attrib=att_t)

    close(ds); close(pds);

    ## Write data end

    @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data for the $(regionfullname(region)) region has been saved into file $(fncr) and moved to the data directory $(folr)."

end
