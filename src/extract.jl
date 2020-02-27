function clisatsubregion(
    productID::AbstractString, date::TimeType;
    path::AbstractString="",
    region::AbstractString="GLB",
    overwrite::Bool=false
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);

    parent = gregionparent(region);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data for the entire $(gregionfullname(region)) region ..."

    if !isdir(clisatregfol(info,parent))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the parent region: $(gregionfullname(parent))")
    end

    folp = clisatrawfol(info,date,parent);
    fncp = joinpath(folp,clisatrawname(info,date,parent));
    if !isdir(folp)
        error("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)) in the parent region: $(gregionfullname(parent))")
    end

    pds  = Dataset(fncp,"r");
    plon = pds["longitude"].var[:]; plat = pds["latitude"].var[:];
    time = pds["time"].var[:];

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

    att_lon = pds["longitude"].attrib;
    att_lat = pds["latitude"].attrib;
    att_t   = pds["time"].attrib;

    folr = clisatrawfol(info,date,region);
    fncr = joinpath(folr,clisatrawname(info,date,region));

    if isfile(fncr)
        @info "$(Dates.now()) - Stale NetCDF file $(fncr) detected.  Overwriting ..."
        rm(fncr);
    end

    ## Write data here

    @debug "$(Dates.now()) - Creating $(info["source"]) $(info["product"]) netCDF file $(fncr) ..."

    ds = Dataset(fncr,"c");
    ds.dim["longitude"] = nlon; ds.dim["latitude"] = nlat; ds.dim["time"] = nt;

    @debug "$(Dates.now()) - Saving $(info["source"]) $(info["product"]) data to netCDF file $(fncr) ..."

    if nvar != 1
        for ii = 1 : nvar
            v = defVar(ds,var_var[ii],Int16,("longitude","latitude","time"),
                       attrib=att_var[ii]);
            v.var[:] = pds[var_var[ii]].var[iWE,iNS,:];
        end
    else;
        v = defVar(ds,var_var[1],Int16,("longitude","latitude","time"),attrib=att_var[1]);
        v.var[:] = pds[var_var[1]].var[iWE,iNS,:];
    end

    defVar(ds,"longitude",rlon,("longitude",),attrib=att_lon)
    defVar(ds,"latitude",rlat,("latitude",),attrib=att_lat)
    defVar(ds,"time",time,("time",),attrib=att_t)

    close(ds); close(pds);

    ## Write data end

    @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data for the $(gregionfullname(region)) region has been saved into file $(fncr) and moved to the data directory $(folr)."

end
