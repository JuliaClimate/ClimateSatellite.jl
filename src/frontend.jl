"""
This file stores all the WIP general functions for ClimateSatellite.jl

"""

function clisatinfo!(productattr::Dict,productID::AbstractString)

    fileinfo = readlines(joinpath(@__DIR__,"../data/info.csv"))
    info = Array{Any,2}(undef,length(fileinfo)-1,10)

    for ii = 1 : length(fileinfo)-1

        row = fileinfo[ii+1];
        if row[1] != '#'; str = split(row,",");
              info[ii,1:5] .= str[1:5]; vinfo = reshape(str[6:end],6,:);
              info[ii,6]  = vinfo[1,:][1]; info[ii,7]  = vinfo[2,:][1];
              info[ii,8]  = vinfo[3,:][1]; info[ii,9]  = vinfo[4,:][1];
              info[ii,10] = vinfo[5,:][1]; info[ii,11] = vinfo[6,:][1];
        else; info[ii,:] .= "N/A";
        end

    end

    ID = (info[:,1] .== productID);

    productattr["short"]   = info[ID,1][1];
    productattr["source"]  = info[ID,2][1];
    productattr["product"] = info[ID,3][1];
    productattr["gridres"] = parse.(Int8,info[ID,4][1]);
    productattr["dayfreq"] = parse.(Int8,info[ID,5][1]);

    productattr["varID"] = info[ID,6]; productattr["standard"] = info[ID,7];
    productattr["units"] = info[ID,9]; productattr["variable"] = info[ID,8];

    productattr["scale"]  = parse.(Float64,info[ID,10]);
    productattr["offset"] = parse.(Float64,info[ID,11]);

    return

end

function clisatdwn(
    date::TimeType;
    productID::AbstractString, email::AbstractString,
    dataroot::AbstractString="", regions::Array{<:AbstractString,1}=["GLB"],
    overwrite::Bool=false
)

    if dataroot == ""; dataroot = clisatroot(productID); end

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

    return info

end

function clisatsave(
    data::Array{<:Real,3}, grid::Vector{Any},
    region::AbstractString, info::Dict, date::TimeType
)

    @info "$(Dates.now()) - Saving $(info["source"]) $(info["product"]) data for the $(regionfullname(region)) region..."

    fnc = clisatncname(info,date,region); lon,lat,t,tunit = grid;
    nlon = size(data,1); nlat = size(data,2); nt = size(data,3);
    nvar = length(info["variable"])

    if nlon != length(lon); error("$(Dates.now()) - nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error("$(Dates.now()) - nlat is $(nlat) but lat contains $(length(lat)) elements") end

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

    var_lon = "longitude"; att_lon = Dict("units"=>"degrees_east","long_name"=>"longitude");
    var_lat = "latitude";  att_lat = Dict("units"=>"degrees_north","long_name"=>"latitude");

    yrstr = @sprintf("%04d",Dates.year(date)); mostr = Dates.month(date);
    var_t = "time"; att_t = Dict("calendar"=>"gregorian","long_name"=>"time",
                                 "units"=>"$(tunit) since $(yrstr)-$(mostr)-1 0:0:0")

    if isfile(fnc)
        @info "$(Dates.now()) - Unfinished netCDF file $(fnc) detected.  Deleting."
        rm(fnc);
    end

    @debug "$(Dates.now()) - Creating $(info["source"]) $(info["product"]) netCDF file $(fnc) ..."

    for ii = 1 : nvar
        nccreate(fnc,var_var[ii],"longitude",nlon,"latitude",nlat,"time",nt,
                 atts=att_var[ii],t=NC_SHORT);
    end

    nccreate(fnc,var_lon,"longitude",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"latitude",nlat,atts=att_lat,t=NC_FLOAT);
    nccreate(fnc,var_t,"time",nlat,atts=att_t,t=NC_INT);

    @debug "$(Dates.now()) - Saving $(info["source"]) $(info["product"]) data to netCDF file $(fnc) ..."

    if nvar != 1
        for ii = 1 : nvar; ncwrite(data[:,:,:,ii],fnc,var_var[ii]); end
    else; ncwrite(data,fnc,var_var[1]);
    end

    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat);
    ncwrite(t,fnc,var_t);

    fol = clisatfol(info,date,region);
    @debug "$(Dates.now()) - Moving $(info["source"]) $(info["product"]) data file $(fnc) to data directory $(fol)"

    if isfile(joinpath(fol,fnc)); @warn "$(Dates.now()) - An older version of $(fnc) exists in the $(fol) directory.  Overwriting." end

    mv(fnc,joinpath(fol,fnc),force=true);

    @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data for the $(regionfullname(region)) region has been saved into file $(fnc) and moved to the data directory $(fol)."

end

# Root Functions
function clisatextractgrid(
    productID::AbstractString, varname::AbstractString,
    start::TimeType, finish::TimeType;
    dataroot::AbstractString="",
    region::AbstractString="GLB"
)

    if dataroot == ""; dataroot = clisatroot(productID); end

    info = Dict{Any,Any}("root"=>dataroot,"email"=>replace(email,"@"=>"%40"));
    clisatinfo!(info,productID);

    if !isdir(clisatfol(info,region))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the $(regionfullname(region))")
    end

    if sum(info["varID"] .== varname) == 0
        error("$(Dates.now()) - There is no varname identifier $(varname) for $(info["source"]) $(info["product"])")
    end

    yrs = Dates.year(start);  mos = Dates.month(start);  dys = Dates.month(start);
    yrf = Dates.year(finish); mof = Dates.month(finish); dyf = Dates.month(finish);
    ndy = Dates.value((finish-start)/Dates.day(1)); nt = info["dayfreq"];
    dvecs = Date(yrs,mos); dvecf = Date(yrf,mof);

    N,S,E,W = regionbounds(reg); res = info["gridres"];
    nlon = length(S:res:N); nlat = length(W:res:E);

    datevec = convert(Array,dvecs:Dates.month(1):dvecf); ndates = length(datevec);
    datavec = zeros(nlon,nlat,ndy*nt);

    for ii = 1 : ndates; dateii = datevec[ii];

        fol = clisatfol(info,dateii,region);
        if !isdir(fol)
            error(("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)).")
        end

        fnc = joinpath(fol,clisatncname(info,dateii,region));

        if     ii == 1 && ii != ndates
            moday = daysinmonth(dateii); ibeg = (dys-1)*nt + 1; iend = (moday+1-dys)*nt;
            datavec[:,:,1:iend] .= ncread(fnc,varname,start=[1,1,ibeg],count=[-1,-1,-1]);
        elseif ii != 1 && ii == ndates
            ibeg = iend+1; iend = dyf*nt;
            datavec[:,:,ibeg:end] .= ncread(fnc,varname,start=[1,1,1],count=[-1,-1,iend]);
        elseif ii == 1 && ii == ndates
            ibeg = (dys-1)*nt + 1; iend = dyf*nt;
            ncread!(fnc,varname,datavec,start=[1,1,ibeg],count=[-1,-1,iend])
        else
            moday = daysinmonth(dateii); ibeg = iend+1; iend = ibeg-1 + moday*nt;
            datavec[:,:,ibeg:iend] .= ncread(fnc,varname);
        end

    end

end

function clisatextractpoint(
    productID::AbstractString, varname::AbstractString,
    start::TimeType, finish::TimeType;
    coord::Array{<:Real,1}=[],
    dataroot::AbstractString="",
    region::AbstractString="GLB"
)

    if dataroot == ""; dataroot = clisatroot(productID); end

    info = Dict{Any,Any}("root"=>dataroot,"email"=>replace(email,"@"=>"%40"));
    clisatinfo!(info,productID);

    if !isdir(clisatfol(info,region))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the $(regionfullname(region))")
    end

    if sum(info["varID"] .== varname) == 0
        error("$(Dates.now()) - There is no varname identifier $(varname) for $(info["source"]) $(info["product"])")
    end

    yrs = Dates.year(start);  mos = Dates.month(start);  dys = Dates.month(start);
    yrf = Dates.year(finish); mof = Dates.month(finish); dyf = Dates.month(finish);
    ndy = Dates.value((finish-start)/Dates.day(1)); nt = info["dayfreq"];
    dvecs = Date(yrs,mos); dvecf = Date(yrf,mof);

    datevec = convert(Array,dvecs:Dates.month(1):dvecf); ndates = length(datevec);
    datavec = zeros(ndy*nt);

    for ii = 1 : ndates; dateii = datevec[ii];

        fol = clisatfol(info,dateii,region);
        if !isdir(fol)
            error(("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)).")
        end

        fnc = joinpath(fol,clisatncname(info,dateii,region));

        if     ii == 1 && ii != ndates
            moday = daysinmonth(dateii); ibeg = (dys-1)*nt + 1; iend = (moday+1-dys)*nt;
            datavec[1:iend] .= ncread(fnc,varname,start=[1,1,ibeg],count=[-1,-1,-1]);
        elseif ii != 1 && ii == ndates
            ibeg = iend+1; iend = dyf*nt;
            datavec[ibeg:end] .= ncread(fnc,varname,start=[1,1,1],count=[-1,-1,iend]);
        elseif ii == 1 && ii == ndates
            ibeg = (dys-1)*nt + 1; iend = dyf*nt;
            ncread!(fnc,varname,datavec,start=[1,1,ibeg],count=[-1,-1,iend])
        else
            moday = daysinmonth(dateii); ibeg = iend+1; iend = ibeg-1 + moday*nt;
            datavec[ibeg:iend] .= ncread(fnc,varname);
        end

    end

end
