"""
This file stores the analysis functions of ClimateSatellite.jl that are important for
climatological analysis of the downloaded data.

Current functionalities stored in this file include:
- Yearly, seasonal and monthly means, maximum, minimum, range and standard deviation:
    - On a general climatological basis
    - For a specified hour
    - Of the diurnal variability
- Zonal and meridional averages of the above
- Saving the results of the abovementioned analysis into yearly NetCDF files
- Retrieval of these data from their respective NetCDF files

Functionalities that are in development include:
- Zonal, meridional and domain averages

Because there is progressively more data every year, analysis will be conducted on a yearly
basis instead of for the entire time-series.  This is to ensure generality and ease-of-use
as the datasets become larger and larger.

"""

function clisatanalysis(
    productID::AbstractString, yr::Integer;
    varname::AbstractString,
    path::AbstractString="", region::AbstractString="GLB"
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);
    clisatvarinfo!(info,productID,varname=varname);

    regfol = clisatregfol(info,region);
    rawfol = clisatrawfol(info,Date(yr),region);
    anafol = clisatanafol(info,Date(yr),region);
    cd(rawfol);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data in $(regionfullname(region)) region during $yr ..."

    lon,lat = clisatlonlat(info); rlon,rlat,rinfo = regiongridvec(region,lon,lat);
    nlon = length(rlon); nlat = length(rlat); nt = info["dayfreq"]+1; grid = [rlon,rlat];

    davg = zeros(Int16,nlon,nlat,nt+1,13); dstd = zeros(Int16,nlon,nlat,nt+1,13);
    dmax = zeros(Int16,nlon,nlat,nt+1,13); dmin = zeros(Int16,nlon,nlat,nt+1,13);

    for mo = 1 : 12; ndy = daysinmonth(yr,mo)

        @info "$(Dates.now()) - Analyzing $(info["source"]) $(info["product"]) data in $(regionfullname(region)) during $(Dates.monthname(mo)) $yr ..."
        ncraw = clisatrawname(productID,Date(yr,mo),region);
        ds = Dataset(ncraw,"r"); vds = ds[varname];
        raw = reshape(vds.var[:],nlon,nlat,(nt-1),ndy);

        @debug "$(Dates.now()) - Extracting monthly diurnal climatological information ..."
        for it = 1 : nt-1, ilat = 1 : nlat, ilon = 1 : nlon

            rawii = raw[ilon,ilat,it,:]; rawii = rawii[rawii.!=-32768];
            davg[ilon,ilat,it,mo] = round(Int16,mean(rawii));
            dstd[ilon,ilat,it,mo] = round(Int16,std(rawii));
            dmax[ilon,ilat,it,mo] = round(Int16,maximum(rawii));
            dmin[ilon,ilat,it,mo] = round(Int16,minimum(rawii));

        end

        @debug "$(Dates.now()) - Extracting monthly averaged climatological information ..."
        davg[:,:,nt,mo] = round.(Int16,mean(davg[:,:,1:nt-1,mo],dims=3));
        dstd[:,:,nt,mo] = round.(Int16,mean(dstd[:,:,1:nt-1,mo],dims=3));
        dmax[:,:,nt,mo] = round.(Int16,maximum(dmax[:,:,1:nt-1,mo],dims=3));
        dmin[:,:,nt,mo] = round.(Int16,minimum(dmin[:,:,1:nt-1,mo],dims=3));

        @debug "$(Dates.now()) - Permuting days and hours dimensions ..."
        raw = permutedims(raw,(1,2,4,3)); tmp = zeros(nlon,nlat,ndy);

        @debug "$(Dates.now()) - Extracting diurnal variability information ..."
        for idy = 1 : ndy, ilat = 1 : nlat, ilon = 1 : nlon

            tmpii = raw[ilon,ilat,idy,:]; tmpii = tmpii[tmpii.!=-32768];
            tmp[ilon,ilat,idy] = maximum(tmpii)/2 - minimum(tmpii)/2;

        end

        @debug "$(Dates.now()) - Extracting monthly diurnal variability information ..."
        davg[:,:,nt+1,mo] = round.(Int16,mean(tmp,dims=3));
        dstd[:,:,nt+1,mo] = round.(Int16,std(tmp,dims=3));
        dmax[:,:,nt+1,mo] = round.(Int16,maximum(tmp,dims=3));
        dmin[:,:,nt+1,mo] = round.(Int16,minimum(tmp,dims=3));

    end

    @info "$(Dates.now()) - Calculating yearly climatology for $(info["source"]) $(info["product"]) in $(regionfullname(region)) during $yr ..."
    davg[:,:,:,end] = round.(Int16,mean(davg[:,:,:,1:12],dims=4));
    dstd[:,:,:,end] = round.(Int16,mean(dstd[:,:,:,1:12],dims=4));
    dmax[:,:,:,end] = round.(Int16,maximum(dmax[:,:,:,1:12],dims=4));
    dmin[:,:,:,end] = round.(Int16,minimum(dmin[:,:,:,1:12],dims=4));

    clisatanasave([davg,dstd,dmax,dmin],grid,productID,varname,yr,path,region,info);

end

function clisatanasave(
    data::Array{Array{Int16,4},1}, grid::Vector{<:Any},
    productID::AbstractString, varname::AbstractString, yr::Integer,
    path::AbstractString, region::AbstractString, info::Dict
)

    @info "$(Dates.now()) - Saving analysed $(info["source"]) $(info["product"]) data in $(regionfullname(region)) for the year $yr ..."

    fol = clisatanafol(productID,Date(yr),region,path=path)
    fnc = joinpath(fol,clisatananame(productID,varname,Date(yr),region));
    rlon,rlat = grid; nlon = length(rlon); nlat = length(rlat); nt = info["dayfreq"];

    if isfile(fnc)
        @info "$(Dates.now()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end

    @debug "$(Dates.now()) - Creating NetCDF file $(fnc) for analyzed $(info["source"]) $(info["product"]) data in $yr ..."

    ds = Dataset(fnc,"c");
    ds.dim["longitude"] = nlon; ds.dim["latitude"] = nlat;
    ds.dim["hour"] = nt; ds.dim["month"] = 12;

    att_var = map(x->Dict(),1:2)
    for ii = 1 : 2
        att_var[ii]["units"]         = info["units"];
        att_var[ii]["standard_name"] = info["standard"];
        att_var[ii]["long_name"]     = info["variable"];
        att_var[ii]["missing_value"] = -32768;
    end
    att_var[1]["scale_factor"] = info["scale"];   att_var[1]["add_offset"] = info["offset"];
    att_var[2]["scale_factor"] = info["scale"]*2; att_var[2]["add_offset"] = 0;

    att_lon = Dict("units"=>"degrees_east","long_name"=>"longitude");
    att_lat = Dict("units"=>"degrees_north","long_name"=>"latitude");

    defVar(ds,"longitude",rlon,("longitude",),attrib=att_lon)
    defVar(ds,"latitude",rlat,("latitude",),attrib=att_lat)

    @debug "$(Dates.now()) - Saving analyzed $(info["source"]) $(info["product"]) data for $yr to NetCDF file $(fnc) ..."

    v = defVar(ds,"domain_yearly_mean_climatology",Int16,
               ("longitude","latitude"),attrib=att_var[1]);
    v.var[:] = data[1][:,:,nt+1,end];

    v = defVar(ds,"domain_yearly_std_climatology",Int16,
               ("longitude","latitude"),attrib=att_var[1]);
    v.var[:] = data[2][:,:,nt+1,end];

    v = defVar(ds,"domain_yearly_maximum_climatology",Int16,
               ("longitude","latitude"),attrib=att_var[1]);
    v.var[:] = data[3][:,:,nt+1,end];

    v = defVar(ds,"domain_yearly_minimum_climatology",Int16,
               ("longitude","latitude"),attrib=att_var[1]);
    v.var[:] = data[4][:,:,nt+1,end];


    v = defVar(ds,"domain_yearly_mean_hourly",Int16,
               ("longitude","latitude","hour"),attrib=att_var[1]);
    v.var[:] = data[1][:,:,1:nt,end];

    v = defVar(ds,"domain_yearly_std_hourly",Int16,
           ("longitude","latitude","hour"),attrib=att_var[1]);
    v.var[:] = data[2][:,:,1:nt,end];

    v = defVar(ds,"domain_yearly_maximum_hourly",Int16,
               ("longitude","latitude","hour"),attrib=att_var[1]);
    v.var[:] = data[3][:,:,1:nt,end];

    v = defVar(ds,"domain_yearly_minimum_hourly",Int16,
               ("longitude","latitude","hour"),attrib=att_var[1]);
    v.var[:] = data[4][:,:,1:nt,end];


    v = defVar(ds,"domain_yearly_mean_diurnalvariance",Int16,
           ("longitude","latitude"),attrib=att_var[2]);
    v.var[:] = data[1][:,:,nt+2,end];

    v = defVar(ds,"domain_yearly_std_diurnalvariance",Int16,
               ("longitude","latitude"),attrib=att_var[2]);
    v.var[:] = data[2][:,:,nt+2,end];

    v = defVar(ds,"domain_yearly_maximum_diurnalvariance",Int16,
               ("longitude","latitude"),attrib=att_var[2]);
    v.var[:] = data[3][:,:,nt+2,end];

    v = defVar(ds,"domain_yearly_minimum_diurnalvariance",Int16,
           ("longitude","latitude"),attrib=att_var[2]);
    v.var[:] = data[4][:,:,nt+2,end];


    v = defVar(ds,"domain_monthly_mean_climatology",Int16,
               ("longitude","latitude","month"),attrib=att_var[1]);
    v.var[:] = data[1][:,:,nt+1,1:12];

    v = defVar(ds,"domain_monthly_std_climatology",Int16,
               ("longitude","latitude","month"),attrib=att_var[1]);
    v.var[:] = data[2][:,:,nt+1,1:12];

    v = defVar(ds,"domain_monthly_maximum_climatology",Int16,
               ("longitude","latitude","month"),attrib=att_var[1]);
    v.var[:] = data[3][:,:,nt+1,1:12];

    v = defVar(ds,"domain_monthly_minimum_climatology",Int16,
               ("longitude","latitude","month"),attrib=att_var[1]);
    v.var[:] = data[4][:,:,nt+1,1:12];


    v = defVar(ds,"domain_monthly_mean_hourly",Int16,
               ("longitude","latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = data[1][:,:,1:nt,1:12];

    v = defVar(ds,"domain_monthly_std_hourly",Int16,
               ("longitude","latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = data[2][:,:,1:nt,1:12];

    v = defVar(ds,"domain_monthly_maximum_hourly",Int16,
               ("longitude","latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = data[3][:,:,1:nt,1:12];

    v = defVar(ds,"domain_monthly_minimum_hourly",Int16,
               ("longitude","latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = data[4][:,:,1:nt,1:12];


    v = defVar(ds,"domain_monthly_mean_diurnalvariance",Int16,
               ("longitude","latitude","month"),attrib=att_var[2]);
    v.var[:] = data[1][:,:,nt+2,1:12];

    v = defVar(ds,"domain_monthly_std_diurnalvariance",Int16,
               ("longitude","latitude","month"),attrib=att_var[2]);
    v.var[:] = data[2][:,:,nt+2,1:12];

    v = defVar(ds,"domain_monthly_maximum_diurnalvariance",Int16,
               ("longitude","latitude","month"),attrib=att_var[2]);
    v.var[:] = data[3][:,:,nt+2,1:12];

    v = defVar(ds,"domain_monthly_minimum_diurnalvariance",Int16,
               ("longitude","latitude","month"),attrib=att_var[2]);
    v.var[:] = data[4][:,:,nt+2,1:12];


    v = defVar(ds,"zonalavg_yearly_mean_climatology",Int16,
               ("latitude",),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[1][:,:,nt+1,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_std_climatology",Int16,
               ("latitude",),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[2][:,:,nt+1,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_maximum_climatology",Int16,
               ("latitude",),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[3][:,:,nt+1,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_minimum_climatology",Int16,
               ("latitude",),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[4][:,:,nt+1,end],dims=1),dims=1);


    v = defVar(ds,"zonalavg_yearly_mean_hourly",Int16,
               ("latitude","hour"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[1][:,:,1:nt,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_std_hourly",Int16,
               ("latitude","hour"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[2][:,:,1:nt,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_maximum_hourly",Int16,
               ("latitude","hour"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[3][:,:,1:nt,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_minimum_hourly",Int16,
               ("latitude","hour"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[4][:,:,1:nt,end],dims=1),dims=1);


    v = defVar(ds,"zonalavg_yearly_mean_diurnalvariance",Int16,
               ("latitude",),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[1][:,:,nt+2,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_std_diurnalvariance",Int16,
               ("latitude",),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[2][:,:,nt+2,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_maximum_diurnalvariance",Int16,
               ("latitude",),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[3][:,:,nt+2,end],dims=1),dims=1);

    v = defVar(ds,"zonalavg_yearly_minimum_diurnalvariance",Int16,
               ("latitude",),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[4][:,:,nt+2,end],dims=1),dims=1);


    v = defVar(ds,"zonalavg_monthly_mean_climatology",Int16,
               ("latitude","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[1][:,:,nt+1,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_std_climatology",Int16,
               ("latitude","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[2][:,:,nt+1,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_maximum_climatology",Int16,
               ("latitude","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[3][:,:,nt+1,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_minimum_climatology",Int16,
               ("latitude","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[4][:,:,nt+1,1:12],dims=1),dims=1);


    v = defVar(ds,"zonalavg_monthly_mean_hourly",Int16,
               ("latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[1][:,:,1:nt,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_std_hourly",Int16,
               ("latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[2][:,:,1:nt,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_maximum_hourly",Int16,
               ("latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[3][:,:,1:nt,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_minimum_hourly",Int16,
               ("latitude","hour","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[4][:,:,1:nt,1:12],dims=1),dims=1);


    v = defVar(ds,"zonalavg_monthly_mean_diurnalvariance",Int16,
               ("latitude","month"),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[1][:,:,nt+2,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_std_diurnalvariance",Int16,
               ("latitude","month"),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[2][:,:,nt+2,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_maximum_diurnalvariance",Int16,
               ("latitude","month"),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[3][:,:,nt+2,1:12],dims=1),dims=1);

    v = defVar(ds,"zonalavg_monthly_minimum_diurnalvariance",Int16,
               ("latitude","month"),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[4][:,:,nt+2,1:12],dims=1),dims=1);


    v = defVar(ds,"meridionalavg_yearly_mean_climatology",Int16,
               ("longitude",),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[1][:,:,nt+1,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_std_climatology",Int16,
               ("longitude",),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[2][:,:,nt+1,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_maximum_climatology",Int16,
               ("longitude",),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[3][:,:,nt+1,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_minimum_climatology",Int16,
               ("longitude",),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[4][:,:,nt+1,end],dims=2),dims=2);


    v = defVar(ds,"meridionalavg_yearly_mean_hourly",Int16,
               ("longitude","hour"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[1][:,:,1:nt,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_std_hourly",Int16,
               ("longitude","hour"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[2][:,:,1:nt,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_maximum_hourly",Int16,
               ("longitude","hour"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[3][:,:,1:nt,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_minimum_hourly",Int16,
               ("longitude","hour"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[4][:,:,1:nt,end],dims=2),dims=2);


    v = defVar(ds,"meridionalavg_yearly_mean_diurnalvariance",Int16,
               ("longitude",),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[1][:,:,nt+2,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_std_diurnalvariance",Int16,
               ("longitude",),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[2][:,:,nt+2,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_maximum_diurnalvariance",Int16,
               ("longitude",),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[3][:,:,nt+2,end],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_yearly_minimum_diurnalvariance",Int16,
               ("longitude",),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[4][:,:,nt+2,end],dims=2),dims=2);


    v = defVar(ds,"meridionalavg_monthly_mean_climatology",Int16,
               ("longitude","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[1][:,:,nt+1,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_std_climatology",Int16,
               ("longitude","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[2][:,:,nt+1,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_maximum_climatology",Int16,
               ("longitude","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[3][:,:,nt+1,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_minimum_climatology",Int16,
               ("longitude","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[4][:,:,nt+1,1:12],dims=2),dims=2);


    v = defVar(ds,"meridionalavg_monthly_mean_hourly",Int16,
               ("longitude","hour","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[1][:,:,1:nt,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_std_hourly",Int16,
               ("longitude","hour","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[2][:,:,1:nt,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_maximum_hourly",Int16,
               ("longitude","hour","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[3][:,:,1:nt,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_minimum_hourly",Int16,
               ("longitude","hour","month"),attrib=att_var[1]);
    v.var[:] = dropdims(mean(data[4][:,:,1:nt,1:12],dims=2),dims=2);


    v = defVar(ds,"meridionalavg_monthly_mean_diurnalvariance",Int16,
               ("longitude","month"),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[1][:,:,nt+2,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_std_diurnalvariance",Int16,
               ("longitude","month"),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[2][:,:,nt+2,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_maximum_diurnalvariance",Int16,
               ("longitude","month"),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[3][:,:,nt+2,1:12],dims=2),dims=2);

    v = defVar(ds,"meridionalavg_monthly_minimum_diurnalvariance",Int16,
               ("longitude","month"),attrib=att_var[2]);
    v.var[:] = dropdims(mean(data[4][:,:,nt+2,1:12],dims=2),dims=2);

    close(ds);

    @info "$(Dates.now()) - Analysed $(info["source"]) $(info["product"]) data for the year $yr in $(regionfullname(region)) has been saved into file $(fnc) and moved to the data directory $(fol)."

end
